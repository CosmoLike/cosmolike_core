#!/usr/bin/python
import numpy as np
import pylab as P
from scipy.integrate import cumtrapz
from scipy.interpolate import interp1d

C = 299792.458 # km/s

omega_m = 0.3156
sigma_8 = 0.831
n_s = 0.9645
w0 = -1.0
wa = 0.0
omega_b = 0.0491685
h0 = 0.6727
MGSigma = 0.0
MGmu = 0.0
alpha_k = 1.0
alpha_b = 0.1
alpha_m = 0.0
alpha_t = 0.0
alpha_M2 = 1.0

def rho_de(a):
    """
    Dimensionless dark energy energy density.
    """
    return (1. - omega_m) * np.exp(3.*wa*(a - 1.)) / a**(3.*(1. + w0 + wa))

def E2(a):
    """
    Dimensionless expansion rate squared, from the Friedmann equation.
    """
    return omega_m / a**3. + rho_de(a)

def omega_x(a):
    """
    Fractional energy density as a function of scale factor.
    """
    _E2 = E2(a)
    oma = omega_m / a**3. / _E2
    odea = rho_de(a) / _E2
    return oma, odea

def interp_planck_mass_sq():
    """
    Interpolation function for the Planck mass-squared as a 
    function of scale factor.
    """
    a = np.logspace(-3., 0., 200)
    oma, odea = omega_x(a)
    integ = alpha_m * odea / (1. - omega_m)
    
    # Cumulative integral
    logM2 = cumtrapz( integ[::-1], 
                      np.log(a[::-1]), 
                      initial=np.log(alpha_M2) )[::-1]
    
    # Interpolation function (in log-space, but takes non-log arguments)
    interp_logM2 = interp1d(np.log(a), logM2, kind='linear', bounds_error=False)
    M2 = lambda aa: np.exp(interp_logM2(np.log(aa)))
    return M2

def mg_functions(a):
    """
    Modified gravity functions G_0 and gamma_0 (leading order terms in G_eff 
    and gamma) for modified growth and slip relations.
    
    Calculates the beta functions, beta_1, beta_2, beta_3, from Eqs. A10-12 of 
    arXiv:1610.09290.
    """
    # Density functions, DE eqn. of state, and time-dep. Planck mass
    ode0 = 1. - omega_m
    oma, odea = omega_x(a)
    wde = w0 + wa*(1. - a)
    Mstar2 = interp_planck_mass_sq()(a)
    H = (100. * h0 / C) * np.sqrt( E2(a) )
    
    # alpha function redshift scaling, alpha_X = alpha_X,0 . f(z)
    fz = odea / ode0
    fzprime = -3.*wde * a*H * oma * odea / ode0
    
    # beta functions
    beta1 = -3. * oma / Mstar2 + alpha_b*fzprime / (a * H) \
          + 1.5*(2. - alpha_b*fz)*(oma + odea*(1. + wde)) 
    beta2 = alpha_b*fz * (1. + alpha_t*fz) + 2. * fz * (alpha_m - alpha_t)
    beta3 = (1. + alpha_t*fz) * beta1 + (1. + alpha_m*fz) * beta2
    
    # Modified growth and slip functions, G_0 and gamma_0
    G0 = 2.*(beta1 + beta2) / Mstar2 / (2.*beta1 + (2. - alpha_b*fz)*beta2)
    gamma0 = beta3 / (beta1 + beta2)
    return G0, gamma0, fz, fzprime / (a*H), H


# Example plot
P.subplot(111)

a = np.logspace(-2., 0., 1000)
G0, gamma0, fz, fzp, H = mg_functions(a)
z = 1./a - 1.

# Numerical deriv
dfz_dtau = np.gradient(fz) / np.gradient(a) * a #*a*H #/aH

P.plot(z, G0, 'r-', lw=1.8, label="$G_0(a)$")
P.plot(z, gamma0, 'b-', lw=1.8, label="$\gamma_0(a)$")
P.plot(z, fz, 'm-', lw=1.8, label="$f(a)$")
P.plot(z, fzp, 'm--', lw=1.8, label="$f^\prime(a)$")
P.plot(z, dfz_dtau, 'g:', lw=1.8, label="$f^\prime(a)$ (numerical)")

P.xlim((0., 3.))
#P.ylim((0., 1.1))

P.legend(loc='upper right', frameon=False)
P.tight_layout()
P.show()
