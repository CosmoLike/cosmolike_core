
// Flag function that throws an error if Horndeski parameters are enabled in 
// the given cosmology. Should be added to all functions that make assumptions 
// that are invalid for Horndeski models
void invalid_for_horndeski(const char* func_name){
    if (cosmology.use_horndeski == 1){
        fprintf(stderr, "Horndeski model is being used, but function '%s' makes assumptions that are invalid for Horndeski.\n", func_name);
        #ifdef EXIT_IF_INVALID_FOR_HORNDESKI
        exit(1);
        #endif
    }
}

double MG_Sigma(double a)
{
//  double aa=a*a;
//  double omegam=cosmology.Omega_m/(aa*a);
  double omegav=omv_vareos(a);
  double hub=hoverh0(a);
  hub = hub*hub;
 
  return cosmology.MGSigma*omegav/hub/cosmology.Omega_v;
}

double horndeski_sigma(double a){
  // Horndeski: Assumes quasi-static limit, as in Eq. 6 of arXiv:1705.04714
  // mu_light = 1 + Sigma (???)
  if (cosmology.use_horndeski != 1) return 0.;
  
  double sigma;
  int status = 0;
  sigma = horndeski_function(a, "mu_light", &status) - 1.; // FIXME: -1?
  
  //printf("Horndeski: sigma = %4.4e\n", sigma);
  //printf("\t%3.3e, %3.3e, %3.3e, %3.3e\n", cosmology.mg_alpha_xk, cosmology.mg_alpha_xb, cosmology.mg_alpha_xm, cosmology.mg_alpha_xt);
  
  if (status != 0){
    fprintf(stderr, "modgrav.c: horndeski_sigma() failed.\n");
    exit(1);
  }
  return sigma;
}

double horndeski_mu(double a){
  // Horndeski: Assumes quasi-static limit, as in Eq. 2 of arXiv:1705.04714
  // horndeski_mu = 1 + mg_mu (???)
  if (cosmology.use_horndeski != 1) return 0.;
  
  double mu;
  int status = 0;
  mu = horndeski_function(a, "mu_eff", &status); // FIXME: -1?
  
  if (status != 0){
    fprintf(stderr, "modgrav.c: horndeski_mu() failed.\n");
    exit(1);
  }
  return mu;
}


void horndeski_quasistatic(double a, double* Geff0, double* gamma0){
    // Quasi-static G_eff,0(a) and gamma_0(a), from Eq. A.1 of 1610.09290.
    double hub, oma, odea, wde, H;
    double Mstar2, fz, fzprime, beta1, beta2, beta3;
    
    // Test whether Horndeski parameters are set to GR values
    int is_gr = (   (cosmology.mg_alpha_xk == 0.) 
                 && (cosmology.mg_alpha_xb == 0.) 
                 && (cosmology.mg_alpha_xm == 0.)
                 && (cosmology.mg_alpha_xt == 0.)
                 && (cosmology.mg_alpha_M2 == 1.)
                );
    if (is_gr == 1){ *Geff0 = 1.; *gamma0 = 1.; return; }
    
    // Background functions
    hub = hoverh0(a); // Dim'less Hubble rate
    oma = cosmology.Omega_m / (a*a*a) / (hub*hub);
    odea = omv_vareos(a) / (hub*hub);
    wde = cosmology.w0 + cosmology.wa*(1. - a);
    H = (100. * cosmology.h0 / constants.lightspeed) * hub;
    
    // alpha function redshift scaling, f(z), where alpha_X = alpha_X,0 * f(z)
    fz = odea / cosmology.Omega_v;
    fzprime = -3.*wde * a*H * oma * odea / cosmology.Omega_v;
    
    // Time dependence of the Planck mass squared
    Mstar2 = horndeski_M2(a);
    
    // beta functions
    beta1 = -3. * oma / Mstar2 + cosmology.mg_alpha_xb*fzprime / (a * H)
          + 1.5*(2. - cosmology.mg_alpha_xb*fz) * (oma + odea*(1. + wde));
    beta2 = cosmology.mg_alpha_xb*fz * (1. + cosmology.mg_alpha_xt*fz) 
          + 2. * fz * (cosmology.mg_alpha_xm - cosmology.mg_alpha_xt);
    beta3 = (1. + cosmology.mg_alpha_xt*fz) * beta1 
          + (1. + cosmology.mg_alpha_xm*fz) * beta2;
    
    // Return modified growth and slip functions, G_0 and gamma_0
    *Geff0 = 2.*(beta1 + beta2) / Mstar2 
           / (2.*beta1 + (2. - cosmology.mg_alpha_xb*fz)*beta2);
    *gamma0 = beta3 / (beta1 + beta2);
}

