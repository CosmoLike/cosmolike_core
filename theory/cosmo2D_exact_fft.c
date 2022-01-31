void f_chi_for_Psi_cl(double* chi_ar, int Nchi, double* f_chi_ar, int ni);
void f_chi_for_Psi_cl_RSD(double* chi_ar, int Nchi, double* f_chi_RSD_ar, int ni);
void f_chi_for_Psi_cl_Mag(double* chi_ar, int Nchi, double* f_chi_Mag_ar, int ni);
void C_cl_mixed(int L, int LMAX, int ni, int nj, double *Cl, double dev, double tolerance);
double w_tomo_nonLimber(int nt, int ni, int nj); //w(theta) including non-Limber+RSD


void f_chi_for_Psi_sh(double* chi_ar, int Nchi, double* f_chi_ar, int nz);
void f_chi_for_Psi_sh_IA(double* chi_ar, int Nchi, double* f_chi_IA_ar, int nz);
void C_gl_mixed(int L, int LMAX, int ni, int nj, double *Cl, double dev, double tolerance);
double w_gamma_t_nonLimber(int nt, int ni, int nj);

//#include "pt.c"

double W_kappa(double a, double fK, double nz);//complete lens efficiency weight
double W_source(double a, double nz); //source redshift distribution (radial weight for IA,source clustering)
double W_gal(double a, double nz); //complete weight for galaxy statistics
double W_HOD(double a, double nz); //galaxy weigth without bias factor (for projecting P_gg instead of P_nl)

double C_cl_tomo(double l, int ni, int nj);  //galaxy clustering power spectrum of galaxies in bins ni, nj
double C_cl_tomo_nointerp(double l, int ni, int nj);
//double C_cl_HOD(double l, int ni);  //galaxy clustering power spectrum of galaxies in bin ni, using HOD model

//double C_gl_tomo(double l, int ni, int nj);  //G-G lensing power spectrum, lens bin ni, source bin nj
//double C_gl_tomo_nointerp(double l, int ni, int nj);
//double C_gl_HOD_tomo(double l, int ni, int nj);  //G-G lensing power spectrum from HOD model, lens bin ni, source bin nj

//double C_shear_tomo(double l, int ni, int nj); //shear tomography power spectra
//double C_shear_tomo_nointerp(double l, int ni, int nj);
/**********************************************************/


double MG_Sigma(double a)
{
//  double aa=a*a;
//  double omegam=cosmology.Omega_m/(aa*a);
  double omegav=omv_vareos(a);
  double hub=hoverh0(a);
  hub = hub*hub;
 
  return cosmology.MGSigma*omegav/hub/cosmology.Omega_v;
}

double dchi_da(double a){
  return 1./(a*a*hoverh0(a));
}
double W_kappa(double a, double fK, double nz){
  double wkappa = 1.5*cosmology.Omega_m*fK/a*g_tomo(a,(int)nz);
  if(cosmology.MGSigma != 0.){
    wkappa *= (1.+MG_Sigma(a));
  }
  return wkappa;
}

double W_mag(double a, double fK, double nz){
  double wmag = 1.5*cosmology.Omega_m*fK/a*g_lens(a,(int)nz);
  if(cosmology.MGSigma != 0.){
    wmag *= (1.+MG_Sigma(a));
  }
  return wmag;
}

double W_gal(double a, double nz){
  double wgal = gbias.b1_function(1./a-1.,(int)nz)*pf_photoz(1./a-1.,(int)nz)*hoverh0(a);
  double wmag = gbias.b_mag[(int)nz]*1.5*cosmology.Omega_m*f_K(chi(a))/a*g_lens(a,(int)nz);
  if(cosmology.MGSigma != 0.){
    wmag *= (1.+MG_Sigma(a));
  }
  return wgal + wmag;
  // return wgal;
}

double W_source(double a, double nz){
  return zdistr_photoz(1./a-1.,(int)nz)*hoverh0(a);
}

double f_rsd (double aa){
  double gamma = 0.55;
  return pow(cosmology.Omega_m /(cosmology.Omega_m +omv_vareos(aa) *aa*aa*aa),gamma);
}
double W_RSD(double l, double a0, double a1, double nz){
  double w = (1+8.*l)/((2*l+1)*(2*l+1))*pf_photoz(1./a0-1.,(int)nz)*hoverh0(a0)*f_rsd(a0);
  w -= 4./(2*l+3)*sqrt((2*l+1.)/(2*l+3.))*pf_photoz(1./a1-1.,(int)nz)*hoverh0(a1)*f_rsd(a1);
  return w;
}

double W_HOD(double a, double nz){
  return pf_photoz(1./a-1.,(int)nz)*hoverh0(a);
}

/*********** Limber approximation integrands for angular power spectra *************/
// use W_i W_j/fK^2 dchi_da(a) as Limber weight for power spectra

double int_for_C_cl_tomo_b2(double a, void *params)
{
  double res,ell, fK, k;
  double *ar = (double *) params;
  double b1 = gbias.b1_function(1./a-1.,(int)ar[0]);
  double b2 = gbias.b2[(int)ar[0]];
  double bs2 = gbias.bs2[(int)ar[0]];
  double g4 = pow(growfac(a)/growfac(1.0),4.);

  if (a >= 1.0) error("a>=1 in int_for_C_cl_tomo");
  
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  
  double s4 = 0.;//PT_sigma4(k);
  double p = Pdelta(k,a);
  double p_c = Pdelta_cluster(k,a);
  res=W_HOD(a,ar[0])*W_HOD(a,ar[1])*dchi_da(a)/fK/fK;
  if(res){
    res= res*(b1*b1*p_c+g4*(b1*b2*PT_d1d2(k)+0.25*b2*b2*(PT_d2d2(k)-2.*s4)+b1*bs2*PT_d1s2(k)+0.5*b2*bs2*(PT_d2s2(k)-4./3.*s4)+.25*bs2*bs2*(PT_s2s2(k)-8./9.*s4)+b1*b3nl_from_b1(b1)*PT_d1d3(k)));
  }
  double w_gal_0 = b1*W_HOD(a, ar[0]) * sqrt(p_c) + gbias.b_mag[(int)ar[0]]*W_mag(a, fK,ar[0])*sqrt(p);
  double w_gal_1 = b1*W_HOD(a, ar[1]) * sqrt(p_c) + gbias.b_mag[(int)ar[1]]*W_mag(a, fK,ar[1])*sqrt(p);
  res += (w_gal_0*W_mag(a, fK,ar[1])+ w_gal_1*W_mag(a, fK,ar[0]))*dchi_da(a)/fK/fK*sqrt(p);

  if (gbias.neutrino_induced_sdb){

	  res=W_HOD(a,ar[0])*W_HOD(a,ar[1])*dchi_da(a)/fK/fK;
	  double f_cb = 1.0-cosmology.Omega_nu/cosmology.Omega_m;
	  double b1_k_0 = gbias.b1_function(1./a-1.,(int)ar[0])* (1.0 + p_lin_cluster(k,a)/p_lin(k,a) * f_cb)/(1.0+f_cb);
	  double b1_k_1 = gbias.b1_function(1./a-1.,(int)ar[1])* (1.0 + p_lin_cluster(k,a)/p_lin(k,a) * f_cb)/(1.0+f_cb);
	  if(res){
	    res= res*(b1_k_0*b1_k_0*p_c+g4*(b1_k_0*b2*PT_d1d2(k)+0.25*b2*b2*(PT_d2d2(k)-2.*s4)+b1_k_0*bs2*PT_d1s2(k)+0.5*b2*bs2*(PT_d2s2(k)-4./3.*s4)+.25*bs2*bs2*(PT_s2s2(k)-8./9.*s4)+b1_k_0*b3nl_from_b1(b1_k_0)*PT_d1d3(k)));
	  }

	  w_gal_0 = b1_k_0*W_HOD(a, ar[0]) * sqrt(p_c) + gbias.b_mag[(int)ar[0]]*W_mag(a, fK,ar[0])*sqrt(p);
	  w_gal_1 = b1_k_1*W_HOD(a, ar[1]) * sqrt(p_c) + gbias.b_mag[(int)ar[1]]*W_mag(a, fK,ar[1])*sqrt(p);

	  res += (w_gal_0*W_mag(a, fK,ar[1])+ w_gal_1*W_mag(a, fK,ar[0]))*dchi_da(a)/fK/fK*sqrt(p);
  }


  return res;
}

double int_for_C_cl_tomo(double a, void *params)
{
  double res,ell, fK, k;
  double *ar = (double *) params;
  if (a >= 1.0) error("a>=1 in int_for_C_cl_tomo");
  
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  double p_c = Pdelta_cluster(k,a);
  double p = Pdelta(k,a);
  res = (gbias.b1_function(1./a-1.,(int)ar[0])*W_HOD(a, ar[0])*sqrt(p_c)+gbias.b_mag[(int)ar[0]]*W_mag(a, fK, ar[0])*sqrt(p));
  res *=(gbias.b1_function(1./a-1.,(int)ar[1])*W_HOD(a, ar[1])*sqrt(p_c)+gbias.b_mag[(int)ar[1]]*W_mag(a, fK, ar[1])*sqrt(p));
  res *= dchi_da(a)/fK/fK;

  if (gbias.neutrino_induced_sdb){


	  double f_cb = 1.0-cosmology.Omega_nu/cosmology.Omega_m;
	  double b1_k_0 = gbias.b1_function(1./a-1.,(int)ar[0])* (1.0 + p_lin_cluster(k,a)/p_lin(k,a) * f_cb)/(1.0+f_cb);//changed this line from ar[0] -> a in p_lin terms
	  double b1_k_1 = gbias.b1_function(1./a-1.,(int)ar[1])* (1.0 + p_lin_cluster(k,a)/p_lin(k,a) * f_cb)/(1.0+f_cb);//changed this line from ar[0] -> a in p_lin terms

	  res = (b1_k_0*W_HOD(a, ar[0])*sqrt(p_c)+gbias.b_mag[(int)ar[0]]*W_mag(a, fK, ar[0])*sqrt(p));
	  res *=(b1_k_1*W_HOD(a, ar[1])*sqrt(p_c)+gbias.b_mag[(int)ar[1]]*W_mag(a, fK, ar[1])*sqrt(p));
	  res *= dchi_da(a)/fK/fK;
  }
  return res;


}


double int_for_C_gl_tomo_b2(double a, void *params)
{
  double res,ell, fK, k;
  double *ar = (double *) params;
  double b1 = gbias.b1_function(1./a-1.,(int)ar[0]);
  double b2 = gbias.b2[(int)ar[0]];
  double bs2 = gbias.bs2[(int)ar[0]];
  double g4 = pow(growfac(a)/growfac(1.0),4.);
  if (a >= 1.0) error("a>=1 in int_for_C_cl_tomo");

  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;

  double p_c = Pdelta_cluster(k,a);
  double p = Pdelta(k,a);
  res= W_HOD(a,ar[0])*W_kappa(a,fK,ar[1])*dchi_da(a)/fK/fK;
  res= res*(b1*sqrt(p)*sqrt(p_c)+g4*(0.5*b2*PT_d1d2(k)+0.5*bs2*PT_d1s2(k)+0.5*b3nl_from_b1(b1)*PT_d1d3(k)));
  res += W_mag(a,fK,ar[0])*W_kappa(a,fK,ar[1])*dchi_da(a)/fK/fK*b1*p;


  if (gbias.neutrino_induced_sdb){
	  double f_cb = 1.0-cosmology.Omega_nu/cosmology.Omega_m;
	  double b1_k = gbias.b1_function(1./a-1.,(int)ar[0])* (1.0 + p_lin_cluster(k,a)/p_lin(k,a) * f_cb)/(1.0+f_cb);
	  res= W_HOD(a,ar[0])*W_kappa(a,fK,ar[1])*dchi_da(a)/fK/fK;
	  res= res*(b1_k*sqrt(p)*sqrt(p_c)+g4*(0.5*b2*PT_d1d2(k)+0.5*bs2*PT_d1s2(k)+0.5*b3nl_from_b1(b1_k)*PT_d1d3(k))); 
	  res += W_mag(a,fK,ar[0])*W_kappa(a,fK,ar[1])*dchi_da(a)/fK/fK*b1*p;
  }
  return res;
}





/*********** angular power spectra - without look-up tables ******************/

double C_cl_tomo_nointerp(double l, int ni, int nj)  //galaxy clustering power spectrum of galaxy bins ni, nj
{ static int init =-1;
  double array[3] = {1.0*ni,1.0*nj,l};
  if (gbias.b2[ni] || gbias.b2[nj]){
    if (ni != nj){
        if (init ==-1){
          printf("\nCalled C_cl(l,z1=%d,z2=%d) with non-linear bias parameters set.\n",ni,nj);
          printf("Cross-clustering beyond linear bias for cross-tomography bins not yet supported.\n");
          printf("Use linear bias only for z1!=z2 clustering\n\n");
          init = 1;}
          return int_gsl_integrate_medium_precision(int_for_C_cl_tomo,(void*)array,fmax(amin_lens(ni),amin_lens(nj)),fmin(amax_lens(ni),amax_lens(nj)),NULL,1000);
        }
    return int_gsl_integrate_medium_precision(int_for_C_cl_tomo_b2,(void*)array,amin_lens(ni),amax_lens(ni),NULL,1000);
  }
  else if (ni == nj){
    // return int_gsl_integrate_medium_precision(int_for_C_cl_tomo,(void*)array,amin_lens(ni),amax_lens(ni),NULL,1000);
    return int_gsl_integrate_medium_precision(int_for_C_cl_tomo,(void*)array,amin_lens(ni),0.999999,NULL,1000);
  }
  // return int_gsl_integrate_medium_precision(int_for_C_cl_tomo,(void*)array,fmax(amin_lens(ni),amin_lens(nj)),fmin(amax_lens(ni),amax_lens(nj)),NULL,1000);
  return int_gsl_integrate_medium_precision(int_for_C_cl_tomo,(void*)array,amin_lens(nj),0.99999,NULL,1000); // zi<=zj
}




/*********** angular power spectra - with look-up tables ******************/

double C_cl_tomo(double l, int ni, int nj)  //galaxy clustering power spectrum of galaxies in bins ni, nj
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  if (ni < 0 || ni >= tomo.clustering_Nbin ||nj < 0 || nj >= tomo.clustering_Nbin){
    printf("C_cl_tomo(l,%d,%d) outside tomo.clustering_Nbin range\nEXIT\n",ni,nj); exit(1);
  }

  if (recompute_clustering(C,G,N,ni,nj))
  {
    if (table==0) {
      table   = create_double_matrix(0, tomo.clustering_Nbin*tomo.clustering_Nbin-1, 0, Ntable.N_ell-1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell);
    }
    for (int i = 0; i < tomo.clustering_Nbin*tomo.clustering_Nbin; i++){table[i][0] = 123456789.0;}
    update_cosmopara(&C); update_nuisance(&N); update_galpara(&G);
  }
  int j = ni*tomo.clustering_Nbin+nj;
  if(table[j][0] > 123456780.0){ //still need to recompute this tomography bin combination
    //printf("recompute C_cl_tomo, %d, %d\n", ni, nj) ;
    double llog = logsmin;
    double result;
    for (int i=0; i<Ntable.N_ell; i++, llog+=ds) {

//      table[j][i]= log(C_cl_RSD_nointerp(exp(llog),ni,nj));
      result = C_cl_tomo_nointerp(exp(llog),ni,nj);
      if(result<=0) table[j][i] = -100; 
      else table[j][i] = log(result);
      table[nj*tomo.clustering_Nbin+ni][i]=table[j][i];
    }
  }
  double f1 = exp(interpol_fitslope(table[j], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1.));
  if (isnan(f1)){f1 = 0.;}
  return f1;
}


/////////////////////
double f_growth(double z){
	double aa = 1./(1+z);
	double gamma = 0.55;
	return pow(cosmology.Omega_m /(cosmology.Omega_m +omv_vareos(aa) *aa*aa*aa),gamma);
}


double int_for_C_cl_lin(double a, void *params)
{
	int status;
	double res,ell, fK, k;
	double *ar = (double *) params;
	ell       = ar[2]+0.5;
	fK     = f_K(chi(a));
	k      = ell/fK;

	double p_c = p_lin_cluster(k,a);
	double p = p_lin(k,a);
	res = (gbias.b1_function(1./a-1.,(int)ar[0])*W_HOD(a, ar[0])*sqrt(p_c)+gbias.b_mag[(int)ar[0]]*W_mag(a, fK, ar[0])*sqrt(p));
	res *=(gbias.b1_function(1./a-1.,(int)ar[1])*W_HOD(a, ar[1])*sqrt(p_c)+gbias.b_mag[(int)ar[1]]*W_mag(a, fK, ar[1])*sqrt(p));
	res *= dchi_da(a)/fK/fK;
	
	if (gbias.neutrino_induced_sdb){



	  	double f_cb = 1.0-cosmology.Omega_nu/cosmology.Omega_m;
	  	double b1_k_0 = gbias.b1_function(1./a-1.,(int)ar[0])* (1.0 + p_c/p * f_cb)/(1.0+f_cb);
	  	double b1_k_1 = gbias.b1_function(1./a-1.,(int)ar[1])* (1.0 + p_c/p * f_cb)/(1.0+f_cb);
	  	res = (b1_k_0*W_HOD(a, ar[0])*sqrt(p_c)+gbias.b_mag[(int)ar[0]]*W_mag(a, fK, ar[0])*sqrt(p));
	  	res *=(b1_k_1*W_HOD(a, ar[1])*sqrt(p_c)+gbias.b_mag[(int)ar[1]]*W_mag(a, fK, ar[1])*sqrt(p));
	  	res *= dchi_da(a)/fK/fK;
  	}
  	//uncomment above lines to implement scale-dependent neutrino bias
	return res;
}

double C_cl_lin_nointerp(double l, int ni, int nj)  //galaxy clustering power spectrum of galaxy bins ni, nj
{
	double array[3] = {1.0*ni,1.0*nj,l};
	// return int_gsl_integrate_medium_precision(int_for_C_cl_lin,(void*)array,fmax(amin_lens(ni),amin_lens(nj)),fmin(amax_lens(ni),amax_lens(nj)),NULL,1000);
	return int_gsl_integrate_medium_precision(int_for_C_cl_lin,(void*)array,fmax(amin_lens(ni),amin_lens(nj)),0.99999,NULL,1000);
}


/////// Integrand for galaxy density
void f_chi_for_Psi_cl(double* chi_ar, int Nchi, double* f_chi_ar, int ni){
	double g0 =1./growfac(1.);//(1+0.2));
	double a, z;
	int i;
	double real_coverH0 = cosmology.coverH0 / cosmology.h0; // unit Mpc
	double pf;
	for(i=0;i<Nchi;i++) {
		a = a_chi(chi_ar[i] / real_coverH0) ; // first convert unit of chi from Mpc to c/H0
		z = 1./a - 1.;
		if( (z<tomo.clustering_zmin[ni]) || (z>tomo.clustering_zmax[ni]) )
		{
			f_chi_ar[i] = 0.;
		}
		else
		{
			pf = (pf_photoz(z,ni)<0.)? 0:pf_photoz(z,ni); // get rid of unphysical negatives
			f_chi_ar[i] = chi_ar[i] * pf*growfac(a)*g0*gbias.b1_function(z,ni)*hoverh0(a)/real_coverH0;
		}
	}
}

// Integrand for galaxy density RSD
void f_chi_for_Psi_cl_RSD(double* chi_ar, int Nchi, double* f_chi_RSD_ar, int ni){
	double g0 =1./growfac(1.);//(1+0.2));
	double a, z;
	int i;
	double real_coverH0 = cosmology.coverH0 / cosmology.h0;
	double pf;
	for(i=0;i<Nchi;i++) {
		a = a_chi(chi_ar[i] / real_coverH0) ; // first convert unit of chi from Mpc to c/H0
		z = 1./a - 1.;
		if( (z<tomo.clustering_zmin[ni]) || (z>tomo.clustering_zmax[ni]) )
		{
			f_chi_RSD_ar[i] = 0.;
		}
		else
		{
			pf = (pf_photoz(z,ni)<0.)? 0:pf_photoz(z,ni); // get rid of unphysical negatives
			f_chi_RSD_ar[i] = -chi_ar[i] * pf*growfac(a)*g0*f_growth(z)*hoverh0(a)/real_coverH0;
		}
	}
}

// Integrand for lensing magnification of galaxy density
void f_chi_for_Psi_cl_Mag(double* chi_ar, int Nchi, double* f_chi_Mag_ar, int ni){
	double g0 =1./growfac(1.);//(1+0.2));
	double a, z, fK;
	int i;
	double real_coverH0 = cosmology.coverH0 / cosmology.h0;
	double window_M, wmag;
	for(i=0;i<Nchi;i++) {
		a = a_chi(chi_ar[i] / real_coverH0) ; // first convert unit of chi from Mpc to c/H0
		z = 1./a - 1.;
		fK = f_K(chi_ar[i]/real_coverH0);

		if( (z>tomo.clustering_zmax[ni]) )
		{
			f_chi_Mag_ar[i] = 0.;
		}
		else
		{
			//printf("Here! a, fK, ni: %lg,%lg,%d\n", a, fK, ni);
			wmag = W_mag(a, fK, (double)ni);
			window_M = wmag/ fK / (real_coverH0*real_coverH0);
			//printf("bmag, wkappa, f_K, real_coverH0, %lg %lg %lg %lg\n", gbias.b_mag[ni], wkappa, fK,real_coverH0);
			// pf = (pf_photoz(z,ni)<0.)? 0:pf_photoz(z,ni); // get rid of unphysical negatives
			// f_chi_Mag_ar[i] = chi_ar[i]/a * window_M*growfac(a)*g0;
			f_chi_Mag_ar[i] = window_M*growfac(a)*g0; // unit [Mpc^-2]
		}
		// printf("%lg\n", f_chi_Mag_ar[i]);
	}
}

// Mixture of non-Limber and Limber of C_cl (galaxy clustering)
void C_cl_mixed(int L, int LMAX, int ni, int nj, double *Cl, double dev, double tolerance) {
	// ni = 4;
	// nj = 4;
	int status;
	int i,j,i_block;
	long l;
	// run 100 ells at a time, and see if switching to Limber is needed.
	// Save runtime for Limber, and save re-creation time of fftw_plan.
	int Nell_block = 100, Nchi = 1000;
	int ell_ar[Nell_block];
	double **k1_ar, **k2_ar, **Fk1_ar, **Fk2_ar, **Fk1_RSD_ar, **Fk2_RSD_ar;
	double **Fk1_Mag_ar, **Fk2_Mag_ar;

	k1_ar = malloc(Nell_block * sizeof(double *));
	k2_ar = malloc(Nell_block * sizeof(double *));
	Fk1_ar = malloc(Nell_block * sizeof(double *));
	Fk2_ar = malloc(Nell_block * sizeof(double *));

	Fk1_RSD_ar = malloc(Nell_block * sizeof(double *));
	Fk2_RSD_ar = malloc(Nell_block * sizeof(double *));

	Fk1_Mag_ar = malloc(Nell_block * sizeof(double *));
	Fk2_Mag_ar = malloc(Nell_block * sizeof(double *));
	for(i=0;i<Nell_block;i++) {
		k1_ar[i] = malloc(Nchi * sizeof(double));
		k2_ar[i] = malloc(Nchi * sizeof(double));
		Fk1_ar[i] = malloc(Nchi * sizeof(double));
		Fk2_ar[i] = malloc(Nchi * sizeof(double));
		Fk1_RSD_ar[i] = malloc(Nchi * sizeof(double));
		Fk2_RSD_ar[i] = malloc(Nchi * sizeof(double));
		Fk1_Mag_ar[i] = malloc(Nchi * sizeof(double));
		Fk2_Mag_ar[i] = malloc(Nchi * sizeof(double));
		for(j=0;j<Nchi;j++) {
			Fk1_ar[i][j] = 0.;
			Fk2_ar[i][j] = 0.;
			Fk1_RSD_ar[i][j] = 0.;
			Fk2_RSD_ar[i][j] = 0.;
			Fk1_Mag_ar[i][j] = 0.;
			Fk2_Mag_ar[i][j] = 0.;
		}
	}

	double chi_ar[Nchi], f1_chi_ar[Nchi], f2_chi_ar[Nchi];
	double f1_chi_RSD_ar[Nchi], f2_chi_RSD_ar[Nchi];
	double f1_chi_Mag_ar[Nchi], f2_chi_Mag_ar[Nchi];

	double real_coverH0 = cosmology.coverH0 / cosmology.h0;
	double chi_min = chi(1./(1.+0.002))*real_coverH0, chi_max = chi(1./(1.+4.))*real_coverH0;
	double dlnchi = log(chi_max/chi_min) / (Nchi - 1.);
	double dlnk = dlnchi;

	for(i=0; i<Nchi; i++) {
		chi_ar[i] = chi_min * exp(dlnchi*i);
	}

	f_chi_for_Psi_cl(chi_ar, Nchi, f1_chi_ar, ni);
	if(ni != nj) {f_chi_for_Psi_cl(chi_ar, Nchi, f2_chi_ar, nj);}

	f_chi_for_Psi_cl_RSD(chi_ar, Nchi, f1_chi_RSD_ar, ni);
	if(ni != nj) {f_chi_for_Psi_cl_RSD(chi_ar, Nchi, f2_chi_RSD_ar, nj);}

	f_chi_for_Psi_cl_Mag(chi_ar, Nchi, f1_chi_Mag_ar, ni);
	if(ni != nj) {f_chi_for_Psi_cl_Mag(chi_ar, Nchi, f2_chi_Mag_ar, nj);}

	i_block = 0;
	double cl_temp;

	config my_config, my_config_RSD, my_config_Mag;
	my_config.nu = 1.;
	my_config.c_window_width = 0.25;
	my_config.derivative = 0;
	my_config.N_pad = 200;
	my_config.N_extrap_low = 0;
	my_config.N_extrap_high = 0;

	my_config_RSD.nu = 1.01;
	my_config_RSD.c_window_width = 0.25;
	my_config_RSD.derivative = 2;
	my_config_RSD.N_pad = 500;
	my_config_RSD.N_extrap_low = 0;
	my_config_RSD.N_extrap_high = 0;

	my_config_Mag.nu = 1.;
	my_config_Mag.c_window_width = 0.25;
	my_config_Mag.derivative = 0;
	my_config_Mag.N_pad = 500;
	my_config_Mag.N_extrap_low = 0;
	my_config_Mag.N_extrap_high = 0;

	double ell_prefactor;
	double k1_cH0;

	while (fabs(dev) > tolerance){
	// while(0){
	// while (L<100){
		//Cl[L] = C_cl_RSD(L,nz,nz);
		for(i=0;i<Nell_block;i++) {ell_ar[i]=i+i_block*Nell_block;}

		cfftlog_ells(chi_ar, f1_chi_ar, Nchi, &my_config, ell_ar, Nell_block, k1_ar, Fk1_ar);
		if(ni != nj) {cfftlog_ells(chi_ar, f2_chi_ar, Nchi, &my_config, ell_ar, Nell_block, k2_ar, Fk2_ar);}
		double f_cb = 1.0-cosmology.Omega_nu/cosmology.Omega_m;

		if (gbias.neutrino_induced_sdb){
			double geff;
			//call to increment in this for-loop after re-weighting of f(\chi)
				//adjust f1_chi_ar for each chi in this for-loop with k-dep factor
			for(j=0;j<Nchi;j++) {
				//need to calculate a g_eff term for each representative redshift to approximate
				//k-dependence in the incremeneted integral. This is the square root of 'growth factors',
				//i.e. p_(k,z) = p_(k,z_bin)*geff(z,z_bin)^2, where geff^2= (G(z)/G(z_bin)*(k_dep_factor(k,z)/k_dep_factor(k,z_bin)))^2

				geff = 1.0;//f_cb * growfac_cluster(a_chi(chi_ar[j] / real_coverH0))/growfac(a_chi(chi_ar[j] / real_coverH0))*growfac(1.)/growfac_cluster(1.);
				
				f1_chi_ar[j] *= geff;
				if(ni != nj) {f2_chi_ar[j] *= geff;}
			}
			//cfftlog_ells_increment(chi_ar, f1_chi_ar, Nchi, &my_config, ell_ar, Nell_block, k1_ar, Fk1_ar);
			//if(ni != nj) {cfftlog_ells_increment(chi_ar, f2_chi_ar, Nchi, &my_config, ell_ar, Nell_block, k2_ar, Fk2_ar);}

		}

		cfftlog_ells(chi_ar, f1_chi_RSD_ar, Nchi, &my_config_RSD, ell_ar, Nell_block, k1_ar, Fk1_RSD_ar);
		if(ni != nj) {cfftlog_ells(chi_ar, f2_chi_RSD_ar, Nchi, &my_config_RSD, ell_ar, Nell_block, k2_ar, Fk2_RSD_ar);}

		// Add in lensing magnification contribution
		cfftlog_ells(chi_ar, f1_chi_Mag_ar, Nchi, &my_config_Mag, ell_ar, Nell_block, k1_ar, Fk1_Mag_ar);
		if(ni != nj) {cfftlog_ells(chi_ar, f2_chi_Mag_ar, Nchi, &my_config_Mag, ell_ar, Nell_block, k2_ar, Fk2_Mag_ar);}
		for(i=0;i<Nell_block;i++) {
			ell_prefactor = ell_ar[i]*(ell_ar[i]+1.);
			for(j=0;j<Nchi;j++) {
				Fk1_Mag_ar[i][j]= (ell_prefactor / (k1_ar[i][j]*k1_ar[i][j]) * (gbias.b_mag[ni]) *  Fk1_Mag_ar[i][j]);
				if(ni != nj) {Fk2_Mag_ar[i][j]= (ell_prefactor / (k2_ar[i][j]*k2_ar[i][j])* (gbias.b_mag[nj]) *  Fk2_Mag_ar[i][j]);}
			}
		}
		//printf("%f\n", zmean(ni));
		for(i=0;i<Nell_block;i++) {
			cl_temp = 0.;
			for(j=0;j<Nchi;j++) {
				// printf("k,Fk: %d,%d, %lf,%lf\n", i,j, k1_ar[i][j], Fk1_ar[i][j]);
				k1_cH0 = k1_ar[i][j] * real_coverH0;
				double plin = p_lin(k1_cH0,1.0);
				double plin_cluster = p_lin_cluster(k1_cH0,1.0);
				//double factor = plin_cluster*(1+2*f_cb*(plin_cluster/plin) + f_cb*f_cb*(plin_cluster/plin)*(plin_cluster/plin))/(1+f_cb)/(1+f_cb);
				double cluster_a = 1.0/(1.0+zmean(ni));
				double factor = plin_cluster * ((1 + f_cb * (p_lin_cluster(k1_cH0, cluster_a)/p_lin(k1_cH0, cluster_a)))/(1+f_cb))*((1 + f_cb * (p_lin_cluster(k1_cH0, cluster_a)/p_lin(k1_cH0, cluster_a)))/(1+f_cb));
				double cluster_a_2  = 1.0/(1.0 + zmean(nj));
				double factor_2 = plin_cluster * ((1 + f_cb * (p_lin_cluster(k1_cH0, cluster_a_2)/p_lin(k1_cH0, cluster_a_2)))/(1+f_cb))*((1 + f_cb * (p_lin_cluster(k1_cH0, cluster_a_2)/p_lin(k1_cH0, cluster_a_2)))/(1+f_cb));
				double fk1;
				double fk2;
				//printf("%f %f\n", factor/plin, cluster_a);
				if(ni == nj) {

					if (gbias.neutrino_induced_sdb){fk1 = (Fk1_RSD_ar[i][j] +Fk1_Mag_ar[i][j])* sqrt(plin) + Fk1_ar[i][j]*sqrt(factor);}
					else{fk1 = (Fk1_RSD_ar[i][j] +Fk1_Mag_ar[i][j])* sqrt(plin) + Fk1_ar[i][j]*sqrt(plin_cluster);}
					cl_temp += fk1*fk1*k1_cH0*k1_cH0*k1_cH0;
					//cl_temp += ((Fk1_RSD_ar[i][j]) * (Fk1_RSD_ar[i][j]) + (Fk1_Mag_ar[i][j])*(Fk1_Mag_ar[i][j]))*k1_cH0*k1_cH0*k1_cH0 *plin;
					//if (gbias.neutrino_induced_sdb){cl_temp += (Fk1_ar[i][j]) * (Fk1_ar[i][j]) *k1_cH0*k1_cH0*k1_cH0 *factor;}
					//else{cl_temp += (Fk1_ar[i][j]) * (Fk1_ar[i][j]) *k1_cH0*k1_cH0*k1_cH0 *plin_cluster;}
				}
				else {
					if (gbias.neutrino_induced_sdb){
						fk1 = (Fk1_RSD_ar[i][j] +Fk1_Mag_ar[i][j])* sqrt(plin) + Fk1_ar[i][j]*sqrt(factor);
						fk2 = (Fk2_RSD_ar[i][j] +Fk2_Mag_ar[i][j])* sqrt(plin) + Fk2_ar[i][j]*sqrt(factor_2);
					}
					else{
						fk1 = (Fk1_RSD_ar[i][j] +Fk1_Mag_ar[i][j])* sqrt(plin) + Fk1_ar[i][j]*sqrt(plin_cluster);
						fk2 = (Fk2_RSD_ar[i][j] +Fk2_Mag_ar[i][j])* sqrt(plin) + Fk2_ar[i][j]*sqrt(plin_cluster);						
					}
					cl_temp += fk1*fk2*k1_cH0*k1_cH0*k1_cH0;
					//cl_temp += ((Fk1_RSD_ar[i][j])*(Fk2_RSD_ar[i][j])+ (Fk1_Mag_ar[i][j])*(Fk2_Mag_ar[i][j])) *k1_cH0*k1_cH0*k1_cH0 *plin;
					//if (gbias.neutrino_induced_sdb){cl_temp += (Fk1_ar[i][j])*(Fk2_ar[i][j]) *k1_cH0*k1_cH0*k1_cH0 *factor;}
					//else{cl_temp += (Fk1_ar[i][j])*(Fk2_ar[i][j]) *k1_cH0*k1_cH0*k1_cH0 *plin_cluster;}
				}
			}
			Cl[ell_ar[i]] = cl_temp * dlnk * 2./M_PI + C_cl_tomo_nointerp(1.*ell_ar[i],ni,nj) - C_cl_lin_nointerp(1.*ell_ar[i],ni,nj);
		}

		i_block++;
		L = i_block*Nell_block -1 ;
		dev = Cl[L]/C_cl_tomo_nointerp((double)L,ni,nj)-1.;
	   // printf("ni,L,Cl[L],dev=%d %d %e %e\n",ni,L,Cl[L],dev);
		// printf("i_block: %d\n", i_block);
	}
	L++;
	//printf("switching to Limber calculation at l = %d %d\n",L, ni);
	// for (l = 1; l < 50; l++){
	// 	Cl[l]=C_cl_tomo_nointerp((double)l,ni,nj);
	// 	// fprintf(OUT, "%d %lg\n", l, Cl[l]);
	// }
	for (l = L; l < LMAX; l++){
		Cl[l]=C_cl_tomo((double)l,ni,nj);
		// fprintf(OUT, "%d %lg\n", l, Cl[l]);
	}
	// printf("finished bin %d\n", ni);
	for(i=0;i<Nell_block;i++) {
		free(k1_ar[i]);free(k2_ar[i]);
		free(Fk1_ar[i]);free(Fk2_ar[i]);
		free(Fk1_RSD_ar[i]);free(Fk2_RSD_ar[i]);
		free(Fk1_Mag_ar[i]);free(Fk2_Mag_ar[i]);
	}
	free(k1_ar);free(k2_ar);
	free(Fk1_ar);free(Fk2_ar);
	free(Fk1_RSD_ar);free(Fk2_RSD_ar);
	free(Fk1_Mag_ar);free(Fk2_Mag_ar);

}

double w_tomo_nonLimber(int nt, int ni, int nj){
	// if(1) return 0.;
	static int LMAX = 100000;
	static int NTHETA = 0;
	static double ** Pl =0;
	static double *Cl =0;
	static double *w_vec =0;
	static cosmopara C;
	static nuisancepara N;
	static galpara G;
	int i,l,nz;
	if (like.Ntheta ==0){
		printf("cosmo2D_real.c:w_tomo_exact: like.Ntheta not initialized\nEXIT\n"); exit(1);
	}
	if (ni != nj){
		printf("cosmo2D_real.c:w_tomo_exact: ni != nj tomography not supported\nEXIT\n"); exit(1);    
	}
	if (Pl ==0){
		Pl =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
		Cl = create_double_vector(0,LMAX-1);
		w_vec = create_double_vector(0,tomo.clustering_Nbin*like.Ntheta-1);
		NTHETA = like.Ntheta;
		double *xmin, *xmax, *Pmin, *Pmax;
		xmin= create_double_vector(0, like.Ntheta-1);
		xmax= create_double_vector(0, like.Ntheta-1);
		double logdt=(log(like.vtmax)-log(like.vtmin))/like.Ntheta;
		Pmin= create_double_vector(0, LMAX+1);
		Pmax= create_double_vector(0, LMAX+1);

		// double *xmid, *Pmid;
		double mythetamin, mythetamax;
		// xmid= create_double_vector(0, like.Ntheta-1);
		// Pmid= create_double_vector(0, LMAX+1);

		for(i=0; i<like.Ntheta ; i++){
			mythetamin = exp(log(like.vtmin)+(i+0.0)*logdt);
			mythetamax = exp(log(like.vtmin)+(i+1.0)*logdt);
			xmin[i]=cos(mythetamin);
			xmax[i]=cos(mythetamax);
			// xmid[i]= cos((2./3.) * (pow(mythetamax,3) - pow(mythetamin,3)) / (mythetamax*mythetamax - mythetamin*mythetamin));
		}

		for (i = 0; i<NTHETA; i ++){
			// printf("Tabulating Legendre coefficients %d/%d\n",i+1, NTHETA);
			gsl_sf_legendre_Pl_array(LMAX, xmin[i],Pmin);
			gsl_sf_legendre_Pl_array(LMAX, xmax[i],Pmax);
			// gsl_sf_legendre_Pl_array(LMAX, xmid[i],Pmid);

			Pl[i][0] = 1.0;
			for (int l = 1; l < LMAX; l ++){
				Pl[i][l] = 1./(4.*M_PI)*(Pmin[l+1]-Pmax[l+1]-Pmin[l-1]+Pmax[l-1])/(xmin[i]-xmax[i]);
				// Pl[i][l] = (2.*l+1)/(4.*M_PI)*Pmid[l];
				// printf("l,%ld\n", l);
			}
		}
		free_double_vector(xmin,0,like.Ntheta-1);
		free_double_vector(xmax,0,like.Ntheta-1);
		free_double_vector(Pmin,0,LMAX+1);
		free_double_vector(Pmax,0,LMAX+1);
		}
		if (recompute_clustering(C,G,N,ni,nj)){
		//required fractional accuracy in C(l)
		double tolerance= 0.01;
		//dev will be the actual difference between exact and Limber calcuation
		double dev;

		for (nz = 0; nz <tomo.clustering_Nbin; nz ++){
			int L = 0;
			// initialize to large value in order to start while loop
			dev=10.*tolerance;
			C_cl_mixed(L, LMAX, nz,nz, Cl, dev, tolerance);
			for (i = 0; i < NTHETA; i++){
				w_vec[nz*like.Ntheta+i] =0;
				for (l = 1; l < LMAX; l++){
					w_vec[nz*like.Ntheta+i]+=Pl[i][l]*Cl[l];
				}
			}
		}

		update_cosmopara(&C);
		update_galpara(&G);
		update_nuisance(&N);
	}
	return w_vec[ni*like.Ntheta+nt];  
}
