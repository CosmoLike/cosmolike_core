double C_cl_non_Limber(int l, int ni, int nj); //includes RSD
double C_cl_RSD(int l, int ni, int nj); //C_cl_Limber + non-Limber RSD terms only
double w_tomo_nonLimber(int nt, int ni, int nj); //w(theta) including non-Limber+RSD
double w_gamma_t_nonLimber(int nt, int ni, int nj);
typedef double (*C_tomo_pointer)(double l, int n1, int n2);

double Psi_Mag(double k, int l,int ni);

double G_taper(double k){
  double s_bao = 5.5/cosmology.coverH0;
  return exp(-k*k*s_bao*s_bao);
}
double f_growth(double z){
	double aa = 1./(1+z);
	double gamma = 0.55;
	return pow(cosmology.Omega_m /(cosmology.Omega_m +omv_vareos(aa) *aa*aa*aa),gamma);
}

double Psi_RSD_z(double z,void *params){
  double *ar = (double *) params;
  int l = (int) ar[0];
  //if (l>50){return 0.;}
  double x = ar[1]*f_K(chi(1./(1.+z)));
  double WRSD =0.;
  if (x < 0.1*sqrt(2.*l) && (2*l+5) < GSL_SF_DOUBLEFACT_NMAX){ //small-x limit for j_l(x) in Eq. 27 of https://arxiv.org/pdf/astro-ph/0605302v2.pdf
    WRSD =(2.*l*l+2.*l-1.)/((2.*l+3.)*(2.*l-1))*pow(x,l)/gsl_sf_doublefact((unsigned int)(2*l+1))*(1-0.5*x*x/(2.*l+3));
    WRSD -= (l+1.)*(l+2.)/((2.*l+1.)*(2.*l+3.))*pow(x,l+2.)/gsl_sf_doublefact((unsigned int)(2*l+5))*(1-0.5*x*x/(2.*l+7));
    if (l>1){WRSD -= l*(l-1.)/((2.*l-1.)*(2.*l+1.))*pow(x,l-2.)/gsl_sf_doublefact((unsigned int)(2*l-3))*(1-0.5*x*x/(2.*l-1));}
  }
  else if (x > 100.*l*l){// large-x limit for j_l(x), to avoid gsl error
    WRSD =(2.*l*l+2.*l-1.)/((2.*l+3.)*(2.*l-1))* sin(x-l*M_PI*0.5)/x;
    WRSD -= (l+1.)*(l+2.)/((2.*l+1.)*(2.*l+3.))* sin(x-(l+2)*M_PI*0.5)/x;
    if (l>1){WRSD -= l*(l-1.)/((2.*l-1.)*(2.*l+1.))*sin(x-(l-2)*M_PI*0.5)/x;}    
  }
  else{ // Eq. 27 of https://arxiv.org/pdf/astro-ph/0605302v2.pdf
    WRSD =(2.*l*l+2.*l-1.)/((2.*l+3.)*(2.*l-1))*gsl_sf_bessel_jl(l,x);
    WRSD -= (l+1.)*(l+2.)/((2.*l+1.)*(2.*l+3.))*gsl_sf_bessel_jl(l+2,x);
    if (l>1){WRSD -= l*(l-1.)/((2.*l-1.)*(2.*l+1.))*gsl_sf_bessel_jl(l-2,x);}
  }
  return WRSD*f_growth(z);
}
double int_for_Psi_RSD(double z, void *params){
  static double g0 = 0.;
  if (g0 ==0){g0 =1./growfac(1.);}
  double *ar = (double *) params;
  return Psi_RSD_z(z,params)*pf_photoz(z,(int)ar[2])*growfac(1./(1+z))*g0;
}
double int_for_Psi_cl (double z, void *params){
  static double g0 = 0.;
  if (g0 ==0){g0 =1./growfac(1.);}
  double *ar = (double *) params;
  int l = (int) ar[0];
  double res, x = ar[1]*f_K(chi(1./(1.+z)));
  if (x < 0.1*sqrt(2.*l) && (2*l+1) < GSL_SF_DOUBLEFACT_NMAX){ // small-x limit for j_l(x), to avoid gsl underflow error
    res =pow(x,l)/gsl_sf_doublefact((unsigned int)(2*l+1))*(1-0.5*x*x/(2.*l+3));
  }
  else if (x > 100.*l*l){// large-x limit for j_l(x), to avoid gsl error
    res = sin(x-l*M_PI*0.5)/x;
  }
  else{
     res = gsl_sf_bessel_jl(l,x);
  }
  return pf_photoz(z,(int)ar[2])*growfac(1./(1+z))*g0*(res*gbias.b1_function(z,(int)ar[2])+Psi_RSD_z(z,params));
}

double int_for_Psi_cl_noRSD (double z, void *params){
  static double g0 = 0.;
  if (g0 ==0){g0 =1./growfac(1.);}
  double *ar = (double *) params;
  int l = (int) ar[0];
  double res, x = ar[1]*f_K(chi(1./(1.+z)));
  if (x < 0.1*sqrt(2.*l) && (2*l+1) < GSL_SF_DOUBLEFACT_NMAX){ // small-x limit for j_l(x), to avoid gsl underflow error
    res =pow(x,l)/gsl_sf_doublefact((unsigned int)(2*l+1))*(1-0.5*x*x/(2.*l+3));
  }
  else if (x > 100.*l*l){// large-x limit for j_l(x), to avoid gsl error
    res = sin(x-l*M_PI*0.5)/x;
  }
  else{
     res = gsl_sf_bessel_jl(l,x);
  }
  return pf_photoz(z,(int)ar[2])*growfac(1./(1+z))*g0*(res*gbias.b1_function(z,(int)ar[2]));
}


double int_for_Psi_cl_withMag (double z, void *params){
  static double g0 = 0.;
  if (g0 ==0){g0 =1./growfac(1.);}
  double *ar = (double *) params;
  int l = (int) ar[0];
  double a = 1./(1.+z);
  double fK = f_K(chi(a));
  double res, x = ar[1]*fK;
  if (x < 0.1*sqrt(2.*l) && (2*l+1) < GSL_SF_DOUBLEFACT_NMAX){ // small-x limit for j_l(x), to avoid gsl underflow error
    res =pow(x,l)/gsl_sf_doublefact((unsigned int)(2*l+1))*(1-0.5*x*x/(2.*l+3));
  }
  else if (x > 100.*l*l){// large-x limit for j_l(x), to avoid gsl error
    res = sin(x-l*M_PI*0.5)/x;
  }
  else{
     res = gsl_sf_bessel_jl(l,x);
  }
  double mag;
  mag = gbias.b_mag[(int)ar[2]] * l*(l+1.)* W_mag(a, fK, ar[2])/x/x * growfac(a)*g0 * dchi_da(a)*a*a;
  return pf_photoz(z,(int)ar[2])*growfac(a)*g0*(res*gbias.b1_function(z,(int)ar[2])+Psi_RSD_z(z,params)) + res*mag;
}

double int_for_Psi_Mag (double z, void *params){
  static double g0 = 0.;
  if (g0 ==0){g0 =1./growfac(1.);}
  double *ar = (double *) params;
  int l = (int) ar[0];
  double a = 1./(1.+z);
  double fK = f_K(chi(a));
  double res, x = ar[1]*fK;
  if (x < 0.1*sqrt(2.*l) && (2*l+1) < GSL_SF_DOUBLEFACT_NMAX){ // small-x limit for j_l(x), to avoid gsl underflow error
    res =pow(x,l)/gsl_sf_doublefact((unsigned int)(2*l+1))*(1-0.5*x*x/(2.*l+3));
  }
  else if (x > 100.*l*l){// large-x limit for j_l(x), to avoid gsl error
    res = sin(x-l*M_PI*0.5)/x;
  }
  else{
     res = gsl_sf_bessel_jl(l,x);
  }
  double mag;
  mag = gbias.b_mag[(int)ar[2]] * l*(l+1.)* W_mag(a, fK, ar[2])/x/x * growfac(a)*g0 * dchi_da(a)*a*a;
  return res*mag;
}

double int_for_Psi_L(double z, void *params){
  static double g0 = 0.;
  if (g0 ==0){g0 =1./growfac(1.);}
  double *ar = (double *) params;
  int l = (int) ar[0];
  double a = 1./(1+z);
  double fK = f_K(chi(a));
  double res, x = ar[1]*fK;
  if (x < 0.1*sqrt(2.*l) && (2*l+1) < GSL_SF_DOUBLEFACT_NMAX){ // small-x limit for j_l(x), to avoid gsl underflow error
    res =pow(x,l)/gsl_sf_doublefact((unsigned int)(2*l+1))*(1-0.5*x*x/(2.*l+3));
  }
  else if (x > 100.*l*l){// large-x limit for j_l(x), to avoid gsl error
    res = sin(x-l*M_PI*0.5)/x;
  }
  else{
     res = gsl_sf_bessel_jl(l,x);
  }
  return growfac(a)*g0*res*W_kappa(a, fK, ar[2])*sqrt((l+2.)*(l+1.)*l*(l-1.))/x/x * dchi_da(a)*a*a;
}

double int_for_Psi_L_withIA(double z, void *params){
  static double g0 = 0.;
  if (g0 ==0){g0 =1./growfac(1.);}
  double *ar = (double *) params;
  int l = (int) ar[0];
  double a = 1./(1+z);
  double fK = f_K(chi(a));
  double res, x = ar[1]*fK;
  if (x < 0.1*sqrt(2.*l) && (2*l+1) < GSL_SF_DOUBLEFACT_NMAX){ // small-x limit for j_l(x), to avoid gsl underflow error
    res =pow(x,l)/gsl_sf_doublefact((unsigned int)(2*l+1))*(1-0.5*x*x/(2.*l+3));
  }
  else if (x > 100.*l*l){// large-x limit for j_l(x), to avoid gsl error
    res = sin(x-l*M_PI*0.5)/x;
  }
  else{
     res = gsl_sf_bessel_jl(l,x);
  }
  double norm = cosmology.Omega_m*nuisance.c1rhocrit_ia*growfac(0.9999)/growfac(a)*nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia),nuisance.eta_ia);
  return growfac(a)*g0*res*(W_kappa(a, fK, ar[2])-norm*W_source(a,ar[2]))*sqrt((l+2.)*(l+1.)*l*(l-1.))/x/x * dchi_da(a)*a*a;
}
double int_for_Psi_IA(double z, void *params){
  static double g0 = 0.;
  if (g0 ==0){g0 =1./growfac(1.);}
  double *ar = (double *) params;
  int l = (int) ar[0];
  double a = 1./(1+z);
  double fK = f_K(chi(a));
  double res, x = ar[1]*fK;
  if (x < 0.1*sqrt(2.*l) && (2*l+1) < GSL_SF_DOUBLEFACT_NMAX){ // small-x limit for j_l(x), to avoid gsl underflow error
    res =pow(x,l)/gsl_sf_doublefact((unsigned int)(2*l+1))*(1-0.5*x*x/(2.*l+3));
  }
  else if (x > 100.*l*l){// large-x limit for j_l(x), to avoid gsl error
    res = sin(x-l*M_PI*0.5)/x;
  }
  else{
     res = gsl_sf_bessel_jl(l,x);
  }
  double norm = cosmology.Omega_m*nuisance.c1rhocrit_ia*growfac(0.9999)/growfac(a)*nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia),nuisance.eta_ia);
  return growfac(a)*g0*res*(-norm*W_source(a,ar[2]))*sqrt((l+2.)*(l+1.)*l*(l-1.))/x/x * dchi_da(a)*a*a;
}

double Psi_cl(double k, int l,int ni){
  double ar[3] = {(double)l,k,(double) ni};
  return int_gsl_integrate_low_precision(int_for_Psi_cl,(void*)ar,tomo.clustering_zmin[ni],tomo.clustering_zmax[ni],NULL,100);
}
double Psi_RSD(double k, int l,int ni){
  double ar[3] = {(double)l,k,(double) ni}; 
  return int_gsl_integrate_low_precision(int_for_Psi_RSD,(void*)ar,tomo.clustering_zmin[ni],tomo.clustering_zmax[ni],NULL,100);
}
double Psi_cl_withMag(double k, int l,int ni){
  double ar[3] = {(double)l,k,(double) ni};
  // printf("chis: %lg, %lg\n", chi(1./(1.+tomo.clustering_zmin[ni]))*cosmology.coverH0 / cosmology.h0,chi(1./(1.+tomo.clustering_zmax[ni]))*cosmology.coverH0 / cosmology.h0);exit(0);
  // return int_gsl_integrate_low_precision(int_for_Psi_cl_withMag,(void*)ar,0.,tomo.clustering_zmax[ni],NULL,100);
  return Psi_Mag(k, l, ni) + int_gsl_integrate_low_precision(int_for_Psi_cl,(void*)ar,tomo.clustering_zmin[ni],tomo.clustering_zmax[ni],NULL,100);
}

double Psi_Mag(double k, int l,int ni){
  double ar[3] = {(double)l,k,(double) ni};
  return int_gsl_integrate_low_precision(int_for_Psi_Mag,(void*)ar,0.,tomo.clustering_zmax[ni],NULL,100);
}

double Psi_cl_noRSD(double k, int l,int ni){
  double ar[3] = {(double)l,k,(double) ni};
  return int_gsl_integrate_low_precision(int_for_Psi_cl_noRSD,(void*)ar,tomo.clustering_zmin[ni],tomo.clustering_zmax[ni],NULL,100);
}
double Psi_L(double k, int l,int ni){
  double ar[3] = {(double)l,k,(double) ni};
  // int i;
  // for(i=0;i<5;i++){printf("i,chimin,chimax:%d,%lg,%lg\n", i, chi(amin_source(ni))*cosmology.coverH0 / cosmology.h0,chi(amax_source(ni))*cosmology.coverH0 / cosmology.h0 );}
  //   exit(0);
  return int_gsl_integrate_low_precision(int_for_Psi_L,(void*)ar,0.,tomo.shear_zmax[ni],NULL,100);
  // return int_gsl_integrate_low_precision(int_for_Psi_L_withIA,(void*)ar,amin_source(ni),amax_source(ni),NULL,100);
}
double Psi_L_withIA(double k, int l,int ni){
  double ar[3] = {(double)l,k,(double) ni};
  // int i;
  // for(i=0;i<5;i++){printf("i,chimin,chimax:%d,%lg,%lg\n", i, chi(amin_source(ni))*cosmology.coverH0 / cosmology.h0,chi(amax_source(ni))*cosmology.coverH0 / cosmology.h0 );}
  //   exit(0);
  return int_gsl_integrate_low_precision(int_for_Psi_L_withIA,(void*)ar,0.,tomo.shear_zmax[ni],NULL,100);
}
double Psi_IA(double k, int l,int ni){
  double ar[3] = {(double)l,k,(double) ni};
  // int i;
  // for(i=0;i<5;i++){printf("i,chimin,chimax:%d,%lg,%lg\n", i, chi(amin_source(ni))*cosmology.coverH0 / cosmology.h0,chi(amax_source(ni))*cosmology.coverH0 / cosmology.h0 );}
  //   exit(0);
  return int_gsl_integrate_low_precision(int_for_Psi_IA,(void*)ar,tomo.shear_zmin[ni],tomo.shear_zmax[ni],NULL,100);
}

double int_for_C_cl_nonLimber (double lk, void *params){
  double k = exp(lk);
  int *ar = (int *) params;
  return k*k*k*pow(Psi_cl_withMag(k,ar[0],ar[1]),2.)*p_lin(k,1.0)*G_taper(k);
  // return k*k*k*pow(Psi_RSD(k,ar[0],ar[1]),2.)*p_lin(k,1.0)*G_taper(k);
}
double int_for_C_cl_RSD (double klog, void *params){
  int *ar = (int *) params;
  double k = exp(klog);
  return k*k*k*pow(Psi_RSD(k,ar[0],ar[1]),2.)*p_lin(k,1.0)*G_taper(k);
}
double int_for_C_cl_nonLimber_tomo (double k, void *params){
  int *ar = (int *) params;
  return k*k*k*Psi_cl(k,ar[0],ar[1])*Psi_cl(k,ar[0],ar[2])*p_lin(k,1.0)*G_taper(k);
}

// for testing Mag
double int_for_C_Mag_nonLimber (double lk, void *params){
  double k = exp(lk);
  int *ar = (int *) params;
  return k*k*k*pow(Psi_Mag(k,ar[0],ar[1]),2.)*p_lin(k,1.0)*G_taper(k);
}

double int_for_C_shear_nonLimber_tomo (double lk, void *params){
  double k = exp(lk);
  int *ar = (int *) params;
  return k*k*k*Psi_L(k,ar[0],ar[1])*Psi_L(k,ar[0],ar[2])*p_lin(k,1.0)*G_taper(k);
}

double int_for_C_cl_lin(double a, void *params)
{
  double res,ell, fK, k;
  double *ar = (double *) params;
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  
  res=W_gal(a,ar[0])*W_gal(a,ar[1])*dchi_da(a)/fK/fK;
  res= res*p_lin(k,a)*G_taper(k);
  return res;
}


double C_cl_lin_nointerp(double l, int ni, int nj)  //galaxy clustering power spectrum of galaxy bins ni, nj
{
  double array[3] = {1.0*ni,1.0*nj,l};
  return int_gsl_integrate_medium_precision(int_for_C_cl_lin,(void*)array,fmax(amin_lens(ni),amin_lens(nj)),0.99999,NULL,1000);
}


double C_cl_non_Limber(int l, int ni, int nj){ //includes RSD too!
  int ar[3] ={l,ni,nj};
  double res;
  if (2*l+5 > GSL_SF_DOUBLEFACT_NMAX){
    printf("cosmo2D_exact.c:C_cl_non_Limber: l = %d too large for stable evaluation of Bessel functions\nSwitching to Limber approximation\n", l);
    return C_cl_tomo_nointerp(1.*l,ni,nj);
  }
  //gsl_set_error_handler_off ();
  if (ni == nj){
    //checked that these boundaries give better than 1% accuracy for l > 2
    // double kmin = fmax(limits.k_min_cH0, 0.25*l/f_K(chi(1./(1.+tomo.clustering_zmax[ni]))));
    // double kmax = fmin(limits.k_max_cH0, 20.*l/f_K(chi(1./(1.+tomo.clustering_zmin[ni]))));
    double kmin = limits.k_min_cH0;
    double kmax = limits.k_max_cH0;
    res = 2.*int_gsl_integrate_low_precision(int_for_C_cl_nonLimber,(void*)ar,log(kmin),log(kmax),NULL,100)/M_PI;
    // res = 2.*int_gsl_integrate_low_precision(int_for_C_Mag_nonLimber,(void*)ar,log(kmin),log(kmax),NULL,100)/M_PI;
  }
  else res = 2.*int_gsl_integrate_low_precision(int_for_C_cl_nonLimber_tomo,(void*)ar,log(limits.k_min_cH0),log(limits.k_max_cH0),NULL,1000)/M_PI;
  //gsl_set_error_handler (NULL);

  // non add non-linear power spectrum contributions in Limber approximation
  res = res + C_cl_tomo_nointerp(1.*l,ni,nj) - C_cl_lin_nointerp(1.*l,ni,nj);
  return res;
}

double C_cl_RSD(int l, int ni, int nj){
  int ar[3] ={l,ni,nj};
  double res;
  if (2*l+5 > GSL_SF_DOUBLEFACT_NMAX){
    printf("cosmo2D_exact.c:C_cl_RSD: l = %d too large for stable evaluation of Bessel functions\nSwitching to Limber approximation\n", l);
    return C_cl_tomo_nointerp(1.*l,ni,nj);
  }
  gsl_set_error_handler_off ();
  if (ni == nj){
    double kmin = fmax(limits.k_min_cH0, 0.25*l/f_K(chi(1./(1.+tomo.clustering_zmax[ni]))));
    double kmax = fmin(limits.k_max_cH0, 20.*l/f_K(chi(1./(1.+tomo.clustering_zmin[ni]))));
    res = 2.*int_gsl_integrate_low_precision(int_for_C_cl_RSD,(void*)ar,log(kmin),log(kmax),NULL,1000)/M_PI;
  }
  else res = 0.;
  gsl_set_error_handler (NULL);
  //add regular galaxy power spectrum contributions in Limber approximation
  res = res + C_cl_tomo_nointerp(1.*l,ni,nj) - C_cl_lin_nointerp(1.*l,ni,nj);
  return res;
}



double w_tomo_nonLimber(int nt, int ni, int nj){
  // if(1){return 0.;}
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
    for(i=0; i<like.Ntheta ; i++){
      xmin[i]=cos(exp(log(like.vtmin)+(i+0.0)*logdt));
      xmax[i]=cos(exp(log(like.vtmin)+(i+1.0)*logdt));
    }

    for (i = 0; i<NTHETA; i ++){
      printf("Tabulating Legendre coefficients %d/%d\n",i+1, NTHETA);
      gsl_sf_legendre_Pl_array(LMAX, xmin[i],Pmin);
      gsl_sf_legendre_Pl_array(LMAX, xmax[i],Pmax);
      for (int l = 1; l < LMAX; l ++){
        //Pl[i][l] = (2*l+1.)/(4.*M_PI)*Pmin[l];
        Pl[i][l] = 1./(4.*M_PI)*(Pmin[l+1]-Pmax[l+1]-Pmin[l-1]+Pmax[l-1])/(xmin[i]-xmax[i]);
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

    char *outfilename = (char*)malloc(40 * sizeof(char));;

    for (nz = 0; nz <tomo.clustering_Nbin; nz ++){
      sprintf(outfilename, "cls2/c_cl_exact_%d_%d_rsd_mag.txt", nz,nz);
      FILE *OUT = fopen(outfilename, "w");

      int L = 1;
      // initialize to large value in order to start while loop
      dev=10.*tolerance;
      // while (fabs(dev) > tolerance){
      while (L<=90){
        //Cl[L] = C_cl_RSD(L,nz,nz);
        //printf("nz: %d\n", nz);
        Cl[L] = C_cl_non_Limber(L,nz,nz);
        dev = Cl[L]/C_cl_tomo_nointerp((double)L,nz,nz)-1.;
        //printf("%d %lg %lg %lg\n", L, Cl[L], C_cl_tomo_nointerp(1.*L,nz,nz), C_cl_lin_nointerp(1.*L,nz,nz));
        fprintf(OUT, "%d %lg %lg %lg\n", L, Cl[L], C_cl_tomo_nointerp(1.*L,nz,nz), C_cl_lin_nointerp(1.*L,nz,nz));
        L = L+1;
      }
      printf("switching to Limber calculation at l = %d\n",L);
      for (l = L; l < LMAX; l++){
        Cl[l]=C_cl_tomo((double)l,nz,nz);
        // fprintf(OUT, "%d %lg %lg %lg\n", l, Cl[l], C_cl_tomo_nointerp(1.*L,nz,nz), C_cl_lin_nointerp(1.*L,nz,nz));
      }
      fclose(OUT);
      // exit(0);
      // if(nz==4){exit(0);}

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


double int_for_C_gl_nonLimber_tomo (double lk, void *params){
  double k = exp(lk);
  int *ar = (int *) params;
  return k*k*k*Psi_cl_withMag(k,ar[0],ar[1])*Psi_L(k,ar[0],ar[2])*p_lin(k,1.0)*G_taper(k);
}

double int_for_C_gl_withIA_nonLimber_tomo (double lk, void *params){
  double k = exp(lk);
  int *ar = (int *) params;
  return k*k*k*Psi_cl_withMag(k,ar[0],ar[1])*Psi_L_withIA(k,ar[0],ar[2])*p_lin(k,1.0)*G_taper(k);
}

double int_for_C_gl_onlyIA_nonLimber_tomo (double lk, void *params){
  double k = exp(lk);
  int *ar = (int *) params;
  return k*k*k*Psi_cl_withMag(k,ar[0],ar[1])*Psi_IA(k,ar[0],ar[2])*p_lin(k,1.0)*G_taper(k);
  // return k*k*k*Psi_cl(k,ar[0],ar[1])*Psi_IA(k,ar[0],ar[2])*p_lin(k,1.0)*G_taper(k);
}

double int_for_C_gl_lin(double a, void *params)
{
  double res,ell, fK, k;
  double *ar = (double *) params;
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  
  double ell_prefactor1 = (ar[2])*(ar[2]+1.);
  double ell_prefactor2 = (ar[2]-1.)*ell_prefactor1*(ar[2]+2.);
  if(ell_prefactor2<=0.) 
    ell_prefactor2=0.;
  else
    ell_prefactor2=sqrt(ell_prefactor2);

  double chi_0,chi_1,a_0,a_1;
  chi_0 = f_K(ell/k);
  chi_1 = f_K((ell+1.)/k);
  if (chi_1 > chi(limits.a_min)){
    return 0;}
  a_0 = a_chi(chi_0);
  a_1 = a_chi(chi_1);

  res=(W_gal(a,ar[0])+W_RSD(ell, a_0, a_1, ar[0]) +W_mag(a,fK,ar[0])*(ell_prefactor1/ell/ell -1.) )*W_kappa(a,fK, ar[1])*dchi_da(a)/fK/fK * ell_prefactor2/ell/ell;
  res= res*p_lin(k,a)*G_taper(k);
  return res;
}

double C_gl_lin_nointerp(double l, int ni, int nj)  //galaxy clustering power spectrum of galaxy bins ni, nj
{
  double array[3] = {1.0*ni,1.0*nj,l};
  return int_gsl_integrate_medium_precision(int_for_C_gl_lin,(void*)array,amin_lens(ni),0.9999,NULL,1000);
}

double int_for_C_gl_IA_lin(double a, void *params)
{
  double res,ell, fK, k;
  double *ar = (double *) params;
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  
  double ell_prefactor1 = (ar[2])*(ar[2]+1.);
  double ell_prefactor2 = (ar[2]-1.)*ell_prefactor1*(ar[2]+2.);
  if(ell_prefactor2<=0.) 
    ell_prefactor2=0.;
  else
    ell_prefactor2=sqrt(ell_prefactor2);

  double chi_0,chi_1,a_0,a_1;
  chi_0 = f_K(ell/k);
  chi_1 = f_K((ell+1.)/k);
  if (chi_1 > chi(limits.a_min)){
    return 0;}
  a_0 = a_chi(chi_0);
  a_1 = a_chi(chi_1);

  double norm = cosmology.Omega_m*nuisance.c1rhocrit_ia*growfac(0.9999)/growfac(a)*nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia),nuisance.eta_ia);

  res=(W_gal(a,ar[0])+W_RSD(ell, a_0, a_1, ar[0]) +W_mag(a,fK,ar[0])*(ell_prefactor1/ell/ell -1.) )*(W_kappa(a,fK, ar[1])-W_source(a,ar[1])*norm)*dchi_da(a)/fK/fK * ell_prefactor2/ell/ell;
  // res=(W_gal(a,ar[0]) )*(W_kappa(a,fK, ar[1])-W_source(a,ar[1])*norm)*dchi_da(a)/fK/fK * ell_prefactor2/ell/ell;
  res= res*p_lin(k,a)*G_taper(k);
  return res;
}

double C_gl_lin_IA_nointerp(double l, int ni, int nj)  //galaxy clustering power spectrum of galaxy bins ni, nj
{
  double array[3] = {1.0*ni,1.0*nj,l};
  // return int_gsl_integrate_medium_precision(int_for_C_gl_IA_lin,(void*)array,amin_lens(ni),amax_lens(ni),NULL,1000);
  return int_gsl_integrate_medium_precision(int_for_C_gl_IA_lin,(void*)array,amin_lens(ni),0.9999,NULL,1000);
}

double C_gl_non_Limber(double l, int ni, int nj){ //includes RSD too!
  int ar[3] ={(int)l,ni,nj};
  double res;
  if (2*l+5 > GSL_SF_DOUBLEFACT_NMAX){
    printf("cosmo2D_exact.c:C_gl_non_Limber: l = %.1f too large for stable evaluation of Bessel functions\nSwitching to Limber approximation\n", l);
    return C_gl_tomo_nointerp(1.*l,ni,nj);
  }
  // gsl_set_error_handler_off ();
  // if (ni == nj){
  //   //checked that these boundaries give better than 1% accuracy for l > 2
  //   double kmin = fmax(limits.k_min_cH0, 0.25*l/f_K(chi(1./(1.+tomo.clustering_zmax[ni]))));
  //   double kmax = fmin(limits.k_max_cH0, 20.*l/f_K(chi(1./(1.+tomo.clustering_zmin[ni]))));
  //   res = 2.*int_gsl_integrate_low_precision(int_for_C_gl_nonLimber,(void*)ar,log(kmin),log(kmax),NULL,100)/M_PI;
  // }
  res = 2.*int_gsl_integrate_low_precision(int_for_C_gl_nonLimber_tomo,(void*)ar,log(limits.k_min_cH0),log(limits.k_max_cH0),NULL,1000)/M_PI;
  //gsl_set_error_handler (NULL);

  // non add non-linear power spectrum contributions in Limber approximation
  res = res + C_gl_tomo_nointerp(1.*l,ni,nj) - C_gl_lin_nointerp(1.*l,ni,nj);
  return res;
}

double C_gl_withIA_non_Limber(double l, int ni, int nj){ //includes RSD too!
  int ar[3] ={(int)l,ni,nj};
  double res;
  if (2*l+5 > GSL_SF_DOUBLEFACT_NMAX){
    printf("cosmo2D_exact.c:C_gl_non_Limber: l = %.1f too large for stable evaluation of Bessel functions\nSwitching to Limber approximation\n", l);
    return C_ggl_IA(1.*l,ni,nj);
  }
  // gsl_set_error_handler_off ();
  // if (ni == nj){
  //   //checked that these boundaries give better than 1% accuracy for l > 2
  //   double kmin = fmax(limits.k_min_cH0, 0.25*l/f_K(chi(1./(1.+tomo.clustering_zmax[ni]))));
  //   double kmax = fmin(limits.k_max_cH0, 20.*l/f_K(chi(1./(1.+tomo.clustering_zmin[ni]))));
  //   res = 2.*int_gsl_integrate_low_precision(int_for_C_gl_nonLimber,(void*)ar,log(kmin),log(kmax),NULL,100)/M_PI;
  // }
  res = 2.*int_gsl_integrate_low_precision(int_for_C_gl_withIA_nonLimber_tomo,(void*)ar,log(limits.k_min_cH0),log(limits.k_max_cH0),NULL,1000)/M_PI;
  // res = 2.*int_gsl_integrate_low_precision(int_for_C_gl_onlyIA_nonLimber_tomo,(void*)ar,log(limits.k_min_cH0),log(limits.k_max_cH0),NULL,1000)/M_PI;
  //gsl_set_error_handler (NULL);

  // non add non-linear power spectrum contributions in Limber approximation
  res = res + C_ggl_IA(1.*l,ni,nj) - C_gl_lin_IA_nointerp(1.*l,ni,nj);
  return res;
}

double w_gamma_t_nonLimber(int nt, int ni, int nj){
  // if(1){return 0.;}
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
    printf("cosmo2D_exact.c:w_gamma_t_tomo: like.Ntheta not initialized\nEXIT\n"); exit(1);
  }
  if (Pl ==0){
    Pl =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
    Cl = create_double_vector(0,LMAX-1);
    w_vec = create_double_vector(0,tomo.ggl_Npowerspectra*like.Ntheta-1);
    NTHETA = like.Ntheta;
    double *xmin, *xmax, *Pmin, *Pmax, *dP;
    xmin= create_double_vector(0, like.Ntheta-1);
    xmax= create_double_vector(0, like.Ntheta-1);
    double logdt=(log(like.vtmax)-log(like.vtmin))/like.Ntheta;
    for(i=0; i<like.Ntheta ; i++){
      xmin[i]=cos(exp(log(like.vtmin)+(i+0.0)*logdt));
      xmax[i]=cos(exp(log(like.vtmin)+(i+1.0)*logdt));
    }
    Pmin= create_double_vector(0, LMAX+1);
    Pmax= create_double_vector(0, LMAX+1);

    for (i = 0; i<NTHETA; i ++){
      printf("Tabulating Legendre coefficients %d/%d\n",i+1, NTHETA);
      gsl_sf_legendre_Pl_array(LMAX, xmin[i],Pmin);
      gsl_sf_legendre_Pl_array(LMAX, xmax[i],Pmax);
      for (int l = 2; l < LMAX; l ++){
        //Pl[i][l] = (2.*l+1)/(4.*M_PI*l*(l+1))*gsl_sf_legendre_Plm(l,2,cos(like.theta[i]));  
        Pl[i][l] = (2.*l+1)/(4.*M_PI*l*(l+1)*(xmin[i]-xmax[i]))
        *((l+2./(2*l+1.))*(Pmin[l-1]-Pmax[l-1])
        +(2-l)*(xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
        -2./(2*l+1.)*(Pmin[l+1]-Pmax[l]));
      }
    }
    free_double_vector(xmin,0,like.Ntheta-1);
    free_double_vector(xmax,0,like.Ntheta-1);
    free_double_vector(Pmin,0,LMAX+1);
    free_double_vector(Pmax,0,LMAX+1);
  }
  if (recompute_ggl(C,G,N,ni)){
    //required fractional accuracy in C(l)
    double tolerance= 0.01;
    //dev will be the actual difference between exact and Limber calcuation
    double dev;

    if (like.IA != 0 && like.IA != 3 && like.IA != 4){printf("cosmo2D_real.c: w_gamma_t_tomo does not support like.IA = %d yet\nEXIT!\n",like.IA); exit(1);}

    C_tomo_pointer C_gl_pointer = &C_gl_tomo;
    C_tomo_pointer C_gl_nointerp_pointer = &C_gl_tomo_nointerp;
    C_tomo_pointer C_gl_nonlimber_pointer = &C_gl_non_Limber;
    C_tomo_pointer C_gl_lin_pointer = &C_gl_lin_nointerp;

    if (like.IA ==3 || like.IA ==4) {
      C_gl_pointer = &C_ggl_IA_tab;
      C_gl_nointerp_pointer = &C_ggl_IA;
      C_gl_nonlimber_pointer = &C_gl_withIA_non_Limber;
      C_gl_lin_pointer = &C_gl_lin_IA_nointerp;
    }

    char *outfilename = (char*)malloc(40 * sizeof(char));;

    for (nz = 0; nz <tomo.ggl_Npowerspectra; nz ++){
      // nz = 0;
      sprintf(outfilename, "cls2/c_gl_exact_%d_%d_mag_IA.txt", ZL(nz),ZS(nz));
      FILE *OUT = fopen(outfilename, "w");

      int L = 2;
      // initialize to large value in order to start while loop
      dev=10.*tolerance;
      while (L<=90){
      // while (fabs(dev) > tolerance){
        //Cl[L] = C_cl_RSD(L,nz,nz);
        printf("nz:%d\n", nz);
        Cl[L] = C_gl_nonlimber_pointer((double)L,ZL(nz),ZS(nz));
        dev = Cl[L]/C_gl_nointerp_pointer((double)L,ZL(nz),ZS(nz))-1.;
       printf("%d %lg %lg %lg\n", L, Cl[L], C_gl_nointerp_pointer(1.*L,ZL(nz),ZS(nz)), C_gl_lin_pointer(1.*L,ZL(nz),ZS(nz)));
        fprintf(OUT, "%d %lg %lg %lg\n", L, Cl[L], C_gl_nointerp_pointer(1.*L,ZL(nz),ZS(nz)), C_gl_lin_pointer(1.*L,ZL(nz),ZS(nz)));
        L = L+1;
      }
      fclose(OUT);
      // exit(0);
      // printf("switching to Limber calculation at l = %d\n",L);
      for (l = L; l < LMAX; l++){
        Cl[l]=C_gl_pointer((double)l,ZL(nz),ZS(nz));
        // fprintf(OUT, "%d %lg 0.0 0.0\n", l, Cl[l]);
      }
      // fclose(OUT);

      for (i = 0; i < NTHETA; i++){
        w_vec[nz*like.Ntheta+i] =0;
        for (l = 1; l < LMAX; l++){
          w_vec[nz*like.Ntheta+i]+=Pl[i][l]*Cl[l];
        }
      }
    }
    // for (nz = 0; nz <tomo.ggl_Npowerspectra; nz ++){
    //  for (l = 1; l < LMAX; l++){
    //    Cl[l]=C_ggl_IA_tab(1.0*l,ZL(nz),ZS(nz));
    //  }
    //  for (i = 0; i < NTHETA; i++){
    //    w_vec[nz*like.Ntheta+i] =0;
    //    for (l = 1; l < LMAX; l++){
    //      w_vec[nz*like.Ntheta+i]+=Pl[i][l]*Cl[l];
    //    }
    //  }
    // }

    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }
  return w_vec[N_ggl(ni,nj)*like.Ntheta+nt];  
}