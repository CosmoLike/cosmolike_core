#ifndef __HALO_FAST__
#define __HALO_FAST__
#include <time.h>
double conc_Ludlow16_fit(double m, double a);
double conc_Ludlow16(double m, double a);

// use this routine in inner I[0-1]j and modify int_rho_nfw if to work with different halo profiles, e.g. to include adiabatic contraction, AGN feedback, etc.
/*halo model building blocks */
double I0j_y (int j, double k1, double k2, double k3, double k4, double a, int is_y[4]);
/*halo model matter power spectrum, bispectrum, trispectrum*/
double p_1h_hmcode(double k, double a);
double p_2h_hmcode(double k, double a);
double Pdelta_hmcode_2020(double k, double a);
double P_my(double k,double a);


/* ********************* Routines for halo model ************* */
/*++++++++++++++++++++++++++++++++++++++++*
*  Variance of the density field         *
*++++++++++++++++++++++++++++++++++++++++*/

double dlogsigma_dlogR(double m)
{
    double dlogm = 0.0001;
    double logm = log(m);
    double m1 = exp(logm + dlogm);
    double dlogsigma2 = log(sigma2(m1) / sigma2(m));
    double dlogR = log(radius(m1) / radius(m));
    return 0.5 * dlogsigma2 / dlogR;
}
double dlogPlin_dlogk(double m)
{
    double k0 = 0.69* 2*M_PI / radius(m);
    double dlogk = 0.001;
    double k1 = exp(log(k0) + dlogk);
    double dlogPlin = log(p_lin(k1,1.) / p_lin(k0,1.));
    return dlogPlin / dlogk;
}


/***************** halo profiles ***************/
/******** mass-concentration relation **********/
double conc(double m, double a) // for matter
{
#ifdef LUDLOW16
  return conc_Ludlow16(m, a);
#endif
#ifdef LUDLOW16_FIT
  return conc_Ludlow16_fit(m, a);
#endif

  double result;
  double M0=pow(10,nuisance.gas_lgM0);
  // result = 9.*pow(nu(m,a),-.29)*pow(growfac(a)/growfac(1.),1.15);// Bhattacharya et al. 2013, Delta = 200 rho_{mean} (Table 2)
  result = 10.14*pow(m/2.e+12,-0.081)*pow(a,1.01); //Duffy et al. 2008 (Delta = 200 mean)
  if(like.feedback_on == 1){
    result *= (1.+nuisance.gas_eps1 + (nuisance.gas_eps2 - nuisance.gas_eps1)/(1.+pow(M0/m, nuisance.gas_beta)) );
  }
    return result;
}

double conc_y(double m, double a) // for y field, with a separate set of gas parameters lgM0,beta,eps1,eps2.
// ONLY USED FOR TEST_CALIB
{
  double result;
  double M0=pow(10,nuisance.gas_lgM0_v2);
  // result = 9.*pow(nu(m,a),-.29)*pow(growfac(a)/growfac(1.),1.15);// Bhattacharya et al. 2013, Delta = 200 rho_{mean} (Table 2)
  result = 10.14*pow(m/2.e+12,-0.081)*pow(a,1.01); //Duffy et al. 2008 (Delta = 200 mean) - If use this, set c_max=50 !!
  if(like.feedback_on == 1){
    result *= (1.+nuisance.gas_eps1_v2 + (nuisance.gas_eps2_v2 - nuisance.gas_eps1_v2)/(1.+pow(M0/m, nuisance.gas_beta_v2)) );
  }
  return result;
}

#ifdef LUDLOW16_FIT
double conc_Ludlow16_fit(double m, double a)
// Ludlow16 fit at Planck cosmology
{
  double c0, beta, gamma1, gamma2, nu0;
  c0 = 3.395 * pow(a, 0.215);
  beta = 0.307 * pow(a, -0.540);
  gamma1 = 0.628 * pow(a, 0.047);
  gamma2 = 0.317 * pow(a, 0.893);
  double a2 = a*a;
  nu0 = (4.135 - 0.564 / a - 0.210 / a2 + 0.0557 / (a2*a) - 0.00348 / (a2*a2)) * growfac(1.) / growfac(a);
  double nu_nu0_ratio = nu(m,a) / nu0;
  return c0 * pow(nu_nu0_ratio, -gamma1) * pow(1+ pow(nu_nu0_ratio, 1./beta), -beta*(gamma2-gamma1));
}
#endif


#ifdef LUDLOW16
double conc_Ludlow16(double m, double a)
{
    static double **table=0;
  static int N_c = 30, N_a = 50;
  static double c_min = 0.608, c_max = 50.; // the method produces cmin ~0.607
  static double a_min, a_max = 0.999;
  static double dc, da;
  static double C_Ludlow = 650., f_Ludlow = 0.02;
  double cc, aa;
  double M2, c3;
  double af, rho_2;
  int i, j;
  printf("here!\n");
  if (table==0) {
    table = create_double_matrix(0, N_c-1, 0, N_a-1);
    a_min = limits.a_min_hm;
    dc = (c_max-c_min) / (N_c - 1);
    da = (a_max - a_min) / (N_a - 1);
    printf("here!-\n");
    for (i=0; i<N_c; i++) {
      cc = c_min + i * dc;
      M2 = (log(2) - 0.5) / (log(1.+cc) - cc / (1.+cc));
      c3 = cc*cc*cc;

      printf("here!--\n");
      for (j=0; j<N_a; j++) {
        aa = a_min + j * da;
        printf("here!---\n");
        rho_2 = delta_Delta(aa) * c3 * M2;
        af = pow(cosmology.Omega_m / (rho_2 / C_Ludlow * (cosmology.Omega_m / (aa*aa*aa) + cosmology.Omega_v) - cosmology.Omega_v), 1./3. );
        
        printf("here!---- [aa,af,M2] = %le, %le, %le \n", aa,af,M2);
        table[i][j] = delta_c(aa) * (1./growfac(af) - 1./growfac(aa)) / erfcinv(M2);
        // af can go beyond 1, Unclear what should we do about it!!!
        printf("here!------\n");
      }
    }
  }

  printf("here!!\n");
  if (a < a_min) {return 0;}
  if (a >= a_max) {a = a_max;}

  int ia = floor((a - a_min) / da);
  double frac_high = (a - (a_min + ia * da)) / da;
  double table_interp[N_c];

  if (ia == N_a-1) {
    ia--;
    frac_high = 0;
  }

  double LHS, sig2rf, sig2r;
  sig2r = sigma2(m);
  sig2rf = sigma2(f_Ludlow * m);
  LHS = sqrt(2.* (sig2rf - sig2r));
  printf("here!!!\n");
  for (i=0; i<N_c; i++) {
    table_interp[i] = table[i][ia] * (1-frac_high) + table[i][ia+1] * frac_high - LHS;
  }
  for (i=1; i<N_c; i++) {
    if (table_interp[i] * table_interp[i-1] <= 0) break;
  }
  double frac_table_low;
  frac_table_low = table_interp[i] / (table_interp[i] - table_interp[i-1]);
  printf("here--\n");
  return c_min + dc * ( (i-1) * frac_table_low + i * (1 - frac_table_low) );
}
#endif
/***********  FT of NFW Profile **************/

double u_nfw_c_interp(double c,double k, double m, double a){
  // static cosmopara C;

  static int N_c = 25, N_x = 80;
  static double cmin = 0.1, cmax = 50.;
  static double xmin = 1e-10, xmax = 5e4; // full range of possible k*R_200/c in the code
  static double logxmin = 0., logxmax=0., dx=0., dc=0.;
  
  static double **table=0;
  
  double cc,xlog, xx, xu;
  double x;
  int i,j;

  x = k * r_Delta(m,a)/c;

  // if (recompute_cosmo3D(C)){ //extend this by halo model parameters if these become part of the model
  //   update_cosmopara(&C);
  if (table==0) {
    table = create_double_matrix(0, N_c-1, 0, N_x-1);
    logxmin = log(xmin); logxmax = log(xmax);
    dx = (logxmax - logxmin)/(N_x-1.); dc = (cmax - cmin)/(N_c-1.);

    cc = cmin;
    for (i=0; i<N_c; i++, cc +=dc) {
      xlog  = logxmin;
      for (j=0; j<N_x; j++, xlog += dx) {
        xx = exp(xlog);
        xu = (1.+cc)*xx;
        table[i][j] = (sin(xx)*(gsl_sf_Si(xu)-gsl_sf_Si(xx))- sinl(cc*xx)/xu +cos(xx)*(gsl_sf_Ci(xu)-gsl_sf_Ci(xx)))*1./(log(1.+cc)-cc/(1.+cc)); 
      }
    }
  }
  if (c < cmin || c >cmax){return 0.;}
  if (log(x) < logxmin || log(x) > logxmax){return 0.;}
  // printf("c, x = %le, %le\n",c,x);
  return interpol2d(table, N_c, cmin, cmax, dc, c, N_x, logxmin, logxmax, dx, log(x), 0.0, 0.0);
}

/***********  FT of bound gas K-S Profile (sinc integral part) **************/

double int_KS(double r, void *params){//Fourier kernel * K-S profile, integrand for u_KS
  double *array = (double*)params;
  double k = array[0];
  double c = array[1];
  double rv = array[2];
  double rs =rv/c;
  double x = r/rs;
  return r*sinl(k*r)/k*pow(log(1.+x)/x, nuisance.gas_Gamma_KS/(nuisance.gas_Gamma_KS-1.));
}

double u_KS(double c, double k, double rv){
  // FFT density of K-S profile, truncated at r_Delta, through direct integration
  // use this routine in inner I[0-1]j and modify int_KS to work with different halo profiles
  
  double array[3] ={k,c,rv};
  return int_gsl_integrate_medium_precision(int_KS, (void*)array, 0,rv,NULL, 5000);
}

double int_F_KS(double x, void *params){
  double *array = (double*)params;
  double y = array[0];
  return x*sinl(x*y)/y * pow(log(1.+x)/x, nuisance.gas_Gamma_KS/(nuisance.gas_Gamma_KS-1.));
}
double F_KS(double c, double kr_s){
  double array[1] ={kr_s};
  return int_gsl_integrate_medium_precision(int_F_KS, (void*)array, 0,c,NULL, 5000);
}

#ifdef RUN_FFT
void fft_F_KS(double *c_arr, int Nc, double *krs_arr, int Nx, double **F) {
  config fftconfig;
  fftconfig.nu = 1.01;
  fftconfig.c_window_width = 0.25;
  fftconfig.derivative = 0;
  fftconfig.N_pad = 0;
  fftconfig.N_extrap_low = 0;
  fftconfig.N_extrap_high = 0;

  int ell=0; // spherical bessel order
  double x[Nx];
  int i,j;
  for(i=0;i<Nx;i++){
    x[i] = (ell+1.)/krs_arr[Nx-1-i]; // determine the integrand sampling x-grid, y[::1] = (ell+1)/x[:] is defined in cfftlog
  }
  double **fx;
  fx = malloc(Nc * sizeof(double *));

  // FILE *output;
  // output = fopen("xfx.txt", "w");
  // fprintf(output, "# c, x, fx\n");

  for(j=0;j<Nc;j++){
    fx[j] = malloc(Nx * sizeof(double));
    for(i=0;i<Nx;i++){
      fx[j][i] = x[i]>c_arr[j] ? 0 : pow(x[i],3)*pow(log(1.+x[i])/x[i], nuisance.gas_Gamma_KS/(nuisance.gas_Gamma_KS-1.));
      // fprintf(output, "%le %le %le\n", c_arr[j], x[i], fx[j][i]);
    }
  }
  // fclose(output);
  // exit(0);

  cfftlog_multiple(x, fx, Nx, Nc, &fftconfig, ell, krs_arr, F);
  for(j=0;j<Nc;j++){ free(fx[j]); }
  free(fx);
}
#endif

double int_F_KS_norm(double x, void *params){
  return x*x * pow(log(1.+x)/x, 1/(nuisance.gas_Gamma_KS-1.));
}
double F_KS_norm(double c){
  double array[1] = {c};
  return int_gsl_integrate_medium_precision(int_F_KS_norm, (void*)array, 0,c,NULL, 5000);
}

double u_KS_normalized_interp(double c,double k, double rv){
  // static cosmopara C;
  static nuisancepara N;

  static int N_c = 20, N_x = 200;
  static double cmin = 0.1, cmax = 50.;
  static double xmin = 1e-10, xmax = 5e3; // full range of possible k*R_200/c in the code
  static double logxmin = 0., logxmax=0., dx=0., dc=0.;
   
  static double **table=0;
  
  double cc,xlog, xx, xu;
  double x;
  int i,j;

  x = k * rv/c;

  double F0;

#ifdef RUN_FFT
  N_x=1024; // for efficient fft
#endif

  if (recompute_halo(N)){ // Only nuisance.gas_Gamma is used in creating the table
    update_nuisance(&N);
    if (table==0) table = create_double_matrix(0, N_c-1, 0, N_x-1);
    logxmin = log(xmin); logxmax = log(xmax);
    dx = (logxmax - logxmin)/(N_x-1.); dc = (cmax - cmin)/(N_c-1.);
    
    cc = cmin;

#ifdef RUN_FFT
    double c_arr[N_c], krs_arr[N_x];

    for (i=0; i<N_c; i++, cc +=dc) { c_arr[i]=cc; }
    xlog  = logxmin;
    for (j=0; j<N_x; j++, xlog += dx) { krs_arr[j]=exp(xlog); }

    fft_F_KS(c_arr, N_c, krs_arr, N_x, table);
    for (i=0; i<N_c; i++) {
      F0 = F_KS_norm(c_arr[i]);
      for (j=0; j<N_x; j++) {
        table[i][j] /= F0;
      }
    }
#else
    for (i=0; i<N_c; i++, cc +=dc) {
      xlog  = logxmin;
      F0 = F_KS_norm(cc);
      for (j=0; j<N_x; j++, xlog += dx) {
        xx = exp(xlog);
        table[i][j] = F_KS(cc,xx)/F0;
        // printf("%le %le %le %le \n", cc, xx, F0, table[i][j]);
      }
    }
#endif
    // exit(0);
  }
  // printf("c, x = %le, %le\n",c,x);
  if (c < cmin || c >cmax){return 0.;}
  if (log(x) < logxmin || log(x) > logxmax){return 0.;}
  return interpol2d(table, N_c, cmin, cmax, dc, c, N_x, logxmin, logxmax, dx, log(x), 0.0, 0.0);
}



double int_KS_norm(double r, void *params){//normalization of u_KS
  double *array = (double*)params;
  double c = array[0];
  double rv = array[1];
  double rs =rv/c;
  double x = r/rs;
  return r*r*pow(log(1.+x)/x, 1./(nuisance.gas_Gamma_KS-1.));
}

double KS_norm(double c, double rv){
  
  double array[2] ={c,rv};
  return int_gsl_integrate_medium_precision(int_KS_norm, (void*)array, 0,rv,NULL, 1000);
}

double frac_bnd(double m){
  double M0=pow(10,nuisance.gas_lgM0);
  return cosmology.omb / cosmology.Omega_m / (1.+ pow(M0/m, nuisance.gas_beta));
}
double frac_ejc(double m){
  double M0=pow(10,nuisance.gas_lgM0);
  double lgm=log10(m);
  double frac_star = nuisance.gas_A_star * exp(-0.5* pow((lgm - nuisance.gas_lgM_star)/nuisance.gas_sigma_star,2));
  if(lgm>nuisance.gas_lgM0 && frac_star<nuisance.gas_A_star/3.) {frac_star = nuisance.gas_A_star/3.;}
  return cosmology.omb / cosmology.Omega_m / (1.+ pow(M0/m, -nuisance.gas_beta)) - frac_star;
}

double u_y_bnd(double c, double k, double m, double a){ //unit: [G(M_solar/h)^2 / (c/H0)]
  double rv = r_Delta(m,a);
  double mu_p = 4./(3.+5*nuisance.gas_f_H);
  double mu_e = 2./(1.+nuisance.gas_f_H);
  // return 2.*nuisance.gas_alpha /(3.*a) * mu_p / mu_e * frac_bnd(m) * m*m/rv * u_KS(c, k, rv)/KS_norm(c, rv);
  // return 2.*nuisance.gas_alpha /(3.*a) * mu_p / mu_e * frac_bnd(m) * m*m/rv * u_KS(c, k, rv)/(F_KS_norm(c)*pow(rv/c,3));
  return 2.*nuisance.gas_alpha /(3.*a) * mu_p / mu_e * frac_bnd(m) * m*m/rv * u_KS_normalized_interp(c, k, rv);
}

double u_y_ejc(double m){
  static double num_p = 1.1892e57; // proton number in 1 solar mass, unit [1/Msun]
   // convert ejected gas temperature (K) to energy eV then to unit [G (Msun/h)^2 / (c/H0) * h]
  // double E_w = pow(10,nuisance.gas_lgT_w) * 8.6173e-5 * 1.81945e-68;
  double E_w = pow(10,nuisance.gas_lgT_w) * 8.6173e-5 * 5.616e-44;
  // printf("Tw: %le, %le\n", nuisance.gas_lgT_w,frac_ejc(m));
  double mu_e = 2./(1.+nuisance.gas_f_H);
  return num_p * m * frac_ejc(m) / mu_e * E_w; // [m]=[Msun/h]; final unit in [G(Msun/h)^2 / (c/H0)]
}

/////////


double inner_I0j_y (double logm, void *para){
  double *array = (double *) para;
  double m = exp(logm);
  long double u = 1.0;
  double a= array[6];
  double c = conc(m,a);

  double c_y;
#ifdef TEST_CALIB
  c_y = conc_y(m,a);
#else
  c_y = c;
#endif

  int l;
  int j = (int)(array[5]);
  double vol = m/(cosmology.rho_crit*cosmology.Omega_m); //rho_crit: [density]=[Msun/h*(H0/c)^3], m: [Msun/h]
  for (l = 0; l< j; l++){
    if(array[7+l]==-2.){ // ni=-2: y-field
      u *= u_y_bnd(c_y,array[l],m,a);
    }else{
      u *= vol*u_nfw_c(c,array[l],m,a);
    } // u_(kappa|g): [volume]=[(c/H0)^3], u_y: [energy]=[G(Msun/h)^2 / (c/H0)]
  } // massfunc*m: [1/volume]=[(c/H0)^-3]
  return massfunc(m,a)*m*u;
}

double I0j_y (int j, double k1, double k2, double k3, double k4,double a, int ni[4]){
  double array[11] = {k1,k2,k3,k4,0.,(double)j,a,(double)ni[0],(double)ni[1],(double)ni[2],(double)ni[3]};
  double result=0., logm;
  double logMmin=log(limits.M_min), logMmax=log(limits.M_max);
  double dlogm = (logMmax-logMmin)/200.;
  for(logm=logMmin;logm<=logMmax;logm+=dlogm){
    result += inner_I0j_y(logm, (void*)array);
  }
  return result*dlogm;
}



double inner_I1j_y (double logm, void *para){
  double *array = (double *) para;
  double m = exp(logm);
  long double u = 1.0;
  double a= array[6];
  double c = conc(m,a);

  double c_y;
#ifdef TEST_CALIB
  c_y = conc_y(m,a);
#else
  c_y = c;
#endif

  int l;
  int j = (int)(array[5]);
  double vol = m/(cosmology.rho_crit*cosmology.Omega_m);
  for (l = 0; l< j; l++){
    if(array[7+l]==-2.){ // ni=-2: y-field
      if(j==1){u *= (u_y_bnd(c_y,array[l],m,a)+u_y_ejc(m));}
      // if(j==1){u *= (u_y_bnd(c,array[l],m,a));}
      else{u *= u_y_bnd(c_y,array[l],m,a);}
    }else{
      u *= vol*u_nfw_c_interp(c,array[l],m,a);
    }
  }
  return massfunc(m,a)*m*u*B1(m,a);
}


double I_11_y (double k,double a){
  double array[11];
  array[0]=k;
  array[5]=1.0;
  array[6]=a;
  array[7]=-2.;
  double result=0., logm;
  double logMmin=log(limits.M_min), logMmax=log(limits.M_max);
  double dlogm = (logMmax-logMmin)/200.;
  for(logm=logMmin;logm<=logMmax;logm+=dlogm){
    result += inner_I1j_y(logm, (void*)array);
  }
  return result*dlogm;
}
// double I_11_y (double k,double a){
//   double array[11];
//   array[0]=k;
//   array[5]=1.0;
//   array[6]=a;
//   array[7]=-2.;

//   return int_gsl_integrate_medium_precision(inner_I1j_y,(void*)array,log(limits.M_min),log(limits.M_max),NULL, 2000);
// }
// double I_11_y (double k,double a){//look-up table for I11_y integral
//   static cosmopara C;
//   static double amin = 0, amax = 0,logkmin = 0., logkmax = 0., dk = 0., da = 0.;
  
//   static double **table_I1=0;
  
//   double aa,klog;
//   int i,j;
  
//   if (recompute_cosmo3D(C)){ //extend this by halo model parameters if these become part of the model
//     update_cosmopara(&C);
//     if (table_I1==0) table_I1 = create_double_matrix(0, Ntable.N_a-1, 0, Ntable.N_k_nlin-1);
//     double array[11];
//     array[5]=1.0;
//     array[7]=-2.;
    
//     amin = limits.a_min;
//     amax = 1.;
//     da = (amax - amin)/(Ntable.N_a);
//     aa = amin;
//     logkmin = log(limits.k_min_cH0);
//     logkmax = log(limits.k_max_cH0);
//     dk = (logkmax - logkmin)/(Ntable.N_k_nlin);
//     for (i=0; i<Ntable.N_a; i++, aa +=da) {
//       array[6] = fmin(aa,0.999);
//       klog  = logkmin;
//       for (j=0; j<Ntable.N_k_nlin; j++, klog += dk) {
//         array[0]= exp(klog);
//         table_I1[i][j] = log(int_gsl_integrate_medium_precision(inner_I1j_y,(void*)array,log(limits.M_min),log(limits.M_max),NULL, 2000));
//       }
//     }
//   }
//   aa = fmin(a,amax-1.1*da);//to avoid interpolation errors near z=0
//   return exp(interpol2d(table_I1, Ntable.N_a, amin, amax, da, aa, Ntable.N_k_nlin, logkmin, logkmax, dk, log(k), 1.0, 1.0));
// }

// double I1j_y (int j, double k1, double k2, double k3,double a, int ni[4]){
//   if (j ==1) {return I_11_y(k1,a);}
//   double array[11] = {k1,k2,k3,0.,0.,(double)j,a,(double)ni[0],(double)ni[1],(double)ni[2],(double)ni[3]};
//   return int_gsl_integrate_medium_precision(inner_I1j_y,(void*)array,log(limits.M_min),log(limits.M_max),NULL, 2000);
// }
double I1j_y (int j, double k1, double k2, double k3,double a, int ni[4]){
  if (j ==1) {return I_11_y(k1,a);}
  double array[11] = {k1,k2,k3,0.,0.,(double)j,a,(double)ni[0],(double)ni[1],(double)ni[2],(double)ni[3]};
  double result=0., logm;
  double logMmin=log(limits.M_min), logMmax=log(limits.M_max);
  double dlogm = (logMmax-logMmin)/200.;
  for(logm=logMmin;logm<=logMmax;logm+=dlogm){
    result += inner_I1j_y(logm, (void*)array);
  }
  return result*dlogm;
}


/*++++++++++++++++++++++++++++++++++++++++*
* 1Halo Terms                             *
*++++++++++++++++++++++++++++++++++++++++*/

double p_1h_hmcode(double k, double a)
{
  double sigma8;
  if (cosmology.sigma_8_cc!=0) sigma8 = cosmology.sigma_8_cc;
  else sigma8 = cosmology.sigma_8
  double k_star = 0.05618/pow(sigma_8*a,1.013)*cosmology.coverH0; // convert to code unit, Table 2, 2009.01858
  double sup = 1./(pow(k/k_star,-4)+1); // suppression at low-k, Eq.17, 2009.01858
  return I0j(2,k,k,0.,0.,a)*sup;
}

double p_1h_y(double k, double a, int n1, int n2)
{
  int ni[]={n1, n2, 0,0};
  // suppress low-k for 1h term same as 1h matter power
  double k_star = 0.05618/pow(cosmology.sigma_8*a,1.013)*cosmology.coverH0; // convert to code unit, Table 2, 2009.01858
  double sup = 1./(pow(k/k_star,-4)+1); // suppression at low-k, Eq.17, 2009.01858
  return I0j_y(2,k,k,0.,0.,a, ni)*sup;
}

/*++++++++++++++++++++++++++++++++++++++++*
* 2Halo Terms                              *
*++++++++++++++++++++++++++++++++++++++++*/

double Tk_EH_nowiggle(double k, double h, double wm, double wb, double T_CMB) {
//astro-ph:9709112
// These only needs to be calculated once
    double rb = wb/wm;     // Baryon ratio
    double e = exp(1.);    // e
    double s = 44.5 * log(9.83/wm) / sqrt(1. + 10. * pow(wb, 0.75));               // Equation (26)
    double alpha = 1. - 0.328 * log(431. * wm) * rb + 0.38 * log(22.3 * wm) * pow(rb, 2); // Equation (31)

// Functions of k
    double Gamma = (wm/h) * (alpha + (1. - alpha) / (1. + pow(0.43 * k[0] * s * h, 4))); // Equation (30)
    double q = k * pow(T_CMB/2.7, 2) / Gamma;   // Equation (28)
    double L = log(2. * e + 1.8 * q);      // Equation (29)
    double C = 14.2 + 731. / (1. + 62.5 * q);   // Equation (29)
    double Tk_nw = L / (L + C * pow(q, 2));    // Equation (29)
return Tk_nw;
}


void smooth_array_Gaussian(double *x, double *f, int n, double sigma) {
    int i, j;
    double *ff, weight, total;
    if (sigma != 0.) {
        ff = (double *)malloc(n * sizeof(double));
        // Save the original input array
        for (i = 0; i < n; i++) {
            ff[i] = f[i];
        }
        // Delete the original array
        for (i = 0; i < n; i++) {
            f[i] = 0.;
        }
        // Apply Gaussian smoothing
        for (i = 0; i < n; i++) {
            total = 0.;
            if (fabs(x[i] - x[0]) < nsig * sigma || fabs(x[i] - x[n-1]) < nsig * sigma) {
                f[i] = ff[i];
            } else {
                for (j = 0; j < n; j++) {
                    weight = exp(-pow(x[i] - x[j], 2) / (2. * pow(sigma, 2)));
                    f[i] = f[i] + ff[j] * weight;
                    total = total + weight;
                }
                f[i] = f[i] / total;
            }
        }

        free(ff);
    }
}




double Pk_smooth(double k) {
    //https://github.com/cmbant/CAMB/blob/master/fortran/halofit.f90
    static double logkmin_wiggle = log(5E-3);
    static double logkmax_wiggle = log(5.0);
    static int nk_wiggle = 512;
    static cosmopara C;
    static double *Pk_nowiggle;
    static double *Pk_ratio;
    static double * Pk_smooth;
    static double * k_array;
    static double dlnk = (logkmax_wiggle-logkmin_wiggle)/(nk_wiggle-1);

    double sigma = 0.25;
    double h = cosmology.h0;
    double omega_m = cosmology.Omega_m;
    double omega_b = cosmology.omb;
    double T_CMB = cosmology.Tcmb;
    double ns = cosmology.n_spec;

    if (recompute_cosmo3D(C)){
        update_cosmopara(&C);
            if (Pk_nowiggle==0){
                Pk_nowiggle  = create_double_vector(0, nk_wiggle-1);
                Pk_ratio = create_double_vector(0, nk_wiggle-1);
                Pk_smooth = create_double_vector(0, nk_wiggle-1);
                k_array = create_double_vector(0, nk_wiggle-1);
            }
        lnki = logkmin_wiggle; 
        for (kin=0; kin<nk_wiggle; kin++, lnki += dlnk) {
            Pk_nowiggle[kin]= pow(exp(lnki), ns) * pow(Tk_EH_nowiggle(exp(lnki), h, omega_m, omega_b, T_CMB), 2);
            Pk_ratio[kin] = p_lin(exp(lnki), 1)/Pk_nowiggle[kin];
            k_array[kin] = lnki
        }
        smooth_array_Gaussian(k_array, pk_ratio, nk_wiggle, sigma);
        for (kin=0; kin<nk_wiggle; kin++) {
            Pk_smooth[kin] = pk_ratio[kin]*Pk_nowiggle[kin]
        }
    }
    return interpol(Pk_smooth,  nk_wiggle, logkmin_wiggle, logkmax_wiggle, dlnk, log(k), 1.0,1.0 );
}



double p_2h_hmcode(double k, double a)
{
  return pow(I_11(k,a),2.0)*p_lin(k,a);
}








double p_2h_yy(double k, double a)
{
  assert(0);//not implemented 
  return pow(I_11_y(k,a),2.0)*p_lin(k,a);
}
double p_2h_my(double k, double a)
{
  assert(0);//not implemented 
  return I_11(k,a)*I_11_y(k,a)*p_lin(k,a);
}

/*********** Look-up table for halomodel matter power spectrum ************/
double Pdelta_hmcode_2020(double k,double a)
{
  static cosmopara C;
  static nuisancepara N;
  static double logkmin = 0., logkmax = 0., dk = 0., da = 0.;

  static int N_a = 20, N_k_nlin = 100;

  static double **table_P_NL=0;
  
  double klog,val,aa,kk;
  int i,j;
  // clock_t t1, t2; double dt;
  if (recompute_cosmo3D(C) || recompute_halo(N)){ //extend this by halo model parameters if these become part of the model
    update_cosmopara(&C); update_nuisance(&N);
    if (table_P_NL!=0) free_double_matrix(table_P_NL,0, N_a-1, 0, N_k_nlin-1);
    table_P_NL = create_double_matrix(0, N_a-1, 0, N_k_nlin-1);
    
    da = (0.999 - limits.a_min_hm)/(N_a*1.0-1);
    logkmin = log(limits.k_min_cH0);
    logkmax = log(limits.k_max_cH0);
    dk = (logkmax - logkmin)/(N_k_nlin*1.0-1);
    aa= limits.a_min_hm;
    for (i=0; i<N_a; i++, aa +=da) {
      if(aa>0.999) aa=.999;
      klog  = logkmin;
      // t1=clock();
      for (j=0; j<N_k_nlin; j++, klog += dk) {
        kk = exp(klog);
        table_P_NL[i][j] = log(p_1h(kk,aa) + p_2h(kk,aa));
        // table_P_NL[i][j] = log(pow(pow(p_1h(kk,aa),0.719) + pow(p_2h(kk,aa),0.719),1/0.719));
        if(isinf(table_P_NL[i][j])) table_P_NL[i][j] = -300.;
        // printf("%le %le %le \n", aa, kk, exp(table_P_NL[i][j]));

        // t1=clock();
        // p_1h(kk,aa);
        // t2=clock();
        // dt = (double)(t2 - t1) / CLOCKS_PER_SEC;
        // printf("time spent 1h %le\n",dt);
        // p_2h(kk,aa);
        // t1=clock();
        // dt = (double)(t1 - t2) / CLOCKS_PER_SEC;
        // printf("time spent 2h %le\n",dt);
      }
      // t2=clock();dt = (double)(t2 - t1) / CLOCKS_PER_SEC;
      // printf("%d, time: %le\n",i, dt);
    }
  }
  klog = log(k);
  if (klog < logkmin || klog >logkmax){return 0.;}
  if (a < limits.a_min_hm){return Pdelta_hmcode_2020(k,limits.a_min_hm)*pow(growfac(a)/growfac(limits.a_min_hm),2);}
  if (a > 0.999){return Pdelta_hmcode_2020(k,0.999)*pow(growfac(a)/growfac(0.999),2);}
  val = interpol2d(table_P_NL, N_a, limits.a_min_hm, 0.999, da, a, N_k_nlin, logkmin, logkmax, dk, klog, 0.0, 0.0);
  return exp(val);
}

double P_yy(double k,double a)
{
  static cosmopara C;
  static nuisancepara N;
  static double logkmin = 0., logkmax = 0., dk = 0., da = 0.;

  static int N_a = 20, N_k_nlin = 100;

  static double **table_P_NL=0;
  
  double klog,val,aa,kk;
  int i,j;
  double p1, p2;
  if (recompute_cosmo3D(C) || recompute_halo(N)){ //extend this by halo model parameters if these become part of the model
    update_cosmopara(&C); update_nuisance(&N);
    if (table_P_NL!=0) free_double_matrix(table_P_NL,0, N_a-1, 0, N_k_nlin-1);
    table_P_NL = create_double_matrix(0, N_a-1, 0, N_k_nlin-1);
    
    da = (0.999 - limits.a_min_hm)/(N_a*1.0-1);
    logkmin = log(limits.k_min_cH0);
    logkmax = log(limits.k_max_cH0);
    dk = (logkmax - logkmin)/(N_k_nlin*1.0-1);
    aa= limits.a_min_hm;
    for (i=0; i<N_a; i++, aa +=da) {
      if(aa>0.999) aa=.999;
      klog  = logkmin;
      for (j=0; j<N_k_nlin; j++, klog += dk) {
        kk = exp(klog);
        p1 = p_1h_y(kk,aa,-2,-2);
        p2 = p_2h_yy(kk,aa);
        table_P_NL[i][j] = log(p1 + p2);
        if(isinf(table_P_NL[i][j])) table_P_NL[i][j] = -300.;
        // printf("%le %le %le %le %le\n", aa, kk, p1, p2, exp(table_P_NL[i][j]));
      }
    }
  }
  klog = log(k);
  if (klog < logkmin || klog >logkmax){return 0.;}
  if (a < limits.a_min_hm){return 0.;}
  if (a > 0.999) {return P_yy(k,0.999);}
  val = interpol2d(table_P_NL, N_a, limits.a_min_hm, 0.999, da, a, N_k_nlin, logkmin, logkmax, dk, klog, 0.0, 0.0);
  return exp(val);
}

double P_my(double k,double a)
{
  static cosmopara C;
  static nuisancepara N;
  static double logkmin = 0., logkmax = 0., dk = 0., da = 0.;
  
  static int N_a = 20, N_k_nlin = 100;

  static double **table_P_NL=0;
  
  double klog,val,aa,kk;
  double p1,p2;
  int i,j;
  if (recompute_cosmo3D(C) || recompute_halo(N)){ //extend this by halo model parameters if these become part of the model
    update_cosmopara(&C); update_nuisance(&N);
    if (table_P_NL!=0) free_double_matrix(table_P_NL,0, N_a-1, 0, N_k_nlin-1);
    table_P_NL = create_double_matrix(0, N_a-1, 0, N_k_nlin-1);
    
    da = (0.999 - limits.a_min_hm)/(N_a*1.0-1);
    logkmin = log(limits.k_min_cH0);
    logkmax = log(limits.k_max_cH0);
    dk = (logkmax - logkmin)/(N_k_nlin*1.0-1);
    aa= limits.a_min_hm;
    for (i=0; i<N_a; i++, aa +=da) {
      if(aa>0.999) aa=.999;
      klog  = logkmin;
      // if(i==N_a-1){
        for (j=0; j<N_k_nlin; j++, klog += dk) {
          kk = exp(klog);
          p1 = p_1h_y(kk,aa,-2,0);
          p2 = p_2h_my(kk,aa);
          table_P_NL[i][j] = p1+p2;
          // table_P_NL[i][j] = log(p_1h_y(kk,aa,-2,0) + p_2h_my(kk,aa));
          if(isinf(table_P_NL[i][j])) table_P_NL[i][j] = -300.;
          // printf("%le %le %le %le %le \n", aa, kk,p1,p2, (table_P_NL[i][j]));
        // }exit(0);
      }
    }
  }
  klog = log(k);
  if (klog < logkmin || klog >logkmax){return 0.;}
  if (a < limits.a_min_hm){return 0.;}
  if (a > 0.999) {return P_my(k,0.999);}
  val = interpol2d(table_P_NL, N_a, limits.a_min_hm, 0.999, da, a, N_k_nlin, logkmin, logkmax, dk, klog, 0.0, 0.0);
  return (val);
}
#endif

/****************** Lookup table for 1-halo term super-sample covariance/halo sample variance term**********/


double I12_SSC_yy (double k,double a){//one-halo term contribution to super sample covariance
  static cosmopara C;
  static nuisancepara N;
  static double amin = 0, amax = 0,logkmin = 0., logkmax = 0., dk = 0., da = 0.;
  
  static double **table_I1=0;
  
  double aa,klog;
  int i,j;

  static int ni[4]={-2,-2,0,0};
  
  if (recompute_cosmo3D(C) || recompute_halo(N)){ //extend this by halo model parameters if these become part of the model
    update_cosmopara(&C); update_nuisance(&N);
    if (table_I1==0) table_I1 = create_double_matrix(0, Ntable.N_a-1, 0, Ntable.N_k_nlin-1);
    // double array[11];
    // array[5]=2.0;
    // array[7]=array[8]=-2.;
    
    amin = limits.a_min_hm;
    amax = 0.999;
    da = (amax - amin)/(Ntable.N_a-1.);
    aa = amin;
    logkmin = log(limits.k_min_cH0);
    logkmax = log(limits.k_max_cH0);
    dk = (logkmax - logkmin)/(Ntable.N_k_nlin-1.);
    for (i=0; i<Ntable.N_a; i++, aa +=da) {
      // array[6] = fmin(aa,0.999);
      klog  = logkmin;
      for (j=0; j<Ntable.N_k_nlin; j++, klog += dk) {
        // array[0]= exp(klog);array[1]= exp(klog);
        // table_I1[i][j] = log(int_gsl_integrate_medium_precision(inner_I1j_y,(void*)array,log(limits.M_min),log(limits.M_max),NULL, 2000));
        table_I1[i][j] = log(I1j_y(2,exp(klog),exp(klog),0.,aa, ni));
      }
    }
  }
  if (a < limits.a_min_hm){return 0.;}
  if (log(k) <logkmin){return I12_SSC_yy(exp(logkmin),a);};
  if (log(k) >logkmax){return 0.0;};
  aa = fmin(a,0.999);
  double res = exp(interpol2d(table_I1, Ntable.N_a, amin, amax, da, aa, Ntable.N_k_nlin, logkmin, logkmax, dk, log(k), 1.0, 1.0));
  if(isnan(res)) {res=0.;}
  return res;
}

double I12_SSC_my (double k,double a){//one-halo term contribution to super sample covariance
  static cosmopara C;
  static nuisancepara N;
  static double amin = 0, amax = 0,logkmin = 0., logkmax = 0., dk = 0., da = 0.;
  
  static double **table_I1=0;
  
  double aa,klog;
  int i,j;

  static int ni[4]={-2,0,0,0};

  if (recompute_cosmo3D(C) || recompute_halo(N)){ //extend this by halo model parameters if these become part of the model
    update_cosmopara(&C); update_nuisance(&N);
    if (table_I1==0) table_I1 = create_double_matrix(0, Ntable.N_a-1, 0, Ntable.N_k_nlin-1);
    // double array[11];
    // array[5]=2.0;
    // array[7]=-2.;
    // array[8]=0.;
    
    amin = limits.a_min_hm;
    amax = 0.999;
    da = (amax - amin)/(Ntable.N_a-1.);
    aa = amin;
    logkmin = log(limits.k_min_cH0);
    logkmax = log(limits.k_max_cH0);
    dk = (logkmax - logkmin)/(Ntable.N_k_nlin-1.);
    for (i=0; i<Ntable.N_a; i++, aa +=da) {
      // array[6] = fmin(aa,0.999);
      klog  = logkmin;
      for (j=0; j<Ntable.N_k_nlin; j++, klog += dk) {
        // array[0]= exp(klog);array[1]= exp(klog);
        // table_I1[i][j] = log(int_gsl_integrate_medium_precision(inner_I1j_y,(void*)array,log(limits.M_min),log(limits.M_max),NULL, 2000));
        table_I1[i][j] = log(I1j_y(2,exp(klog),exp(klog),0.,aa, ni));
        if(isinf(table_I1[i][j])) table_I1[i][j] = -300.;
      }
    }
  }
  if (a < limits.a_min_hm){return 0.;}
  if (log(k) <logkmin){return I12_SSC_my(exp(logkmin),a);};
  if (log(k) >logkmax){return 0.0;};
  aa = fmin(a,0.999);
  double res = exp(interpol2d(table_I1, Ntable.N_a, amin, amax, da, aa, Ntable.N_k_nlin, logkmin, logkmax, dk, log(k), 1.0, 1.0));
  if(isnan(res)) {res=0.;}
  return res;
}





