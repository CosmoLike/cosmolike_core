#ifndef __HALO_FAST__
#define __HALO_FAST__
#include <time.h>
/*relations for converting mass -> radius, using \Delta = 200 \rho_m for consistency with mass & bias function */
double m_Delta(double r, double a);
double r_Delta(double m, double a);

double bias_norm(double a);
/* fitting functions for mass and bias */
double massfunc(double m, double a);
double bias_norm(double a);
double mass_norm(double a);
double B1 (double m, double a);
/* halo density profiles */
double conc(double m, double a);//mass-concentration relation
double u_nfw_c(double c,double k, double m, double a); //analytic expression for Fourier transform of NFW density profile
double u_nfw(double c,double k, double m, double a);   // Fourier transform of NFW density profile, truncated at r_Delta, through direct integration
// use this routine in inner I[0-1]j and modify int_rho_nfw if to work with different halo profiles, e.g. to include adiabatic contraction, AGN feedback, etc.
/*halo model building blocks */
double I0j (int j, double k1, double k2, double k3, double k4, double a);
double I0j_y (int j, double k1, double k2, double k3, double k4, double a, int is_y[4]);
double I_11 (double k,double a);
double I1j (int j, double k1, double k2, double k3, double a);
/*halo model matter power spectrum, bispectrum, trispectrum*/
double p_1h(double k, double a);
double p_2h(double k, double a);
double Pdelta_halo(double k, double a);
double P_my(double k,double a);
/*look up table for 1-h halo sample variance term */
double I12_SSC (double k,double a);
/* simple abundance matching routines for source galaxies */
double n_s_cmv (double a); //comoving source galaxy density based on n_gal +n(z)
double b_ngmatched(double a, double n_cmv); //simple abundance matching routine to get bias for source galaxies from redshift distribution
double b_source(double a); //lookup table for b1 of source galaxies

/*==============================================================*/

/* ********************* Routines for halo model ************* */

/************** begin halo properties *****************/

/********* routines to convert scales/radii to halo masses, and vice versa *********/

double delta_c(double a) /*set to 1.686 for consistency with Tinker mass & bias definition */
{
  return 1.686;
}
double delta_Delta(double a)
{
  return 200.0;//using Tinker et al mass & bias functions with \Delta = 200 rho_{mean}
}

double rho_Delta(double a) //virial density in solar masses/h/(H0/c)^3
{
  //using Tinker et al mass & bias functions with \Delta = 200 rho_{mean}
  return (delta_Delta(a) *cosmology.rho_crit*cosmology.Omega_m);
}

double m_Delta(double r, double a){
  return 4.*M_PI/3.0*pow(r,3.0)*rho_Delta(a);
}
double r_Delta(double m, double a) //calculate r_Delta in c/H0 given m in (solar masses/h)
{
  return pow(3./(4.*M_PI)*(m/rho_Delta(a)),1./3.);
}

// double r_s(double m, double a)
// {
//   return r_Delta(m,a)/conc(m,a);
// }

double radius(double m)
{
  return pow(3./4.*m/(M_PI*cosmology.rho_crit*cosmology.Omega_m),1./3.);
}

/*++++++++++++++++++++++++++++++++++++++++*
*  Variance of the density field         *
*++++++++++++++++++++++++++++++++++++++++*/

double sigma2_integrand(double x, void * params)   // inner integral
{
  double *array = (double*)params;
  double k= x/array[0];
  //refactored FT of spherical top-hat to avoid numerica divergence of 1/x
  return p_lin(k,1.0)*pow(3.*gsl_sf_bessel_j1(x)/array[0],2.)/(array[0]*2.*M_PI*M_PI);
}
double sigma2(double m)
{
  static cosmopara C;
  
  static double *table_S2;
  static double dm = .0, logmmin =1., logmmax = 1.;
  double mlog,result,array[1],x1;
  int i,j;
  if (recompute_cosmo3D(C)){
    update_cosmopara(&C);
    if (table_S2==0){
      table_S2  = create_double_vector(0, Ntable.N_S2-1);
      logmmin = log(limits.M_min/2.0);
      logmmax = log(limits.M_max*2.0);
      dm = (logmmax - logmmin)/(Ntable.N_S2);
    }
    mlog = logmmin;
    
    for (i=0; i<Ntable.N_S2; i++, mlog += dm) {
      array[0] = radius(exp(mlog));
      result  = int_gsl_integrate_medium_precision(sigma2_integrand,(void*)array,0.,14.1,NULL,1000);
      table_S2[i]=log(result);
    }
  }
  return exp(interpol(table_S2, Ntable.N_S2, logmmin, logmmax, dm,log(m), 1.0,1.0 ));
}

double nu(double m, double a){
  static int init = 1;
  if ((cosmology.Omega_nu > 0 || cosmology.M_nu >0) && init){
    fprintf(stderr,"halo.c does not support cosmologies with massive neutrinos\n");
//    exit(1);
    init =0;
  }
	return delta_c(a)/(sqrt(sigma2(m))*growfac(a)/growfac(1.));
}

double nulogm_a1(double lgm,void * params)
{
  (void)(params);
  return log(nu(exp(lgm),1.));
}


double dlognudlogm(double m)
{
  static cosmopara C;
  
  static double *table_DS;
  static double dm = .0, logmmin = 1.0,logmmax = 1.0;
  double mlog;
  int i;
  if (recompute_cosmo3D(C)){
    update_cosmopara(&C);

    if (table_DS==0){
      table_DS  = create_double_vector(0, Ntable.N_DS-1);
      logmmin = log(limits.M_min/2.0);
      logmmax = log(limits.M_max*2.0);
      dm = (logmmax - logmmin)/(Ntable.N_DS);
    }
    gsl_function F;
    double result, abserr;
  
    F.function = &nulogm_a1;
    F.params = 0;
    mlog = logmmin;
    for (i=0; i<Ntable.N_DS; i++, mlog += dm) {
      gsl_deriv_central (&F,mlog, 0.1*mlog, &result, &abserr);
      table_DS[i]=result;
    }
  }
  return interpol(table_DS, Ntable.N_DS, logmmin, logmmax, dm,log(m), 1.0,1.0);
}

/*************** mass function & halo bias (Tinker et al.) **********/
double f_tinker(double n, double a_in)
{  //Eqs. (8-12) + Table 4 from Tinker et al. 2010
  //aa = alpha, b = beta, c = gamma, f = phi, e = eta
  double a = fmax(0.25,a_in); //limit fit range of mass function evolution to z <= 3, as discussed after Eq. 12 of http://arxiv.org/pdf/1001.3162v2.pdf
	double aa,b,c,f,e;
  aa = 0.368;
  b = 0.589*pow(a,-0.2); c = 0.864*pow(a,0.01); f = -0.729*pow(a,.08); e = -0.243*pow(a,-0.27);
	return aa*(1.0+pow(b*n,-2.0*f))*pow(n,2.0*e)*exp(-c*n*n/2.0);
}
double fnu_tinker(double n, double a)
{
  return f_tinker(n,a)*n;
}

double B1_nu (double n,double a){
  // Eq (6) + Table 2 from Tinker et al. 2010
	double A,aa,B,b,C,c,y;
	y= log10(200.0);
	A = 1.0+0.24*y*exp(-pow(4.0/y,4.0));
	aa = 0.44*y-0.88; B = 0.183; b = 1.5;
	C = 0.019+0.107*y+0.19*exp(-pow(4.0/y,4.0)); c =2.4;
	return 1.0-A*pow(n,aa)/(pow(n,aa)+pow(delta_c(a),aa)) + B*pow(n,b)+C*pow(n,c);
}

// normalization for redshift evolution of mass function (Eq. (8) in Tinker et al. 2010), enforcing mass to be unbiased on average
double mass_norm_integrand(double n, void * params)   // inner integral
{
  double *array = (double*)params;
  double a = array[0];
  return f_tinker(n,a);
}

double mass_norm(double a)
{
  static cosmopara C;
  static double *table_BN;
  static double da = 0., amin = 0., amax =0.;
  double aa,result,array[1];
  int i;
  if (recompute_cosmo3D(C)){
    update_cosmopara(&C);
    if (table_BN==0){
      table_BN  = create_double_vector(0, Ntable.N_a-1);
       amin = limits.a_min;
      amax = 1.;
      da = (amax - amin)/(Ntable.N_a-1.);
    }
    aa= amin;
    for (i=0; i<Ntable.N_a; i++, aa += da) {
      array[0] = aa;
      result = int_gsl_integrate_medium_precision(mass_norm_integrand,(void*)array,0.0,10,NULL,5000);
      table_BN[i]=result;
      }
    }
  return interpol(table_BN, Ntable.N_a, amin, amax, da,fmin(a,amax-da), 1.0,1.0 );
}

//also correct bias for to halo mass cuts (so that large-scale 2-h term matches PT results at all redshifts)

double bias_norm_integrand (double n,void * params){
  double *array = (double*)params;
  double a = array[0];  
  return B1_nu(n,a)*f_tinker(n,a);
}

double bias_norm(double a)
{
	static cosmopara C;
	static double *table_BN;
	static double da = .0, amin = 0., amax =0.;
	double aa,result,array[1];
	int i;
  if (recompute_cosmo3D(C)){
    update_cosmopara(&C);
  	if (table_BN==0){
      table_BN  = create_double_vector(0, Ntable.N_a-1);
      amin = limits.a_min; //1./(1+redshift.shear_zdistrpar_zmax);
      amax = 1.;
      da = (amax - amin)/(Ntable.N_a-1.);
    }
		aa= amin;		
		for (i=0; i<Ntable.N_a-1; i++, aa += da) {
			array[0] = aa;
			result = int_gsl_integrate_medium_precision(bias_norm_integrand,(void*)array,nu(limits.M_min,aa),nu(limits.M_max,aa),NULL,5000);
			table_BN[i]=result;
		}
    table_BN[Ntable.N_a-1] = 1.;
	}
	return interpol(table_BN, Ntable.N_a, amin,amax, da,fmin(a,amax-da), 1.0,1.0 );

//   printf("yup2\n");
//	double res = interpol(table_BN, Ntable.N_a, amin,amax, da,fmin(a,amax-da), 1.0,1.0 );
//   printf("yup3\n");
//   return res;

}

double massfunc(double m, double a){
	return fnu_tinker(nu(m,a),a)*cosmology.rho_crit*cosmology.Omega_m/m/m*dlognudlogm(m);
}


double B1 (double m,double a){ //b(m,a) based on redshift evolution fits in Tinker et al. paper, no additional normalization
	return B1_nu(nu(m,a),a);
}
double B1_normalized (double m,double a){ //divide by bias norm only in matter spectra, not in HOD modeling/cluster analyses
	return B1_nu(nu(m,a),a)/bias_norm(a);

//   printf("yok0\n");
//   double res = B1_nu(nu(m,a),a);
//   printf("yok1\n");
//   res *= bias_norm(a);
//   printf("yok2\n");
//   return res;
}
/***************** halo profiles ***************/
/******** mass-concentration relation **********/
double conc(double m, double a) // for matter
{
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

/***********  FT of NFW Profile **************/

double int_rho_nfw(double r, void *params){//Fourier kernel * NFW profile, integrand for u_nfw 
  double *array = (double*)params;
  double k = array[0];
  double c = array[3];
  double rv = array[4];
  double rs =rv/c;
  return r*sinl(k*r)/k*pow(rs,-3.)*1./(log(1.+c)-c/(1.+c))/(r/rs)*pow(1+r/rs,-2.);
}

double u_nfw(double c, double k, double m, double a){
  // FFT density of NFW profile, truncated at r_Delta, through direct integration
  // use this routine in inner I[0-1]j and modify int_rho_nfw to work with different halo profiles
  
  double array[5] ={k,m,a,c,r_Delta(m,a)};
  return int_gsl_integrate_medium_precision(int_rho_nfw, (void*)array, 0,array[4],NULL, 5000);
}


double u_nfw_c(double c,double k, double m, double aa){// analytic FT of NFW profile, from Cooray & Sheth 01
  double x, xu;
  x = k * r_Delta(m,aa)/c;
  xu = (1.+c)*x;
  return (sin(x)*(gsl_sf_Si(xu)-gsl_sf_Si(x))- sinl(c*x)/xu +cos(x)*(gsl_sf_Ci(xu)-gsl_sf_Ci(x)))*1./(log(1.+c)-c/(1.+c));
}

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

double inner_I0j (double logm, void *para){
  double *array = (double *) para;
  double m = exp(logm);
  long double u = 1.0;
  double a= array[6];
  double c = conc(m,a);
  int l;
  int j = (int)(array[5]);
  for (l = 0; l< j; l++){
#ifdef SLOW
    u = u*u_nfw_c(c,array[l],m,a);
#else
    u = u*u_nfw_c_interp(c,array[l],m,a);
#endif
  }
  return massfunc(m,a)*m*pow(m/(cosmology.rho_crit*cosmology.Omega_m),(double)j)*u;
}

#ifdef SLOW
double I0j (int j, double k1, double k2, double k3, double k4,double a){
  double array[7] = {k1,k2,k3,k4,0.,(double)j,a};
  return int_gsl_integrate_medium_precision(inner_I0j,(void*)array,log(limits.M_min),log(limits.M_max),NULL, 2000);
}
#else
double I0j (int j, double k1, double k2, double k3, double k4,double a){
  double array[7] = {k1,k2,k3,k4,0.,(double)j,a};
  double result=0., logm;
  double logMmin=log(limits.M_min), logMmax=log(limits.M_max);
  double dlogm = (logMmax-logMmin)/200.;
  for(logm=logMmin;logm<=logMmax;logm+=dlogm){
    result += inner_I0j(logm, (void*)array);
  }
  return result*dlogm;
}
#endif

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
// double I0j_y (int j, double k1, double k2, double k3, double k4,double a, int ni[4]){
//   double array[11] = {k1,k2,k3,k4,0.,(double)j,a,(double)ni[0],(double)ni[1],(double)ni[2],(double)ni[3]};
//   return int_gsl_integrate_medium_precision(inner_I0j_y,(void*)array,log(limits.M_min),log(limits.M_max),NULL, 2000);
// }


double inner_I1j (double logm, void *para){
  double *array = (double *) para;
  double m = exp(logm);
  long double u = 1.0;
  double a= array[6];
  double c = conc(m,a);
  int l;
  int j = (int)(array[5]);
  for (l = 0; l< j; l++){
#ifdef SLOW
    u = u*u_nfw_c(c,array[l],m,a);
#else
    u = u*u_nfw_c_interp(c,array[l],m,a);
#endif
  }
  return massfunc(m,a)*m*pow(m/(cosmology.rho_crit*cosmology.Omega_m),(double)j)*u*B1_normalized(m,a);
}

// double I_11 (double k,double a){//look-up table for I11 integral
//   double array[7];
//   array[0]=k;
//   array[5]=1.0;
//   array[6]=a;
//   return int_gsl_integrate_medium_precision(inner_I1j,(void*)array,log(limits.M_min),log(limits.M_max),NULL, 2000);
// }
#ifdef SLOW
double I_11 (double k,double a){//look-up table for I11 integral
  static cosmopara C;
  static nuisancepara N;
  static double amin = 0, amax = 0,logkmin = 0., logkmax = 0., dk = 0., da = 0.;
  
  static double **table_I1=0;
  
  double aa,klog;
  int i,j;
  
  if (recompute_cosmo3D(C) || recompute_halo(N)){ //extend this by halo model parameters if these become part of the model
    update_cosmopara(&C); update_nuisance(&N);
    if (table_I1==0) table_I1 = create_double_matrix(0, Ntable.N_a-1, 0, Ntable.N_k_nlin-1);
    double array[7];
    array[5]=1.0;
    
    amin = limits.a_min;
    amax = 1.;
    da = (amax - amin)/(Ntable.N_a);
    aa = amin;
    logkmin = log(limits.k_min_cH0);
    logkmax = log(limits.k_max_cH0);
    dk = (logkmax - logkmin)/(Ntable.N_k_nlin);
    for (i=0; i<Ntable.N_a; i++, aa +=da) {
      array[6] = fmin(aa,0.999);
      klog  = logkmin;
      for (j=0; j<Ntable.N_k_nlin; j++, klog += dk) {
        array[0]= exp(klog);
        table_I1[i][j] = log(int_gsl_integrate_medium_precision(inner_I1j,(void*)array,log(limits.M_min),log(limits.M_max),NULL, 2000));
      }
    }
  }
  aa = fmin(a,amax-1.1*da);//to avoid interpolation errors near z=0
  return exp(interpol2d(table_I1, Ntable.N_a, amin, amax, da, aa, Ntable.N_k_nlin, logkmin, logkmax, dk, log(k), 1.0, 1.0));
}
#else
double I_11 (double k,double a){
  double array[7];
  array[0]=k;
  array[5]=1.0;
  array[6]=a;
  double result=0., logm;
  double logMmin=log(limits.M_min), logMmax=log(limits.M_max);
  double dlogm = (logMmax-logMmin)/200.;
  for(logm=logMmin;logm<=logMmax;logm+=dlogm){
    result += inner_I1j(logm, (void*)array);
  }
  return result*dlogm;
}
#endif

#ifdef SLOW
double I1j (int j, double k1, double k2, double k3,double a){
  if (j ==1) {return I_11(k1,a);}
  double array[7] = {k1,k2,k3,0.,0.,(double)j,a};
  return int_gsl_integrate_medium_precision(inner_I1j,(void*)array,log(limits.M_min),log(limits.M_max),NULL, 2000);
}
#else
double I1j (int j, double k1, double k2, double k3,double a){
  if (j ==1) {return I_11(k1,a);}
  double array[7] = {k1,k2,k3,0.,0.,(double)j,a};
  double result=0., logm;
  double logMmin=log(limits.M_min), logMmax=log(limits.M_max);
  double dlogm = (logMmax-logMmin)/200.;
  for(logm=logMmin;logm<=logMmax;logm+=dlogm){
    result += inner_I1j(logm, (void*)array);
  }
  return result*dlogm;
}
#endif

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

double p_1h(double k, double a)
{
  double k_star = 0.05618/pow(cosmology.sigma_8*a,1.013)*cosmology.coverH0; // convert to code unit, Table 2, 2009.01858
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
* 2Halo Terms                         	   *
*++++++++++++++++++++++++++++++++++++++++*/
double p_2h(double k, double a)
{
  return pow(I_11(k,a),2.0)*p_lin(k,a);
}

double p_2h_yy(double k, double a)
{
  return pow(I_11_y(k,a),2.0)*p_lin(k,a);
}
double p_2h_my(double k, double a)
{
  return I_11(k,a)*I_11_y(k,a)*p_lin(k,a);
}

/******* halomodel power without tabulation********/
#ifdef NOTAB
double Pdelta_halo(double k,double a){
  return p_1h(k,a) + p_2h(k,a);
}
double P_yy(double k,double a){
  return p_1h_y(k,a,-2,-2) + p_2h_yy(k,a);
}
double P_my(double k,double a){
  return p_1h_y(k,a,-2,0) + p_2h_my(k,a);
}

#else
/*********** Look-up table for halomodel matter power spectrum ************/
double Pdelta_halo(double k,double a)
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
  if (a < limits.a_min_hm){return Pdelta_halo(k,limits.a_min_hm)*pow(growfac(a)/growfac(limits.a_min_hm),2);}
  if (a > 0.999){return Pdelta_halo(k,0.999)*pow(growfac(a)/growfac(0.999),2);}
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
double I12_SSC (double k,double a){//one-halo term contribution to super sample covariance
  static cosmopara C;
  static nuisancepara N;
  static double amin = 0, amax = 0,logkmin = 0., logkmax = 0., dk = 0., da = 0.;
  
  static double **table_I1=0;
  
  double aa,klog;
  int i,j;

  if (recompute_cosmo3D(C) || recompute_halo(N)){ //extend this by halo model parameters if these become part of the model
    update_cosmopara(&C); update_nuisance(&N);
    if (table_I1==0) table_I1 = create_double_matrix(0, Ntable.N_a-1, 0, Ntable.N_k_nlin-1);
    // double array[7];
    // array[5]=2.0;
    
    amin = limits.a_min;
    // amax = 1.-1.e-5;
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
        // table_I1[i][j] = log(int_gsl_integrate_medium_precision(inner_I1j,(void*)array,log(limits.M_min),log(limits.M_max),NULL, 2000));
        table_I1[i][j] = log(I1j(2,exp(klog),exp(klog),0.,aa));
      }
    }
  }
  if (log(k) <logkmin){return I12_SSC(exp(logkmin),a);};
  if (log(k) >logkmax){return 0.0;};
  aa = fmin(a,0.999);
  double res = exp(interpol2d(table_I1, Ntable.N_a, amin, amax, da, aa, Ntable.N_k_nlin, logkmin, logkmax, dk, log(k), 1.0, 1.0));
  if(isnan(res)) {res=0.;}
  return res;
}

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




/****************** simplistic abundance matching to get approximate b(z) of source galaxies **************/
double n_s_cmv (double a){
  double dV_dz = pow(f_K(chi(a)),2.0)/hoverh0(a);//comoving dV/dz(z = 1./a-1) per radian^2
  return zdistr_photoz(1./a-1.,-1)*survey.n_gal*survey.n_gal_conversion_factor/dV_dz;//dN/dz/radian^2/(dV/dz/radian^2)
}

double int_for_ngmatched (double logm, void *para){
  double *array = (double *) para;
  return massfunc(exp(logm),array[0])*exp(logm);
}
double int_for_b_ngmatched (double logm, void *para){
  double *array = (double *) para;
  return massfunc(exp(logm),array[0])*B1(exp(logm),array[0])*exp(logm);
}

double b_ngmatched(double a, double n_cmv){
  double array[1] = {a};
  double mlim = log(5.e+14);
  double dm = 0.1;
  double n = int_gsl_integrate_medium_precision(int_for_ngmatched, (void*)array, mlim,log(limits.M_max),NULL, 5000);
  while (n < n_cmv){
    mlim -=dm;
    n += int_gsl_integrate_medium_precision(int_for_ngmatched, (void*)array, mlim,mlim+dm,NULL, 5000);
  }
  return int_gsl_integrate_medium_precision(int_for_b_ngmatched, (void*)array, mlim,log(limits.M_max),NULL, 5000)/n;
}

double b_source(double a){ //lookup table for b1 of source galaxies
  
  static double *table_b;
  static double ngal = 0,mag_lim = 0.;
  static double da =0., amin = 0., amax =1.;
  if (mag_lim != survey.m_lim || ngal != survey.n_gal){//if survey/redshift distribution changes, hopefully one of these will change too
    //uncertainty in abundance matching is larger than change with cosmology, so save time and don't recompute if only cosmology changes
    double aa;
    int i;
    
    ngal = survey.n_gal;
    mag_lim = survey.m_lim;
     amin = limits.a_min;// 1./(1.+redshift.shear_zdistrpar_zmax);
    amax= 1./(1.+redshift.shear_zdistrpar_zmin);
    da = (amax - amin)/(1.*Ntable.N_a);
    
    if (table_b==0){
      table_b  = create_double_vector(0, Ntable.N_a-1);
    }
    
    aa= amin;
    for (i=0; i<Ntable.N_a; i++, aa += da) {
      table_b[i] = b_ngmatched(aa,n_s_cmv(aa));
    }
  }
  return interpol(table_b, Ntable.N_a, amin, amax, da,a, 1.0,1.0 );
}
#endif