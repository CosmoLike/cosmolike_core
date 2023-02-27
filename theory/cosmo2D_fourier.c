//#include "pt.c"
#include "cosmo2D_exact_fft.c"
#include "cosmo2D_fullsky_TATT.c"



//double W_kappa(double a, double fK, double nz);//complete lens efficiency weight
//double W_source(double a, double nz); //source redshift distribution (radial weight for IA,source clustering)
//double W_gal(double a, double nz); //complete weight for galaxy statistics
//double W_HOD(double a, double nz); //galaxy weigth without bias factor (for projecting P_gg instead of P_nl)

//double C_cl_tomo(double l, int ni, int nj);  //galaxy clustering power spectrum of galaxies in bins ni, nj
//double C_cl_tomo_nointerp(double l, int ni, int nj);
//double C_cl_HOD(double l, int ni);  //galaxy clustering power spectrum of galaxies in bin ni, using HOD model

double C_gl_tomo(double l, int ni, int nj);  //G-G lensing power spectrum, lens bin ni, source bin nj
double C_gl_tomo_nointerp(double l, int ni, int nj);
double C_gl_HOD_tomo(double l, int ni, int nj);  //G-G lensing power spectrum from HOD model, lens bin ni, source bin nj

double C_shear_tomo(double l, int ni, int nj); //shear tomography power spectra
double C_shear_tomo_nointerp(double l, int ni, int nj);
/**********************************************************/



/*********** Limber approximation integrands for angular power spectra *************/
// use W_i W_j/fK^2 dchi_da(a) as Limber weight for power spectra



double int_for_C_cl_tomo_RSD(double k, void *params)
{


  double res,ell, chi_0, chi_1, a_0, a_1;
//  k = exp(lgk);
  double *ar = (double *) params;
  ell       = ar[2]+0.5;
  chi_0 = f_K(ell/k);
  chi_1 = f_K((ell+1.)/k);
  if (chi_1 > chi(limits.a_min)){
    return 0;}
  a_0 = a_chi(chi_0);
  a_1 = a_chi(chi_1);
//  res=(W_gal(a_0,ar[0]))*(W_gal(a_0,ar[1]));
  //res=(W_gal(a_0,ar[0])+W_RSD(ell,a_0,a_1,ar[0]))*(W_gal(a_0,ar[1])+W_RSD(ell,a_0,a_1,ar[1]));
  //res= res*Pdelta(k,a_0);
  double p_0 = Pdelta(k,a_0);
  double p_c_0 = p_0;
  double p_1 = Pdelta(k,a_1);
  double p_c_1 = p_1;
  double f_cb = 1.0-cosmology.Omega_nu/cosmology.Omega_m;

  double b1_0 = gbias.b1_function(1./a_0-1.,(int)ar[0]);
  double b1_1 = gbias.b1_function(1./a_1-1.,(int)ar[1]);

  if (gbias.cluster_tracer>0.0){
    if (gbias.neutrino_induced_sdb>0.0){
      b1_0 *= ((1.0 + p_lin_cluster(k,a_0)/p_lin(k,a_0) * f_cb)/(1.0+f_cb));
      b1_1 *= ((1.0 + p_lin_cluster(k,a_1)/p_lin(k,a_1) * f_cb)/(1.0+f_cb));

    }
  }
  else{
    if (gbias.neutrino_induced_sdb>0.0){
      b1_0 *= bias_factor_mult((int)ar[0],k);
      b1_1 *= bias_factor_mult((int)ar[1],k);
    }
    p_c_0 = Pdelta_cluster(k,a_0);
    p_c_1 = Pdelta_cluster(k,a_1);
  }



  //double b1_k_0 = gbias.b1_function(1./a_0-1.,(int)ar[0])* (1.0 + p_lin_cluster(k,a_0)/p_lin(k,a_0) * f_cb)/(1.0+f_cb);
  res = (b1_0*W_HOD(a_0, ar[0])*sqrt(p_c_0)+gbias.b_mag[(int)ar[0]]*W_mag(a_0, chi_0, ar[0])*sqrt(p_0)+W_RSD(ell,a_0,a_1,ar[0])*sqrt(p_c_0));
  res *=(b1_1*W_HOD(a_1, ar[1])*sqrt(p_c_1)+gbias.b_mag[(int)ar[1]]*W_mag(a_1, chi_1, ar[1])*sqrt(p_1)+W_RSD(ell,a_0,a_1,ar[1])*sqrt(p_c_1));
  
  
  //uncomment above lines to implement scale-dependent neutrino bias
  

  return res;
}

double int_for_C_cl_HOD(double a, void *params)
{
  double res,ell, fK, k;
  double *ar = (double *) params;
  if (a >= 1.0) error("a>=1 in int_for_p_2");
  
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  
  res= W_HOD(a,ar[0])*W_HOD(a,ar[0])*dchi_da(a)/fK/fK;
  if (res !=0){res= res*P_gg(k,a,(int)ar[0]);}
  return res;
}


double int_for_C_gl_tomo(double a, void *params) // Add RSD
{
  return int_for_C_ggl_IA_TATT(a, params);
  double *ar = (double *) params;
  double res,ell, fK, k;
  if (a >= 1.0) error("a>=1 in int_for_C_gl_tomo");

  double ell_prefactor1 = (ar[2])*(ar[2]+1.);
  double ell_prefactor2 = (ar[2]-1.)*ell_prefactor1*(ar[2]+2.);
  if(ell_prefactor2<=0.) 
    ell_prefactor2=0.;
  else
    ell_prefactor2=sqrt(ell_prefactor2);

  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  double chi_0,chi_1,a_0,a_1;
  chi_0 = f_K(ell/k);
  chi_1 = f_K((ell+1.)/k);
  if (chi_1 > chi(limits.a_min)){
    return 0;}
  a_0 = a_chi(chi_0);
  a_1 = a_chi(chi_1);

  double wgal = W_gal(a,ar[0]);
  wgal += W_mag(a,fK,ar[0])*(ell_prefactor1/ell/ell -1.) ;
  // res= (wgal+W_RSD(ell, a_0, a_1, ar[0]))*W_kappa(a,fK,ar[1])*dchi_da(a)/fK/fK  * ell_prefactor2/ell/ell;
  res= (wgal)*W_kappa(a,fK,ar[1])*dchi_da(a)/fK/fK  * ell_prefactor2/ell/ell;
  // res= (W_gal(a,ar[0]) + W_mag(a,fK,ar[0])*gbias.b_mag[(int)ar[0]] *(ell_prefactor1/ell/ell-1.) )*W_kappa(a,fK,ar[1])*dchi_da(a)/fK/fK /fK/fK * ell_prefactor2/k/k;
  res= res*Pdelta(k,a);
  return res;
}

double int_for_C_gl_HOD_tomo(double a, void *params)
{
  double *ar = (double *) params;
  double res,ell, fK, k;
  if (a >= 1.0) error("a>=1 in int_for_p_2");
  
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  
  res= W_HOD(a,ar[0])*W_kappa(a,fK,ar[1])*dchi_da(a)/fK/fK;
  if (res !=0){res= res*P_gm(k,a, (int) ar[0]);}
  return res;
}

double int_for_C_shear_tomo(double a, void *params)
{
  double *ar = (double *) params;
  double res,ell, fK, k;
  if (a >= 1.0) error("a>=1 in int_for_C_shear_tomo");
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  res= W_kappa(a,fK,ar[0])*W_kappa(a,fK,ar[1])*dchi_da(a)/fK/fK;
  res= res*Pdelta(k,a); 
  return res;
}

/*********** angular power spectra - without look-up tables ******************/
double C_cl_RSD_nointerp(double l, int ni, int nj)  //galaxy clustering power spectrum of galaxy bins ni, nj
{ 
  double array[3] = {1.0*ni,1.0*nj,l};
  if (gbias.b2[ni] || gbias.b2[nj]){
          printf("\nCalled C_cl_RSD_nointerp(l,z1=%d,z2=%d) with non-linear bias parameters set.\n",ni,nj);
          printf("RSD beyond linear bias not yet supported.\n");
          printf("Use linear bias only for clustering + RSD\n\n");
          exit(1);
  }
  double chi_max = chi(fmax(amin_lens(ni),amin_lens(nj)));
  double chi_min = chi(fmin(amax_lens(ni),amax_lens(nj)));
  double k_min = (l+0.5)/f_K(chi_max);
  double k_max = (l+0.5)/f_K(chi_min);
  return int_gsl_integrate_medium_precision(int_for_C_cl_tomo_RSD,(void*)array,k_min,k_max,NULL,1000)*2./(2*l+1.);
}



double C_gl_tomo_nointerp(double l, int ni, int nj)  //G-G lensing power spectrum, lens bin ni, source bin nj
{
  // l=1.0*2;
  if(l==1.) {return 0.;}
  double array[3] = {(double)ni,(double)nj,l};
  if (gbias.b2[ni] || gbias.b2[nj]){
    return int_gsl_integrate_low_precision(int_for_C_gl_tomo_b2,(void*)array,amin_lens(ni),amax_lens(ni),NULL,1000);
  }
  if (gbias.hod[ni][0] > 10 && gbias.hod[ni][0] < 16) {return int_gsl_integrate_low_precision(int_for_C_gl_HOD_tomo,(void*)array,amin_lens(ni),amax_lens(ni),NULL,1000);}

  // double res = int_gsl_integrate_medium_precision(int_for_C_gl_tomo,(void*)array,amin_lens(ni),amax_lens(ni),NULL,1000);
  double res = int_gsl_integrate_medium_precision(int_for_C_gl_tomo,(void*)array,amin_lens(ni),0.99999,NULL,1000);
  // printf("C_gl_tomo_nointerp(l=%lg,ni=%d,nj=%d):%lg\n",l,ni,nj,res);
  // exit(0);
  return res;
  // return int_gsl_integrate_medium_precision(int_for_C_gl_tomo,(void*)array,amin_lens(ni),amax_lens(ni),NULL,1000);
}

double C_shear_tomo_nointerp(double l, int ni, int nj) //shear tomography power spectra of source galaxy bins ni, nj
{
  double res;
  double array[3] = {(double)ni, (double)nj,l};
  // printf("ni,nj:%d,%d\n", ni,nj );
  int j,k;
  if (ni <= nj){j =nj; k = ni;}
  else{j = ni; k = nj;}
  res=int_gsl_integrate_medium_precision(int_for_C_shear_tomo,(void*)array,amin_source(j),amax_source(k),NULL,1000);
  // printf("res:%lg\n", res);
  return res;
  
}

/*********** angular power spectra - with look-up tables ******************/

double C_cl_HOD(double l, int ni)  //galaxy clustering power spectrum of galaxies in bin ni, using HOD model
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  if (ni < 0 || ni >= tomo.clustering_Nbin){
    printf("Invalid tomography option in C_gl_HOD(l,%d)\nEXIT\n",ni);
    exit(EXIT_FAILURE);
  }
  if (recompute_clustering(C,G,N,ni,ni)){
    if (table==0) {
      table   = create_double_matrix(0, tomo.clustering_Nbin, 0, Ntable.N_ell-1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell);
    }
    double llog,array[3];
    int i,k;

    for (k=0; k<tomo.clustering_Nbin; k++){
      array[0]=(double) k;
      llog = logsmin;
      for (i=0; i<Ntable.N_ell; i++, llog+=ds){
        array[2] = exp(llog);
        table[k][i] =log(int_gsl_integrate_low_precision(int_for_C_cl_HOD,(void*)array,amin_lens(k),amax_lens(k),NULL,1000));
      }
    }
    update_cosmopara(&C); update_nuisance(&N); update_galpara(&G);
  }
  
  double f1 = exp(interpol(table[ni], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1., 1.));
  if (isnan(f1)){f1 = 0.;}
  return f1;
}
double C_gl_tomo(double l, int ni, int nj)  //G-G lensing power spectrum, lens bin ni, source bin nj
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  double f1 = 0.;
  
  if (ni < 0 || ni >= tomo.clustering_Nbin ||nj < 0 || nj >= tomo.shear_Nbin){
    printf("C_gl_tomo(l,%d,%d) outside tomo.X_Nbin range\nEXIT\n",ni,nj); exit(1);
  }
  
  if (recompute_ggl(C,G,N,ni)){
    if (table==0){
      table   = create_double_matrix(0, tomo.ggl_Npowerspectra-1, 0, Ntable.N_ell-1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell-1.);
    }
    int i,k;
    double llog;
    
    for (k=0; k<tomo.ggl_Npowerspectra; k++) {
      llog = logsmin;
      for (i=0; i<Ntable.N_ell; i++, llog+=ds) {
          table[k][i]= log(C_gl_tomo_nointerp(exp(llog),ZL(k),ZS(k)));
      }
    }
    
    update_cosmopara(&C); update_nuisance(&N); update_galpara(&G);
    
  }
  if(test_zoverlap(ni,nj)){f1 = exp(interpol_fitslope(table[N_ggl(ni,nj)], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1.));}
  if (isnan(f1)){f1 = 0;}
  return f1;
}

double C_gl_HOD_tomo(double l, int ni, int nj)  //G-G lensing power spectrum, lens bin ni, source bin nj
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  double f1 = 0.;
  
  if (ni < 0 || ni >= tomo.clustering_Nbin ||nj < 0 || nj >= tomo.shear_Nbin){
    printf("C_gl_tomo(l,%d,%d) outside tomo.X_Nbin range\nEXIT\n",ni,nj); exit(1);
  }
  
  if (recompute_ggl(C,G,N,ni)){
    if (table==0){
      table   = create_double_matrix(0, tomo.ggl_Npowerspectra-1, 0, Ntable.N_ell-1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell);
    }
    int i,k;
    double llog,array[3];
    
    for (k=0; k<tomo.ggl_Npowerspectra; k++) {
      array[0]=(double) ZL(k); array[1]=(double) ZS(k);
      llog = logsmin;
      for (i=0; i<Ntable.N_ell; i++, llog+=ds) {
        table[k][i]= log(int_gsl_integrate_low_precision(int_for_C_gl_HOD_tomo,(void*)array,amin_source(ZS(k)),amax_lens(ZL(k)),NULL,1000));
      }
    }
    
    update_cosmopara(&C); update_nuisance(&N); update_galpara(&G);
    
  }
  if(test_zoverlap(ni,nj)){f1 = exp(interpol(table[N_ggl(ni,nj)], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1.,1.));}
  if (isnan(f1)){f1 = 0;}
  return f1;
}

double C_shear_tomo(double l, int ni, int nj)  //shear power spectrum of source galaxies in bins ni, nj
{
  // printf("hallo\n");
  static cosmopara C;
  static nuisancepara N;
  
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  if (ni < 0 || ni >= tomo.shear_Nbin ||nj < 0 || nj >= tomo.shear_Nbin){
    printf("C_shear_tomo(l,%d,%d) outside tomo.shear_Nbin range\nEXIT\n",ni,nj); exit(1);
  }
  
  if (recompute_shear(C,N)){
    if (table==0) {
      table   = create_double_matrix(0, tomo.shear_Npowerspectra-1, 0, Ntable.N_ell-1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell-1.);
    }
    
    double llog;
    int i,k;
    
    for (k=0; k<tomo.shear_Npowerspectra; k++) {
      llog = logsmin;
      for (i=0; i<Ntable.N_ell; i++, llog+=ds) {
        // printf("Hier!\n");
        // printf("C_shear_tomo, l, Z1(k),Z2(k): %lg,%lg,%d,%d:\n",log(C_shear_tomo_nointerp(exp(llog),Z1(k),Z2(k))),exp(llog),Z1(k),Z2(k));
        table[k][i]= log(C_shear_tomo_nointerp(exp(llog),Z1(k),Z2(k)));
       // printf("table[k][i]: %lg:\n",table[k][i]);
      }
    }
    update_cosmopara(&C); update_nuisance(&N); 
  }
  double f1 = exp(interpol_fitslope(table[N_shear(ni,nj)], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1.));
  if (isnan(f1)){f1 = 0.;}
  // printf("%le %d %d %le\n", l, ni, nj, f1);
  return f1;
}
