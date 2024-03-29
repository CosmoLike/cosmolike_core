/* This new file aims to rewrite the bin-averaged real-space fullsky 3x2pt covariances
// in covariances_real_binned_fullsky.c
// but also provide bin-averaged mixed/real(fullsky)/fourier space 6x2pt cov
// in a unified way
// Current status: binned-Fourier-space still needs testing;
// real-space: good
// mixed-space: not sure!
// -- Xiao Fang */

// fourier space 3x2pt Gaussian cov, without and with pure noise term
double func_for_cov_G_shear_noNN(double l, int *ar);
double func_for_cov_G_cl_noNN(double l, int *ar);
double func_for_cov_G_gl_noNN(double l, int *ar);
double func_for_cov_G_shear(double l, int *ar);
double func_for_cov_G_cl(double l, int *ar);
double func_for_cov_G_gl(double l, int *ar);
double func_for_cov_G_cl_shear(double l, int *ar);
double func_for_cov_G_cl_gl(double l, int *ar);
double func_for_cov_G_gl_shear(double l, int *ar);
// fourier space Gaussian cov with CMB lensing, pure noise included
double func_for_cov_G_gk(double l, int *ar);
double func_for_cov_G_ks(double l, int *ar);
double func_for_cov_G_gk_shear(double l, int *ar);
double func_for_cov_G_gk_gl(double l, int *ar);
double func_for_cov_G_gk_cl(double l, int *ar);
double func_for_cov_G_ks_shear(double l, int *ar);
double func_for_cov_G_ks_gl(double l, int *ar);
double func_for_cov_G_ks_cl(double l, int *ar);
double func_for_cov_G_ks_gk(double l, int *ar);
double func_for_cov_G_kk_shear(double l, int *ar);
double func_for_cov_G_kk_gl(double l, int *ar);
double func_for_cov_G_kk_cl(double l, int *ar);
double func_for_cov_G_kk_gk(double l, int *ar);
double func_for_cov_G_kk_ks(double l, int *ar);
//
double func_for_cov_G_kk(double l, int *ar);
// fourier space NG cov look-up tables
double bin_cov_NG_shear_shear(double l1,double l2, int *z_ar);
double bin_cov_NG_gl_gl(double l1,double l2, int *z_ar);
double bin_cov_NG_cl_cl(double l1,double l2, int *z_ar);
double bin_cov_NG_cl_shear(double l1,double l2, int *z_ar);
double bin_cov_NG_cl_gl(double l1,double l2, int *z_ar);
double bin_cov_NG_gl_shear(double l1,double l2, int *z_ar);
// bridge functions to unpack z_ar array to call the fourier cov routines
double cov_NG_shear_shear_bridge(double ll1, double ll2, int *z_ar);
double cov_NG_gl_gl_bridge(double ll1, double ll2, int *z_ar);
double cov_NG_cl_cl_bridge(double ll1, double ll2, int *z_ar);
double cov_NG_cl_gl_bridge(double ll1, double ll2, int *z_ar);
double cov_NG_cl_shear_bridge(double ll1, double ll2, int *z_ar);
double cov_NG_gl_shear_bridge(double ll1, double ll2, int *z_ar);
// tabulate fourier covNG
void tabulate_covNG(char *cov_type, int *z_ar, int Ntab, double **table);
// look-up tables for fullsky transform legendre polys
double Pl_tab(int itheta, int ell);
double Pl2_tab(int itheta, int ell);
double Glplus_tab(int itheta, int ell);
double Glminus_tab(int itheta, int ell);
double CMB_lensing_mat(int ibp, int L); // CMB lensing binning+correction mat
// CMB beam smoothing kernel
double GaussianBeam(double fwhm, int ell, double ell_left, double ell_right);
double GaussianBeam_test(double theta_fwhm, int ell_1, int ell_2, double ell_min, double ell_max, 
						int power_index_1, int power_index_2);
// 3x2pt Gaussian cov pure noise term on the diagonal, without masking effect
void pure_noise_xipm_xipm(int *z_ar, double *theta, double *dtheta, double *N);
void pure_noise_gl_gl(int *z_ar, double *theta, double *dtheta, double *N);
void pure_noise_cl_cl(int *z_ar, double *theta, double *dtheta, double *N);
// fullsky real-space cov, template routine
void cov_real_binned_fullsky(double **cov, double **covNG, char *realcov_type, int *z_ar, int FLAG_NG, double *theta, double *dtheta);
// fullsky real-space covs
void cov_shear_shear_real_binned_fullsky(double **cov, double **covNG, int z1,int z2,int z3,int z4,int pm1,int pm2, int FLAG_NG, double *theta, double *dtheta);
void cov_gl_shear_real_binned_fullsky(double **cov, double **covNG, int z1,int z2,int z3,int z4,int pm, int FLAG_NG, double *theta, double *dtheta);
void cov_cl_shear_real_binned_fullsky(double **cov, double **covNG, int z1,int z2,int z3,int z4,int pm, int FLAG_NG, double *theta, double *dtheta);
void cov_cl_gl_real_binned_fullsky(double **cov, double **covNG, int z1,int z2,int z3,int z4, int FLAG_NG, double *theta, double *dtheta);
void cov_gl_gl_real_binned_fullsky(double **cov, double **covNG, int z1,int z2,int z3,int z4, int FLAG_NG, double *theta, double *dtheta);
void cov_cl_cl_real_binned_fullsky(double **cov, double **covNG, int z1,int z2,int z3,int z4, int FLAG_NG, double *theta, double *dtheta);
//
void cov_kk_shear_real_binned_fullsky(double **cov, double **covNG, int z3,int z4,int pm, int FLAG_NG, double *theta, double *dtheta);
void cov_kk_gl_real_binned_fullsky(double **cov, double **covNG, int z3,int z4, int FLAG_NG, double *theta, double *dtheta);
void cov_kk_cl_real_binned_fullsky(double **cov, double **covNG, int z3,int z4, int FLAG_NG, double *theta, double *dtheta);
void cov_kk_gk_real_binned_fullsky(double **cov, double **covNG, int zl, int FLAG_NG, double *theta, double *dtheta);
void cov_kk_ks_real_binned_fullsky(double **cov, double **covNG, int zs, int FLAG_NG, double *theta, double *dtheta);
void cov_kk_kk_real_binned_fullsky(double **cov, double **covNG, int FLAG_NG, double *theta, double *dtheta);
// mixed cov 
void cov_kk_shear_mix_binned_fullsky(double **cov, double **covNG, int z3,int z4,int pm, int FLAG_NG, double *theta, double *dtheta, double *ell);
void cov_kk_gl_mix_binned_fullsky(double **cov, double **covNG, int z3,int z4, int FLAG_NG, double *theta, double *dtheta, double *ell);
void cov_kk_cl_mix_binned_fullsky(double **cov, double **covNG, int z3,int z4, int FLAG_NG, double *theta, double *dtheta, double *ell);
void cov_kk_gk_mix_binned_fullsky(double **cov, double **covNG, int zl, int FLAG_NG, double *theta, double *dtheta, double *ell);
void cov_kk_ks_mix_binned_fullsky(double **cov, double **covNG, int zs, int FLAG_NG, double *theta, double *dtheta, double *ell);
//
// Fourier-space cov, template routine
void cov_fourier_binned(double **cov, double **covNG, char *cov_type, int *z_ar, int FLAG_NG, double *ell);
// Fourier-space covs band averaged
void cov_shear_shear_fourier_binned(double **cov, double **covNG, int z1,int z2,int z3,int z4, int FLAG_NG, double *ell);
void cov_gl_shear_fourier_binned(double **cov, double **covNG, int z1,int z2,int z3,int z4, int FLAG_NG, double *ell);
void cov_cl_shear_fourier_binned(double **cov, double **covNG, int z1,int z2,int z3,int z4, int FLAG_NG, double *ell);
void cov_cl_gl_fourier_binned(double **cov, double **covNG, int z1,int z2,int z3,int z4, int FLAG_NG, double *ell);
void cov_gl_gl_fourier_binned(double **cov, double **covNG, int z1,int z2,int z3,int z4, int FLAG_NG, double *ell);
void cov_cl_cl_fourier_binned(double **cov, double **covNG, int z1,int z2,int z3,int z4, int FLAG_NG, double *ell);
//
void cov_kk_shear_fourier_binned(double **cov, double **covNG, int z3,int z4, int FLAG_NG, double *ell);
void cov_kk_gl_fourier_binned(double **cov, double **covNG, int z3,int z4, int FLAG_NG, double *ell);
void cov_kk_cl_fourier_binned(double **cov, double **covNG, int z3,int z4, int FLAG_NG, double *ell);
void cov_kk_gk_fourier_binned(double **cov, double **covNG, int zl, int FLAG_NG, double *ell);
void cov_kk_ks_fourier_binned(double **cov, double **covNG, int zs, int FLAG_NG, double *ell);
void cov_kk_kk_fourier_binned(double **cov, double **covNG, int FLAG_NG, double *ell);

//
// double cmb.fsky = 0.673;
/////////////////////////

// bridge functions to unpack z_ar array to call the fourier cov routines
double cov_NG_shear_shear_bridge(double ll1, double ll2, int *z_ar){
  return cov_NG_shear_shear_tomo(ll1,ll2,z_ar[0],z_ar[1],z_ar[2],z_ar[3]);
}
double cov_NG_gl_gl_bridge(double ll1, double ll2, int *z_ar){
  return cov_NG_gl_gl_tomo(ll1,ll2,z_ar[0],z_ar[1],z_ar[2],z_ar[3]);
}
double cov_NG_cl_cl_bridge(double ll1, double ll2, int *z_ar){
  return cov_NG_cl_cl_tomo(ll1,ll2,z_ar[0],z_ar[1],z_ar[2],z_ar[3]);
}
double cov_NG_cl_gl_bridge(double ll1, double ll2, int *z_ar){
  return cov_NG_cl_gl_tomo(ll1,ll2,z_ar[0],z_ar[1],z_ar[2],z_ar[3]);
}
double cov_NG_cl_shear_bridge(double ll1, double ll2, int *z_ar){
  return cov_NG_cl_shear_tomo(ll1,ll2,z_ar[0],z_ar[1],z_ar[2],z_ar[3]);
}
double cov_NG_gl_shear_bridge(double ll1, double ll2, int *z_ar){
  return cov_NG_gl_shear_tomo(ll1,ll2,z_ar[0],z_ar[1],z_ar[2],z_ar[3]);
}
double cov_NG_gk_shear_bridge(double ll1, double ll2, int *z_ar){
  return cov_NG_gk_ss(ll1,ll2,z_ar[0],z_ar[1],z_ar[2]);
}
double cov_NG_ks_shear_bridge(double ll1, double ll2, int *z_ar){
  return cov_NG_ks_ss(ll1,ll2,z_ar[0],z_ar[1],z_ar[2]);
}
double cov_NG_gk_gl_bridge(double ll1, double ll2, int *z_ar){
  return cov_NG_gk_gs(ll1,ll2,z_ar[0],z_ar[1],z_ar[2]);
}
double cov_NG_ks_gl_bridge(double ll1, double ll2, int *z_ar){
  return cov_NG_gs_ks(ll1,ll2,z_ar[1],z_ar[2],z_ar[0]); // Note: switched z_ar order to match previous fourier implementation
}
double cov_NG_gk_cl_bridge(double ll1, double ll2, int *z_ar){
  return cov_NG_gg_gk(ll1,ll2,z_ar[1],z_ar[2],z_ar[0]); // Note: switched z_ar order to match previous fourier implementation
}
double cov_NG_ks_cl_bridge(double ll1, double ll2, int *z_ar){
  return cov_NG_gg_ks(ll1,ll2,z_ar[1],z_ar[2],z_ar[0]); // Note: switched z_ar order to match previous fourier implementation
}
double cov_NG_gk_gk_bridge(double ll1, double ll2, int *z_ar){
  return cov_NG_gk_gk(ll1,ll2,z_ar[0],z_ar[1]);
}
double cov_NG_ks_gk_bridge(double ll1, double ll2, int *z_ar){
  return cov_NG_gk_ks(ll1,ll2,z_ar[1],z_ar[0]); // Note: switched z_ar order to match previous fourier implementation
}
double cov_NG_ks_ks_bridge(double ll1, double ll2, int *z_ar){
  return cov_NG_ks_ks(ll1,ll2,z_ar[0],z_ar[1]);
}
double cov_NG_kk_shear_bridge(double ll1, double ll2, int *z_ar){
  return cov_NG_kk_ss(ll1,ll2,z_ar[0],z_ar[1]);
}
double cov_NG_kk_gl_bridge(double ll1, double ll2, int *z_ar){
  return cov_NG_gs_kk(ll1,ll2,z_ar[0],z_ar[1]);
}
double cov_NG_kk_cl_bridge(double ll1, double ll2, int *z_ar){
  return cov_NG_gg_kk(ll1,ll2,z_ar[0],z_ar[1]);
}
double cov_NG_kk_gk_bridge(double ll1, double ll2, int *z_ar){
  return cov_NG_gk_kk(ll1,ll2,z_ar[0]);
}
double cov_NG_kk_ks_bridge(double ll1, double ll2, int *z_ar){
  return cov_NG_kk_ks(ll1,ll2,z_ar[0]);
}

// masking effect on pure noise
double w_mask(double theta_min){
  static int NTHETA = 0;
  static double *w_vec =0;
  int i,l;
  if (like.theta ==NULL || like.Ntheta < 1){
    printf("covariances_real_binned.c:w_mask: like.theta or like.Ntheta not initialized\nEXIT\n");
    exit(1);
  }
  if (w_vec ==0){
    w_vec = create_double_vector(0,like.Ntheta-1);
    NTHETA = like.Ntheta;
    FILE *F1;
    F1 = fopen(covparams.C_FOOTPRINT_FILE,"r");
    if (F1 != NULL) { //covparams.C_FOOTPRINT_FILE exists, use healpix C_mask(l) to compute mask correlation function
      fclose(F1);
      int lbins = line_count(covparams.C_FOOTPRINT_FILE);
      double *Cl;
      Cl = create_double_vector(0,lbins-1);
      F1=fopen(covparams.C_FOOTPRINT_FILE,"r");
      for (int i = 0; i < lbins; i++){
        int tmp;
        double tmp2;
        fscanf(F1,"%d %le\n",&tmp, &tmp2);
        Cl[i] = tmp2;
      }
      fclose(F1);

      printf("\nTabulating w_mask(theta) from mask power spectrum %s\n",covparams.C_FOOTPRINT_FILE);
      for (i = 0; i < NTHETA; i++){
        w_vec[i] =0.;
        for (l = 0; l < lbins; l++){
          w_vec[i]+=Cl[l]*(2.*l+1)/(4.*M_PI)*gsl_sf_legendre_Pl(l,cos(like.theta[i]));
        }
      }
      free_double_vector(Cl,0,lbins-1);
    }
    else{ //covparams.C_FOOTPRINT_FILE does not exit, ignore boundary effects
      printf("covparams.C_FOOTPRINT_FILE = %s not found\nNo boundary effect correction applied\n",covparams.C_FOOTPRINT_FILE);
      for (i = 0; i<NTHETA; i ++){
        w_vec[i] = 1.0;
        printf("w_mask[%d] = %e\n",i, w_vec[i]);
      }      
    }
  }
  i = 0;
  while(like.theta[i]< theta_min){
    i ++;
  }
  return w_vec[i];  
}

// tabulate fourier covNG
void tabulate_covNG(char *cov_type, int *z_ar, int Ntab, double **table){
  double logsmin = log(1.);
  double logsmax = log(5.e+4);

  double (*func_for_cov_NG)(double, double, int*);
  if(strcmp(cov_type, "shear_shear")==0) {
    func_for_cov_NG = &cov_NG_shear_shear_bridge;
  } else if(strcmp(cov_type, "gl_shear")==0) {
    func_for_cov_NG = &cov_NG_gl_shear_bridge;
  } else if(strcmp(cov_type, "cl_shear")==0) {
    func_for_cov_NG = &cov_NG_cl_shear_bridge;
  } else if(strcmp(cov_type, "gl_gl")==0) {
    func_for_cov_NG = &cov_NG_gl_gl_bridge;
  } else if(strcmp(cov_type, "cl_gl")==0) {
    func_for_cov_NG = &cov_NG_cl_gl_bridge;
  } else if(strcmp(cov_type, "cl_cl")==0) {
    func_for_cov_NG = &cov_NG_cl_cl_bridge;
  } else if(strcmp(cov_type, "gk_shear")==0) {
    func_for_cov_NG = &cov_NG_gk_shear_bridge;
  } else if(strcmp(cov_type, "ks_shear")==0) {
    func_for_cov_NG = &cov_NG_ks_shear_bridge;
  } else if(strcmp(cov_type, "gk_gl")==0) {
    func_for_cov_NG = &cov_NG_gk_gl_bridge;
  } else if(strcmp(cov_type, "ks_gl")==0) {
    func_for_cov_NG = &cov_NG_ks_gl_bridge;
  } else if(strcmp(cov_type, "gk_cl")==0) {
    func_for_cov_NG = &cov_NG_gk_cl_bridge;
  } else if(strcmp(cov_type, "ks_cl")==0) {
    func_for_cov_NG = &cov_NG_ks_cl_bridge;
  } else if(strcmp(cov_type, "gk_gk")==0) {
    func_for_cov_NG = &cov_NG_gk_gk_bridge;
  } else if(strcmp(cov_type, "ks_gk")==0) {
    func_for_cov_NG = &cov_NG_ks_gk_bridge;
  } else if(strcmp(cov_type, "ks_ks")==0) {
    func_for_cov_NG = &cov_NG_ks_ks_bridge;
  } else if(strcmp(cov_type, "kk_shear")==0) {
    func_for_cov_NG = &cov_NG_kk_shear_bridge;
  } else if(strcmp(cov_type, "kk_gl")==0) {
    func_for_cov_NG = &cov_NG_kk_gl_bridge;
  } else if(strcmp(cov_type, "kk_cl")==0) {
    func_for_cov_NG = &cov_NG_kk_cl_bridge;
  } else if(strcmp(cov_type, "kk_gk")==0) {
    func_for_cov_NG = &cov_NG_kk_gk_bridge;
  } else if(strcmp(cov_type, "kk_ks")==0) {
    func_for_cov_NG = &cov_NG_kk_ks_bridge;
  } else {
    printf("tabulate_covNG: cov_type \"%s\" not defined!\n", cov_type); exit(1);
  }

  int i,j;
  double llog1,llog2,ll1,ll2;
  double ds = (logsmax - logsmin)/(Ntab - 1.);
  llog1 = logsmin;
  for (i=0; i<Ntab; i++, llog1+=ds) {
    ll1 = exp(llog1);
    llog2 = logsmin;
    for (j=0; j<Ntab; j++, llog2+=ds) {
      ll2 = exp(llog2);
      table[i][j]=func_for_cov_NG(ll1,ll2,z_ar);
    }
  }
}

/**************** look-up tables for covariance  *********************/
double bin_cov_NG_shear_shear(double l1,double l2, int *z_ar){
  static int Z1 = -42, Z2 = -42, Z3 = -42, Z4 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  static char cov_type[] = "shear_shear";
  int z1,z2,z3,z4;
  z1=z_ar[0]; z2=z_ar[1]; z3=z_ar[2]; z4=z_ar[3];
  logsmin = log(1.); logsmax = log(5.e+4);
  ds = (logsmax - logsmin)/(Ntab - 1.);
  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 || Z4!=z4 ){
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    tabulate_covNG(cov_type, z_ar, Ntab, table);
    Z1=z1; Z2=z2; Z3=z3; Z4=z4;
  }
  res = 0.;
  if(l1<=1e-10 || l2<=1e-10) {return 0;}
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
    // res *= filter_cov_fourier(llog1, llog2, logsmax, log(4.e4));
  return res;
}
double bin_cov_NG_gl_gl(double l1,double l2, int *z_ar){
  static int Z1 = -42, Z2 = -42, Z3 = -42, Z4 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  static char cov_type[] = "gl_gl";
  int z1,z2,z3,z4;
  z1=z_ar[0]; z2=z_ar[1]; z3=z_ar[2]; z4=z_ar[3];
  logsmin = log(1.); logsmax = log(5.e+4);
  ds = (logsmax - logsmin)/(Ntab - 1.);

  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 || Z4!=z4 ){
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    tabulate_covNG(cov_type, z_ar, Ntab, table);
    Z1=z1; Z2=z2; Z3=z3; Z4=z4;
  }
  res = 0.;
  if(l1<=1e-10 || l2<=1e-10) {return 0;}
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
  return res;
}
double bin_cov_NG_cl_cl(double l1,double l2, int *z_ar){
  static int Z1 = -42, Z2 = -42, Z3 = -42, Z4 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  static char cov_type[] = "cl_cl";
  int z1,z2,z3,z4;
  z1=z_ar[0]; z2=z_ar[1]; z3=z_ar[2]; z4=z_ar[3];
  logsmin = log(1.); logsmax = log(5.e+4);
  ds = (logsmax - logsmin)/(Ntab - 1.);

  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 || Z4!=z4 ){
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    tabulate_covNG(cov_type, z_ar, Ntab, table);
    Z1=z1; Z2=z2; Z3=z3; Z4=z4;
  }
  res = 0.;
  if(l1<=1e-10 || l2<=1e-10) {return 0;}
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
  return res;
}
double bin_cov_NG_cl_shear(double l1,double l2, int *z_ar){
  static int Z1 = -42, Z2 = -42, Z3 = -42, Z4 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  static char cov_type[] = "cl_shear";
  int z1,z2,z3,z4;
  z1=z_ar[0]; z2=z_ar[1]; z3=z_ar[2]; z4=z_ar[3];
  logsmin = log(1.); logsmax = log(5.e+4);
  ds = (logsmax - logsmin)/(Ntab - 1.);

  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 || Z4!=z4 ){
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    tabulate_covNG(cov_type, z_ar, Ntab, table);
    Z1=z1; Z2=z2; Z3=z3; Z4=z4;
  }
  res = 0.;
  if(l1<=1e-10 || l2<=1e-10) {return 0;}
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
  return res;
}
double bin_cov_NG_cl_gl(double l1,double l2, int *z_ar){
  static int Z1 = -42, Z2 = -42, Z3 = -42, Z4 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  static char cov_type[] = "cl_gl";
  int z1,z2,z3,z4;
  z1=z_ar[0]; z2=z_ar[1]; z3=z_ar[2]; z4=z_ar[3];
  logsmin = log(1.); logsmax = log(5.e+4);
  ds = (logsmax - logsmin)/(Ntab - 1.);

  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 || Z4!=z4 )
    {
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    tabulate_covNG(cov_type, z_ar, Ntab, table);
    Z1=z1; Z2=z2; Z3=z3; Z4=z4;
    }
  res = 0.;
  if(l1<=1e-10 || l2<=1e-10) {return 0;}
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
  return res;
}
double bin_cov_NG_gl_shear(double l1,double l2, int *z_ar){
  static int Z1 = -42, Z2 = -42, Z3 = -42, Z4 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  static char cov_type[] = "gl_shear";
  int z1,z2,z3,z4;
  z1=z_ar[0]; z2=z_ar[1]; z3=z_ar[2]; z4=z_ar[3];
  logsmin = log(1.); logsmax = log(5.e+4);
  ds = (logsmax - logsmin)/(Ntab - 1.);

  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 || Z4!=z4 ){
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    tabulate_covNG(cov_type, z_ar, Ntab, table);
    Z1=z1; Z2=z2; Z3=z3; Z4=z4;
  }
  res = 0.;
  if(l1<=1e-10 || l2<=1e-10) {return 0;}
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
  return res;
}

double bin_cov_NG_gk_shear(double l1,double l2, int *z_ar){
  static int Z1 = -42, Z2 = -42, Z3 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  static char cov_type[] = "gk_shear";
  int z1,z2,z3;
  z1=z_ar[0]; z2=z_ar[1]; z3=z_ar[2];
  logsmin = log(1.); logsmax = log(5.e+4);
  ds = (logsmax - logsmin)/(Ntab - 1.);

  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 ){
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    tabulate_covNG(cov_type, z_ar, Ntab, table);
    Z1=z1; Z2=z2; Z3=z3;
  }
  res = 0.;
  if(l1<=1e-10 || l2<=1e-10) {return 0;}
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
  return res;
}

double bin_cov_NG_ks_shear(double l1,double l2, int *z_ar){
  static int Z1 = -42, Z2 = -42, Z3 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  static char cov_type[] = "ks_shear";
  int z1,z2,z3;
  z1=z_ar[0]; z2=z_ar[1]; z3=z_ar[2];
  logsmin = log(1.); logsmax = log(5.e+4);
  ds = (logsmax - logsmin)/(Ntab - 1.);

  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 ){
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    tabulate_covNG(cov_type, z_ar, Ntab, table);
    Z1=z1; Z2=z2; Z3=z3;
  }
  res = 0.;
  if(l1<=1e-10 || l2<=1e-10) {return 0;}
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
  return res;
}

double bin_cov_NG_gk_gl(double l1,double l2, int *z_ar){
  static int Z1 = -42, Z2 = -42, Z3 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  static char cov_type[] = "gk_gl";
  int z1,z2,z3;
  z1=z_ar[0]; z2=z_ar[1]; z3=z_ar[2];
  logsmin = log(1.); logsmax = log(5.e+4);
  ds = (logsmax - logsmin)/(Ntab - 1.);

  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 ){
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    tabulate_covNG(cov_type, z_ar, Ntab, table);
    Z1=z1; Z2=z2; Z3=z3;
  }
  res = 0.;
  if(l1<=1e-10 || l2<=1e-10) {return 0;}
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
  return res;
}

double bin_cov_NG_ks_gl(double l1,double l2, int *z_ar){
  static int Z1 = -42, Z2 = -42, Z3 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  static char cov_type[] = "ks_gl";
  int z1,z2,z3;
  z1=z_ar[0]; z2=z_ar[1]; z3=z_ar[2];
  logsmin = log(1.); logsmax = log(5.e+4);
  ds = (logsmax - logsmin)/(Ntab - 1.);

  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 ){
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    tabulate_covNG(cov_type, z_ar, Ntab, table);
    Z1=z1; Z2=z2; Z3=z3;
  }
  res = 0.;
  if(l1<=1e-10 || l2<=1e-10) {return 0;}
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
  return res;
}

double bin_cov_NG_gk_cl(double l1,double l2, int *z_ar){
  static int Z1 = -42, Z2 = -42, Z3 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  static char cov_type[] = "gk_cl";
  int z1,z2,z3;
  z1=z_ar[0]; z2=z_ar[1]; z3=z_ar[2];
  logsmin = log(1.); logsmax = log(5.e+4);
  ds = (logsmax - logsmin)/(Ntab - 1.);

  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 ){
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    tabulate_covNG(cov_type, z_ar, Ntab, table);
    Z1=z1; Z2=z2; Z3=z3;
  }
  res = 0.;
  if(l1<=1e-10 || l2<=1e-10) {return 0;}
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
  return res;
}

double bin_cov_NG_ks_cl(double l1,double l2, int *z_ar){
  static int Z1 = -42, Z2 = -42, Z3 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  static char cov_type[] = "ks_cl";
  int z1,z2,z3;
  z1=z_ar[0]; z2=z_ar[1]; z3=z_ar[2];
  logsmin = log(1.); logsmax = log(5.e+4);
  ds = (logsmax - logsmin)/(Ntab - 1.);

  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 ){
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    tabulate_covNG(cov_type, z_ar, Ntab, table);
    Z1=z1; Z2=z2; Z3=z3;
  }
  res = 0.;
  if(l1<=1e-10 || l2<=1e-10) {return 0;}
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
  return res;
}

double bin_cov_NG_gk_gk(double l1,double l2, int *z_ar){
  static int Z1 = -42, Z2 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  static char cov_type[] = "gk_gk";
  int z1,z2;
  z1=z_ar[0]; z2=z_ar[1];
  logsmin = log(1.); logsmax = log(5.e+4);
  ds = (logsmax - logsmin)/(Ntab - 1.);

  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2){
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    tabulate_covNG(cov_type, z_ar, Ntab, table);
    Z1=z1; Z2=z2;
  }
  res = 0.;
  if(l1<=1e-10 || l2<=1e-10) {return 0;}
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
  return res;
}

double bin_cov_NG_ks_gk(double l1,double l2, int *z_ar){
  static int Z1 = -42, Z2 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  static char cov_type[] = "ks_gk";
  int z1,z2;
  z1=z_ar[0]; z2=z_ar[1];
  logsmin = log(1.); logsmax = log(5.e+4);
  ds = (logsmax - logsmin)/(Ntab - 1.);

  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2){
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    tabulate_covNG(cov_type, z_ar, Ntab, table);
    Z1=z1; Z2=z2;
  }
  res = 0.;
  if(l1<=1e-10 || l2<=1e-10) {return 0;}
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
  return res;
}

double bin_cov_NG_ks_ks(double l1,double l2, int *z_ar){
  static int Z1 = -42, Z2 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  static char cov_type[] = "ks_ks";
  int z1,z2;
  z1=z_ar[0]; z2=z_ar[1];
  logsmin = log(1.); logsmax = log(5.e+4);
  ds = (logsmax - logsmin)/(Ntab - 1.);

  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2){
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    tabulate_covNG(cov_type, z_ar, Ntab, table);
    Z1=z1; Z2=z2;
  }
  res = 0.;
  if(l1<=1e-10 || l2<=1e-10) {return 0;}
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
  return res;
}

double bin_cov_NG_kk_shear(double l1,double l2, int *z_ar){
  static int Z1 = -42, Z2 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  static char cov_type[] = "kk_shear";
  int z1,z2;
  z1=z_ar[0]; z2=z_ar[1];
  logsmin = log(1.); logsmax = log(5.e+4);
  ds = (logsmax - logsmin)/(Ntab - 1.);

  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 ){
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    tabulate_covNG(cov_type, z_ar, Ntab, table);
    Z1=z1; Z2=z2;
  }
  res = 0.;
  if(l1<=1e-10 || l2<=1e-10) {return 0;}
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
  return res;
}

double bin_cov_NG_kk_gl(double l1,double l2, int *z_ar){
  static int Z1 = -42, Z2 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  static char cov_type[] = "kk_gl";
  int z1,z2;
  z1=z_ar[0]; z2=z_ar[1];
  logsmin = log(1.); logsmax = log(5.e+4);
  ds = (logsmax - logsmin)/(Ntab - 1.);

  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 ){
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    tabulate_covNG(cov_type, z_ar, Ntab, table);
    Z1=z1; Z2=z2;
  }
  res = 0.;
  if(l1<=1e-10 || l2<=1e-10) {return 0;}
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
  return res;
}

double bin_cov_NG_kk_cl(double l1,double l2, int *z_ar){
  static int Z1 = -42, Z2 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  static char cov_type[] = "kk_cl";
  int z1,z2;
  z1=z_ar[0]; z2=z_ar[1];
  logsmin = log(1.); logsmax = log(5.e+4);
  ds = (logsmax - logsmin)/(Ntab - 1.);

  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 ){
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    tabulate_covNG(cov_type, z_ar, Ntab, table);
    Z1=z1; Z2=z2;
  }
  res = 0.;
  if(l1<=1e-10 || l2<=1e-10) {return 0;}
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
  return res;
}

double bin_cov_NG_kk_gk(double l1,double l2, int *z_ar){
  static int Z1 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  static char cov_type[] = "kk_gk";
  int z1;
  z1=z_ar[0];
  logsmin = log(1.); logsmax = log(5.e+4);
  ds = (logsmax - logsmin)/(Ntab - 1.);

  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 ){
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    tabulate_covNG(cov_type, z_ar, Ntab, table);
    Z1=z1;
  }
  res = 0.;
  if(l1<=1e-10 || l2<=1e-10) {return 0;}
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
  return res;
}

double bin_cov_NG_kk_ks(double l1,double l2, int *z_ar){
  static int Z1 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  static char cov_type[] = "kk_ks";
  int z1;
  z1=z_ar[0];
  logsmin = log(1.); logsmax = log(5.e+4);
  ds = (logsmax - logsmin)/(Ntab - 1.);

  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1){
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    tabulate_covNG(cov_type, z_ar, Ntab, table);
    Z1=z1;
  }
  res = 0.;
  if(l1<=1e-10 || l2<=1e-10) {return 0;}
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
  return res;
}

double bin_cov_NG_kk_kk(double l1,double l2, int *z_ar){
  static int Z1 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  static char cov_type[] = "kk_kk";
  int z1;
  z1=z_ar[0];
  logsmin = log(1.); logsmax = log(5.e+4);
  ds = (logsmax - logsmin)/(Ntab - 1.);

  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1){
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    tabulate_covNG(cov_type, z_ar, Ntab, table);
    Z1=z1;
  }
  res = 0.;
  if(l1<=1e-10 || l2<=1e-10) {return 0;}
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
  return res;
}

double Glplus_tab(int itheta, int ell) {
  
  int i;
  static int LMAX = 50000;
  static double **Glplus =0;
  if (Glplus ==0){
    Glplus =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);

    double *xmin, *xmax, *Pmin, *Pmax, *dPmin, *dPmax;
    xmin= create_double_vector(0, like.Ntheta-1);
    xmax= create_double_vector(0, like.Ntheta-1);
    double logdt=(log(like.vtmax)-log(like.vtmin))/like.Ntheta;
    // printf("like.vtmax,like.vtmin,like.Ntheta,%lg,%lg,%lg\n",like.vtmax,like.vtmin,like.Ntheta);
    for(i=0; i<like.Ntheta ; i++){
      xmin[i]=cos(exp(log(like.vtmin)+(i+0.0)*logdt));
      xmax[i]=cos(exp(log(like.vtmin)+(i+1.0)*logdt));
    }
    Pmin= create_double_vector(0, LMAX+1);
    Pmax= create_double_vector(0, LMAX+1);
    dPmin= create_double_vector(0, LMAX+1);
    dPmax= create_double_vector(0, LMAX+1);
    for (i = 0; i<like.Ntheta; i ++){
      gsl_sf_legendre_Pl_deriv_array(LMAX, xmin[i],Pmin,dPmin);
      gsl_sf_legendre_Pl_deriv_array(LMAX, xmax[i],Pmax,dPmax);
      // printf("xmax[%d]:%lg\n", i,xmax[i]);
      for (int l = 2; l < LMAX; l ++){
        Glplus[i][l] =(2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(

        -l*(l-1.)/2*(l+2./(2*l+1)) * (Pmin[l-1]-Pmax[l-1])
        -l*(l-1.)*(2.-l)/2         * (xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
        +l*(l-1.)/(2.*l+1)           * (Pmin[l+1]-Pmax[l+1])

        +(4-l)   * (dPmin[l]-dPmax[l])
        +(l+2)   * (xmin[i]*dPmin[l-1] - xmax[i]*dPmax[l-1] - Pmin[l-1] + Pmax[l-1])

        +2*(l-1) * (xmin[i]*dPmin[l]   - xmax[i]*dPmax[l]   - Pmin[l] + Pmax[l])
        -2*(l+2) * (dPmin[l-1]-dPmax[l-1])

        )/(xmin[i]-xmax[i]);           
      }
      Glplus[i][0] = 0.; Glplus[i][1] = 0.;
    }
    free_double_vector(xmin,0,like.Ntheta-1);
    free_double_vector(xmax,0,like.Ntheta-1);
    free_double_vector(Pmin,0,LMAX+1);
    free_double_vector(Pmax,0,LMAX+1);
    free_double_vector(dPmin,0,LMAX+1);
    free_double_vector(dPmax,0,LMAX+1);
  }
  return Glplus[itheta][ell];
}

double Glminus_tab(int itheta, int ell) {
  
  int i;
  static int LMAX = 50000;
  static double **Glminus =0;
  if (Glminus ==0){
    Glminus =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);

    double *xmin, *xmax, *Pmin, *Pmax, *dPmin, *dPmax;
    xmin= create_double_vector(0, like.Ntheta-1);
    xmax= create_double_vector(0, like.Ntheta-1);
    double logdt=(log(like.vtmax)-log(like.vtmin))/like.Ntheta;
    // printf("like.vtmax,like.vtmin,like.Ntheta,%lg,%lg,%lg\n",like.vtmax,like.vtmin,like.Ntheta);
    for(i=0; i<like.Ntheta ; i++){
      xmin[i]=cos(exp(log(like.vtmin)+(i+0.0)*logdt));
      xmax[i]=cos(exp(log(like.vtmin)+(i+1.0)*logdt));
    }
    Pmin= create_double_vector(0, LMAX+1);
    Pmax= create_double_vector(0, LMAX+1);
    dPmin= create_double_vector(0, LMAX+1);
    dPmax= create_double_vector(0, LMAX+1);
    for (i = 0; i<like.Ntheta; i ++){
      gsl_sf_legendre_Pl_deriv_array(LMAX, xmin[i],Pmin,dPmin);
      gsl_sf_legendre_Pl_deriv_array(LMAX, xmax[i],Pmax,dPmax);
      // printf("xmax[%d]:%lg\n", i,xmax[i]);
      for (int l = 2; l < LMAX; l ++){
        Glminus[i][l] =(2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(

        -l*(l-1.)/2*(l+2./(2*l+1)) * (Pmin[l-1]-Pmax[l-1])
        -l*(l-1.)*(2.-l)/2         * (xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
        +l*(l-1.)/(2.*l+1)           * (Pmin[l+1]-Pmax[l+1])

        +(4-l)   * (dPmin[l]-dPmax[l])
        +(l+2)   * (xmin[i]*dPmin[l-1] - xmax[i]*dPmax[l-1] - Pmin[l-1] + Pmax[l-1])

        -2*(l-1) * (xmin[i]*dPmin[l]   - xmax[i]*dPmax[l]   - Pmin[l] + Pmax[l])
        +2*(l+2) * (dPmin[l-1]-dPmax[l-1])

        )/(xmin[i]-xmax[i]);        
      }
      Glminus[i][0] = 0.; Glminus[i][1] = 0.;
    }
    free_double_vector(xmin,0,like.Ntheta-1);
    free_double_vector(xmax,0,like.Ntheta-1);
    free_double_vector(Pmin,0,LMAX+1);
    free_double_vector(Pmax,0,LMAX+1);
    free_double_vector(dPmin,0,LMAX+1);
    free_double_vector(dPmax,0,LMAX+1);
  }
  return Glminus[itheta][ell];
}

double Pl_tab(int itheta, int ell) {
  
  int i;
  static int LMAX = 50000;
  static double **Pl =0;
  if (Pl ==0){
    Pl =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);

    double *xmin, *xmax, *Pmin, *Pmax, *dPmin, *dPmax;
    xmin= create_double_vector(0, like.Ntheta-1);
    xmax= create_double_vector(0, like.Ntheta-1);
    double logdt=(log(like.vtmax)-log(like.vtmin))/like.Ntheta;
    // printf("like.vtmax,like.vtmin,like.Ntheta,%lg,%lg,%lg\n",like.vtmax,like.vtmin,like.Ntheta);
    for(i=0; i<like.Ntheta ; i++){
      xmin[i]=cos(exp(log(like.vtmin)+(i+0.0)*logdt));
      xmax[i]=cos(exp(log(like.vtmin)+(i+1.0)*logdt));
    }
    Pmin= create_double_vector(0, LMAX+1);
    Pmax= create_double_vector(0, LMAX+1);
    dPmin= create_double_vector(0, LMAX+1);
    dPmax= create_double_vector(0, LMAX+1);
    for (i = 0; i<like.Ntheta; i ++){
      gsl_sf_legendre_Pl_deriv_array(LMAX, xmin[i],Pmin,dPmin);
      gsl_sf_legendre_Pl_deriv_array(LMAX, xmax[i],Pmax,dPmax);
      // printf("xmax[%d]:%lg\n", i,xmax[i]);
      for (int l = 1; l < LMAX; l ++){
        Pl[i][l] = 1./(4.*M_PI)*(Pmin[l+1]-Pmax[l+1]-Pmin[l-1]+Pmax[l-1])/(xmin[i]-xmax[i]);     
      }
      Pl[i][0] = 1./(4.*M_PI);
    }
    free_double_vector(xmin,0,like.Ntheta-1);
    free_double_vector(xmax,0,like.Ntheta-1);
    free_double_vector(Pmin,0,LMAX+1);
    free_double_vector(Pmax,0,LMAX+1);
    free_double_vector(dPmin,0,LMAX+1);
    free_double_vector(dPmax,0,LMAX+1);
  }
  return Pl[itheta][ell];
}

double Pl2_tab(int itheta, int ell) {
  
  int i;
  static int LMAX = 50000;
  static double **Pl2 =0;
  if (Pl2 ==0){
    Pl2 =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);

    double *xmin, *xmax, *Pmin, *Pmax, *dPmin, *dPmax;
    xmin= create_double_vector(0, like.Ntheta-1);
    xmax= create_double_vector(0, like.Ntheta-1);
    double logdt=(log(like.vtmax)-log(like.vtmin))/like.Ntheta;
    // printf("like.vtmax,like.vtmin,like.Ntheta,%lg,%lg,%lg\n",like.vtmax,like.vtmin,like.Ntheta);
    for(i=0; i<like.Ntheta ; i++){
      xmin[i]=cos(exp(log(like.vtmin)+(i+0.0)*logdt));
      xmax[i]=cos(exp(log(like.vtmin)+(i+1.0)*logdt));
    }
    Pmin= create_double_vector(0, LMAX+1);
    Pmax= create_double_vector(0, LMAX+1);
    dPmin= create_double_vector(0, LMAX+1);
    dPmax= create_double_vector(0, LMAX+1);
    for (i = 0; i<like.Ntheta; i ++){
      gsl_sf_legendre_Pl_deriv_array(LMAX, xmin[i],Pmin,dPmin);
      gsl_sf_legendre_Pl_deriv_array(LMAX, xmax[i],Pmax,dPmax);
      // printf("xmax[%d]:%lg\n", i,xmax[i]);
      for (int l = 1; l < LMAX; l ++){
        Pl2[i][l] = (2.*l+1)/(4.*M_PI*l*(l+1)*(xmin[i]-xmax[i]))
        *((l+2./(2*l+1.))*(Pmin[l-1]-Pmax[l-1])
        +(2-l)*(xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
        -2./(2*l+1.)*(Pmin[l+1]-Pmax[l+1]));
      }
      Pl2[i][0] = 0.;
    }
    free_double_vector(xmin,0,like.Ntheta-1);
    free_double_vector(xmax,0,like.Ntheta-1);
    free_double_vector(Pmin,0,LMAX+1);
    free_double_vector(Pmax,0,LMAX+1);
    free_double_vector(dPmin,0,LMAX+1);
    free_double_vector(dPmax,0,LMAX+1);
  }
  return Pl2[itheta][ell];
}

double CMB_lensing_mat(int ibp, int L){
  int iL;
  static double **binning_correction_matrix = NULL;
  if (binning_correction_matrix == NULL){
    binning_correction_matrix = create_double_matrix(0, like.Nbp-1, 
      0, like.Ncl-1);
    FILE *F1;
    F1 = fopen(covparams.BINMAT_FILE,"r");
    if (F1 != NULL) {
      fclose(F1);
      int nbp = line_count(covparams.BINMAT_FILE);
      
      F1=fopen(covparams.BINMAT_FILE, "r");
      for (int i = 0; i < like.Nbp; i++){
        for (int j = 0; j < like.Ncl; j++){
          double tmp2;
          fscanf(F1,"%le ",&tmp2);
          binning_correction_matrix[i][j] = tmp2;
        }
      }
      fclose(F1);
    }
    else{
      printf("ERROR: Cant't open file %s\n", covparams.BINMAT_FILE);
      exit(-1);
    }  
  }
  iL = L - like.lmin;
  return binning_correction_matrix[ibp][iL];
}

/*  Calculate Gaussian Beam Kernel
    
    F(ell) = B(ell)H(ell - ell_min)H(ell_max - ell)
    where B(ell) is the Gaussian kernel
    B(ell) = exp( -ell*(ell+1)/ell_beam^2 )
    and ell_beam = sqrt(16*ln(2)) / theta_FWHM
    H(ell) is the Heaviside step function

    Inputs:
    -------
    theta_fwhm: FWHM of the CMB Gaussian beam, in units of rad
    ell:        angular mode ell where the beam kernel is evaluated
    ell_min:    large-scale cut
    ell_max:    small-scale cut

    Outputs:
    --------
    BeamKernel: The Gaussian kernel
*/
double GaussianBeam(double theta_fwhm, int ell, double ell_left, double ell_right){
  double ell_min = ell_left - 10;
  double ell_max = ell_right + 50;
  double ell_beam = sqrt(16.0 * log(2.0)) / theta_fwhm;
  double BeamKernel, window_func, u=0;
  // Turn-On Gaussian Beam Smoothing
  if(theta_fwhm > 0.0){
	// smooth kernel
    BeamKernel = exp(-1.*ell*(ell+1.) / (ell_beam * ell_beam));
    // window function
    if((ell<ell_min) || (ell>ell_max)){
      window_func = 0.;
    }
    // tapering the edges
    else if((ell_min <= ell) && (ell<ell_left)){
      u = (ell - ell_min)/(ell_left - ell_min);
      window_func = u - sin( constants.twopi * u) / constants.twopi;
    }
    else if ((ell_right < ell) && (ell <= ell_max)){
      u = (ell_max - ell)/(ell_max - ell_right);
      window_func = u - sin( constants.twopi * u) / constants.twopi; 
    }
    else{
      window_func = 1.0;
    }
    BeamKernel *= window_func;
  }
  // Turn-Off Gaussian Beam Smoothing
  else{BeamKernel = 1.0;}
  return BeamKernel;
}

double GaussianBeam_test(double theta_fwhm, int ell_1, int ell_2, double ell_min, double ell_max, 
						int power_index_1, int power_index_2){
    double ell_beam = sqrt(16.0 * log(2.0)) / theta_fwhm;
    double BeamKernel;
    
    if( ((power_index_1==0) && (power_index_2==0)) || (theta_fwhm<=0.0) ){BeamKernel = 1.0;}
    else{
		double ell2 = (double)power_index_1*ell_1*(ell_1+1.0) + power_index_2*ell_2*(ell_2+1.0);
        BeamKernel = exp( -1.*ell2 / (ell_beam * ell_beam) );
        if(ell_1<ell_min || ell_1>ell_max || ell_2<ell_min || ell_2>ell_max){BeamKernel = 0.;}
    }
    return BeamKernel;
}

// Pure noise term in Gaussian cov
void pure_noise_xipm_xipm(int *z_ar, double *theta, double *dtheta, double *N) {
  int z1,z2,z3,z4;
  int i;
  z1=z_ar[0]; z2=z_ar[1]; z3=z_ar[2]; z4=z_ar[3];
  if (z1 ==z3 && z2 ==z4){
    for(i=0;i<like.Ntheta;i++) {
      N[i] = 2.* pow(survey.sigma_e,4.0)/(M_PI*(2.*theta[i]+dtheta[i])*dtheta[i]*4.*nsource(z1)*nsource(z2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor);
    } //(number of galaxy pairs in the survey contributing to annulus of width Dtheta centered at theta1)^-1
  }
  if (z1 ==z4 && z2 ==z3){
    for(i=0;i<like.Ntheta;i++) {
      N[i] += 2.* pow(survey.sigma_e,4.0)/(M_PI*(2.*theta[i]+dtheta[i])*dtheta[i]*4.*nsource(z1)*nsource(z2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor);
    }
  } // factor of 2 already included
}
void pure_noise_gl_gl(int *z_ar, double *theta, double *dtheta, double *N) {
  int zl1, zs1, zl2, zs2;
  int i;
  zl1 = z_ar[0]; zs1 = z_ar[1]; zl2 = z_ar[2]; zs2 = z_ar[3];
  if (zl1 ==zl2 && zs1 ==zs2){
    for(i=0;i<like.Ntheta;i++) {
      N[i] = pow(survey.sigma_e,2.0)/(2.0*M_PI*(2.*theta[i]+dtheta[i])*dtheta[i]*nlens(zl1)*nsource(zs2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor);
    }
  }
}
void pure_noise_cl_cl(int *z_ar, double *theta, double *dtheta, double *N) {
  int z1,z2,z3,z4;
  int i;
  z1=z_ar[0]; z2=z_ar[1]; z3=z_ar[2]; z4=z_ar[3];
  if (z1 ==z3 && z2 ==z4){
    for(i=0;i<like.Ntheta;i++) {
      N[i] = 1./(M_PI*(2.*theta[i]+dtheta[i])*dtheta[i]*nlens(z1)*nlens(z2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor); //(number of galaxy pairs in the survey contributing to annulus of width Dtheta centered at theta1)^-1
    }
  }
  if (z1 ==z4 && z2 ==z3){
    for(i=0;i<like.Ntheta;i++) {
      N[i] += 1./(M_PI*(2.*theta[i]+dtheta[i])*dtheta[i]*nlens(z1)*nlens(z2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor);
    }
  }
}

/////// full sky covs
// Template routine
void cov_real_binned_fullsky(double **cov, double **covNG, char *realcov_type, int *z_ar, int FLAG_NG, double *theta, double *dtheta){

  int i,j;
  static int LMAX = 50000;
  int CMB_smooth_1 = 0, CMB_smooth_2 = 0; // flag for CMB smoothing kernel

  double N[like.Ntheta];
  for(i=0;i<like.Ntheta;i++) {N[i] = 0.;}

  double (*func_for_cov_G)(double, int*);
  double (*func_bin_cov_NG)(double, double, int*);
  double (*func_P1)(int, int);
  double (*func_P2)(int, int);

  if(strcmp(realcov_type, "xi+_xi+")==0) {
    func_for_cov_G  = &func_for_cov_G_shear_noNN;
    func_bin_cov_NG = &bin_cov_NG_shear_shear;
    func_P1 = &Glplus_tab;
    func_P2 = &Glplus_tab;
    pure_noise_xipm_xipm(z_ar, theta, dtheta, N); //C+- doesn't have the diagonal shot noise term
  } else if(strcmp(realcov_type, "xi-_xi-")==0) {
    func_for_cov_G  = &func_for_cov_G_shear_noNN;
    func_bin_cov_NG = &bin_cov_NG_shear_shear;
    func_P1 = &Glminus_tab;
    func_P2 = &Glminus_tab;
    pure_noise_xipm_xipm(z_ar, theta, dtheta, N); //C+- doesn't have the diagonal shot noise term
  } else if(strcmp(realcov_type, "xi+_xi-")==0) {
    func_for_cov_G  = &func_for_cov_G_shear_noNN;
    func_bin_cov_NG = &bin_cov_NG_shear_shear;
    func_P1 = &Glplus_tab;
    func_P2 = &Glminus_tab;
  } else if(strcmp(realcov_type, "xi-_xi+")==0) {
    func_for_cov_G  = &func_for_cov_G_shear_noNN;
    func_bin_cov_NG = &bin_cov_NG_shear_shear;
    func_P1 = &Glminus_tab;
    func_P2 = &Glplus_tab;
  } else if(strcmp(realcov_type, "gl_xi+")==0) {
    func_for_cov_G  = &func_for_cov_G_gl_shear;
    func_bin_cov_NG = &bin_cov_NG_gl_shear;
    func_P1 = &Pl2_tab;
    func_P2 = &Glplus_tab;
  } else if(strcmp(realcov_type, "gl_xi-")==0) {
    func_for_cov_G  = &func_for_cov_G_gl_shear;
    func_bin_cov_NG = &bin_cov_NG_gl_shear;
    func_P1 = &Pl2_tab;
    func_P2 = &Glminus_tab;
  } else if(strcmp(realcov_type, "cl_xi+")==0) {
    func_for_cov_G  = &func_for_cov_G_cl_shear;
    func_bin_cov_NG = &bin_cov_NG_cl_shear;
    func_P1 = &Pl_tab;
    func_P2 = &Glplus_tab;
  } else if(strcmp(realcov_type, "cl_xi-")==0) {
    func_for_cov_G  = &func_for_cov_G_cl_shear;
    func_bin_cov_NG = &bin_cov_NG_cl_shear;
    func_P1 = &Pl_tab;
    func_P2 = &Glminus_tab;
  } else if(strcmp(realcov_type, "gl_gl")==0) {
    func_for_cov_G  = &func_for_cov_G_gl_noNN;
    func_bin_cov_NG = &bin_cov_NG_gl_gl;
    func_P1 = &Pl2_tab;
    func_P2 = &Pl2_tab;
    pure_noise_gl_gl(z_ar, theta, dtheta, N);
  } else if(strcmp(realcov_type, "cl_gl")==0) {
    func_for_cov_G  = &func_for_cov_G_cl_gl;
    func_bin_cov_NG = &bin_cov_NG_cl_gl;
    func_P1 = &Pl_tab;
    func_P2 = &Pl2_tab;
  } else if(strcmp(realcov_type, "cl_cl")==0) {
    func_for_cov_G  = &func_for_cov_G_cl_noNN;
    func_bin_cov_NG = &bin_cov_NG_cl_cl;
    func_P1 = &Pl_tab;
    func_P2 = &Pl_tab;
    pure_noise_cl_cl(z_ar, theta, dtheta, N);
  } else if(strcmp(realcov_type, "gk_xi+")==0) {
    func_for_cov_G  = &func_for_cov_G_gk_shear;
    func_bin_cov_NG = &bin_cov_NG_gk_shear;
    func_P1 = &Pl_tab;
    func_P2 = &Glplus_tab;
    CMB_smooth_1 = 1;
    CMB_smooth_2 = 0;
  } else if(strcmp(realcov_type, "gk_xi-")==0) {
    func_for_cov_G  = &func_for_cov_G_gk_shear;
    func_bin_cov_NG = &bin_cov_NG_gk_shear;
    func_P1 = &Pl_tab;
    func_P2 = &Glminus_tab;
    CMB_smooth_1 = 1;
    CMB_smooth_2 = 0;
  } else if(strcmp(realcov_type, "ks_xi+")==0) {
    func_for_cov_G  = &func_for_cov_G_ks_shear;
    func_bin_cov_NG = &bin_cov_NG_ks_shear;
    func_P1 = &Pl2_tab;
    func_P2 = &Glplus_tab;
    CMB_smooth_1 = 1;
    CMB_smooth_2 = 0;
  } else if(strcmp(realcov_type, "ks_xi-")==0) {
    func_for_cov_G  = &func_for_cov_G_ks_shear;
    func_bin_cov_NG = &bin_cov_NG_ks_shear;
    func_P1 = &Pl2_tab;
    func_P2 = &Glminus_tab;
    CMB_smooth_1 = 1;
    CMB_smooth_2 = 0;
  } else if(strcmp(realcov_type, "gk_gl")==0) {
    func_for_cov_G  = &func_for_cov_G_gk_gl;
    func_bin_cov_NG = &bin_cov_NG_gk_gl;
    func_P1 = &Pl_tab;
    func_P2 = &Pl2_tab;
    CMB_smooth_1 = 1;
    CMB_smooth_2 = 0;
  } else if(strcmp(realcov_type, "ks_gl")==0) {
    func_for_cov_G  = &func_for_cov_G_ks_gl;
    func_bin_cov_NG = &bin_cov_NG_ks_gl;
    func_P1 = &Pl2_tab;
    func_P2 = &Pl2_tab;
    CMB_smooth_1 = 1;
    CMB_smooth_2 = 0;
  } else if(strcmp(realcov_type, "gk_cl")==0) {
    func_for_cov_G  = &func_for_cov_G_gk_cl;
    func_bin_cov_NG = &bin_cov_NG_gk_cl;
    func_P1 = &Pl_tab;
    func_P2 = &Pl_tab;
    CMB_smooth_1 = 1;
    CMB_smooth_2 = 0;
  } else if(strcmp(realcov_type, "ks_cl")==0) {
    func_for_cov_G  = &func_for_cov_G_ks_cl;
    func_bin_cov_NG = &bin_cov_NG_ks_cl;
    func_P1 = &Pl2_tab;
    func_P2 = &Pl_tab;
    CMB_smooth_1 = 1;
    CMB_smooth_2 = 0;
  } else if(strcmp(realcov_type, "gk_gk")==0) {
    func_for_cov_G  = &func_for_cov_G_gk;
    func_bin_cov_NG = &bin_cov_NG_gk_gk;
    func_P1 = &Pl_tab;
    func_P2 = &Pl_tab;
    CMB_smooth_1 = 1;
    CMB_smooth_2 = 1;
  } else if(strcmp(realcov_type, "ks_gk")==0) {
    func_for_cov_G  = &func_for_cov_G_ks_gk;
    func_bin_cov_NG = &bin_cov_NG_ks_gk;
    func_P1 = &Pl2_tab;
    func_P2 = &Pl_tab;
    CMB_smooth_1 = 1;
    CMB_smooth_2 = 1;
  } else if(strcmp(realcov_type, "ks_ks")==0) {
    func_for_cov_G  = &func_for_cov_G_ks;
    func_bin_cov_NG = &bin_cov_NG_ks_ks;
    func_P1 = &Pl2_tab;
    func_P2 = &Pl2_tab;
    CMB_smooth_1 = 1;
    CMB_smooth_2 = 1;
  } else if(strcmp(realcov_type, "kk_xi+")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_shear;
    func_bin_cov_NG = &bin_cov_NG_kk_shear;
    func_P1 = &Pl_tab;
    func_P2 = &Glplus_tab;
    CMB_smooth_1 = 0;
    CMB_smooth_2 = 0;
  } else if(strcmp(realcov_type, "kk_xi-")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_shear;
    func_bin_cov_NG = &bin_cov_NG_kk_shear;
    func_P1 = &Pl_tab;
    func_P2 = &Glminus_tab;
    CMB_smooth_1 = 0;
    CMB_smooth_2 = 0;
  } else if(strcmp(realcov_type, "kk_gl")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_gl;
    func_bin_cov_NG = &bin_cov_NG_kk_gl;
    func_P1 = &Pl_tab;
    func_P2 = &Pl2_tab;
    CMB_smooth_1 = 0;
    CMB_smooth_2 = 0;
  } else if(strcmp(realcov_type, "kk_cl")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_cl;
    func_bin_cov_NG = &bin_cov_NG_kk_cl;
    func_P1 = &Pl_tab;
    func_P2 = &Pl_tab;
    CMB_smooth_1 = 0;
    CMB_smooth_2 = 0;
  } else if(strcmp(realcov_type, "kk_gk")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_gk;
    func_bin_cov_NG = &bin_cov_NG_kk_gk;
    func_P1 = &Pl_tab;
    func_P2 = &Pl_tab;
    CMB_smooth_1 = 0;
    CMB_smooth_2 = 1;
  } else if(strcmp(realcov_type, "kk_ks")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_ks;
    func_bin_cov_NG = &bin_cov_NG_kk_ks;
    func_P1 = &Pl_tab;
    func_P2 = &Pl2_tab;
    CMB_smooth_1 = 0;
    CMB_smooth_2 = 1;
  } else if(strcmp(realcov_type, "kk_kk")==0) {
    func_for_cov_G  = &func_for_cov_G_kk;
    func_bin_cov_NG = &bin_cov_NG_kk_kk;
    func_P1 = &Pl_tab;
    func_P2 = &Pl_tab;
    CMB_smooth_2 = 0;
    CMB_smooth_2 = 0;
  } else {
    printf("cov_real_binned_fullsky: realcov_type \"%s\" not defined!\n", realcov_type); exit(1);
  }

  double l1_double,tri;
  int l1,l2;
  for(i=0; i<like.Ntheta ; i++){
    for(j=0; j<like.Ntheta ; j++){
      cov[i][j] = 0.;
      covNG[i][j] = 0.;
    }
  }

  double triP;
  double covGl1, cov_g_l;
  double beam1, beam2;
  //double  beam;

  for (l1 = 0; l1 < LMAX; l1++){
    l1_double = (double)l1;
	
	beam1 = pow( GaussianBeam(cmb.fwhm, l1, covparams.lmin, covparams.lmax), CMB_smooth_1);
	//beam1 = GaussianBeam_test(cmb.fwhm, l1, covparams.lmin, covparams.lmax, CMB_smooth_1);

    // printf("l1,%d\n", l1);
    cov_g_l = func_for_cov_G(l1_double, z_ar);
    for(i=0; i<like.Ntheta ; i++){
      covGl1 = cov_g_l * func_P1(i,l1);
      // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
      for(j=0; j<like.Ntheta ; j++){
		beam2 = pow( GaussianBeam(cmb.fwhm, l1, covparams.lmin, covparams.lmax), CMB_smooth_2);
		//beam = GaussianBeam_test(cmb.fwhm, l1, l1, covparams.lmin, covparams.lmax, CMB_smooth_1, CMB_smooth_2);
        cov[i][j] += covGl1 * func_P2(j,l1) * beam1 * beam2;
      }
    }

    if(FLAG_NG){
      for (l2 = 0; l2 < LMAX; l2++){
        tri = func_bin_cov_NG(l1_double,(double)l2,z_ar);
		beam2 = pow( GaussianBeam(cmb.fwhm, l2, covparams.lmin, covparams.lmax), CMB_smooth_2);
		//beam = GaussianBeam_test(cmb.fwhm, l1, l2, covparams.lmin, covparams.lmax, CMB_smooth_1, CMB_smooth_2);

        for(i=0; i<like.Ntheta ; i++){
          triP = tri * func_P1(i,l1);
          // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
          for(j=0; j<like.Ntheta ; j++){
            covNG[i][j] += triP * func_P2(j,l2) * beam1 * beam2;
          }
        }
      }
    }
  }
  for(i=0; i<like.Ntheta; i++){if(N[i]){cov[i][i] += N[i]/w_mask(theta[i]);}}
}


// Template routine
void cov_mix_binned_fullsky(double **cov, double **covNG, char *mixcov_type, int *z_ar, int FLAG_NG, double *theta, double *dtheta, double *ell){

  int i,j;
  static int LMAX = 50000;
  int CMB_smooth_1 = 0, CMB_smooth_2 = 0; // flag for CMB smoothing kernel

  // double N[like.Ntheta];
  // for(i=0;i<like.Ntheta;i++) {N[i] = 0.;}

  double (*func_for_cov_G)(double, int*);
  double (*func_bin_cov_NG)(double, double, int*);
  double (*func_P2)(int, int);

  if(strcmp(mixcov_type, "kk_xi+")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_shear;
    func_bin_cov_NG = &bin_cov_NG_kk_shear;
    func_P1 = &CMB_lensing_mat;
    func_P2 = &Glplus_tab;
    CMB_smooth_1 = 0;
    CMB_smooth_2 = 0;
  } else if(strcmp(mixcov_type, "kk_xi-")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_shear;
    func_bin_cov_NG = &bin_cov_NG_kk_shear;
    func_P1 = &CMB_lensing_mat;
    func_P2 = &Glminus_tab;
    CMB_smooth_1 = 0;
    CMB_smooth_2 = 0;
  } else if(strcmp(mixcov_type, "kk_gl")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_gl;
    func_bin_cov_NG = &bin_cov_NG_kk_gl;
    func_P1 = &CMB_lensing_mat;
    func_P2 = &Pl2_tab;
    CMB_smooth_1 = 0;
    CMB_smooth_2 = 0;
  } else if(strcmp(mixcov_type, "kk_cl")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_cl;
    func_bin_cov_NG = &bin_cov_NG_kk_cl;
    func_P1 = &CMB_lensing_mat;
    func_P2 = &Pl_tab;
    CMB_smooth_1 = 0;
    CMB_smooth_2 = 0;
  } else if(strcmp(mixcov_type, "kk_gk")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_gk;
    func_bin_cov_NG = &bin_cov_NG_kk_gk;
    func_P1 = &CMB_lensing_mat;
    func_P2 = &Pl_tab;
    CMB_smooth_1 = 0;
    CMB_smooth_2 = 1;
  } else if(strcmp(mixcov_type, "kk_ks")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_ks;
    func_bin_cov_NG = &bin_cov_NG_kk_ks;
    func_P1 = &CMB_lensing_mat;
    func_P2 = &Pl2_tab;
    CMB_smooth_1 = 0;
    CMB_smooth_2 = 1;
  } else {
    printf("cov_mix_binned_fullsky: mixcov_type \"%s\" not defined!\n", mixcov_type); exit(1);
  }

  double l1_double,tri;
  int l1,l2;
  double beam1, beam2;
  //double beam;
  for(i=0; i<like.Ncl ; i++){
    for(j=0; j<like.Ntheta ; j++){
      cov[i][j] = 0.;
      covNG[i][j] = 0.;
    }
  }

  double triP;
  double covGl1;

  double cov_g_l[LMAX];
  //printf("lmin,lmax, %d, %d\n", (int)floor(like.lmin),(int)ceil(like.lmax));
  for(l1=(int)ceil(like.lmin); l1<(int)ceil(like.lmax)+1; l1++){
    cov_g_l[l1] = func_for_cov_G((double)l1, z_ar);
    //printf("cov_g_l: %d, %le\n", l1,cov_g_l[l1]);
  }
  int l1_min, l1_max;
  int N_l1;
  
  // Loop 1: Harmonic space
  /* TODO
      - i = [0, Ncl-1] -> [0, Nbp-1]
      - decide l1_min and l1_max
  */
  for(i=0; i<like.Nbp; i++){
    l1_min = (int)ceil(ell[i]);
    l1_max = (int)ceil(ell[i+1])-1;
    N_l1 = l1_max - l1_min + 1;

    // Loop 2: Configuration space
    for(j=0; j<like.Ntheta ; j++){
      
      // sum over ell
      for(l1=l1_min; l1<=l1_max; l1++){
		    
        // Gaussian
        beam1 = pow( GaussianBeam(cmb.fwhm, l1, 
          covparams.lmin, covparams.lmax), CMB_smooth_1 + CMB_smooth_2);

        cov[i][j] += cov_g_l[l1]*func_P1(i,l1)*func_P2(j,l1)/N_l1*beam1;
        
        // Non-Gaussian
        if(FLAG_NG){
          l1_double = (double)l1;
          for (l2 = 0; l2 < LMAX; l2++){
            beam1 = pow( GaussianBeam(cmb.fwhm, l1, 
              covparams.lmin, covparams.lmax), CMB_smooth_1);
            beam2 = pow( GaussianBeam(cmb.fwhm, l2, 
              covparams.lmin, covparams.lmax), CMB_smooth_2);
            tri = func_bin_cov_NG(l1_double,(double)l2,z_ar);
            
            covNG[i][j] += tri * func_P2(j,l2) / N_l1 * beam1 * beam2;
          }
        }
      }// end sum over ell
    }// end Loop 2
  }// End Loop 1
}

// Template routine: Fourier averaged band power
void cov_fourier_binned(double **cov, double **covNG, char *cov_type, int *z_ar, int FLAG_NG, double *ell){

  int i,j;
  static int LMAX = 50000;
  // Did not implement Gaussian smoothing since this one is a theoretical study.

  // double N[like.Ntheta];
  // for(i=0;i<like.Ntheta;i++) {N[i] = 0.;}

  double (*func_for_cov_G)(double, int*);
  double (*func_bin_cov_NG)(double, double, int*);


  if(strcmp(cov_type, "ss_ss")==0) {
    func_for_cov_G  = &func_for_cov_G_shear;
    func_bin_cov_NG = &bin_cov_NG_shear_shear;
  } else if(strcmp(cov_type, "gl_ss")==0) {
    func_for_cov_G  = &func_for_cov_G_gl_shear;
    func_bin_cov_NG = &bin_cov_NG_gl_shear;
  } else if(strcmp(cov_type, "cl_ss")==0) {
    func_for_cov_G  = &func_for_cov_G_cl_shear;
    func_bin_cov_NG = &bin_cov_NG_cl_shear;
  } else if(strcmp(cov_type, "gl_gl")==0) {
    func_for_cov_G  = &func_for_cov_G_gl;
    func_bin_cov_NG = &bin_cov_NG_gl_gl;
  } else if(strcmp(cov_type, "cl_gl")==0) {
    func_for_cov_G  = &func_for_cov_G_cl_gl;
    func_bin_cov_NG = &bin_cov_NG_cl_gl;
  } else if(strcmp(cov_type, "cl_cl")==0) {
    func_for_cov_G  = &func_for_cov_G_cl;
    func_bin_cov_NG = &bin_cov_NG_cl_cl;
  } else if(strcmp(cov_type, "gk_ss")==0) {
    func_for_cov_G  = &func_for_cov_G_gk_shear;
    func_bin_cov_NG = &bin_cov_NG_gk_shear;
  } else if(strcmp(cov_type, "ks_ss")==0) {
    func_for_cov_G  = &func_for_cov_G_ks_shear;
    func_bin_cov_NG = &bin_cov_NG_ks_shear;
  } else if(strcmp(cov_type, "gk_gl")==0) {
    func_for_cov_G  = &func_for_cov_G_gk_gl;
    func_bin_cov_NG = &bin_cov_NG_gk_gl;
  } else if(strcmp(cov_type, "ks_gl")==0) {
    func_for_cov_G  = &func_for_cov_G_ks_gl;
    func_bin_cov_NG = &bin_cov_NG_ks_gl;
  } else if(strcmp(cov_type, "gk_cl")==0) {
    func_for_cov_G  = &func_for_cov_G_gk_cl;
    func_bin_cov_NG = &bin_cov_NG_gk_cl;
  } else if(strcmp(cov_type, "ks_cl")==0) {
    func_for_cov_G  = &func_for_cov_G_ks_cl;
    func_bin_cov_NG = &bin_cov_NG_ks_cl;
  } else if(strcmp(cov_type, "gk_gk")==0) {
    func_for_cov_G  = &func_for_cov_G_gk;
    func_bin_cov_NG = &bin_cov_NG_gk_gk;
  } else if(strcmp(cov_type, "ks_gk")==0) {
    func_for_cov_G  = &func_for_cov_G_ks_gk;
    func_bin_cov_NG = &bin_cov_NG_ks_gk;
  } else if(strcmp(cov_type, "ks_ks")==0) {
    func_for_cov_G  = &func_for_cov_G_ks;
    func_bin_cov_NG = &bin_cov_NG_ks_ks;

  } else if(strcmp(cov_type, "kk_ss")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_shear;
    func_bin_cov_NG = &bin_cov_NG_kk_shear;
  } else if(strcmp(cov_type, "kk_gl")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_gl;
    func_bin_cov_NG = &bin_cov_NG_kk_gl;
  } else if(strcmp(cov_type, "kk_cl")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_cl;
    func_bin_cov_NG = &bin_cov_NG_kk_cl;
  } else if(strcmp(cov_type, "kk_gk")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_gk;
    func_bin_cov_NG = &bin_cov_NG_kk_gk;
  } else if(strcmp(cov_type, "kk_ks")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_ks;
    func_bin_cov_NG = &bin_cov_NG_kk_ks;

  } else if(strcmp(cov_type, "kk_kk")==0) {
    func_for_cov_G  = &func_for_cov_G_kk;
    func_bin_cov_NG = &bin_cov_NG_kk_kk;

  } else {
    printf("cov_fourier_binned: cov_type \"%s\" not defined!\n", cov_type); exit(1);
  }

  double l1_double,tri;
  int l1,l2;
  for(i=0; i<like.Ncl ; i++){
    for(j=0; j<like.Ncl ; j++){
      cov[i][j] = 0.;
      covNG[i][j] = 0.;
    }
  }

  double triP;
  double covGl1;

  double cov_g_l[LMAX];
    // printf("lmin,lmax, %d, %d\n", (int)floor(like.lmin),(int)ceil(like.lmax));
  for(l1=(int)ceil(like.lmin); l1<(int)ceil(like.lmax)+1; l1++){
    cov_g_l[l1] = func_for_cov_G((double)l1, z_ar);
    // printf("cov_g_l: %d, %le\n", l1,cov_g_l[l1]);
  }

  int l1_min, l1_max, l2_min, l2_max;
  int N_l1, N_l2;
  for(i=0; i<like.Ncl; i++){
    l1_min = (int)ceil(ell[i]);
    l1_max = (int)ceil(ell[i+1])-1;
    N_l1 = l1_max - l1_min + 1;
    for(l1=l1_min; l1<=l1_max; l1++){
      cov[i][i] += cov_g_l[l1];

      if(FLAG_NG){
        l1_double = (double)l1;
        for(j=0; j<like.Ncl ; j++){
          l2_min = (int)ceil(ell[j]);
          l2_max = (int)ceil(ell[j+1])-1;
          N_l2 = l2_max - l2_min + 1;
          for (l2 = l2_min; l2 <= l2_max; l2++){
            tri = func_bin_cov_NG(l1_double,(double)l2,z_ar);
            covNG[i][j] += tri/(N_l2*N_l1);
            // printf("here!!!%d, %d, %le\n", z_ar[0], z_ar[1], bin_cov_NG_shear_shear(3.366055e+01, 3.366055e+01,z_ar));
          }
          // covNG[i][j] /= (float)(N_l2);
        }
      }
    }
    cov[i][i] /= (float)(N_l1*N_l1);
  }
}

// Real cov
void cov_shear_shear_real_binned_fullsky(double **cov, double **covNG, int z1,int z2,int z3,int z4,int pm1,int pm2, int FLAG_NG, double *theta, double *dtheta){
  static char xipp[]="xi+_xi+", xipm[]="xi+_xi-", ximp[]="xi-_xi+", ximm[]="xi-_xi-";
  char realcov_type[8];
  if(pm1==1) {
    if(pm2==1) { //++
      strcpy(realcov_type, xipp);
    }else { //+-
      strcpy(realcov_type, xipm);
    }
  }else {
    if(pm2==1) { //-+
      strcpy(realcov_type, ximp);
    }else { //--
      strcpy(realcov_type, ximm);
    }
  }
  int z_ar[4];
  z_ar[0]=z1; z_ar[1]=z2; z_ar[2]=z3; z_ar[3]=z4;
  cov_real_binned_fullsky(cov, covNG, realcov_type, z_ar, FLAG_NG, theta, dtheta);
}

void cov_gl_shear_real_binned_fullsky(double **cov, double **covNG, int z1,int z2,int z3,int z4,int pm, int FLAG_NG, double *theta, double *dtheta){
  static char xip[]="gl_xi+", xim[]="gl_xi-";
  char realcov_type[8];
  if(pm==1) {strcpy(realcov_type, xip);}
  else {strcpy(realcov_type, xim);}
  int z_ar[4];
  z_ar[0]=z1; z_ar[1]=z2; z_ar[2]=z3; z_ar[3]=z4;
  cov_real_binned_fullsky(cov, covNG, realcov_type, z_ar, FLAG_NG, theta, dtheta);
}

void cov_cl_shear_real_binned_fullsky(double **cov, double **covNG, int z1,int z2,int z3,int z4,int pm, int FLAG_NG, double *theta, double *dtheta){
  static char xip[]="cl_xi+", xim[]="cl_xi-";
  char realcov_type[8];
  if(pm==1) {strcpy(realcov_type, xip);}
  else {strcpy(realcov_type, xim);}
  int z_ar[4];
  z_ar[0]=z1; z_ar[1]=z2; z_ar[2]=z3; z_ar[3]=z4;
  cov_real_binned_fullsky(cov, covNG, realcov_type, z_ar, FLAG_NG, theta, dtheta);
}

void cov_cl_gl_real_binned_fullsky(double **cov, double **covNG, int z1,int z2,int z3,int z4, int FLAG_NG, double *theta, double *dtheta){
  char realcov_type[] = "cl_gl";
  int z_ar[4];
  z_ar[0]=z1; z_ar[1]=z2; z_ar[2]=z3; z_ar[3]=z4;
  cov_real_binned_fullsky(cov, covNG, realcov_type, z_ar, FLAG_NG, theta, dtheta);
}

void cov_gl_gl_real_binned_fullsky(double **cov, double **covNG, int z1,int z2,int z3,int z4, int FLAG_NG, double *theta, double *dtheta){
  char realcov_type[] = "gl_gl";
  int z_ar[4];
  z_ar[0]=z1; z_ar[1]=z2; z_ar[2]=z3; z_ar[3]=z4;
  cov_real_binned_fullsky(cov, covNG, realcov_type, z_ar, FLAG_NG, theta, dtheta);
}

void cov_cl_cl_real_binned_fullsky(double **cov, double **covNG, int z1,int z2,int z3,int z4, int FLAG_NG, double *theta, double *dtheta){
  char realcov_type[] = "cl_cl";
  int z_ar[4];
  z_ar[0]=z1; z_ar[1]=z2; z_ar[2]=z3; z_ar[3]=z4;
  cov_real_binned_fullsky(cov, covNG, realcov_type, z_ar, FLAG_NG, theta, dtheta);
}

// CMB lensing
void cov_ks_shear_real_binned_fullsky(double **cov, double **covNG, int z1,int z3,int z4,int pm, int FLAG_NG, double *theta, double *dtheta){
  static char xip[]="ks_xi+", xim[]="ks_xi-";
  char realcov_type[8];
  if(pm==1) {strcpy(realcov_type, xip);}
  else {strcpy(realcov_type, xim);}
  int z_ar[3];
  z_ar[0]=z1; z_ar[1]=z3; z_ar[2]=z4;
  cov_real_binned_fullsky(cov, covNG, realcov_type, z_ar, FLAG_NG, theta, dtheta);
}

void cov_gk_shear_real_binned_fullsky(double **cov, double **covNG, int z1,int z3,int z4,int pm, int FLAG_NG, double *theta, double *dtheta){
  static char xip[]="gk_xi+", xim[]="gk_xi-";
  char realcov_type[8];
  if(pm==1) {strcpy(realcov_type, xip);}
  else {strcpy(realcov_type, xim);}
  int z_ar[3];
  z_ar[0]=z1; z_ar[1]=z3; z_ar[2]=z4;
  cov_real_binned_fullsky(cov, covNG, realcov_type, z_ar, FLAG_NG, theta, dtheta);
}

void cov_gk_gl_real_binned_fullsky(double **cov, double **covNG, int z1,int z3,int z4, int FLAG_NG, double *theta, double *dtheta){
  char realcov_type[] = "gk_gl";
  int z_ar[3];
  z_ar[0]=z1; z_ar[1]=z3; z_ar[2]=z4;
  cov_real_binned_fullsky(cov, covNG, realcov_type, z_ar, FLAG_NG, theta, dtheta);
}

void cov_ks_gl_real_binned_fullsky(double **cov, double **covNG, int z1,int z3,int z4, int FLAG_NG, double *theta, double *dtheta){
  char realcov_type[] = "ks_gl";
  int z_ar[3];
  z_ar[0]=z1; z_ar[1]=z3; z_ar[2]=z4;
  cov_real_binned_fullsky(cov, covNG, realcov_type, z_ar, FLAG_NG, theta, dtheta);
}

void cov_gk_cl_real_binned_fullsky(double **cov, double **covNG, int z1,int z3,int z4, int FLAG_NG, double *theta, double *dtheta){
  char realcov_type[] = "gk_cl";
  int z_ar[3];
  z_ar[0]=z1; z_ar[1]=z3; z_ar[2]=z4;
  cov_real_binned_fullsky(cov, covNG, realcov_type, z_ar, FLAG_NG, theta, dtheta);
}

void cov_ks_cl_real_binned_fullsky(double **cov, double **covNG, int z1,int z3,int z4, int FLAG_NG, double *theta, double *dtheta){
  char realcov_type[] = "ks_cl";
  int z_ar[3];
  z_ar[0]=z1; z_ar[1]=z3; z_ar[2]=z4;
  cov_real_binned_fullsky(cov, covNG, realcov_type, z_ar, FLAG_NG, theta, dtheta);
}

void cov_gk_gk_real_binned_fullsky(double **cov, double **covNG, int z1,int z2, int FLAG_NG, double *theta, double *dtheta){
  char realcov_type[] = "gk_gk";
  int z_ar[2];
  z_ar[0]=z1; z_ar[1]=z2;
  cov_real_binned_fullsky(cov, covNG, realcov_type, z_ar, FLAG_NG, theta, dtheta);
}

void cov_ks_gk_real_binned_fullsky(double **cov, double **covNG, int z1,int z2, int FLAG_NG, double *theta, double *dtheta){
  char realcov_type[] = "ks_gk";
  int z_ar[2];
  z_ar[0]=z1; z_ar[1]=z2;
  cov_real_binned_fullsky(cov, covNG, realcov_type, z_ar, FLAG_NG, theta, dtheta);
}

void cov_ks_ks_real_binned_fullsky(double **cov, double **covNG, int z1,int z2, int FLAG_NG, double *theta, double *dtheta){
  char realcov_type[] = "ks_ks";
  int z_ar[2];
  z_ar[0]=z1; z_ar[1]=z2;
  cov_real_binned_fullsky(cov, covNG, realcov_type, z_ar, FLAG_NG, theta, dtheta);
}
///// kk real
// real cov 
void cov_kk_shear_real_binned_fullsky(double **cov, double **covNG, int z3,int z4,int pm, int FLAG_NG, double *theta, double *dtheta){
  static char xip[]="kk_xi+", xim[]="kk_xi-";
  char mixcov_type[8];
  if(pm==1) {strcpy(mixcov_type, xip);}
  else {strcpy(mixcov_type, xim);}
  int z_ar[2];
  z_ar[0]=z3; z_ar[1]=z4;
  cov_real_binned_fullsky(cov, covNG, mixcov_type, z_ar, FLAG_NG, theta, dtheta);
}

void cov_kk_gl_real_binned_fullsky(double **cov, double **covNG, int z3,int z4, int FLAG_NG, double *theta, double *dtheta){
  char mixcov_type[] = "kk_gl";
  int z_ar[2];
  z_ar[0]=z3; z_ar[1]=z4;
  cov_real_binned_fullsky(cov, covNG, mixcov_type, z_ar, FLAG_NG, theta, dtheta);
}

void cov_kk_cl_real_binned_fullsky(double **cov, double **covNG, int z3,int z4, int FLAG_NG, double *theta, double *dtheta){
  char mixcov_type[] = "kk_cl";
  int z_ar[2];
  z_ar[0]=z3; z_ar[1]=z4;
  cov_real_binned_fullsky(cov, covNG, mixcov_type, z_ar, FLAG_NG, theta, dtheta);
}

void cov_kk_gk_real_binned_fullsky(double **cov, double **covNG, int zl, int FLAG_NG, double *theta, double *dtheta){
  char mixcov_type[] = "kk_gk";
  int z_ar[1];
  z_ar[0]=zl;
  cov_real_binned_fullsky(cov, covNG, mixcov_type, z_ar, FLAG_NG, theta, dtheta);
}

void cov_kk_ks_real_binned_fullsky(double **cov, double **covNG, int zs, int FLAG_NG, double *theta, double *dtheta){
  char mixcov_type[] = "kk_ks";
  int z_ar[1];
  z_ar[0]=zs;
  cov_real_binned_fullsky(cov, covNG, mixcov_type, z_ar, FLAG_NG, theta, dtheta);
}

void cov_kk_kk_real_binned_fullsky(double **cov, double **covNG, int FLAG_NG, double *theta, double *dtheta){
  char mixcov_type[] = "kk_kk";
  int z_ar[1];
  z_ar[0]=-1;
  cov_real_binned_fullsky(cov, covNG, mixcov_type, z_ar, FLAG_NG, theta, dtheta);
}


////


// mixed cov 
void cov_kk_shear_mix_binned_fullsky(double **cov, double **covNG, int z3,int z4,int pm, int FLAG_NG, double *theta, double *dtheta, double *ell){
  static char xip[]="kk_xi+", xim[]="kk_xi-";
  char mixcov_type[8];
  if(pm==1) {strcpy(mixcov_type, xip);}
  else {strcpy(mixcov_type, xim);}
  int z_ar[2];
  z_ar[0]=z3; z_ar[1]=z4;
  cov_mix_binned_fullsky(cov, covNG, mixcov_type, z_ar, FLAG_NG, theta, dtheta, ell);
}

void cov_kk_gl_mix_binned_fullsky(double **cov, double **covNG, int z3,int z4, int FLAG_NG, double *theta, double *dtheta, double *ell){
  char mixcov_type[] = "kk_gl";
  int z_ar[2];
  z_ar[0]=z3; z_ar[1]=z4;
  cov_mix_binned_fullsky(cov, covNG, mixcov_type, z_ar, FLAG_NG, theta, dtheta, ell);
}

void cov_kk_cl_mix_binned_fullsky(double **cov, double **covNG, int z3,int z4, int FLAG_NG, double *theta, double *dtheta, double *ell){
  char mixcov_type[] = "kk_cl";
  int z_ar[2];
  z_ar[0]=z3; z_ar[1]=z4;
  cov_mix_binned_fullsky(cov, covNG, mixcov_type, z_ar, FLAG_NG, theta, dtheta, ell);
}

void cov_kk_gk_mix_binned_fullsky(double **cov, double **covNG, int zl, int FLAG_NG, double *theta, double *dtheta, double *ell){
  char mixcov_type[] = "kk_gk";
  int z_ar[1];
  z_ar[0]=zl;
  cov_mix_binned_fullsky(cov, covNG, mixcov_type, z_ar, FLAG_NG, theta, dtheta, ell);
}

void cov_kk_ks_mix_binned_fullsky(double **cov, double **covNG, int zs, int FLAG_NG, double *theta, double *dtheta, double *ell){
  char mixcov_type[] = "kk_ks";
  int z_ar[1];
  z_ar[0]=zs;
  cov_mix_binned_fullsky(cov, covNG, mixcov_type, z_ar, FLAG_NG, theta, dtheta, ell);
}

/// Full Fourier 6x2pt cov - band averaged


// Real cov
void cov_shear_shear_fourier_binned(double **cov, double **covNG, int z1,int z2,int z3,int z4, int FLAG_NG, double *ell){
  char cov_type[]= "ss_ss";
  int z_ar[4];
  z_ar[0]=z1; z_ar[1]=z2; z_ar[2]=z3; z_ar[3]=z4;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_gl_shear_fourier_binned(double **cov, double **covNG, int z1,int z2,int z3,int z4, int FLAG_NG, double *ell){
  char cov_type[]= "gl_ss";
  int z_ar[4];
  z_ar[0]=z1; z_ar[1]=z2; z_ar[2]=z3; z_ar[3]=z4;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_cl_shear_fourier_binned(double **cov, double **covNG, int z1,int z2,int z3,int z4, int FLAG_NG, double *ell){
  char cov_type[]= "cl_ss";
  int z_ar[4];
  z_ar[0]=z1; z_ar[1]=z2; z_ar[2]=z3; z_ar[3]=z4;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_cl_gl_fourier_binned(double **cov, double **covNG, int z1,int z2,int z3,int z4, int FLAG_NG, double *ell){
  char cov_type[] = "cl_gl";
  int z_ar[4];
  z_ar[0]=z1; z_ar[1]=z2; z_ar[2]=z3; z_ar[3]=z4;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_gl_gl_fourier_binned(double **cov, double **covNG, int z1,int z2,int z3,int z4, int FLAG_NG, double *ell){
  char cov_type[] = "gl_gl";
  int z_ar[4];
  z_ar[0]=z1; z_ar[1]=z2; z_ar[2]=z3; z_ar[3]=z4;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_cl_cl_fourier_binned(double **cov, double **covNG, int z1,int z2,int z3,int z4, int FLAG_NG, double *ell){
  char cov_type[] = "cl_cl";
  int z_ar[4];
  z_ar[0]=z1; z_ar[1]=z2; z_ar[2]=z3; z_ar[3]=z4;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

// CMB lensing
void cov_ks_shear_fourier_binned(double **cov, double **covNG, int z1,int z3,int z4, int FLAG_NG, double *ell){
  char cov_type[] = "ks_ss";
  int z_ar[3];
  z_ar[0]=z1; z_ar[1]=z3; z_ar[2]=z4;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_gk_shear_fourier_binned(double **cov, double **covNG, int z1,int z3,int z4, int FLAG_NG, double *ell){
  char cov_type[] = "gk_ss";
  int z_ar[3];
  z_ar[0]=z1; z_ar[1]=z3; z_ar[2]=z4;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_gk_gl_fourier_binned(double **cov, double **covNG, int z1,int z3,int z4, int FLAG_NG, double *ell){
  char cov_type[] = "gk_gl";
  int z_ar[3];
  z_ar[0]=z1; z_ar[1]=z3; z_ar[2]=z4;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_ks_gl_fourier_binned(double **cov, double **covNG, int z1,int z3,int z4, int FLAG_NG, double *ell){
  char cov_type[] = "ks_gl";
  int z_ar[3];
  z_ar[0]=z1; z_ar[1]=z3; z_ar[2]=z4;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_gk_cl_fourier_binned(double **cov, double **covNG, int z1,int z3,int z4, int FLAG_NG, double *ell){
  char cov_type[] = "gk_cl";
  int z_ar[3];
  z_ar[0]=z1; z_ar[1]=z3; z_ar[2]=z4;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_ks_cl_fourier_binned(double **cov, double **covNG, int z1,int z3,int z4, int FLAG_NG, double *ell){
  char cov_type[] = "ks_cl";
  int z_ar[3];
  z_ar[0]=z1; z_ar[1]=z3; z_ar[2]=z4;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_gk_gk_fourier_binned(double **cov, double **covNG, int z1,int z2, int FLAG_NG, double *ell){
  char cov_type[] = "gk_gk";
  int z_ar[2];
  z_ar[0]=z1; z_ar[1]=z2;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_ks_gk_fourier_binned(double **cov, double **covNG, int z1,int z2, int FLAG_NG, double *ell){
  char cov_type[] = "ks_gk";
  int z_ar[2];
  z_ar[0]=z1; z_ar[1]=z2;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_ks_ks_fourier_binned(double **cov, double **covNG, int z1,int z2, int FLAG_NG, double *ell){
  char cov_type[] = "ks_ks";
  int z_ar[2];
  z_ar[0]=z1; z_ar[1]=z2;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}
///// kk fourier
void cov_kk_shear_fourier_binned(double **cov, double **covNG, int z3,int z4, int FLAG_NG, double *ell){
  char cov_type[] = "kk_ss";
  int z_ar[2];
  z_ar[0]=z3; z_ar[1]=z4;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_kk_gl_fourier_binned(double **cov, double **covNG, int z3,int z4, int FLAG_NG, double *ell){
  char cov_type[] = "kk_gl";
  int z_ar[2];
  z_ar[0]=z3; z_ar[1]=z4;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_kk_cl_fourier_binned(double **cov, double **covNG, int z3,int z4, int FLAG_NG, double *ell){
  char cov_type[] = "kk_cl";
  int z_ar[2];
  z_ar[0]=z3; z_ar[1]=z4;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_kk_gk_fourier_binned(double **cov, double **covNG, int zl, int FLAG_NG, double *ell){
  char cov_type[] = "kk_gk";
  int z_ar[1];
  z_ar[0]=zl;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_kk_ks_fourier_binned(double **cov, double **covNG, int zs, int FLAG_NG, double *ell){
  char cov_type[] = "kk_ks";
  int z_ar[1];
  z_ar[0]=zs;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_kk_kk_fourier_binned(double **cov, double **covNG, int FLAG_NG, double *ell){
  char cov_type[] = "kk_kk";
  int z_ar[1];
  z_ar[0]=-1;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}




/********** Functions for differnt covariances ***************************************/
double func_for_cov_G_shear_noNN(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[0]; n2 = ar[1]; n3 = ar[2]; n4 = ar[3];
  // printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_shear_tomo(l,n1,n3);
  C24 = C_shear_tomo(l,n2,n4);
  C14 = C_shear_tomo(l,n1,n4);
  C23 = C_shear_tomo(l,n2,n3);
  
  if (n1 == n3){N13= pow(survey.sigma_e,2.0)/(2.*nsource(n1)*survey.n_gal_conversion_factor);}
  if (n1 == n4){N14= pow(survey.sigma_e,2.0)/(2.*nsource(n1)*survey.n_gal_conversion_factor);}
  if (n2 == n3){N23= pow(survey.sigma_e,2.0)/(2.*nsource(n2)*survey.n_gal_conversion_factor);}
  if (n2 == n4){N24= pow(survey.sigma_e,2.0)/(2.*nsource(n2)*survey.n_gal_conversion_factor);}
  
  return (C13*C24+C13*N24+N13*C24 + C14*C23+C14*N23+N14*C23)*4.*M_PI/(survey.area*survey.area_conversion_factor* (2.*l+1.));
}

double func_for_cov_G_shear(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[0]; n2 = ar[1]; n3 = ar[2]; n4 = ar[3];
  // printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_shear_tomo(l,n1,n3);
  C24 = C_shear_tomo(l,n2,n4);
  C14 = C_shear_tomo(l,n1,n4);
  C23 = C_shear_tomo(l,n2,n3);
  
  if (n1 == n3){N13= pow(survey.sigma_e,2.0)/(2.*nsource(n1)*survey.n_gal_conversion_factor);}
  if (n1 == n4){N14= pow(survey.sigma_e,2.0)/(2.*nsource(n1)*survey.n_gal_conversion_factor);}
  if (n2 == n3){N23= pow(survey.sigma_e,2.0)/(2.*nsource(n2)*survey.n_gal_conversion_factor);}
  if (n2 == n4){N24= pow(survey.sigma_e,2.0)/(2.*nsource(n2)*survey.n_gal_conversion_factor);}
  
  return ((C13+N13)*(C24+N24) + (C14+N14)*(C23+N23))*4.*M_PI/(survey.area*survey.area_conversion_factor* (2.*l+1.));
}

double func_for_cov_G_cl_noNN(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[0]; n2 = ar[1]; n3 = ar[2]; n4 = ar[3];
  //printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_cl_tomo(l,n1,n3);
  C24 = C_cl_tomo(l,n2,n4);
  C14 = C_cl_tomo(l,n1,n4);
  C23 = C_cl_tomo(l,n2,n3);
  
  if (n1 == n3){N13= 1./(nlens(n1)*survey.n_gal_conversion_factor);}
  if (n1 == n4){N14= 1./(nlens(n1)*survey.n_gal_conversion_factor);}
  if (n2 == n3){N23= 1./(nlens(n2)*survey.n_gal_conversion_factor);}
  if (n2 == n4){N24= 1./(nlens(n2)*survey.n_gal_conversion_factor);}

  return (C13*C24+C13*N24+N13*C24 + C14*C23+C14*N23+N14*C23)*4.*M_PI/(survey.area*survey.area_conversion_factor* (2.*l+1.));
}
double func_for_cov_G_cl(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[0]; n2 = ar[1]; n3 = ar[2]; n4 = ar[3];
  //printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_cl_tomo(l,n1,n3);
  C24 = C_cl_tomo(l,n2,n4);
  C14 = C_cl_tomo(l,n1,n4);
  C23 = C_cl_tomo(l,n2,n3);
  
  if (n1 == n3){N13= 1./(nlens(n1)*survey.n_gal_conversion_factor);}
  if (n1 == n4){N14= 1./(nlens(n1)*survey.n_gal_conversion_factor);}
  if (n2 == n3){N23= 1./(nlens(n2)*survey.n_gal_conversion_factor);}
  if (n2 == n4){N24= 1./(nlens(n2)*survey.n_gal_conversion_factor);}

  return ((C13+N13)*(C24+N24) + (C14+N14)*(C23+N23))*4.*M_PI/(survey.area*survey.area_conversion_factor* (2.*l+1.));
}

double func_for_cov_G_cl_gl(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[0]; n2 = ar[1]; n3 = ar[2]; n4 = ar[3];
  //printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_cl_tomo(l,n1,n3);
  C24 = C_gl_tomo(l,n2,n4);
  C14 = C_gl_tomo(l,n1,n4);
  C23 = C_cl_tomo(l,n2,n3);
  if (n1 == n3){N13= 1./(nlens(n1)*survey.n_gal_conversion_factor);}
  if (n2 == n3){N23= 1./(nlens(n1)*survey.n_gal_conversion_factor);}
  
  return (C13*C24+C13*N24+C24*N13+C14*C23+C14*N23+C23*N14)*4.*M_PI/(survey.area*survey.area_conversion_factor* (2.*l+1.));
}


double func_for_cov_G_cl_shear(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[0]; n2 = ar[1]; n3 = ar[2]; n4 = ar[3];
  //printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_gl_tomo(l,n1,n3);
  C24 = C_gl_tomo(l,n2,n4);
  C14 = C_gl_tomo(l,n1,n4);
  C23 = C_gl_tomo(l,n2,n3);

  return (C13*C24+ C14*C23)*4.*M_PI/(survey.area*survey.area_conversion_factor* (2.*l+1.));
}

double func_for_cov_G_gl_noNN(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[0]; n2 = ar[1]; n3 = ar[2]; n4 = ar[3];
  //printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_cl_tomo(l,n1,n3);
  C24 = C_shear_tomo(l,n2,n4);
  C14 = C_gl_tomo(l,n1,n4);
  C23 = C_gl_tomo(l,n3,n2);
  
  if (n1 == n3){N13= 1./(nlens(n1)*survey.n_gal_conversion_factor);}
  if (n2 == n4){N24= pow(survey.sigma_e,2.0)/(2.0*nsource(n2)*survey.n_gal_conversion_factor);}
  
  return (C13*C24+C13*N24+N13*C24 + C14*C23+C14*N23+N14*C23)*4.*M_PI/(survey.area*survey.area_conversion_factor* (2.*l+1.));
}
double func_for_cov_G_gl(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[0]; n2 = ar[1]; n3 = ar[2]; n4 = ar[3];
  //printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_cl_tomo(l,n1,n3);
  C24 = C_shear_tomo(l,n2,n4);
  C14 = C_gl_tomo(l,n1,n4);
  C23 = C_gl_tomo(l,n3,n2);
  
  if (n1 == n3){N13= 1./(nlens(n1)*survey.n_gal_conversion_factor);}
  if (n2 == n4){N24= pow(survey.sigma_e,2.0)/(2.0*nsource(n2)*survey.n_gal_conversion_factor);}
  
  return ((C13+N13)*(C24+N24) + (C14+N14)*(C23+N23))*4.*M_PI/(survey.area*survey.area_conversion_factor* (2.*l+1.));
}

double func_for_cov_G_gl_shear(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[0]; n2 = ar[1]; n3 = ar[2]; n4 = ar[3];
  //printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_gl_tomo(l,n1,n3);
  C24 = C_shear_tomo(l,n2,n4);
  C14 = C_gl_tomo(l,n1,n4);
  C23 = C_shear_tomo(l,n2,n3);
  if (n2 == n4){N24= pow(survey.sigma_e,2.0)/(2.0*nsource(n2)*survey.n_gal_conversion_factor);}
  if (n2 == n3){N23= pow(survey.sigma_e,2.0)/(2.0*nsource(n2)*survey.n_gal_conversion_factor);}
 
  return (C13*C24+C13*N24+N13*C24 + C14*C23+C14*N23+N14*C23)*4.*M_PI/(survey.area*survey.area_conversion_factor* (2.*l+1.));
}

// CMB lensing
double func_for_cov_G_gk_shear(double l, int *ar){
  double fsky = survey.area*survey.area_conversion_factor/(4.*M_PI);
  double C13, C14, C23, C24;
  int zl,zs1,zs2;
  zl = ar[0]; zs1 = ar[1]; zs2 = ar[2];
  C13 = C_gl_tomo(l,zl, zs1);
  C24 = C_ks(l, zs2);
  C14 = C_gl_tomo(l, zl, zs2);
  C23 = C_ks(l, zs1);
  return (C13*C24+C14*C23)/((2.*l+1.)*fsky);
}

double func_for_cov_G_ks_shear(double l, int *ar){
  double fsky = survey.area*survey.area_conversion_factor/(4.*M_PI);
  double C13, C14, C23, C24;
  int z1,z2,z3;
  z1 = ar[0]; z2 = ar[1]; z3 = ar[2];
  C13 = C_ks(l,z2);
  C24 = C_shear_tomo(l, z1, z3);
  C14 = C_ks(l,z3);
  C23 = C_shear_tomo(l, z1, z2);
  double N24 = 0, N23 = 0;
  if (z1 == z3){
    N24= pow(survey.sigma_e,2.0)/(2.0*nsource(z1)*survey.n_gal_conversion_factor);
  }
  if (z1 == z2){
    N23 = pow(survey.sigma_e,2.0)/(2.0*nsource(z1)*survey.n_gal_conversion_factor);
  }
  return (C13*(C24+N24)+C14*(C23+N23))/((2.*l+1.)*fsky);
}

double func_for_cov_G_gk_gl(double l, int *ar){
  double fsky = survey.area*survey.area_conversion_factor/(4.*M_PI);
  double C13, C14, C23, C24;
  int zl1,zl2,zs;
  zl1 = ar[0]; zl2 = ar[1]; zs = ar[2];
  C13 = C_cl_tomo(l,zl1, zl2);
  C24 = C_ks(l, zs);
  C14 = C_gl_tomo(l, zl1, zs);
  C23 = C_gk(l, zl2);
  double N = 0;
  if (zl1==zl2) N = 1./(nlens(zl1)*survey.n_gal_conversion_factor);
  return ((C13+N)*C24+C14*C23)/((2.*l+1.)*fsky);
}

double func_for_cov_G_ks_gl(double l, int *ar){
  double fsky = survey.area*survey.area_conversion_factor/(4.*M_PI);
  double C13, C14, C23, C24;
  int zl,zs1,zs2;
  zs1 = ar[0]; zl = ar[1]; zs2 = ar[2];
  C13 = C_gk(l, zl);
  C24 = C_shear_tomo(l, zs1, zs2);
  C14 = C_ks(l, zs2);
  C23 = C_gl_tomo(l, zl, zs1);
  double N = 0;
  if (zs1==zs2) N = pow(survey.sigma_e,2.0)/(2.0*nsource(zs1)*survey.n_gal_conversion_factor);
  return (C13*(C24+N)+C14*C23)/((2.*l+1.)*fsky);
}

double func_for_cov_G_gk_cl(double l, int *ar){
  double C13, C14, C23, C24, N13=0, N14=0, N23=0, N24=0;
  int n1,n3,n4;
  n1 = ar[0]; n3 = ar[1]; n4 = ar[2];
  double fsky = survey.area*survey.area_conversion_factor/(4.*M_PI);
  C13 = C_cl_tomo(l,n1,n3);
  C24 = C_gk(l,n4);
  C23 = C_gk(l,n3);
  C14 = C_cl_tomo(l,n1,n4);
  if (n1==n3){
    N13 = 1./(nlens(n1)*survey.n_gal_conversion_factor);
  }
  if (n1==n4){
    N14 = 1./(nlens(n1)*survey.n_gal_conversion_factor);
  }
  return ((C13+N13)*C24+(C14+N14)*C23)/((2.*l+1.)*fsky);
}

double func_for_cov_G_ks_cl(double l, int *ar){
  double fsky = survey.area*survey.area_conversion_factor/(4.*M_PI);
  double C13, C14, C23, C24;
  int n1,n2,zs;
  zs = ar[0]; n1 = ar[1]; n2 = ar[2];
  C13 = C_gk(l, n1);
  C24 = C_gl_tomo(l, n2, zs);
  C23 = C_gl_tomo(l, n1, zs);
  C14 = C_gk(l, n2);
  return (C13*C24+C14*C23)/((2.*l+1.)*fsky);
}

double func_for_cov_G_gk(double l, int *ar){
  double C13, C14, C23, C24, N13=0, N14=0, N23=0, N24=0;
  double fsky = survey.area*survey.area_conversion_factor/(4.*M_PI);
  int n1,n3;
  n1 = ar[0]; n3 = ar[1];
  C13 = C_cl_tomo(l,n1,n3);
  C24 = C_kk(l);
  C14 = C_gk(l,n1);
  C23 = C_gk(l,n3);
  if (n1 == n3){
    N13 = 1./(nlens(n1)*survey.n_gal_conversion_factor);
  }
  N24=kappa_reconstruction_noise(l);
  return ((C13+N13)*(C24+N24)+C14*C23)/((2.*l+1.)*fsky);
}

double func_for_cov_G_ks_gk(double l, int *ar){
  double fsky = survey.area*survey.area_conversion_factor/(4.*M_PI);
  double C13, C14, C23, C24;
  int zl,zs;
  zs = ar[0]; zl = ar[1];
  C13 = C_gk(l,zl);
  C24 = C_ks(l, zs);
  C23 = C_gl_tomo(l, zl, zs);
  C14 = C_kk(l);
  double N14 = kappa_reconstruction_noise(l);
  return (C13*C24+(C14+N14)*C23)/((2.*l+1.)*fsky);
}

double func_for_cov_G_ks(double l, int *ar){
  double C13, C14, C23, C24, N13=0, N14=0, N23=0, N24=0;
  int z1,z3;
  z1 = ar[0]; z3 = ar[1];
  double fsky = survey.area*survey.area_conversion_factor/(4.*M_PI);
  C13 = C_shear_tomo(l,z1,z3);
  C24 = C_kk(l);
  C14 = C_ks(l,z1);
  C23 = C_ks(l,z3);
  if (z1 == z3){
    N13= pow(survey.sigma_e,2.0)/(2.0*nsource(z1)*survey.n_gal_conversion_factor);
  }
  N24=kappa_reconstruction_noise(l);
  // for(l=30.;l<3000.;l++){
  //     C13 = C_shear_tomo(l,z1,z3);
  // C24 = C_kk(l);
  // C14 = C_ks(l,z1);
  // C23 = C_ks(l,z3);
  // if (z1 == z3){
  //   N13= pow(survey.sigma_e,2.0)/(2.0*nsource(z1)*survey.n_gal_conversion_factor);
  // }
  // N24=kappa_reconstruction_noise(l);
  // printf("%d,%d, %le,%le,%le,%le ,%le,%le,%le,%le,\n",z1,z3,C13,C24,C14,C23,N13,N24, l,((C13+N13)*(C24+N24)+C14*C23)/((2.*l+1.)*fsky));}
  // exit(0);
  return ((C13+N13)*(C24+N24)+C14*C23)/((2.*l+1.)*fsky);
}

double func_for_cov_G_kk_shear(double l, int *ar){
  double C13, C14, C23, C24;
  int z1,z2;
  z1 = ar[0]; z2 = ar[1];
  C13 = C_ks(l,z1);
  C23 = C13;
  C24 = C_ks(l, z2);
  C14 = C24;
  return (C13*C24+C14*C23)/((2.*l+1.)*cmb.fsky);
}

double func_for_cov_G_kk_gl(double l, int *ar){
  double C13, C14, C23, C24;
  int zl,zs;
  zl = ar[0]; zs = ar[1];
  C13 = C_gk(l, zl);
  C14 = C13;
  C24 = C_ks(l, zs);
  C23 = C24;
  return (C13*C24+C14*C23)/((2.*l+1.)*cmb.fsky);
}

double func_for_cov_G_kk_cl(double l, int *ar){
  double C13, C14, C23, C24;
  int n1, n2;
  n1=ar[0]; n2=ar[1];
  C13 = C_gk(l,n1);
  C24 = C_gk(l,n2);
  C14 = C13;
  C23 = C24;
  return (C13*C24+C14*C23)/((2.*l+1.)*cmb.fsky);
}

double func_for_cov_G_kk_gk(double l, int *ar){
  double C13, C14, C23, C24;
  int zl;
  zl=ar[0];
  C13 = C_gk(l,zl);
  C14 = C13;
  C24 = C_kk(l);
  C23 = C24;
  double N = kappa_reconstruction_noise(l);
  return (C13*(C24+N)+C14*(C23+N))/((2.*l+1.)*cmb.fsky);
}

double func_for_cov_G_kk_ks(double l, int *ar){
  double C13, C24, C14, C23;
  int zs;
  zs = ar[0];
  C13 = C_kk(l);
  C23 = C13;
  C14 = C_ks(l, zs);
  C24 = C14;
  double N = kappa_reconstruction_noise(l);
  double ans = ( (C13+N)*C24 + C14*(C23+N) )/( (2.*l+1.)*cmb.fsky );
  if(isfinite(ans)){return ans;}
  else{
    printf("cov_G_kk_ks non finite!\n");
    printf("C13=%e C23=%e C14=%e C24=%e N=%e fsky=%e\n",C13, C23, C14, C24, N, cmb.fsky);
	return 0.0;
  }
  //return ( (C13+N)*C24 + C14*(C23+N) )/( (2.*l+1.)*cmb.fsky );
}


double func_for_cov_G_kk(double l, int *ar){
   double C, N;
   C = C_kk(l);
   N = kappa_reconstruction_noise(l);
   return 2.*pow((C+N), 2) / ((2.*l+1.)*cmb.fsky);
}
