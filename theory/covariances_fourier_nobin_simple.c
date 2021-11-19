/* 
// but also provide non-binned fourier space-only 10x2pt cov
// in a unified way
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

// 3x2pt Gaussian cov pure noise term on the diagonal, without masking effect
// void pure_noise_xipm_xipm(int *z_ar, double *theta, double *dtheta, double *N);
// void pure_noise_gl_gl(int *z_ar, double *theta, double *dtheta, double *N);
// void pure_noise_cl_cl(int *z_ar, double *theta, double *dtheta, double *N);

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

  } else if(strcmp(cov_type, "yy_shear")==0) {
    func_for_cov_NG = &cov_NG_yy_shear_bridge; // ADD MORE
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



int func_for_Cl_Nl(double *Cij, double *Nij, double l, int ni, int nj, int is_lsi, int is_lsj){
  if(is_lsi==1 && is_lsj==1){// case: ll
    *Cij = C_cl_tomo(l,ni,nj);
    if(ni==nj) {*Nij = 1./(nlens(ni)*survey.n_gal_conversion_factor);}
    return 0;
  }
  if(is_lsi*ls_lsj==-1){ // case: ls, need to determine lens and src bins
    if(is_lsi==1){*Cij = C_gl_tomo(l,ni,nj);}
    else{*Cij = C_gl_tomo(l,nj,ni);}
    return 0;
  }
  if(is_lsi==-1 && is_lsj==-1){ // case: ss
    *Cij = C_shear_tomo(l,ni,nj);
    if(ni==nj) {*Nij= pow(survey.sigma_e,2.0)/(2.0*nsource(ni)*survey.n_gal_conversion_factor);}
    return 0;
  }
  if(ni==-1 || nj==-1){ // at least one kcmb field
    if(is_lsi==1){*Cij = C_gk(l,ni);return 0;}
    else if(is_lsj==1){*Cij = C_gk(l,nj);return 0;}
    else if(is_lsi==-1){*Cij = C_ks(l,ni);return 0;}
    else if(is_lsj==-1){*Cij = C_ks(l,nj);return 0;}
    else if(ni==-2 || nj==-2){*Cij = C_ky(l);return 0;}
    else if(ni==-1 && nj==-1){*Cij = C_kk(l); *Nij = kappa_reconstruction_noise(l); return 0;}
  }
  if(ni==-2 || nj==-2){ // at least one y field, but ky not included
    if(is_lsi==1){*Cij = C_gy(l,ni);return 0;}
    else if(is_lsj==1){*Cij = C_gy(l,nj);return 0;}
    else if(is_lsi==-1){*Cij = C_sy(l,ni);return 0;}
    else if(is_lsj==-1){*Cij = C_sy(l,nj);return 0;}
    else if(ni==-2 && nj==-2){*Cij = C_yy(l); *Nij = y_reconstruction_noise(l);return 0;}
  }
  printf("func_for_Cl_Nl: field combination not supported!\n");
  exit(1);
}

double cov_G_AB_CD(char ABCD[2][4], double l, double delta_l, int z_ar[4], int is_ls[4]){
  int i,j;
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  int is_ls1, is_ls2, is_ls3, is_ls4;
  double fsky = survey.area*survey.area_conversion_factor/(4.*M_PI);

  n1 = z_ar[0]; n2 = z_ar[1]; n3 = z_ar[2]; n4 = z_ar[3];
  is_ls1 = is_ls[0]; is_ls2 = is_ls[1]; is_ls3 = is_ls[2]; is_ls4 = is_ls[3];

  func_for_Cl_Nl(&C13, &N13, l,n1,n3,is_ls1,is_ls3);
  func_for_Cl_Nl(&C14, &N14, l,n1,n4,is_ls1,is_ls4);
  func_for_Cl_Nl(&C23, &N23, l,n2,n3,is_ls2,is_ls3);
  func_for_Cl_Nl(&C24, &N24, l,n2,n4,is_ls2,is_ls4);
  
  char cmbprobes[3][4] = {"kk", "ky", "yy"};
  for(i=0;i<2;i++){
    for(j=0;j<3;j++){
      if(strcmp(ABCD[i],cmbprobes[j])==0){
        fsky = (fsky<cmb.fsky) ? cmb.fsky : fsky; // if one of the probes is from CMB, then use the cmb area (if it's larger)
        i=2; j=3; break; // break out double for-loop
      }
    }
  }

  return ((C13+N13)*(C24+N24)+(C14+N14)*(C23+N23))/(fsky*(2.*l+1.)*delta_l);
}
