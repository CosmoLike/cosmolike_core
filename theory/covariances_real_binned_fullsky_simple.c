// fourier space 3x2pt Gaussian cov, without pure noise term
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
// double fsky_planck = 0.673;
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

  double N[like.Ntheta];
  for(i=0;i<like.Ntheta;i++) {N[i] = 0.;}

  double (*func_for_cov_G)(double, int*);
  double (*func_bin_cov_NG)(double, double, int*);
  double (*func_P1)(int, int);
  double (*func_P2)(int, int);

  if(strcmp(realcov_type, "xi+_xi+")==0) {
    func_for_cov_G  = &func_for_cov_G_shear;
    func_bin_cov_NG = &bin_cov_NG_shear_shear;
    func_P1 = &Glplus_tab;
    func_P2 = &Glplus_tab;
    pure_noise_xipm_xipm(z_ar, theta, dtheta, N); //C+- doesn't have the diagonal shot noise term
  } else if(strcmp(realcov_type, "xi-_xi-")==0) {
    func_for_cov_G  = &func_for_cov_G_shear;
    func_bin_cov_NG = &bin_cov_NG_shear_shear;
    func_P1 = &Glminus_tab;
    func_P2 = &Glminus_tab;
    pure_noise_xipm_xipm(z_ar, theta, dtheta, N); //C+- doesn't have the diagonal shot noise term
  } else if(strcmp(realcov_type, "xi+_xi-")==0) {
    func_for_cov_G  = &func_for_cov_G_shear;
    func_bin_cov_NG = &bin_cov_NG_shear_shear;
    func_P1 = &Glplus_tab;
    func_P2 = &Glminus_tab;
  } else if(strcmp(realcov_type, "xi-_xi+")==0) {
    func_for_cov_G  = &func_for_cov_G_shear;
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
    func_for_cov_G  = &func_for_cov_G_gl;
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
    func_for_cov_G  = &func_for_cov_G_cl;
    func_bin_cov_NG = &bin_cov_NG_cl_cl;
    func_P1 = &Pl_tab;
    func_P2 = &Pl_tab;
    pure_noise_cl_cl(z_ar, theta, dtheta, N);
  } else if(strcmp(realcov_type, "gk_xi+")==0) {
    func_for_cov_G  = &func_for_cov_G_gk_shear;
    func_bin_cov_NG = &bin_cov_NG_gk_shear;
    func_P1 = &Pl_tab;
    func_P2 = &Glplus_tab;
  } else if(strcmp(realcov_type, "gk_xi-")==0) {
    func_for_cov_G  = &func_for_cov_G_gk_shear;
    func_bin_cov_NG = &bin_cov_NG_gk_shear;
    func_P1 = &Pl_tab;
    func_P2 = &Glminus_tab;
  } else if(strcmp(realcov_type, "ks_xi+")==0) {
    func_for_cov_G  = &func_for_cov_G_ks_shear;
    func_bin_cov_NG = &bin_cov_NG_ks_shear;
    func_P1 = &Pl2_tab;
    func_P2 = &Glplus_tab;
  } else if(strcmp(realcov_type, "ks_xi-")==0) {
    func_for_cov_G  = &func_for_cov_G_ks_shear;
    func_bin_cov_NG = &bin_cov_NG_ks_shear;
    func_P1 = &Pl2_tab;
    func_P2 = &Glminus_tab;
  } else if(strcmp(realcov_type, "gk_gl")==0) {
    func_for_cov_G  = &func_for_cov_G_gk_gl;
    func_bin_cov_NG = &bin_cov_NG_gk_gl;
    func_P1 = &Pl_tab;
    func_P2 = &Pl2_tab;
  } else if(strcmp(realcov_type, "ks_gl")==0) {
    func_for_cov_G  = &func_for_cov_G_ks_gl;
    func_bin_cov_NG = &bin_cov_NG_ks_gl;
    func_P1 = &Pl2_tab;
    func_P2 = &Pl2_tab;
  } else if(strcmp(realcov_type, "gk_cl")==0) {
    func_for_cov_G  = &func_for_cov_G_gk_cl;
    func_bin_cov_NG = &bin_cov_NG_gk_cl;
    func_P1 = &Pl_tab;
    func_P2 = &Pl_tab;
  } else if(strcmp(realcov_type, "ks_cl")==0) {
    func_for_cov_G  = &func_for_cov_G_ks_cl;
    func_bin_cov_NG = &bin_cov_NG_ks_cl;
    func_P1 = &Pl2_tab;
    func_P2 = &Pl_tab;
  } else if(strcmp(realcov_type, "gk_gk")==0) {
    func_for_cov_G  = &func_for_cov_G_gk;
    func_bin_cov_NG = &bin_cov_NG_gk_gk;
    func_P1 = &Pl_tab;
    func_P2 = &Pl_tab;
  } else if(strcmp(realcov_type, "ks_gk")==0) {
    func_for_cov_G  = &func_for_cov_G_ks_gk;
    func_bin_cov_NG = &bin_cov_NG_ks_gk;
    func_P1 = &Pl2_tab;
    func_P2 = &Pl_tab;
  } else if(strcmp(realcov_type, "ks_ks")==0) {
    func_for_cov_G  = &func_for_cov_G_ks;
    func_bin_cov_NG = &bin_cov_NG_ks_ks;
    func_P1 = &Pl2_tab;
    func_P2 = &Pl2_tab;

  } else if(strcmp(realcov_type, "kk_xi+")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_shear;
    func_bin_cov_NG = &bin_cov_NG_kk_shear;
    func_P1 = &Pl_tab;
    func_P2 = &Glplus_tab;
  } else if(strcmp(realcov_type, "kk_xi-")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_shear;
    func_bin_cov_NG = &bin_cov_NG_kk_shear;
    func_P1 = &Pl_tab;
    func_P2 = &Glminus_tab;
  } else if(strcmp(realcov_type, "kk_gl")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_gl;
    func_bin_cov_NG = &bin_cov_NG_kk_gl;
    func_P1 = &Pl_tab;
    func_P2 = &Pl2_tab;
  } else if(strcmp(realcov_type, "kk_cl")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_cl;
    func_bin_cov_NG = &bin_cov_NG_kk_cl;
    func_P1 = &Pl_tab;
    func_P2 = &Pl_tab;
  } else if(strcmp(realcov_type, "kk_gk")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_gk;
    func_bin_cov_NG = &bin_cov_NG_kk_gk;
    func_P1 = &Pl_tab;
    func_P2 = &Pl_tab;
  } else if(strcmp(realcov_type, "kk_ks")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_ks;
    func_bin_cov_NG = &bin_cov_NG_kk_ks;
    func_P1 = &Pl_tab;
    func_P2 = &Pl2_tab;
  } else if(strcmp(realcov_type, "kk_kk")==0) {
    func_for_cov_G  = &func_for_cov_G_kk;
    func_bin_cov_NG = &bin_cov_NG_kk_kk;
    func_P1 = &Pl_tab;
    func_P2 = &Pl_tab;

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

  for (l1 = 0; l1 < LMAX; l1++){
    l1_double = (double)l1;
    // printf("l1,%d\n", l1);
    cov_g_l = func_for_cov_G(l1_double, z_ar);
    for(i=0; i<like.Ntheta ; i++){
      covGl1 = cov_g_l * func_P1(i,l1);
      // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
      for(j=0; j<like.Ntheta ; j++){
        cov[i][j] += covGl1 * func_P2(j,l1);
      }
    }

    if(FLAG_NG){
      for (l2 = 0; l2 < LMAX; l2++){
        tri = func_bin_cov_NG(l1_double,(double)l2,z_ar);
        for(i=0; i<like.Ntheta ; i++){
          triP = tri * func_P1(i,l1);
          // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
          for(j=0; j<like.Ntheta ; j++){
            covNG[i][j] += triP * func_P2(j,l2);
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

  // double N[like.Ntheta];
  // for(i=0;i<like.Ntheta;i++) {N[i] = 0.;}

  double (*func_for_cov_G)(double, int*);
  double (*func_bin_cov_NG)(double, double, int*);
  double (*func_P2)(int, int);

  if(strcmp(mixcov_type, "kk_xi+")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_shear;
    func_bin_cov_NG = &bin_cov_NG_kk_shear;
    func_P2 = &Glplus_tab;
  } else if(strcmp(mixcov_type, "kk_xi-")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_shear;
    func_bin_cov_NG = &bin_cov_NG_kk_shear;
    func_P2 = &Glminus_tab;
  } else if(strcmp(mixcov_type, "kk_gl")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_gl;
    func_bin_cov_NG = &bin_cov_NG_kk_gl;
    func_P2 = &Pl2_tab;
  } else if(strcmp(mixcov_type, "kk_cl")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_cl;
    func_bin_cov_NG = &bin_cov_NG_kk_cl;
    func_P2 = &Pl_tab;
  } else if(strcmp(mixcov_type, "kk_gk")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_gk;
    func_bin_cov_NG = &bin_cov_NG_kk_gk;
    func_P2 = &Pl_tab;
  } else if(strcmp(mixcov_type, "kk_ks")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_ks;
    func_bin_cov_NG = &bin_cov_NG_kk_ks;
    func_P2 = &Pl2_tab;
  } else {
    printf("cov_mix_binned_fullsky: mixcov_type \"%s\" not defined!\n", mixcov_type); exit(1);
  }

  double l1_double,tri;
  int l1,l2;
  for(i=0; i<like.Ncl ; i++){
    for(j=0; j<like.Ntheta ; j++){
      cov[i][j] = 0.;
      covNG[i][j] = 0.;
    }
  }

  double triP;
  double covGl1;

  double cov_g_l[LMAX];
    // printf("lmin,lmax, %d, %d\n", (int)floor(like.lmin),(int)ceil(like.lmax));
  for(l1=(int)floor(like.lmin); l1<(int)ceil(like.lmax); l1++){
    cov_g_l[l1] = func_for_cov_G((double)l1, z_ar);
    // printf("cov_g_l: %d, %le\n", l1,cov_g_l[l1]);
  }

  for(i=0; i<like.Ncl; i++){
    for(j=0; j<like.Ntheta ; j++){
      for(l1=(int)ceil(ell[i]); l1<ell[i+1]; l1++){
        cov[i][j] += cov_g_l[l1] * func_P2(j,l1) / ((int)floor(ell[i+1])-(int)ceil(ell[i])+1); // rewrite as window function
        if(FLAG_NG){
          l1_double = (double)l1;
          for (l2 = 0; l2 < LMAX; l2++){
            tri = func_bin_cov_NG(l1_double,(double)l2,z_ar);
            covNG[i][j] += tri * func_P2(j,l2);
          }
        }
      }
    }
  }
}


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
// mixed cov 
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


/********** Functions for differnt covariances ***************************************/
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

  return (C13*C24+C13*N24+N13*C24 + C14*C23+C14*N23+N14*C23)*4.*M_PI/(survey.area*survey.area_conversion_factor* (2.*l+1.));
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
  
  return (C13*C24+C13*N24+N13*C24 + C14*C23+C14*N23+N14*C23)*4.*M_PI/(survey.area*survey.area_conversion_factor* (2.*l+1.));
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
  return (C13*C24+C14*C23)/((2.*l+1.)*fsky_planck);
}

double func_for_cov_G_kk_gl(double l, int *ar){
  double C13, C14, C23, C24;
  int zl,zs;
  zl = ar[0]; zs = ar[1];
  C13 = C_gk(l, zl);
  C14 = C13;
  C24 = C_ks(l, zs);
  C23 = C24;
  return (C13*C24+C14*C23)/((2.*l+1.)*fsky_planck);
}

double func_for_cov_G_kk_cl(double l, int *ar){
  double C13, C14, C23, C24;
  int n1, n2;
  n1=ar[0]; n2=ar[1];
  C13 = C_gk(l,n1);
  C24 = C_gk(l,n2);
  C14 = C13;
  C23 = C24;
  return (C13*C24+C14*C23)/((2.*l+1.)*fsky_planck);
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
  return (C13*(C24+N)+C14*(C23+N))/((2.*l+1.)*fsky_planck);
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
  return ((C13+N)*C24+C14*(C23+N))/((2.*l+1.)*fsky_planck);
}


double func_for_cov_G_kk(double l, int *ar){
   double C, N;
   C = C_kk(l);
   N = kappa_reconstruction_noise(l);
   return 2.*pow((C+N), 2) / ((2.*l+1.)*fsky_planck);
}
