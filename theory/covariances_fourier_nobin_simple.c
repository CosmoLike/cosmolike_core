/* 
// non-binned fourier space-only 10x2pt cov
// in a unified way
// -- Xiao Fang */


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
    *Cij = C_cl_tomo_nointerp(l,ni,nj);
    if(ni==nj) {*Nij = 1./(nlens(ni)*survey.n_gal_conversion_factor);}
    return 0;
  }
  if(is_lsi*ls_lsj==-1){ // case: ls, need to determine lens and src bins
    if(is_lsi==1){*Cij = C_gl_tomo_nointerp(l,ni,nj);}
    else{*Cij = C_gl_tomo_nointerp(l,nj,ni);}
    return 0;
  }
  if(is_lsi==-1 && is_lsj==-1){ // case: ss
    *Cij = C_shear_tomo_nointerp(l,ni,nj);
    if(ni==nj) {*Nij= pow(survey.sigma_e,2.0)/(2.0*nsource(ni)*survey.n_gal_conversion_factor);}
    return 0;
  }
  if(ni==-1 || nj==-1){ // at least one kcmb field
    if(is_lsi==1){*Cij = C_gk_nointerp(l,ni);return 0;}
    else if(is_lsj==1){*Cij = C_gk_nointerp(l,nj);return 0;}
    else if(is_lsi==-1){*Cij = C_ks_nointerp(l,ni);return 0;}
    else if(is_lsj==-1){*Cij = C_ks_nointerp(l,nj);return 0;}
    else if(ni==-2 || nj==-2){*Cij = C_ky_nointerp(l);return 0;}
    else if(ni==-1 && nj==-1){*Cij = C_kk_nointerp(l); *Nij = kappa_reconstruction_noise(l); return 0;}
  }
  if(ni==-2 || nj==-2){ // at least one y field, but ky not included
    if(is_lsi==1){*Cij = C_gy_nointerp(l,ni);return 0;}
    else if(is_lsj==1){*Cij = C_gy_nointerp(l,nj);return 0;}
    else if(is_lsi==-1){*Cij = C_sy_nointerp(l,ni);return 0;}
    else if(is_lsj==-1){*Cij = C_sy_nointerp(l,nj);return 0;}
    else if(ni==-2 && nj==-2){*Cij = C_yy_nointerp(l); *Nij = y_reconstruction_noise(l);return 0;}
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


double tab_cov_NG_AB_CD(char ABCD[2][4], double l1, double l2, int z_ar[4], int is_ls[4]){
  return 0;
}

double bgal_a(double a, double nz){
  return gbias.b1_function(1./a-1.,(int)nz);
}
double inner_project_tri_cov_AB_CD(double a,void *params)
{
  double k[2],fK,weights,res = 0.;
  double *ar = (double *) params;
  int i,j;
  fK = f_K(chi(a));
  k[0] = (ar[0]+0.5)/fK;
  k[1] = (ar[1]+0.5)/fK;
  
  weights = dchi_da(a);
  for(i=2;i<6;i++){
    if(ar[i]==-1.){
      weights *= W_k(a,fK);
    }else if(ar[i]==-2.){
      weights *= W_y(a,fK);
      printf("NG cov for y is not yet supported!\n"); exit(1);
    }else if(is_ls[4+i]==1.){
      weights *= W_gal(a,ar[i]);
    }else{
      weights *= W_kappa(a,fK,ar[i]);
    }
  }

  double ssc[2];
  double sig_b, fsky1, fsky2, fsky_larger;
  fsky1 = ar[10]; fsky2 = ar[11];

  for(i=0;i<2;i++){
    ssc[i]=delP_SSC(k[i],a); // no galaxy density field
    for(j=0;j<2;j++){// if galaxy density field, rescale wrt mean, Eq(A12) in CosmoLike paper
      if(ar[6+2*i+j]==1){ssc[i] -= bgal_a(a,ar[2+2*i+j])*Pdelta(k[i],a);}
    }
  }

  if(fsky1==fsky2){
    sig_b = survey_variance(a,fsky1);
  }else{
    sig_b = cross_survey_variance(a,fsky1,fsky2);
  }
  fsky_larger = (fsky1>fsky2 ? fsky1 : fsky2);

  if(weights >0.){
    if (covparams.cng) {res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(fsky_larger*41253.0*survey.area_conversion_factor);}
    // printf("res cNG:%lg , ", res);
    res += ssc[0]*ssc[1]*sig_b*pow(fK,-4.); //SSC
  }
  // printf("res +SSC, weights:%lg, %lg\n", res, weights);
  res *= weights;
  return res;
}

double cov_NG_AB_CD(char ABCD[2][4], double l1,double l2, int z_ar[4], int is_ls[4]){
  double amin[2], amax[2], a1,a2,array[12];
  int zmin, zmax, zlen, zs;
  int i;
  int fsky[2];

  fsky_gal = survey.area/41253.0;
  for(i=0;i<2;i++){
    if(strcmp(ABCD[i], "ss")==0){
      zmin = (z_ar[2*i] < z_ar[2*i+1] ? z_ar[2*i] : z_ar[2*i+1]);
      amin[i] = amin_source(zmin); amax[i] = amax_source(zmin); fsky[i]=fsky_gal;
    }else if(strcmp(ABCD[i], "ls")==0 || strcmp(ABCD[i], "lk")==0 || strcmp(ABCD[i], "ly")==0){
      if(is_ls[2*i]==1){zlen = z_ar[2*i];}
      else{zlen = z_ar[2*i+1];}
      amin[i] = amin_lens(zlen); a2 = amax_lens(zlen); fsky[i]=fsky_gal;
    }else if(strcmp(ABCD[i], "ll")==0){
      zmin = (z_ar[2*i] < z_ar[2*i+1] ? z_ar[2*i] : z_ar[2*i+1]);
      zmax = (z_ar[2*i] > z_ar[2*i+1] ? z_ar[2*i] : z_ar[2*i+1]);
      amin[i] = amin_lens(zmax); amax[i] = amax_lens(zmin); fsky[i]=fsky_gal;
    }else if(strcmp(ABCD[i], "ks")==0){
      if(is_ls[2*i]==-1){zs = z_ar[2*i];}
      else{zs = z_ar[2*i+1];}
      amin[i] = amin_source(zs); amax[i] = amax_source(zs); fsky[i]=fsky_gal;
    }else if(strcmp(ABCD[i], "kk")==0 || strcmp(ABCD[i], "ky")==0 || strcmp(ABCD[i], "yy")==0){
      amin[i] = limits.a_min*(1.+1.e-5); amax[i] = 1.-1.e-5; fsky[i]=cmb.fsky;
    }
  }
  a1 = (amin[0]>amin[1] ? amin[0] : amin[1]);
  a2 = (amax[0]<amax[1] ? amax[0] : amax[1]);

  array[0] = l1;
  array[1] = l2;
  for(i=0;i<4;i++){
    array[2+i] = (double) z_ar[i];
    array[6+i] = (double) is_ls[i];
  }
  for(i=0;i<2;i++){array[10+i] = fsky[i];}
  return int_gsl_integrate_low_precision(inner_project_tri_cov_AB_CD,(void*)array,a1,a2,NULL,1000);
}