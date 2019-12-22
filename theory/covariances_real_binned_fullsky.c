/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
// Note the look-up tables for power spectrum covariances are recomputed if one the redshift bins changes
//      Loop over theta1, theta2 first, before computing the next combination of redshift bins
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
double w_mask(double theta_min);// angular correlation funtion of survey mask, computed as sum P_l(cos(like.theta[nt])*C_mask(l)
double cov_G_shear_binned(double theta1_min, double theta1_max,double theta2_min, double theta2_max, int z1,int z2,int z3,int z4,int pm1,int pm2); //Version of Gaussian cov calculation for wide bins
double cov_NG_shear_binned(double theta1_min, double theta1_max,double theta2_min, double theta2_max, int z1,int z2,int z3,int z4,int pm1,int pm2);

double cov_NG_cl_cl_real_binned(double theta1_min, double theta1_max,double theta2_min, double theta2_max, int z1,int z2, int z3, int z4);
double cov_G_cl_cl_real_binned(double theta1_min, double theta1_max,double theta2_min, double theta2_max, int z1,int z2, int z3, int z4);

double cov_NG_gl_gl_real_binned(double theta1_min, double theta1_max,double theta2_min, double theta2_max, int z1,int z2, int z3, int z4);
double cov_G_gl_gl_real_binned(double theta1_min, double theta1_max,double theta2_min, double theta2_max, int z1l,int z1s, int z2l, int z2s);

double cov_NG_cl_gl_real_binned(double theta1_min, double theta1_max,double theta2_min, double theta2_max, int z1,int z2, int zl, int zs);//z1,z2 clustering bins; zl,zs g-g lensing bins
double cov_G_cl_gl_real_binned(double theta1_min, double theta1_max,double theta2_min, double theta2_max, int z1,int z2, int zl, int zs);

double cov_NG_gl_shear_real_binned(double theta1_min, double theta1_max,double theta2_min, double theta2_max, int zl,int zs, int z3, int z4, int pm);
double cov_G_gl_shear_real_binned(double theta1_min, double theta1_max,double theta2_min, double theta2_max,int zl,int zs, int z3, int z4,int pm);

double cov_NG_cl_shear_real_binned(double theta1_min, double theta1_max,double theta2_min, double theta2_max, int z1,int z2, int z3, int z4, int pm); //z1,z2 clustering bins; z3,z4 shear bins
double cov_G_cl_shear_real_binned(double theta1_min, double theta1_max,double theta2_min, double theta2_max, int z1,int z2, int z3, int z4, int pm);

double func_for_cov_G_shear(double l, int *ar);
double func_for_cov_G_cl(double l, int *ar);
double func_for_cov_G_gl(double l, int *ar);
double func_for_cov_G_cl_shear(double l, int *ar);
double func_for_cov_G_cl_gl(double l, int *ar);
double func_for_cov_G_gl_shear(double l, int *ar);

double filter_cov_fourier(double l1, double l2, double lmax, double lpivot) {
  long i,j;
  double W;
  double l_interval = lmax - lpivot;
  if(l_interval<=0) {return 1.;}

  if(l1<=lpivot){
    W = 1.;
  }
  else if(l1>=lmax){
    W = 0.;
  }
  else{
    W = (lmax - l1) / l_interval - 1./(2.*M_PI) * sin(2.*(lmax - l1)*M_PI/l_interval);
  }

  if(l2<=lpivot){
    W *= 1.;
  }
  else if(l2>=lmax){
    W *= 0.;
  }
  else{
    W *= (lmax - l2) / l_interval - 1./(2.*M_PI) * sin(2.*(lmax - l2)*M_PI/l_interval);
  }
  return W;
}


/************************* covariance routines for angular correlation functions *********************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
// Note the look-up tables for power spectrum covariances are recomputed if one the redshift bins changes
//      Loop over theta1, theta2 first, before computing the next combination of redshift bins
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
double J0_binned(double l, double tmin, double tmax){
  return 2./(l*tmax*tmax-l*tmin*tmin)*(tmax*gsl_sf_bessel_J1(l*tmax)-tmin*gsl_sf_bessel_J1(l*tmin));
}
double xJ2 (double x, void *params){
  double *ar = (double *) params;
  return 2.*x*gsl_sf_bessel_Jn(2,x*ar[0]);
}
double J2_binned(double l, double tmin, double tmax){
  if (tmin*l < 1.0){
    double array[1] ={l};
    return int_gsl_integrate_high_precision(xJ2,(void*)array,tmin,tmax,NULL,1000)/(tmax*tmax-tmin*tmin);
  }
  return 1./(l*l*tmax*tmax-l*l*tmin*tmin)*(4.*gsl_sf_bessel_J0(l*tmin)-4.*gsl_sf_bessel_J0(l*tmax)+2.*l*tmin*gsl_sf_bessel_J1(l*tmin)-2.*l*tmax*gsl_sf_bessel_J1(l*tmax));
}

double xJ4 (double x, void *params){
  double *ar = (double *) params;
  return 2.*x*gsl_sf_bessel_Jn(4,x*ar[0]);
}
double J4_binned(double l, double tmin, double tmax){
  if (tmin*l < 1.0){
  	double array[1] ={l};
  	return int_gsl_integrate_high_precision(xJ4,(void*)array,tmin,tmax,NULL,1000)/(tmax*tmax-tmin*tmin);
  }
  double j1_max = gsl_sf_bessel_J1(l*tmax);
  double j1_min = gsl_sf_bessel_J1(l*tmin);

  return 2./(l*l*l*(tmax*tmax-tmin*tmin))*(((8.-l*l*tmin*tmin)*j1_min/tmin+(l*l*tmax*tmax-8.)*j1_max/tmax)-8.*l*(gsl_sf_bessel_Jn(2,l*tmax)-gsl_sf_bessel_Jn(2,l*tmin)));
}

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
      double **Pl, *xmin, *xmax, *Pmin, *Pmax;
      Pl =create_double_matrix(0, like.Ntheta-1, 0, lbins);
      xmin= create_double_vector(0, like.Ntheta-1);
      xmax= create_double_vector(0, like.Ntheta-1);
      double logdt=(log(like.vtmax)-log(like.vtmin))/like.Ntheta;
      Pmin= create_double_vector(0, lbins);
      Pmax= create_double_vector(0, lbins);
      for(i=0; i<like.Ntheta ; i++){
        xmin[i]=cos(exp(log(like.vtmin)+(i+0.0)*logdt));
        xmax[i]=cos(exp(log(like.vtmin)+(i+1.0)*logdt));
      }
      for (i = 0; i < like.Ntheta; i++){
        gsl_sf_legendre_Pl_array(lbins, xmin[i],Pmin);
        gsl_sf_legendre_Pl_array(lbins, xmax[i],Pmax);
        Pl[i][0] = 1.0/(4.*M_PI);
        for (int l = 1; l < lbins; l ++){
          Pl[i][l] = 1./(4.*M_PI)*(Pmin[l+1]-Pmax[l+1]-Pmin[l-1]+Pmax[l-1])/(xmin[i]-xmax[i]);
        }
      }
      free_double_vector(xmin,0,like.Ntheta-1);
      free_double_vector(xmax,0,like.Ntheta-1);
      free_double_vector(Pmin,0,lbins);
      free_double_vector(Pmax,0,lbins);
      for (i = 0; i < NTHETA; i++){
        w_vec[i] =0.;
        for (l = 0; l < lbins; l++){
          w_vec[i]+=Cl[l]*Pl[i][l];
        }
        printf("w_mask[%d] = %e\n",i, w_vec[i]);
      }
      free_double_vector(Cl,0,lbins-1);
      free_double_matrix(Pl, 0, like.Ntheta-1, 0, lbins);
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

double C_gl_tomo_all(double l, int ni, int nj)  //slower version of G-G lensing power spectrum, lens bin ni, source bin nj - tabulated for all lens-source combinations without overlap criterion
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  double f1 = 0.;
  
  if (ni < 0 || ni >= tomo.clustering_Nbin ||nj < 0 || nj >= tomo.shear_Nbin){
    printf("C_gl_tomo_all(l,%d,%d) outside tomo.X_Nbin range\nEXIT\n",ni,nj); exit(1);
  }
  
  if (recompute_ggl(C,G,N,ni)){
    if (table==0){
      if (tomo.shear_Nbin> 10){printf("tomo.shear_Nbin too large for look-up table in C_gl_tomo_all\nEXIT\n");exit(1);}
      table   = create_double_matrix(0, tomo.clustering_Nbin*10-1, 0, Ntable.N_ell-1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell);
    }
    int i,j,k;
    double llog;
    
    for (k=0; k<tomo.clustering_Nbin; k++) {
      for (j=0; j<tomo.shear_Nbin; j++) {
        llog = logsmin;
        for (i=0; i<Ntable.N_ell; i++, llog+=ds) {
          table[k*10+j][i]= log(C_gl_tomo_nointerp(exp(llog),k,j));
        }
      }
    }
    
    update_cosmopara(&C); update_nuisance(&N); update_galpara(&G);
  }
  f1 = exp(interpol(table[10*ni+nj], Ntable.N_ell, logsmin, logsmax, ds, log(l), 0.,0.));
  if (isnan(f1)){f1 = 0;}
  return f1;
  // return C_gl_tomo(l,ni,nj);
}
/**************** look-up tables for covariance  *********************/
double bin_cov_NG_shear_shear_tomo(double l1,double l2, int z1, int z2, int z3, int z4){
  static int Z1 = -42;
  static int Z2 = -42;
  static int Z3 = -42;
  static int Z4 = -42;
  static int Ntab = 20;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 || Z4!=z4 )
    {
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    logsmin = log(1.);
    logsmax = log(5.e+4);
    ds = (logsmax - logsmin)/(Ntab - 1.);
    llog1 = logsmin;
    for (i=0; i<Ntab; i++, llog1+=ds) {
      ll1 = exp(llog1);
      llog2 = logsmin;
      for (j=0; j<Ntab; j++, llog2+=ds) {
        ll2 = exp(llog2);
        table[i][j]=log(cov_NG_shear_shear_tomo(ll1,ll2,z1,z2,z3,z4));
      }
    }
    Z1=z1; Z2=z2; Z3=z3; Z4=z4;
    }
  res = 0.;
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = exp(interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0));}
    res *= filter_cov_fourier(llog1, llog2, logsmax, log(4.e4));
  return res;
}
double bin_cov_NG_gl_gl_tomo(double l1,double l2, int z1, int z2, int z3, int z4){
  static int Z1 = -42;
  static int Z2 = -42;
  static int Z3 = -42;
  static int Z4 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 || Z4!=z4 )
    {
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    logsmin = log(1.);
    logsmax = log(5.e+4);
    ds = (logsmax - logsmin)/(Ntab - 1.);
    llog1 = logsmin;
    for (i=0; i<Ntab; i++, llog1+=ds) {
      ll1 = exp(llog1);
      llog2 = logsmin;
      for (j=0; j<Ntab; j++, llog2+=ds) {
        ll2 = exp(llog2);
        table[i][j]= cov_NG_gl_gl_tomo(ll1,ll2,z1,z2,z3,z4);
      }
    }
    Z1=z1; Z2=z2; Z3=z3; Z4=z4;
    }
  res = 0.;
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
  return res;
}
double bin_cov_NG_cl_cl_tomo(double l1,double l2, int z1, int z2, int z3, int z4){
  static int Z1 = -42;
  static int Z2 = -42;
  static int Z3 = -42;
  static int Z4 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 || Z4!=z4 )
    {
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    logsmin = log(1.);
    logsmax = log(5.e+4);
    ds = (logsmax - logsmin)/(Ntab - 1.);
    llog1 = logsmin;
    for (i=0; i<Ntab; i++, llog1+=ds) {
      ll1 = exp(llog1);
      llog2 = logsmin;
      for (j=0; j<Ntab; j++, llog2+=ds) {
        ll2 = exp(llog2);
        table[i][j]=cov_NG_cl_cl_tomo(ll1,ll2,z1,z2,z3,z4);
      }
    }
    Z1=z1; Z2=z2; Z3=z3; Z4=z4;
  }
  res = 0.;
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
  return res;
}
double bin_cov_NG_cl_shear_tomo(double l1,double l2, int z1, int z2, int z3, int z4){
  static int Z1 = -42;
  static int Z2 = -42;
  static int Z3 = -42;
  static int Z4 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 || Z4!=z4 )
    {
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    logsmin = log(1.);
    logsmax = log(5.e+4);
    ds = (logsmax - logsmin)/(Ntab - 1.);
    llog1 = logsmin;
    for (i=0; i<Ntab; i++, llog1+=ds) {
      ll1 = exp(llog1);
      llog2 = logsmin;
      for (j=0; j<Ntab; j++, llog2+=ds) {
        ll2 = exp(llog2);
        table[i][j]=cov_NG_cl_shear_tomo(ll1,ll2,z1,z2,z3,z4);
        // printf("cov_NG_cl_shear_tomo(%lg,%lg),%lg\n",ll1,ll2, table[i][j]); 
      }
    }
    Z1=z1; Z2=z2; Z3=z3; Z4=z4;
    }
  res = 0.;
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
  return res;
}
double bin_cov_NG_cl_gl_tomo(double l1,double l2, int z1, int z2, int z3, int z4){
  static int Z1 = -42;
  static int Z2 = -42;
  static int Z3 = -42;
  static int Z4 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 || Z4!=z4 )
    {
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    logsmin = log(1.);
    logsmax = log(5.e+4);
    ds = (logsmax - logsmin)/(Ntab - 1.);
    llog1 = logsmin;
    for (i=0; i<Ntab; i++, llog1+=ds) {
      ll1 = exp(llog1);
      llog2 = logsmin;
      for (j=0; j<Ntab; j++, llog2+=ds) {
        ll2 = exp(llog2);
        table[i][j]=cov_NG_cl_gl_tomo(ll1,ll2,z1,z2,z3,z4);
      }
    }
    Z1=z1; Z2=z2; Z3=z3; Z4=z4;
    }
  res = 0.;
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
  return res;
}
double bin_cov_NG_gl_shear_tomo(double l1,double l2, int z1, int z2, int z3, int z4){
  static int Z1 = -42;
  static int Z2 = -42;
  static int Z3 = -42;
  static int Z4 = -42;
  static int Ntab = 40;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  int i,j;
  double res, llog1,llog2,ll1,ll2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 || Z4!=z4 )
    {
    if (table==0) { table = create_double_matrix(0, Ntab-1, 0, Ntab-1);}
    logsmin = log(1.);
    logsmax = log(5.e+4);
    ds = (logsmax - logsmin)/(Ntab - 1.);
    llog1 = logsmin;
    for (i=0; i<Ntab; i++, llog1+=ds) {
      ll1 = exp(llog1);
      llog2 = logsmin;
      for (j=0; j<Ntab; j++, llog2+=ds) {
        ll2 = exp(llog2);
        table[i][j]=cov_NG_gl_shear_tomo(ll1,ll2,z1,z2,z3,z4);
        // printf("cov_NG_gl_shear_tomo(%lg,%lg),%lg\n",ll1,ll2, cov_NG_gl_shear_tomo(ll1,ll2,z1,z2,z3,z4)); 
      }
    }
    Z1=z1; Z2=z2; Z3=z3; Z4=z4;
    }
  res = 0.;
  llog1=log(l1);
  llog2=log(l2);
  if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
    res = interpol2d(table, Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
  return res;
}


/////// full sky covs

void cov_shear_shear_real_binned_fullsky(double **cov, double **covNG, int z1,int z2,int z3,int z4,int pm1,int pm2, int FLAG_NG, double *theta, double *dtheta){
  
  int i,j;
  static int LMAX = 50000;
  static double **Glplus =0;
  static double **Glminus =0;
  if (Glplus ==0){
    Glplus =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
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
      double x = cos(like.theta[i]);
      gsl_sf_legendre_Pl_deriv_array(LMAX, xmin[i],Pmin,dPmin);
      gsl_sf_legendre_Pl_deriv_array(LMAX, xmax[i],Pmax,dPmax);
      // printf("xmax[%d]:%lg\n", i,xmax[i]);
      for (int l = 3; l < LMAX; l ++){
        /*double plm = gsl_sf_legendre_Plm(l,2,x);
        double plm_1 = gsl_sf_legendre_Plm(l-1,2,x);
        Glplus[i][l] = (2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))
        *(plm*((4-l+2.*x*(l-1))/(1-x*x)-l*(l+1)/2)
        +plm_1*(l-1,2,x)*(l+2)*(x-2)/(1-x*x));


        Glminus[i][l] = (2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))
        *(plm*(l,2,x)*((4-l-2.*x*(l-1))/(1-x*x)-l*(l+1)/2)
        +plm_1*(l-1,2,x)*(l+2)*(x+2)/(1-x*x));*/

        Glplus[i][l] =(2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(

        -l*(l-1.)/2*(l+2./(2*l+1)) * (Pmin[l-1]-Pmax[l-1])
        -l*(l-1.)*(2.-l)/2         * (xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
        +l*(l-1)/(2*l+1)           * (Pmin[l+1]-Pmax[l+1])

        +(4-l)   * (dPmin[l]-dPmax[l])
        +(l+2)   * (xmin[i]*dPmin[l-1] - xmax[i]*dPmax[l-1] - Pmin[l-1] + Pmax[l-1])

        +2*(l-1) * (xmin[i]*dPmin[l]   - xmax[i]*dPmax[l]   - Pmin[l] + Pmax[l])
        -2*(l+2) * (dPmin[l-1]-dPmax[l-1])

        )/(xmin[i]-xmax[i]);           

        Glminus[i][l] =(2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(

        -l*(l-1.)/2*(l+2./(2*l+1)) * (Pmin[l-1]-Pmax[l-1])
        -l*(l-1.)*(2.-l)/2         * (xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
        +l*(l-1)/(2*l+1)           * (Pmin[l+1]-Pmax[l+1])

        +(4-l)   * (dPmin[l]-dPmax[l])
        +(l+2)   * (xmin[i]*dPmin[l-1] - xmax[i]*dPmax[l-1] - Pmin[l-1] + Pmax[l-1])

        -2*(l-1) * (xmin[i]*dPmin[l]   - xmax[i]*dPmax[l]   - Pmin[l] + Pmax[l])
        +2*(l+2) * (dPmin[l-1]-dPmax[l-1])

        )/(xmin[i]-xmax[i]);
        
        // printf("Pmin[%d], Pmax[%d]:%lg, %lg\n", l,l, Pmin[l], Pmax[l]);
        // printf("Glplus[%d][%d],%lg\n", i,l, Glplus[i][l]);
      }
    }
    free_double_vector(xmin,0,like.Ntheta-1);
    free_double_vector(xmax,0,like.Ntheta-1);
    free_double_vector(Pmin,0,LMAX+1);
    free_double_vector(Pmax,0,LMAX+1);
    free_double_vector(dPmin,0,LMAX+1);
    free_double_vector(dPmax,0,LMAX+1);
  }

  double l1_double,tri;
  int l1,l2;
  for(i=0; i<like.Ntheta ; i++){
    for(j=0; j<like.Ntheta ; j++){
      cov[i][j] = 0.;
      covNG[i][j] = 0.;
    }
  }

  // printf("i,j:%d,%d\n", i,j);
  double triGl1;
  double covGl1, cov_g_l;
  int ar[4];
  ar[0] = z1; ar[1] = z2; ar[2] = z3; ar[3] = z4;

  double N[like.Ntheta];
  for(i=0;i<like.Ntheta;i++) {N[i] = 0.;}
  if (z1 ==z3 && z2 ==z4 && pm1 == pm2){ //&& pm1 == pm2 required as C+- doesn't have the diagonal shot noise term
    for(i=0;i<like.Ntheta;i++) {
      N[i] += pow(survey.sigma_e,4.0)/(M_PI*(2.*theta[i]+dtheta[i])*dtheta[i]*4.*nsource(z1)*nsource(z2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor);
    } //(number of galaxy pairs in the survey contributing to annulus of width Dtheta centered at theta1)^-1
  }
  if (z1 ==z4 && z2 ==z3 && pm1 == pm2){ //&& pm1 == pm2 required as C+- doesn't have the diagonal shot noise term
    for(i=0;i<like.Ntheta;i++) {
      N[i] += pow(survey.sigma_e,4.0)/(M_PI*(2.*theta[i]+dtheta[i])*dtheta[i]*4.*nsource(z1)*nsource(z2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor);
    }
  }
  for(i=0;i<like.Ntheta;i++) {if(N[i]) N[i] /= w_mask(theta[i]); printf("N[i]: %lg\n", N[i]);}

  if(pm1>0 && pm2>0){
    for (l1 = 3; l1 < LMAX; l1++){
      l1_double = (double)l1;
      cov_g_l = func_for_cov_G_shear(l1_double, ar);
      for(i=0; i<like.Ntheta ; i++){
        covGl1 = cov_g_l * Glplus[i][l1];
        // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
        for(j=0; j<like.Ntheta ; j++){
          cov[i][j] += covGl1 * Glplus[j][l1];
        }
      }
      // printf("l1,%d\n", l1);
      if(FLAG_NG){
        for (l2 = 3; l2 < LMAX; l2++){
          tri = bin_cov_NG_shear_shear_tomo(l1_double,(double)l2,z1,z2,z3,z4);
          for(i=0; i<like.Ntheta ; i++){
            triGl1 = tri * Glplus[i][l1];
            // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
            for(j=0; j<like.Ntheta ; j++){
              covNG[i][j] += triGl1 * Glplus[j][l2];
            }
          }
        }
      }
    }
  }
  else if(pm1>0 && pm2==0){
    for (l1 = 3; l1 < LMAX; l1++){
      l1_double = (double)l1;
      cov_g_l = func_for_cov_G_shear(l1_double, ar);
      for(i=0; i<like.Ntheta ; i++){
        covGl1 = cov_g_l * Glplus[i][l1];
        // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
        for(j=0; j<like.Ntheta ; j++){
          cov[i][j] += covGl1 * Glminus[j][l1];
        }
      }
      if(FLAG_NG){
        for (l2 = 3; l2 < LMAX; l2++){
          tri = bin_cov_NG_shear_shear_tomo(l1_double,(double)l2,z1,z2,z3,z4);
          for(i=0; i<like.Ntheta ; i++){
            triGl1 = tri * Glplus[i][l1];
            for(j=0; j<like.Ntheta ; j++){
              covNG[i][j] += triGl1 * Glminus[j][l2];
            }
          }
        }
      }
    }
  }
  else if(pm1==0 && pm2>0){
    for (l1 = 3; l1 < LMAX; l1++){
      l1_double = (double)l1;
      cov_g_l = func_for_cov_G_shear(l1_double, ar);
      for(i=0; i<like.Ntheta ; i++){
        covGl1 = cov_g_l * Glminus[i][l1];
        // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
        for(j=0; j<like.Ntheta ; j++){
          cov[i][j] += covGl1 * Glplus[j][l1];
        }
      }
      if(FLAG_NG){
        for (l2 = 3; l2 < LMAX; l2++){
          tri = bin_cov_NG_shear_shear_tomo(l1_double,(double)l2,z1,z2,z3,z4);
          for(i=0; i<like.Ntheta ; i++){
            triGl1 = tri * Glminus[i][l1];
            for(j=0; j<like.Ntheta ; j++){
              covNG[i][j] += triGl1 * Glplus[j][l2];
            }
          }
        }
      }
    }
  }
  else{
    for (l1 = 3; l1 < LMAX; l1++){
      l1_double = (double)l1;
      cov_g_l = func_for_cov_G_shear(l1_double, ar);
      for(i=0; i<like.Ntheta ; i++){
        covGl1 = cov_g_l * Glminus[i][l1];
        // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
        for(j=0; j<like.Ntheta ; j++){
          cov[i][j] += covGl1 * Glminus[j][l1];
        }
      }
      if(FLAG_NG){
        for (l2 = 3; l2 < LMAX; l2++){
          tri = bin_cov_NG_shear_shear_tomo(l1_double,(double)l2,z1,z2,z3,z4);
          for(i=0; i<like.Ntheta ; i++){
            triGl1 = tri * Glminus[i][l1];
            for(j=0; j<like.Ntheta ; j++){
              covNG[i][j] += triGl1 * Glminus[j][l2];
            }
          }
        }
      }
    }
  }
  for(i=0; i<like.Ntheta; i++){cov[i][i] += 2.* N[i];}
}

void cov_gl_shear_real_binned_fullsky(double **cov, double **covNG, int z1,int z2,int z3,int z4,int pm, int FLAG_NG, double *theta, double *dtheta){
  
  int i,j;
  static int LMAX = 50000;
  static double **Glplus =0;
  static double **Glminus =0;
  static double **Pl =0;
  if (Glplus ==0){
    Glplus =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
    Glminus =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
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
      double x = cos(like.theta[i]);
      gsl_sf_legendre_Pl_deriv_array(LMAX, xmin[i],Pmin,dPmin);
      gsl_sf_legendre_Pl_deriv_array(LMAX, xmax[i],Pmax,dPmax);
      // printf("xmax[%d]:%lg\n", i,xmax[i]);
      for (int l = 2; l < LMAX; l ++){
        /*double plm = gsl_sf_legendre_Plm(l,2,x);
        double plm_1 = gsl_sf_legendre_Plm(l-1,2,x);
        Glplus[i][l] = (2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))
        *(plm*((4-l+2.*x*(l-1))/(1-x*x)-l*(l+1)/2)
        +plm_1*(l-1,2,x)*(l+2)*(x-2)/(1-x*x));


        Glminus[i][l] = (2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))
        *(plm*(l,2,x)*((4-l-2.*x*(l-1))/(1-x*x)-l*(l+1)/2)
        +plm_1*(l-1,2,x)*(l+2)*(x+2)/(1-x*x));*/

        Glplus[i][l] =(2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(

        -l*(l-1.)/2*(l+2./(2*l+1)) * (Pmin[l-1]-Pmax[l-1])
        -l*(l-1.)*(2.-l)/2         * (xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
        +l*(l-1)/(2*l+1)           * (Pmin[l+1]-Pmax[l+1])

        +(4-l)   * (dPmin[l]-dPmax[l])
        +(l+2)   * (xmin[i]*dPmin[l-1] - xmax[i]*dPmax[l-1] - Pmin[l-1] + Pmax[l-1])

        +2*(l-1) * (xmin[i]*dPmin[l]   - xmax[i]*dPmax[l]   - Pmin[l] + Pmax[l])
        -2*(l+2) * (dPmin[l-1]-dPmax[l-1])

        )/(xmin[i]-xmax[i]);           

        Glminus[i][l] =(2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(

        -l*(l-1.)/2*(l+2./(2*l+1)) * (Pmin[l-1]-Pmax[l-1])
        -l*(l-1.)*(2.-l)/2         * (xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
        +l*(l-1)/(2*l+1)           * (Pmin[l+1]-Pmax[l+1])

        +(4-l)   * (dPmin[l]-dPmax[l])
        +(l+2)   * (xmin[i]*dPmin[l-1] - xmax[i]*dPmax[l-1] - Pmin[l-1] + Pmax[l-1])

        -2*(l-1) * (xmin[i]*dPmin[l]   - xmax[i]*dPmax[l]   - Pmin[l] + Pmax[l])
        +2*(l+2) * (dPmin[l-1]-dPmax[l-1])

        )/(xmin[i]-xmax[i]);
        
        // printf("Pmin[%d], Pmax[%d]:%lg, %lg\n", l,l, Pmin[l], Pmax[l]);
        // printf("Glplus[%d][%d],%lg\n", i,l, Glplus[i][l]);
        Pl[i][l] = (2.*l+1)/(4.*M_PI*l*(l+1)*(xmin[i]-xmax[i]))
        *((l+2./(2*l+1.))*(Pmin[l-1]-Pmax[l-1])
        +(2-l)*(xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
        -2./(2*l+1.)*(Pmin[l+1]-Pmax[l]));
      }
      Glplus[i][2] = 0.; Glminus[i][2] = 0.;
    }
    free_double_vector(xmin,0,like.Ntheta-1);
    free_double_vector(xmax,0,like.Ntheta-1);
    free_double_vector(Pmin,0,LMAX+1);
    free_double_vector(Pmax,0,LMAX+1);
    free_double_vector(dPmin,0,LMAX+1);
    free_double_vector(dPmax,0,LMAX+1);
  }

  double l1_double,tri;
  int l1,l2;
  for(i=0; i<like.Ntheta ; i++){
    for(j=0; j<like.Ntheta ; j++){
      cov[i][j] = 0.;
      covNG[i][j] = 0.;
    }
  }

  // printf("i,j:%d,%d\n", i,j);
  double triP;
  double covGl1, cov_g_l;
  int ar[4];
  ar[0] = z1; ar[1] = z2; ar[2] = z3; ar[3] = z4;

  if(pm>0){
    for (l1 = 2; l1 < LMAX; l1++){
      l1_double = (double)l1;
      // printf("l1,%d\n", l1);
      cov_g_l = func_for_cov_G_gl_shear(l1_double, ar);
      for(i=0; i<like.Ntheta ; i++){
        covGl1 = cov_g_l * Pl[i][l1];
        // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
        for(j=0; j<like.Ntheta ; j++){
          cov[i][j] += covGl1 * Glplus[j][l1];
        }
      }

      if(FLAG_NG){
        for (l2 = 3; l2 < LMAX; l2++){
          tri = bin_cov_NG_gl_shear_tomo(l1_double,(double)l2,z1,z2,z3,z4);
          for(i=0; i<like.Ntheta ; i++){
            triP = tri * Pl[i][l1];
            // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
            for(j=0; j<like.Ntheta ; j++){
              covNG[i][j] += triP * Glplus[j][l2];
            }
          }
        }
      }
    }
  }
  else{
    for (l1 = 2; l1 < LMAX; l1++){
      l1_double = (double)l1;
      cov_g_l = func_for_cov_G_gl_shear(l1_double, ar);
      for(i=0; i<like.Ntheta ; i++){
        covGl1 = cov_g_l * Pl[i][l1];
        // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
        for(j=0; j<like.Ntheta ; j++){
          cov[i][j] += covGl1 * Glminus[j][l1];
        }
      }

      if(FLAG_NG){
        for (l2 = 3; l2 < LMAX; l2++){
          tri = bin_cov_NG_gl_shear_tomo(l1_double,(double)l2,z1,z2,z3,z4);
          for(i=0; i<like.Ntheta ; i++){
            triP = tri * Pl[i][l1];
            for(j=0; j<like.Ntheta ; j++){
              covNG[i][j] += triP * Glminus[j][l2];
            }
          }
        }
      }
    }
  }
}

void cov_cl_shear_real_binned_fullsky(double **cov, double **covNG, int z1,int z2,int z3,int z4,int pm, int FLAG_NG, double *theta, double *dtheta){
  
  int i,j;
  static int LMAX = 10000;
  static double **Glplus =0;
  static double **Glminus =0;
  static double **Pl =0;
  if (Glplus ==0){
    Glplus =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
    Glminus =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
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
      double x = cos(like.theta[i]);
      gsl_sf_legendre_Pl_deriv_array(LMAX, xmin[i],Pmin,dPmin);
      gsl_sf_legendre_Pl_deriv_array(LMAX, xmax[i],Pmax,dPmax);
      // printf("xmax[%d]:%lg\n", i,xmax[i]);
      for (int l = 3; l < LMAX; l ++){
        /*double plm = gsl_sf_legendre_Plm(l,2,x);
        double plm_1 = gsl_sf_legendre_Plm(l-1,2,x);
        Glplus[i][l] = (2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))
        *(plm*((4-l+2.*x*(l-1))/(1-x*x)-l*(l+1)/2)
        +plm_1*(l-1,2,x)*(l+2)*(x-2)/(1-x*x));


        Glminus[i][l] = (2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))
        *(plm*(l,2,x)*((4-l-2.*x*(l-1))/(1-x*x)-l*(l+1)/2)
        +plm_1*(l-1,2,x)*(l+2)*(x+2)/(1-x*x));*/

        Glplus[i][l] =(2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(

        -l*(l-1.)/2*(l+2./(2*l+1)) * (Pmin[l-1]-Pmax[l-1])
        -l*(l-1.)*(2.-l)/2         * (xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
        +l*(l-1)/(2*l+1)           * (Pmin[l+1]-Pmax[l+1])

        +(4-l)   * (dPmin[l]-dPmax[l])
        +(l+2)   * (xmin[i]*dPmin[l-1] - xmax[i]*dPmax[l-1] - Pmin[l-1] + Pmax[l-1])

        +2*(l-1) * (xmin[i]*dPmin[l]   - xmax[i]*dPmax[l]   - Pmin[l] + Pmax[l])
        -2*(l+2) * (dPmin[l-1]-dPmax[l-1])

        )/(xmin[i]-xmax[i]);           

        Glminus[i][l] =(2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(

        -l*(l-1.)/2*(l+2./(2*l+1)) * (Pmin[l-1]-Pmax[l-1])
        -l*(l-1.)*(2.-l)/2         * (xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
        +l*(l-1)/(2*l+1)           * (Pmin[l+1]-Pmax[l+1])

        +(4-l)   * (dPmin[l]-dPmax[l])
        +(l+2)   * (xmin[i]*dPmin[l-1] - xmax[i]*dPmax[l-1] - Pmin[l-1] + Pmax[l-1])

        -2*(l-1) * (xmin[i]*dPmin[l]   - xmax[i]*dPmax[l]   - Pmin[l] + Pmax[l])
        +2*(l+2) * (dPmin[l-1]-dPmax[l-1])

        )/(xmin[i]-xmax[i]);
        
        Pl[i][l] = 1./(4.*M_PI)*(Pmin[l+1]-Pmax[l+1]-Pmin[l-1]+Pmax[l-1])/(xmin[i]-xmax[i]);
        // printf("Pmin[%d], Pmax[%d]:%lg, %lg\n", l,l, Pmin[l], Pmax[l]);
        // printf("Glplus[%d][%d],%lg\n", i,l, Glplus[i][l]);
      }
      Pl[i][1] = 1./(4.*M_PI)*(Pmin[2]-Pmax[2]-Pmin[0]+Pmax[0])/(xmin[i]-xmax[i]);
      Pl[i][2] = 1./(4.*M_PI)*(Pmin[3]-Pmax[3]-Pmin[1]+Pmax[1])/(xmin[i]-xmax[i]);
      Glplus[i][1] = 0.; Glplus[i][2] = 0.;
      Glminus[i][1] = 0.; Glminus[i][2] = 0.;
    }
    free_double_vector(xmin,0,like.Ntheta-1);
    free_double_vector(xmax,0,like.Ntheta-1);
    free_double_vector(Pmin,0,LMAX+1);
    free_double_vector(Pmax,0,LMAX+1);
    free_double_vector(dPmin,0,LMAX+1);
    free_double_vector(dPmax,0,LMAX+1);
  }

  double l1_double,tri;
  int l1,l2;
  for(i=0; i<like.Ntheta ; i++){
    for(j=0; j<like.Ntheta ; j++){
      cov[i][j] = 0.;
      covNG[i][j] = 0.;
    }
  }

  // printf("i,j:%d,%d\n", i,j);
  double triP;
  double covGl1, cov_g_l;
  int ar[4];
  ar[0] = z1; ar[1] = z2; ar[2] = z3; ar[3] = z4;

  if(pm>0){
    for (l1 = 1; l1 < LMAX; l1++){
      l1_double = (double)l1;
      // printf("l1,%d\n", l1);
      cov_g_l = func_for_cov_G_cl_shear(l1_double, ar);
      for(i=0; i<like.Ntheta ; i++){
        covGl1 = cov_g_l * Pl[i][l1];
        // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
        for(j=0; j<like.Ntheta ; j++){
          cov[i][j] += covGl1 * Glplus[j][l1];
        }
      }

      if(FLAG_NG){
        for (l2 = 3; l2 < LMAX; l2++){
          tri = bin_cov_NG_cl_shear_tomo(l1_double,(double)l2,z1,z2,z3,z4);
          for(i=0; i<like.Ntheta ; i++){
            triP = tri * Pl[i][l1];
            // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
            for(j=0; j<like.Ntheta ; j++){
              covNG[i][j] += triP * Glplus[j][l2];
            }
          }
        }
      }
    }
  }
  else{
    for (l1 = 1; l1 < LMAX; l1++){
      l1_double = (double)l1;
      cov_g_l = func_for_cov_G_cl_shear(l1_double, ar);
      for(i=0; i<like.Ntheta ; i++){
        covGl1 = cov_g_l * Pl[i][l1];
        // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
        for(j=0; j<like.Ntheta ; j++){
          cov[i][j] += covGl1 * Glminus[j][l1];
        }
      }

      if(FLAG_NG){
        for (l2 = 3; l2 < LMAX; l2++){
          tri = bin_cov_NG_cl_shear_tomo(l1_double,(double)l2,z1,z2,z3,z4);
          // printf("tri:%lg\n", tri);
          for(i=0; i<like.Ntheta ; i++){
            triP = tri * Pl[i][l1];
            for(j=0; j<like.Ntheta ; j++){
              covNG[i][j] += triP * Glminus[j][l2];
            }
          }
        }
      }
    }
  }
}

void cov_cl_gl_real_binned_fullsky(double **cov, double **covNG, int z1,int z2,int z3,int z4, int FLAG_NG, double *theta, double *dtheta){
  
  int i,j;
  static int LMAX = 50000;
  static double **Pl =0;
  static double **Pl2=0;
  if (Pl ==0){
    Pl =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
    Pl2=create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);

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
      double x = cos(like.theta[i]);
      gsl_sf_legendre_Pl_array(LMAX, xmin[i],Pmin);
      gsl_sf_legendre_Pl_array(LMAX, xmax[i],Pmax);
      Pl[i][1] = 1./(4.*M_PI)*(Pmin[2]-Pmax[2]-Pmin[0]+Pmax[0])/(xmin[i]-xmax[i]);
      Pl2[i][1]= 0.;
      for (int l = 2; l < LMAX; l ++){
        Pl[i][l] = 1./(4.*M_PI)*(Pmin[l+1]-Pmax[l+1]-Pmin[l-1]+Pmax[l-1])/(xmin[i]-xmax[i]);
        Pl2[i][l] = (2.*l+1)/(4.*M_PI*l*(l+1)*(xmin[i]-xmax[i]))
        *((l+2./(2*l+1.))*(Pmin[l-1]-Pmax[l-1])
        +(2-l)*(xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
        -2./(2*l+1.)*(Pmin[l+1]-Pmax[l]));
      }
    }
    free_double_vector(xmin,0,like.Ntheta-1);
    free_double_vector(xmax,0,like.Ntheta-1);
    free_double_vector(Pmin,0,LMAX+1);
    free_double_vector(Pmax,0,LMAX+1);
    free_double_vector(dPmin,0,LMAX+1);
    free_double_vector(dPmax,0,LMAX+1);
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
  int ar[4];
  ar[0] = z1; ar[1] = z2; ar[2] = z3; ar[3] = z4;

  for (l1 = 1; l1 < LMAX; l1++){
    l1_double = (double)l1;
    cov_g_l = func_for_cov_G_cl_gl(l1_double, ar);
    for(i=0; i<like.Ntheta ; i++){
      covGl1 = cov_g_l * Pl[i][l1];
      // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
      for(j=0; j<like.Ntheta ; j++){
        cov[i][j] += covGl1 * Pl2[j][l1];
      }
    }

    if(FLAG_NG){
      for (l2 = 2; l2 < LMAX; l2++){
        tri = bin_cov_NG_cl_gl_tomo(l1_double,(double)l2,z1,z2,z3,z4);
        for(i=0; i<like.Ntheta ; i++){
          triP = tri * Pl[i][l1];
          for(j=0; j<like.Ntheta ; j++){
            covNG[i][j] += triP * Pl2[j][l2];
          }
        }
      }
    }
  }
}

void cov_gl_gl_real_binned_fullsky(double **cov, double **covNG, int z1,int z2,int z3,int z4, int FLAG_NG, double *theta, double *dtheta){
  
  int i,j;
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
      double x = cos(like.theta[i]);
      gsl_sf_legendre_Pl_array(LMAX, xmin[i],Pmin);
      gsl_sf_legendre_Pl_array(LMAX, xmax[i],Pmax);
      for (int l = 2; l < LMAX; l ++){
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
    free_double_vector(dPmin,0,LMAX+1);
    free_double_vector(dPmax,0,LMAX+1);
  }

  double l1_double,tri;
  int l1,l2;
  for(i=0; i<like.Ntheta ; i++){
    for(j=0; j<like.Ntheta ; j++){
      cov[i][j] = 0.;
      covNG[i][j] = 0.;
    }
  }

  double N[like.Ntheta];
  for(i=0;i<like.Ntheta;i++) {N[i] = 0.;}

  int zl1, zs1, zl2, zs2;
  zl1 = z1; zs1 = z2; zl2 = z3; zs2 = z4;
  if (zl1 ==zl2 && zs1 ==zs2){
    for(i=0;i<like.Ntheta;i++) {
      N[i] += pow(survey.sigma_e,2.0)/(2.0*M_PI*(2.*theta[i]+dtheta[i])*dtheta[i]*nlens(zl1)*nsource(zs2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor);
    }
  }
  for(i=0;i<like.Ntheta;i++) {if(N[i]) N[i] /= w_mask(theta[i]);}

  double triP;
  double covGl1, cov_g_l;
  int ar[4];
  ar[0] = z1; ar[1] = z2; ar[2] = z3; ar[3] = z4;

  for (l1 = 2; l1 < LMAX; l1++){
    l1_double = (double)l1;
    cov_g_l = func_for_cov_G_gl(l1_double, ar);
    for(i=0; i<like.Ntheta ; i++){
      covGl1 = cov_g_l * Pl[i][l1];
      // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
      for(j=0; j<like.Ntheta ; j++){
        cov[i][j] += covGl1 * Pl[j][l1];
      }
    }

    if(FLAG_NG){
      for (l2 = 2; l2 < LMAX; l2++){
        tri = bin_cov_NG_gl_gl_tomo(l1_double,(double)l2,z1,z2,z3,z4);
        for(i=0; i<like.Ntheta ; i++){
          triP = tri * Pl[i][l1];
          for(j=0; j<like.Ntheta ; j++){
            covNG[i][j] += triP * Pl[j][l2];
          }
        }
      }
    }
  }
  for(i=0; i<like.Ntheta; i++){cov[i][i] += N[i];}
}

void cov_cl_cl_real_binned_fullsky(double **cov, double **covNG, int z1,int z2,int z3,int z4, int FLAG_NG, double *theta, double *dtheta){
  
  int i,j;
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
      double x = cos(like.theta[i]);
      gsl_sf_legendre_Pl_array(LMAX, xmin[i],Pmin);
      gsl_sf_legendre_Pl_array(LMAX, xmax[i],Pmax);
      for (int l = 1; l < LMAX; l ++){
        Pl[i][l] = 1./(4.*M_PI)*(Pmin[l+1]-Pmax[l+1]-Pmin[l-1]+Pmax[l-1])/(xmin[i]-xmax[i]);
      }
    }
    free_double_vector(xmin,0,like.Ntheta-1);
    free_double_vector(xmax,0,like.Ntheta-1);
    free_double_vector(Pmin,0,LMAX+1);
    free_double_vector(Pmax,0,LMAX+1);
    free_double_vector(dPmin,0,LMAX+1);
    free_double_vector(dPmax,0,LMAX+1);
  }

  double l1_double,tri;
  int l1,l2;
  for(i=0; i<like.Ntheta ; i++){
    for(j=0; j<like.Ntheta ; j++){
      cov[i][j] = 0.;
      covNG[i][j] = 0.;
    }
  }

  double N[like.Ntheta];
  for(i=0;i<like.Ntheta;i++) {N[i] = 0.;}

  if (z1 ==z3 && z2 ==z4){
    for(i=0;i<like.Ntheta;i++) {
      N[i] += 1./(M_PI*(2.*theta[i]+dtheta[i])*dtheta[i]*nlens(z1)*nlens(z2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor); //(number of galaxy pairs in the survey contributing to annulus of width Dtheta centered at theta1)^-1
    }
  }
  if (z1 ==z4 && z2 ==z3){
    for(i=0;i<like.Ntheta;i++) {
      N[i] += 1./(M_PI*(2.*theta[i]+dtheta[i])*dtheta[i]*nlens(z1)*nlens(z2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor);
    }
  } 
  for(i=0;i<like.Ntheta;i++) {if(N[i]) N[i] /= w_mask(theta[i]);}

  double triP;
  double covGl1, cov_g_l;
  int ar[4];
  ar[0] = z1; ar[1] = z2; ar[2] = z3; ar[3] = z4;

  for (l1 = 1; l1 < LMAX; l1++){
    l1_double = (double)l1;
    cov_g_l = func_for_cov_G_cl(l1_double, ar);
    for(i=0; i<like.Ntheta ; i++){
      covGl1 = cov_g_l * Pl[i][l1];
      // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
      for(j=0; j<like.Ntheta ; j++){
        cov[i][j] += covGl1 * Pl[j][l1];
      }
    }

    if(FLAG_NG){
      for (l2 = 1; l2 < LMAX; l2++){
        tri = bin_cov_NG_cl_cl_tomo(l1_double,(double)l2,z1,z2,z3,z4);
        for(i=0; i<like.Ntheta ; i++){
          triP = tri * Pl[i][l1];
          for(j=0; j<like.Ntheta ; j++){
            covNG[i][j] += triP * Pl[j][l2];
          }
        }
      }
    }
  }
  for(i=0; i<like.Ntheta; i++){cov[i][i] += N[i];}
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
  C24 = C_gl_tomo_all(l,n2,n4);
  C14 = C_gl_tomo_all(l,n1,n4);
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
  C14 = C_gl_tomo_all(l,n1,n4);
  C23 = C_gl_tomo_all(l,n3,n2);
  
  if (n1 == n3){N13= 1./(nlens(n1)*survey.n_gal_conversion_factor);}
  if (n2 == n4){N24= pow(survey.sigma_e,2.0)/(2.0*nsource(n2)*survey.n_gal_conversion_factor);}
  
  return (C13*C24+C13*N24+N13*C24 + C14*C23+C14*N23+N14*C23)*4.*M_PI/(survey.area*survey.area_conversion_factor* (2.*l+1.));
}

double func_for_cov_G_gl_shear(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[0]; n2 = ar[1]; n3 = ar[2]; n4 = ar[3];
  //printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_gl_tomo_all(l,n1,n3);
  C24 = C_shear_tomo(l,n2,n4);
  C14 = C_gl_tomo_all(l,n1,n4);
  C23 = C_shear_tomo(l,n2,n3);
  if (n2 == n4){N24= pow(survey.sigma_e,2.0)/(2.0*nsource(n2)*survey.n_gal_conversion_factor);}
  if (n2 == n3){N23= pow(survey.sigma_e,2.0)/(2.0*nsource(n2)*survey.n_gal_conversion_factor);}
 
  return (C13*C24+C13*N24+N13*C24 + C14*C23+C14*N23+N14*C23)*4.*M_PI/(survey.area*survey.area_conversion_factor* (2.*l+1.));
}
