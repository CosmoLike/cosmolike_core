/************************* covariance routines for shear shear angular correlation functions *********************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
// Note the look-up tables for power spectrum covariances are recomputed if one the redshift bins changes
//      Loop over theta1, theta2 first, before computing the next combination of redshift bins
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
double HSV_shear_shear_tomo(double l1, double l2,int z1, int z2, int z3, int z4);

/************************* routines for angular trispectrum terms ************************/
double project_tri_2h_cov_shear_shear_tomo(double l1,double l2, int z1, int z2, int z3, int z4);
double project_tri_1h_cov_shear_shear_tomo(double l1,double l2, int z1, int z2, int z3, int z4);
double project_tri_lin_cov_shear_shear_tomo(double l1,double l2, int z1, int z2, int z3, int z4);

/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
// Note the look-up tables for power spectrum covariances are recomputed if one the redshift bins changes
//      Loop over theta1, theta2 first, before computing the next combination of redshift bins
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/


/************************* covariance routines for angular power spectra ************************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       in order to speed up calculations, these routines restrict different terms to the following regimes:
       tri_lin tabulated in the range exp(limits.halo_s_lin_min)< l1,l2 <exp(limits.halo_s_lin_max)
       tri_1h  tabulated in the range exp(limits.halo_s_1h_min)< l1,l2 <exp(limits.halo_s_1h_max)
       tri_2h  tabulated in the range exp(limits.halo_s_2h_min)< l1,l2 <exp(limits.halo_s_2h_max)
       HSV is only calculated is sqrt(l1*l2) > 100.
       outside these ranges the routines return 0., NOT the actual value of the different terms!
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
*/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
// Note the look-up tables for power spectrum covariances are recomputed if one the redshift bins changes
//      Loop over theta1, theta2 first, before computing the next combination of redshift bins
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/

/********************* cosmic shear covariance ******************************/

double inner_project_tri_lin_cov_shear_shear_tomo(double a,void *params)
{
  double k1,k2,res,wa,weight1,m04;
  double *ar = (double *) params;
  wa = chi(a);
  k1 = ar[0]/wa;
  k2 = ar[1]/wa;
  //printf("inner_project_tri_lin_cov_tomo %le %le %le\n",wa,k1,k2);
  m04 = tri_lin_cov(k1,k2,a);
  weight1 = g_tomo(a,(int) ar[2])*g_tomo(a, (int) ar[3])*g_tomo(a, (int) ar[4])*g_tomo(a, (int) ar[5])/a/a/a/a/a/a/hoverh0(a)/wa/wa;
  //printf("inner_project_tri_lin_cov_tomo %le %le %le %le\n",ar[2],ar[3],ar[4],ar[5]);
  //printf("inner_project_tri_lin_cov_tomo %le %le %le %le %le %le %le\n",1./a-1.,g_tomo(a,(int) ar[2]),g_tomo(a, (int) ar[3]),g_tomo(a, (int) ar[4]),g_tomo(a, (int) ar[5]),weight1,m04);
  res = weight1*m04;
  return res;
}

double inner_project_tri_1h_cov_shear_shear_tomo(double a, void *params)
{
  double k1,k2,res,wa,weight1,m04;
  double *ar = (double *) params;
  wa = chi(a);
  k1 = ar[0]/wa;
  k2 = ar[1]/wa;
  m04 = tri_1h_cov(k1,k2,a);
  weight1 = g_tomo(a,(int) ar[2])*g_tomo(a, (int) ar[3])*g_tomo(a, (int) ar[4])*g_tomo(a, (int) ar[5])/a/a/a/a/a/a/hoverh0(a)/wa/wa;
  res = weight1*m04;
  return res;
}

double inner_project_tri_2h_cov_shear_shear_tomo(double a,void *params)
{
  double k1,k2,res,wa,weight1,m04;
  double *ar = (double *) params;
  wa = chi(a);
  k1 = ar[0]/wa;
  k2 = ar[1]/wa;
  m04 = tri_2h_cov(k1,k2,a);
  weight1 = g_tomo(a,(int) ar[2])*g_tomo(a, (int) ar[3])*g_tomo(a, (int) ar[4])*g_tomo(a, (int) ar[5])/a/a/a/a/a/a/hoverh0(a)/wa/wa;
  res = weight1*m04;
  return res;
}

double project_tri_2h_cov_shear_shear_tomo(double l1,double l2, int z1, int z2, int z3, int z4){
  static int Z1 = -42;
  static int Z2 = -42;
  static int Z3 = -42;
  static int Z4 = -42;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  int i,j;
  double res,as,val,slog1,slog2;
  
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 || Z4!=z4 )
    {
    if (table!=0) free_double_matrix(table,0, Ntable.N_halo_s2-1, 0, Ntable.N_halo_s2-1);
    table = create_double_matrix(0, Ntable.N_halo_s2-1, 0, Ntable.N_halo_s2-1);
    double array[6];
    logsmin = log(limits.halo_s_2h_min);
    logsmax = log(limits.halo_s_2h_max);
    ds = (logsmax - logsmin)/(Ntable.N_halo_s2 - 1.);
    slog1 = logsmin;
    as = 1./(redshift.shear_zdistrpar_zmax +1.);
    array[2] = (double) z1;
    array[3] = (double) z2;
    array[4] = (double) z3;
    array[5] = (double) z4;
    
    for (i=0; i<Ntable.N_halo_s2; i++, slog1+=ds) {
      array[0] = exp(slog1);
      slog2 = slog1;
      for (j=i; j<Ntable.N_halo_s2; j++, slog2+=ds) {
        array[1] = exp(slog2);
        res = pow(3./2.*cosmology.Omega_m,4.)*int_gsl_integrate_low_precision(inner_project_tri_2h_cov_shear_shear_tomo,(void*)array,as,0.999999,NULL,1000);
        table[i][j]=log(res);
        table[j][i]=log(res);        
      }
    }
    Z1=z1;
    Z2=z2;
    Z3=z3;
    Z4=z4;
    }
  val = 0.;
  slog1=log(l1);
  slog2=log(l2);
  if (slog1 > logsmin && slog2 > logsmin && slog1 < logsmax && slog2 < logsmax){
    val = exp(interpol2d(table, Ntable.N_halo_s2, logsmin, logsmax, ds, slog1, Ntable.N_halo_s2, logsmin, logsmax, ds, slog2,1.0,1.0));}
  return val;
}


double project_tri_1h_cov_shear_shear_tomo(double l1,double l2, int z1, int z2, int z3, int z4){
  static int Z1 = -42;
  static int Z2 = -42;
  static int Z3 = -42;
  static int Z4 = -42;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  int i,j;
  double res,as,val,slog1,slog2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 || Z4!=z4 )
    {
    if (table!=0) free_double_matrix(table,0, Ntable.N_halo_s1-1, 0, Ntable.N_halo_s1-1);
    table = create_double_matrix(0, Ntable.N_halo_s1-1, 0, Ntable.N_halo_s1-1);
    double array[6];
    logsmin = log(limits.halo_s_1h_min);
    logsmax = log(limits.halo_s_1h_max);
    ds = (logsmax - logsmin)/(Ntable.N_halo_s1 - 1.);
    slog1 = logsmin;
    as = 1./(redshift.shear_zdistrpar_zmax +1.);
    array[2] = (double) z1;
    array[3] = (double) z2;
    array[4] = (double) z3;
    array[5] = (double) z4;
    for (i=0; i<Ntable.N_halo_s1; i++, slog1+=ds) {
      array[0] = exp(slog1);
      slog2 = slog1;
      for (j=i; j<Ntable.N_halo_s1; j++, slog2+=ds) {
        array[1] = exp(slog2);
        res =pow(3./2.*cosmology.Omega_m,4.)*int_gsl_integrate_low_precision(inner_project_tri_1h_cov_shear_shear_tomo,(void*)array,as,0.999999,NULL,1000);
        table[i][j]=log(res);
        table[j][i]=log(res);
      }
    }
    Z1=z1;
    Z2=z2;
    Z3=z3;
    Z4=z4;
    }
  val = 0.;
  slog1=log(l1);
  slog2=log(l2);
  if (slog1 > logsmin && slog2 > logsmin && slog1 < logsmax && slog2 < logsmax){
    val = exp(interpol2d(table, Ntable.N_halo_s1, logsmin, logsmax, ds, slog1, Ntable.N_halo_s1, logsmin, logsmax, ds, slog2,1.0,1.0));}
  return val;
}


double project_tri_lin_cov_shear_shear_tomo(double l1,double l2, int z1, int z2, int z3, int z4){
  static int Z1 = -42;
  static int Z2 = -42;
  static int Z3 = -42;
  static int Z4 = -42;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  int i,j;
  double res,as,val,slog1,slog2;
  if (Z1!=z1 || Z2!=z2 || Z3!=z3 || Z4!=z4 )
    {
    if (table!=0) free_double_matrix(table,0, Ntable.N_halo_lin-1, 0, Ntable.N_halo_lin-1);
    table = create_double_matrix(0, Ntable.N_halo_lin-1, 0, Ntable.N_halo_lin-1);
    double array[6];
    
    logsmin = log(limits.halo_s_lin_min);
    logsmax = log(limits.halo_s_lin_max);
    ds = (logsmax - logsmin)/(Ntable.N_halo_lin - 1.);
    slog1 = logsmin;
    as = 1./(redshift.shear_zdistrpar_zmax +1.);
    array[2] = (double) z1;
    array[3] = (double) z2;
    array[4] = (double) z3;
    array[5] = (double) z4;
    for (i=0; i<Ntable.N_halo_lin; i++, slog1+=ds) {
      array[0] = exp(slog1);
      slog2 = slog1;
      for (j=i; j<Ntable.N_halo_lin; j++, slog2+=ds) {
        array[1] = exp(slog2);
        res =pow(3./2.*cosmology.Omega_m,4.)*int_gsl_integrate_low_precision(inner_project_tri_lin_cov_shear_shear_tomo,(void*)array,as,0.999999,NULL,1000);
        table[i][j]=log(res);
        table[j][i]=log(res);
      }
    }
    Z1=z1;
    Z2=z2;
    Z3=z3;
    Z4=z4;
    }
  val = 0.;
  slog1=log(l1);
  slog2=log(l2);
  if (slog1 > logsmin && slog2 > logsmin && slog1 < logsmax && slog2 < logsmax){
    val = exp(interpol2d(table, Ntable.N_halo_lin, logsmin, logsmax, ds, slog1, Ntable.N_halo_lin, logsmin, logsmax, ds, slog2,1.0,1.0));}
  return val;
}

double int_for_HSV_shear_shear_tomo(double a, void *params){
  double *ar = (double *) params;
  double weight1,weight2,wa,res =0.;
  wa=chi(a);
  weight1 = (g_tomo(a,(int) ar[2])*wa/a)*(g_tomo(a,(int) ar[3])*wa/a)*(g_tomo(a,(int) ar[4])*wa/a)*(g_tomo(a,(int) ar[5])*wa/a);
  weight2 = pow((a*a*hoverh0(a)),-1.);
  res= weight1*weight2/pow(wa,6.)*wa*wa;
  //res=res*survey_variance(a,ar[6])*bI02(ar[0]/wa,a)*bI02(ar[1]/wa,a); // beat coupling term unreliably in flat sky approximation for large surveys. aka LSST
  res=res*survey_variance(a,ar[6])*(bI02(ar[0]/wa,a)*bI02(ar[1]/wa,a)+68./21.*68./21.*p_lin(ar[0]/wa,a)*p_lin(ar[1]/wa,a));
	return res;
}

double HSV_shear_shear_tomo(double l1, double l2,int z1, int z2, int z3, int z4){
  static int Z1 = -42;
  static int Z2 = -42;
  static int Z3 = -42;
  static int Z4 = -42;
  static double FSKY = -42.;
  static double **table=0;
  static double ds = .0, logsmin = .0, logsmax = .0;
  int i,j;
  double res,val,slog1,slog2;

  if (Z1!=z1 || Z2!=z2 || Z3!=z3 || Z4!=z4 || FSKY != survey.area/41253.0)
    {
    if (table!=0) free_double_matrix(table,0, Ntable.N_halo_hsv-1, 0, Ntable.N_halo_hsv-1);
    table = create_double_matrix(0, Ntable.N_halo_hsv-1, 0, Ntable.N_halo_hsv-1);
    double array[7];
//    printf("CALCULATING HSV TABLE\n");
    logsmin = log(limits.halo_s_hsv_min);
    logsmax = log(limits.halo_s_hsv_max);
    ds = (logsmax - logsmin)/(Ntable.N_halo_hsv - 1.);
    slog1 = logsmin;
    array[2] = (double) z1;
    array[3] = (double) z2;
    array[4] = (double) z3;
    array[5] = (double) z4;
    array[6] = survey.area/41253.0;
    
    for (i=0; i<Ntable.N_halo_hsv; i++, slog1+=ds) {
      array[0] = exp(slog1);
      slog2 = slog1;
      for (j=i; j<Ntable.N_halo_hsv; j++, slog2+=ds) {
        array[1] = exp(slog2);
        res =int_gsl_integrate_low_precision(int_for_HSV_shear_shear_tomo,(void*)array,1./(1.0+redshift.shear_zdistrpar_zmax),0.9999,NULL,1000);
        table[i][j]=log(res);
        table[j][i]=log(res);
     //   printf("Element %d %d finished\n",i,j);
      }
    }
    Z1=z1;
    Z2=z2;
    Z3=z3;
    Z4=z4;
    FSKY = array[6];
    printf("CALCULATING HSV TABLE FINISHED\n");
    }
  val = 0.;
  slog1=log(l1);
  slog2=log(l2);
  if (slog1 > logsmin && slog2 > logsmin && slog1 < logsmax && slog2 < logsmax){
    val = exp(interpol2d(table, Ntable.N_halo_hsv, logsmin, logsmax, ds, slog1, Ntable.N_halo_hsv, logsmin, logsmax, ds, slog2,1.0,1.0));}
  return pow(3./2.*cosmology.Omega_m,4.)*val;
}


/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
// Note the look-up tables for power spectrum covariances are recomputed if one the redshift bins changes
//      Loop over theta1, theta2 first, before computing the next combination of redshift bins
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/********************* end covariance routines for angular correlation functions *********************/


