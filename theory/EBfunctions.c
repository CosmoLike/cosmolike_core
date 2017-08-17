
/**************** COSEBIS filter functions ********************/
double integrand_j0_Tfunc(double vartheta, void *params);
double W_calc_j0(double ell,int ORDER,double max,double min);
double integrand_j2_Tfunc(double vartheta, void *params);
double W_calc_j2(double ell,int ORDER,double max,double min);
void create_Ztable_lin(double vt_max,double vt_min,double thetarad,char *filename,int table_Nring_z);
void create_Ztable_log(double vt_max,double vt_min,double thetarad,char *filename,int table_Nring_z);
double T_lin_plus_COSEBIS(double vartheta,int ORDER,double vt_max,double vt_min);
double T_lin_minus_integrand_COSEBIS(double t, void *params);
double T_lin_minus_COSEBIS(double theta,int ORDER,double vt_max,double vt_min);
double T_log_plus_COSEBIS(double vartheta,int ORDER,double vt_max,double vt_min);
double T_log_minus_integrand_COSEBIS(double theta, void *params);
double T_log_minus_COSEBIS(double vartheta,int ORDER,double vt_max,double vt_min);
void create_Tm_COSEBIs_table(double vt_max,double vt_min,char *filename2,int table_Tm_COSEBIs,int N_ORDER);
void normalization(int vartheta_max,int vt_min,double *store);
void roots(int vartheta_max,int vt_min,double store[20][21]);
/**************** multiprobe COSEBIS ***************/
double integrand_Tm_log_shear_shear(double vartheta, void *p);
double integrand_Tp_log_shear_shear(double vartheta, void *p);
double integrand_Tp_log_shear_position(double vartheta, void *p);
double integrand_Tp_log_position_position(double vartheta, void *p);
double integrand_Tp_log_shear_magnification(double vartheta, void *p);
double integrand_Tp_log_position_magnification(double vartheta, void *p);
void COS_shear_shear_tomo(double *E,double *B,int ORDER,double vt_max, double vt_min,int ni, int nj);
void COS_pos_pos_tomo(double *E,int ORDER,double vt_max, double vt_min,int ni, int nj);
void COS_shear_mag_tomo(double *E,int ORDER,double vt_max, double vt_min,int ni, int nj);
void COS_shear_pos_tomo(double *E,int ORDER,double vt_max, double vt_min,int ni, int nj);
void COS_mag_pos_tomo(double *E,int ORDER,double vt_max, double vt_min,int ni, int nj);
void COS_mag_mag_tomo(double *E,int ORDER,double vt_max, double vt_min,int ni, int nj);
/*************** ring statistic filter functions *******************/
void Y_minus(double x, double eta, double *ym);
void Y_plus(double x, double eta, double *yp);
double integrand2_minus(double theta2, void *params);
double integrand1_minus(double theta1, void *params);
double integrand2_plus(double theta2, void *params);
double integrand1_plus(double theta1, void *params);
double Z_minus(double vartheta,double *array);
double Z_plus(double vartheta, double *array);
double W1(double theta_1, double zeta_1, double zeta_2);
double W2(double theta_2, double zeta_3, double zeta_4);
/*************** MAP filter functions ******************/
double T_plus_MAP(double x);
double T_minus_MAP(double x);


double integrand_j0_Tfunc(double vartheta, void *params)
{
  double ell,j0,f,Tp;
  int ORDER;
  
  double *array = (double *) params;
  ell = array[0];
  ORDER = (int) array[1]; 
  double vt_max= array[2];
  double vt_min= array[3];
  //printf("%le %le %le %d\n",vartheta,vt_max/constants.arcmin,vt_min/constants.arcmin,ORDER);
  j0=gsl_sf_bessel_J0(vartheta*ell);
  Tp=T_log_plus_COSEBIS(vartheta,ORDER,vt_max,vt_min);
  f= vartheta*j0*Tp;
  //printf("%le %le %le\n",vartheta,j0,Tp);
  return f;
}

double W_calc_j0(double ell,int ORDER,double max,double min)
{
  double array[4],res,error;
  gsl_function Fpp;
  Fpp.function=&integrand_j0_Tfunc;
  Fpp.params = array;
  array[0]=ell;
  array[1]=ORDER*1.0;
  array[2]=max*constants.arcmin;
  array[3]=min*constants.arcmin;
  double upper=max*constants.arcmin;
  double lower=min*constants.arcmin;
  
  gsl_integration_cquad_workspace *w2;
  w2=gsl_integration_cquad_workspace_alloc(10000);
  gsl_integration_cquad(&Fpp,lower,upper,0,1.e-7,w2,&res,&error,0);
  //printf("%le %le %le\n",lower, upper,res);
  gsl_integration_cquad_workspace_free(w2);
  //printf("%le %le %le\n",lower, upper,res);
  return res;
}

double integrand_j2_Tfunc(double vartheta, void *params)
{
  double ell,j2,f,Tp;
  int ORDER;
  
  double *array = (double *) params;
  ell = array[0];
  ORDER = (int) array[1]; 
  double vt_max= array[2];
  double vt_min= array[3];
  //printf("%le %le %le %d\n",vartheta,vt_max,vt_min,ORDER);
  j2=gsl_sf_bessel_Jn(2,vartheta*ell);
  Tp=T_log_plus_COSEBIS(vartheta,ORDER,vt_max,vt_min);
  f= vartheta*j2*Tp;
  //	printf("%le %le %le\n",vartheta,j0,Tp);
  return f;
}


double W_calc_j2(double ell,int ORDER,double max,double min)
{
  double array[4],res,error;
  gsl_function Fpp;
  Fpp.function=&integrand_j2_Tfunc;
  Fpp.params = array;
  array[0]=ell;
  array[1]=ORDER*1.0;
  array[2]=max*constants.arcmin;
  array[3]=min*constants.arcmin;
  double upper=max*constants.arcmin;
  double lower=min*constants.arcmin;
  
  gsl_integration_cquad_workspace *w2;
  w2=gsl_integration_cquad_workspace_alloc(10000);
  gsl_integration_cquad(&Fpp,lower,upper,0,1.e-7,w2,&res,&error,0);
  //printf("%le %le %le\n",lower, upper,res);
  gsl_integration_cquad_workspace_free(w2);
  //printf("%le %le %le\n",lower, upper,res);
  return res;
}

///////////////////////////////////////////////////////////////////
// Calculation of z_table
///////////////////////////////////////////////////////////////////


void create_Ztable_lin(double vt_max,double vt_min,double thetarad,char *filename,int table_Nring_z)
{
  double zeta1,zeta2,zeta3,zeta4,psi,eta,array[6],zplus,zminus,vartheta,dvartheta;
  
  int i;
  FILE *F;    
  
  psi=thetarad; // psi=vt in the notation of SK07
  eta=vt_min/psi;
  printf("CALCULATING NEW Z_TABLE psi=%le eta=%le\n",psi,eta);
  
  zeta1=(1.-eta)/8.*psi;
  zeta2=3.*(1.-eta)/8.*psi;
  zeta3=(5.*eta +3.)/8.*psi;
  zeta4=(3.*eta+5.)/8.*psi;
  
  array[0]=zeta1; //zeta_1
  array[1]=zeta2;  //zeta_2
  array[2]=zeta3; //zeta_3
  array[3]=zeta4;  //zeta_4
  array[4]=0.;
  array[5]=0.;
  //printf("arraycheck Z_table %le %le %le %le %le %le\n",array[0], array[1],array[2], array[3],array[4], array[5]);
  
  F=fopen(filename,"w");
  // printf("%s\n",filename);
  dvartheta=(psi - vt_min)/(table_Nring_z-1);
  for (i=0;i<table_Nring_z;i++){
    vartheta=vt_min+i*dvartheta;
    zminus=Z_minus(vartheta,array);
    zplus=Z_plus(vartheta,array);
    fprintf(F,"%le %le %le\n",vartheta,zplus,zminus);
  }
  fclose(F);
}

void create_Ztable_log(double vt_max,double vt_min,double thetarad,char *filename,int table_Nring_z)
{
  double zeta1,zeta2,zeta3,zeta4,psi,eta,array[6],zplus,zminus,vartheta,dvartheta;
  
  int i;
  FILE *F;    
  
  psi=thetarad; // psi=vt in the notation of SK07
  eta=vt_min/psi;
  printf("CALCULATING NEW Z_TABLE psi=%le eta=%le\n",psi,eta);
  
  zeta1=(1.-eta)/8.*psi;
  zeta2=3.*(1.-eta)/8.*psi;
  zeta3=(5.*eta +3.)/8.*psi;
  zeta4=(3.*eta+5.)/8.*psi;
  
  array[0]=zeta1; //zeta_1
  array[1]=zeta2;  //zeta_2
  array[2]=zeta3; //zeta_3
  array[3]=zeta4;  //zeta_4
  array[4]=0.;
  array[5]=0.;
  //printf("arraycheck Z_table %le %le %le %le %le %le\n",array[0], array[1],array[2], array[3],array[4], array[5]);
  
  F=fopen(filename,"w");
  // printf("%s\n",filename);
  dvartheta=(log(psi) - log(vt_min))/(table_Nring_z-1);
  for (i=0;i<table_Nring_z;i++){
    vartheta=exp(log(vt_min)+i*dvartheta);
    zminus=Z_minus(vartheta,array);
    zplus=Z_plus(vartheta,array);
    fprintf(F,"%le %le %le\n",vartheta,zplus,zminus);
  }
  fclose(F);
}


///////////////////////////////////////////////////////////////////
// COSEBIS LINEAR FILTER FUNCTIONS
///////////////////////////////////////////////////////////////////

double T_lin_plus_COSEBIS(double vartheta,int ORDER,double vt_max,double vt_min)
{
  double B,Y,F=0.0,x,P;
  
  B=(vt_max - vt_min)/(vt_max + vt_min);
  x=(2.0*vartheta-vt_min-vt_max)/(vt_max-vt_min);
  //printf("%le\n",x);
  if((x<(-1.0)) || (x>1.0)){
    printf("ERROR vt_min=%le vt=%le vt_max=%le\n",vt_min,vartheta,vt_max);
    return 0.0;
    
  }
  if (ORDER==1){
    Y=8.0*(25.0+5.0*B*B+6.0*pow(B,4.0))/5.0;
    F=1.0/sqrt(Y)*(3.0*B*B-5.0-6.0*B*x+3.0*(5.0-B*B)*x*x);
  }
  
  if (ORDER==2){
    Y=8.0*(25.0+5.0*B*B+6.0*pow(B,4.0))*(175.0+35.0*B*B+45.0*pow(B,4.0)+pow(B,6.0));
    F=1.0/sqrt(Y)*(B*B*B*(25.0+3.0*B*B)-15.0*(35.0+9.0*B*B+8.0*pow(B,4.0))*x-15.0*B*B*B*(3.0+B*B)*x*x+35.0*(25.0+5.0*B*B+6.0*pow(B,4.0))*x*x*x);
  }
  
  if (ORDER>2){
    P=gsl_sf_legendre_Pl(ORDER+1,x);
    F=sqrt((2.0*ORDER+3.0)*0.5)*P;
  }
  return F;
}

double T_lin_minus_integrand_COSEBIS(double t, void *params)
{
  double *array = (double *) params;
  int ORDER = (int) array[0];
  double theta = array[1];
  double vt_max = array[2];
  double vt_min = array[3];
  return 4.0/pow(theta,2.0)*t*T_lin_plus_COSEBIS(t,ORDER,vt_max,vt_min)*(1.0-3.0*pow(t/theta,2.0));
  
}

double T_lin_minus_COSEBIS(double theta,int ORDER,double vt_max,double vt_min)
{
  double result,error,array[4];
  gsl_function H;
  gsl_integration_workspace *w;
  w = gsl_integration_workspace_alloc (1000);
  H.function = &T_lin_minus_integrand_COSEBIS;
  array[0] = 1.0*ORDER;
  array[1] = theta;
  array[2]=vt_max;
  array[3]=vt_min;
  //printf("%le %le\n",array[1],array[2]);
  H.params = &array;
  gsl_integration_qag (&H,vt_min,theta, 1e-8, 1e-6, 1000, GSL_INTEG_GAUSS41,w, &result, &error);
  gsl_integration_workspace_free(w);
  return result+ T_lin_plus_COSEBIS(theta,ORDER,vt_max,vt_min);
}


////////////////////////////////////////////////////////////////////
// COSEBIS LOGARITHMIC FILTER FUNCTIONS
///////////////////////////////////////////////////////////////////


double T_log_plus_COSEBIS(double vartheta,int ORDER,double vt_max,double vt_min)
{
  int i,j,k;
  double F,x,store_root[20][21],store_norm[20];
  static int MAXI=-42;
  static double **root_array=0;
  static double *norm_vector=0;
  x = log(vartheta/vt_min);
  int max=round(vt_max/constants.arcmin);	
  int min=round(vt_min/constants.arcmin);	
  if (MAXI != max){
    printf("Enter read coefficients max=%d\n",max);
    
    roots(max,min,store_root);
    normalization(max,min,store_norm);
    root_array= create_double_matrix(0, 19, 0,20);
    norm_vector= create_double_vector(0, 19);
    for (j=0;j<20;j++){
      norm_vector[j]=store_norm[j];
      for (k=0;k<21;k++){
	root_array[j][k]=store_root[j][k];
	//printf("%d %d %le %le\n",j,k,store_norm[j],store_root[j][k]);
      }
    }
    MAXI=max;		
  }
  F = norm_vector[ORDER-1];
  for (i = 0; i <=ORDER; i ++){ 
    F =F*(x-root_array[ORDER-1][i]); 
  }
  return F;
}



double T_log_minus_integrand_COSEBIS(double theta, void *params)
{
  double *array = (double *) params;
  int ORDER = (int) array[0];
  double vartheta = array[1];
  double vt_max = array[2];
  double vt_min = array[3];
  //	 printf("%le %le %le\n",array[0],array[1],array[2]);
  //double t = exp(x);
  return 4.0/pow(vartheta,2.0)*theta*T_log_plus_COSEBIS(theta,ORDER,vt_max,vt_min)*(1.0-3.0*pow(theta/vartheta,2.0));
}

double T_log_minus_COSEBIS(double vartheta,int ORDER,double vt_max,double vt_min)
{
  double result,error,array[4];
  gsl_function H;
  
  gsl_integration_workspace *w;
  w = gsl_integration_workspace_alloc (5000);
  H.function = &T_log_minus_integrand_COSEBIS;
  array[0] = 1.0*ORDER;
  array[1] = vartheta;
  array[2]=vt_max;
  array[3]=vt_min;
  //printf("%le %le %le\n",array[0],array[1],array[2]);
  H.params = &array;
  gsl_integration_qag (&H,vt_min,vartheta, 1e-8, 1e-4, 5000, GSL_INTEG_GAUSS41,w, &result, &error);
  gsl_integration_workspace_free(w);
  
  return result+ T_log_plus_COSEBIS(vartheta,ORDER,vt_max,vt_min);
}

void create_Tm_COSEBIs_table(double vt_max,double vt_min,char *filename2,int table_Tm_COSEBIs,int N_ORDER)
{
  printf("CALCULATING TM_table for Interpolation\n");
  int i,j;
  double dvt,log_thetamin,log_thetamax,theta;
  double *tm;
  tm= create_double_vector(0, N_ORDER-1);
  FILE *F2;
  F2 = fopen(filename2,"w");
  log_thetamin = log(vt_min);
  log_thetamax = log(vt_max);
  dvt = (log_thetamax - log_thetamin)/(table_Tm_COSEBIs-1);
  for (i=0; i<table_Tm_COSEBIs; i++) {
    theta=exp(log_thetamin+i*dvt);
    for (j=0;j<N_ORDER;j++){
      tm[j]=T_log_minus_COSEBIS(theta,j+1,vt_max,vt_min);
    }
    fprintf(F2,"%le %le %le %le %le %le %le %le %le %le %le \n",theta,tm[0],tm[1],tm[2],tm[3],tm[4],tm[5],tm[6],tm[7],tm[8],tm[9]);
  }
  fclose(F2);
}

void normalization(int vartheta_max,int vt_min,double *store)
{
  int j;
  
  if ((vartheta_max==299)&&(vt_min==15)){printf("Entering DES BCC 15-300 arcmin norm max=%d\n",vartheta_max);
  double store2[20]={
    1.073888485, 
    2.004396917, 
    2.929895961, 
    3.956203031, 
    5.241826328, 
    6.927608704, 
    9.157453425, 
    12.10979405, 
    16.01638711, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0   
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }
  
  if ((vartheta_max==390)&&(vt_min==5)){printf("Entering DES BCC>5arcmin norm max=%d\n",vartheta_max);
  double store2[20]={
    0.57112, 
    0.736016,
    0.773603,
    0.742261,
    0.688064,
    0.63182,
    0.579248,
    0.531086,
    0.480688,
    0.239549,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0, 
    0.0, 
    0.0,
    0.0, 
    0.0
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }
  if ((vartheta_max==390)&&(vt_min==5)){printf("Entering DES BCC>5arcmin norm max=%d\n",vartheta_max);
  double store2[20]={
    0.57112, 
    0.736016,
    0.773603,
    0.742261,
    0.688064,
    0.63182,
    0.579248,
    0.531086,
    0.480688,
    0.239549,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0, 
    0.0, 
    0.0,
    0.0, 
    0.0
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }
  
  if ((vartheta_max==390)&&(vt_min==2)){printf("Entering DES BCC norm max=%d\n",vartheta_max);
  double store2[20]={
    0.451318,
    0.468736,
    0.411671,
    0.331715,
    0.257048,
    0.196482,
    0.149677,
    0.113973,
    0.0864252,
    0.0179121,    
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0, 
    0.0, 
    0.0,
    0.0, 
    0.0
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }
  
  if ((vartheta_max==398)&&(vt_min==0)){printf("Entering Sato norm max=%d\n",vartheta_max);
  double store2[20]={
    0.349933,
    0.252056,
    0.160831,
    0.0961824,
    0.0553504,
    0.0312688,
    0.0175412,
    0.00982081,
    0.00510995,
    0.00034739,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0, 
    0.0, 
    0.0,
    0.0, 
    0.0
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }
  
  if ((vartheta_max==200)&&(vt_min==0)){printf("Entering Sato norm max=%d\n",vartheta_max);
  double store2[20]={
    0.373799,
    0.303171,
    0.21559,
    0.142491,
    0.0905009,
    0.0564889,
    0.035054,
    0.0217229,
    0.013455,
    0.00160248,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0, 
    0.0, 
    0.0,
    0.0, 
    0.0
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }
  
  if ((vartheta_max==100)&&(vt_min==0)){printf("Entering Sato norm max=%d\n",vartheta_max);
  double store2[20]={
    0.410694,
    0.381528,
    0.305128,
    0.225015,
    0.159429,
    0.111216,
    0.0772371,
    0.0535959,
    0.0370586,
    0.0130154,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0, 
    0.0, 
    0.0,
    0.0, 
    0.0
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }
  
  if ((vartheta_max==195)&&(vt_min==1)){printf("Entering norm max=%d\n",vartheta_max);
  double store2[20]={
    0.451319, 
    0.468736, 
    0.411672,
    0.331716, 
    0.257048, 
    0.196483, 
    0.149677, 
    0.113974, 
    0.086627, 
    0.0311375,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0, 
    0.0, 
    0.0,
    0.0, 
    0.0
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }
  
  
  if ((vartheta_max==388)&&(vt_min==1)){printf("Entering norm max=%d\n",vartheta_max);
  double store2[20]={
    0.399721,
    0.358188,
    0.2778,
    0.199048,
    0.137017,
    0.0928108,
    0.0625637,
    0.0421317,
    0.025781,
    0.00212743,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0, 
    0.0, 
    0.0,
    0.0, 
    0.0
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }
  
  // For 2PCF_270711
  if ((vartheta_max==117)&&(vt_min==1)){printf("Entering norm max=%d\n",vartheta_max);
  double store2[20]={
    0.509234,
    0.595968, 
    0.578008, 
    0.51246, 
    0.437852, 
    0.369864, 
    0.3117, 
    0.262655, 
    0.220575, 
    0.0510563,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0, 
    0.0, 
    0.0,
    0.0, 
    0.0
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }
  
  // For xi_tree...max 117.195 min 1.02452
  /*  if ((vartheta_max==117)&&(vt_min==1)){printf("Entering norm max=%d\n",vartheta_max);
   *  double store2[20]={
   *    0.509238,
   *    0.595977,
   *    0.57802,
   *    0.512474,
   *    0.437866,
   *    0.369878,
   *    0.311714,
   *    0.262668,
   *    0.220539,
   *    0.164363,
   *    0.0,
   *    0.0,
   *    0.0,
   *    0.0,
   *    0.0,
   *    0.0, 
   *    0.0, 
   *    0.0,
   *    0.0, 
   *    0.0
};
for (j=0;j<20;j++){
  store[j]=store2[j];
}
}*/
  
  
   if ((vartheta_max==500)&&(vt_min==3)){printf("Entering norm max=%d\n",vartheta_max);
  double store2[20]={
    0.4641339843, 
    0.4965759991, 
    0.4470336511, 
    0.368822439, 
    0.2927571965, 
    0.2293480907, 
    0.1791114698, 
    0.1398313535, 
    0.1087635522, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0   
  };
  for (j=0;j<20;j++){
  store[j]=store2[j];
  printf("%le\n",store[j]);
  }
}

 if ((vartheta_max==488)&&(vt_min==3)){printf("Entering norm max=%d\n",vartheta_max);
  double store2[20]={
    0.4694546521, 
    0.508185571, 
    0.4619592912, 
    0.3847121943, 
    0.308290813, 
    0.24388075, 
    0.1923451578, 
    0.1516532618, 
    0.1191605172,
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0, 
    0.0   
  };
  for (j=0;j<20;j++){
  store[j]=store2[j];
  //printf("%le\n",store[j]);
  }
}





  if ((vartheta_max==20)&&(vt_min==1)){printf("Entering norm max=%d\n",vartheta_max);
  double store2[20]={
    1.0912236457227,
    2.0515514155032,
    3.0200045713541,
    4.1089725504215, 
    5.4877342836006,
    7.3113138257473,
    9.7430797473534,
    12.988788739557, 
    17.320989628122,
    23.103307540560,
    30.821184089026,
    41.122696633507,
    54.873038382641,
    73.227255938616,
    97.727362675356,
    130.43199593575, 
    174.08952374005,
    232.36916188038,
    310.16954243516,
    414.03069048616
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }	
  if ((vartheta_max==40)&&(vt_min==1)){printf("Entering norm max=%d\n",vartheta_max);
  double store2[20]={
    0.74182574373416,
    1.1427470175461,
    1.3995666240945,
    1.5676550138487,
    1.7075209543777,
    1.8489954901318,
    2.0011836632941,
    2.1665367244194,
    2.3462493898256, 
    2.5414457374284, 
    2.7533548317110, 
    2.9833299238703, 
    3.2328531843402,
    3.5035419048401,
    3.7971568425455,
    4.1156123119701,
    4.4609877776564, 
    4.8355408394539, 
    5.2417215855707, 
    5.6821883399040
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }
  
  if ((vartheta_max==60)&&(vt_min==1)){printf("Entering norm max=%d\n",vartheta_max);
  double store2[20]={
    0.62618860754059,
    0.86408858423908,
    0.96221141411689,
    0.97808681153931,
    0.96268108560172,
    0.93991266153018, 
    0.91664705575816, 
    0.89411744314843,
    0.87238917617592, 
    0.85138263780436, 
    0.83102597640681, 
    0.81126476023476,
    0.79205715850229, 
    0.77336992988749, 
    0.75517581974460, 
    0.73745186833324, 
    0.72017828068489,
    0.70333765116739,
    0.68691441859013,
    0.67089447509482
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }
  if ((vartheta_max==80)&&(vt_min==1)){printf("Entering norm max=%d\n",vartheta_max);
  double store2[20]={
    0.56626142634976,
    0.72487039842188,
    0.75760770936181,
    0.72286484459930,
    0.66621862407434,
    0.60814839419719,
    0.55422723209528,
    0.50512071593879,
    0.46048851335780,
    0.41989531485430, 
    0.38294726493836,
    0.34929748405514, 
    0.31863837394303,
    0.29069509674596,
    0.26522074291435, 
    0.24199264322169, 
    0.22080944025790, 
    0.20148868631831, 
    0.18386482329605, 
    0.16778745141276
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }
  if ((vartheta_max==100)&&(vt_min==1)){printf("Entering norm max=%d\n",vartheta_max);
  double store2[20]={
    0.52878252965174,
    0.63975031401998,
    0.63785462623257, 
    0.58098657301416,
    0.51037835659613,
    0.44356168105533,
    0.38469256317696,
    0.33362383889802,
    0.28940711757330, 
    0.25110735825538, 
    0.21791428062071, 
    0.18913450292727, 
    0.16417316105271, 
    0.14251832270774, 
    0.12372842709774, 
    0.10742200308332, 
    0.093269113711252, 
    0.080984153565525, 
    0.070319735669244, 
    0.061061473507429
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }
  if ((vartheta_max==120)&&(vt_min==1)){printf("Entering norm max=%d\n",vartheta_max);
  double store2[20]={
    0.50276036838747, 
    0.58156235743233,
    0.55859316476501,
    0.49059835295600, 
    0.41513314202853,
    0.34721293364857,
    0.28969604756657,
    0.24167441885499,
    0.20166014252390, 
    0.16830881397076,
    0.14049791903637, 
    0.11729836885785, 
    0.097940100750002, 
    0.081783629332177,
    0.068297143791217, 
    0.057037926410621, 
    0.047637154147253, 
    0.039787395978856,
    0.033232287310470, 
    0.027757978814574
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }
  if ((vartheta_max==140)&&(vt_min==1)){printf("Entering norm max=%d\n",vartheta_max);
  double store2[20]={
    0.48344728165833,
    0.53886395525006,
    0.50188899910901, 
    0.42785139918169,
    0.35114724020862,
    0.28465527086006, 
    0.23011712131113,
    0.18598640175380,
    0.15035130683575,
    0.12157094827480, 
    0.098317041577689, 
    0.079521954706800, 
    0.064326797194670,
    0.052039626511560,
    0.042102396276127,
    0.034064697373383, 
    0.027562790728743, 
    0.022302805493337, 
    0.018047241862012,
    0.014604109686068
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }
  if ((vartheta_max==160)&&(vt_min==1)){printf("Entering norm max=%d\n",vartheta_max);
  double store2[20]={
    0.46843547273784,
    0.50595938643927, 
    0.45908923131898, 
    0.38164654436586, 
    0.30528297177407, 
    0.24105607229218, 
    0.18976306473082,
    0.14933851964669, 
    0.11754884805712, 
    0.092546625889411,
    0.072875124999194,
    0.057392805692747, 
    0.045204572349726, 
    0.035607765803491, 
    0.028050302574079, 
    0.022098128734794, 
    0.017409825375348,
    0.013716743371921, 
    0.010807435803052, 
    0.0085154428136202
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }
  if ((vartheta_max==180)&&(vt_min==1)){printf("Entering norm max=%d\n",vartheta_max);
  double store2[20]={
    0.45636444390726,
    0.47967723071642,
    0.42549534821079, 
    0.34612909151039,
    0.27082170506591,
    0.20906636927620, 
    0.16086254562106,
    0.12372537450580,
    0.095179146603768, 
    0.073235145685632, 
    0.056360375331525,
    0.043379823093303, 
    0.033392471054205, 
    0.025706732151380, 
    0.019791363281447, 
    0.015238061239278,
    0.011732882374609, 
    0.0090343606843047, 
    0.0069567308998912, 
    0.0053570526996466
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }
  if ((vartheta_max==200)&&(vt_min==1)){printf("Entering norm max=%d\n",vartheta_max);
  double store2[20]={
    0.44640304554138,
    0.45810301495360,
    0.39832983425476, 
    0.31792067632480,
    0.24398554164518, 
    0.18466181945520,
    0.13927181201200,
    0.10499079083366,
    0.079161022095734, 
    0.059698757409474,
    0.045029364316117,
    0.033969257156831, 
    0.025628504882317, 
    0.019337400845341, 
    0.014591619820024,
    0.011011186444115, 
    0.0083097078198266, 
    0.0062712645112597, 
    0.0047330334253004, 
    0.0035722099203729
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }
  if ((vartheta_max==220)&&(vt_min==1)){printf("Entering norm max=%d\n",vartheta_max);
  double store2[20]={
    0.43801257147586,
    0.44000818187654,
    0.37584196745600, 
    0.29493541721038,
    0.22249302107249, 
    0.16546670963520,
    0.12259994087143,
    0.090790945134080,
    0.067245038408870,
    0.049816205127788, 
    0.036911158751909, 
    0.027352973780891,
    0.020272088989055, 
    0.015025541606172,
    0.011137619919591, 
    0.0082561940390852, 
    0.0061205217104915, 
    0.0045374801912399, 
    0.0033640012149185, 
    0.0024940807239186
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }
  if ((vartheta_max==240)&&(vt_min==1)){printf("Entering norm max=%d\n",vartheta_max);
  double store2[20]={
    0.43082726853786,
    0.42456563848438, 
    0.35687139809238,
    0.27581552214413,
    0.20488747859690,
    0.14999400582798,
    0.10937999071121, 
    0.079716287583726,
    0.058105215871721,
    0.042361840112644, 
    0.030889519635248,
    0.022527211417000, 
    0.016430504175724, 
    0.011984835740584, 
    0.0087426673504077, 
    0.0063779497901702, 
    0.0046530671033214, 
    0.0033948086628911, 
    0.0024768882023462,
    0.0018072177289358
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }
  
  if ((vartheta_max==260)&&(vt_min==1)){printf("Entering norm max=%d\n",vartheta_max);
  double store2[20]={
    0.42458933795335, 
    0.41119665397929, 
    0.34061735680836, 
    0.25963879330725,
    0.19019636139118, 
    0.13726819929304, 
    0.098666546152480,
    0.070874378986115,
    0.050916802391828,
    0.036586727870032, 
    0.026294362051095, 
    0.018899997633525, 
    0.013586516479499, 
    0.0097676995128339,
    0.0070227494844915, 
    0.0050494882623604,
    0.0036308529346908, 
    0.0026108849426756, 
    0.0018775086165015, 
    0.0013501720869739
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }
  if ((vartheta_max==280)&&(vt_min==1)){printf("Entering norm max=%d\n",vartheta_max);
  double store2[20]={
    0.41911138239619,
    0.39948321127190,
    0.32650855963886, 
    0.24575664895965,
    0.17774609008824, 
    0.12662424678200, 
    0.089825255992900,
    0.063675883098726,
    0.045143772331248,
    0.032011811856103, 
    0.022703836468210,
    0.016104542020660, 
    0.011424701060780,
    0.0081054880246688,
    0.0057510116438607, 
    0.0040806999176473,
    0.0028956508473178, 
    0.0020548282562317, 
    0.0014582097797411, 
    0.0010348501691530
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }
  if ((vartheta_max==300)&&(vt_min==1)){printf("Entering norm max=%d\n",vartheta_max);
  double store2[20]={
    0.41425361441223, 
    0.38911523687435, 
    0.31412595430654, 
    0.23369944289826,
    0.16705564367617,
    0.11759396611014, 
    0.082415938803209,
    0.057717719987604,
    0.040424765888775, 
    0.028318767258058, 
    0.019841646760076, 
    0.013904049391605,
    0.0097443371303574, 
    0.0068296934457931,
    0.0047871927388211, 
    0.0033557220153272,
    0.0023524054892629, 
    0.0016491349909581, 
    0.0011561531937024, 
    0.00081056446407846
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }
  if ((vartheta_max==320)&&(vt_min==1)){printf("Entering norm max=%d\n",vartheta_max);
  double store2[20]={
    0.40990944415918,
    0.37985746834025,
    0.30315481308290, 
    0.22311857845032,
    0.15777247645705, 
    0.10983867180043,
    0.076124331949080, 
    0.052715888118344, 
    0.036508514691081, 
    0.025289150125277,
    0.017520673074454,
    0.012140250095296, 
    0.0084130207783686,
    0.0058306134505022, 
    0.0040411723023437, 
    0.0028010823237015, 
    0.0019416258278669, 
    0.0013459317493912, 
    0.00093303019929416, 
    0.00064681692042618
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }
  if ((vartheta_max==340)&&(vt_min==1)){printf("Entering norm max=%d\n",vartheta_max);
  double store2[20]={
    0.40599604345134,
    0.37152789943633, 
    0.29335391520001,
    0.21374974888997,
    0.14963237908404,
    0.10310756932397,
    0.070720527677440, 
    0.048465124781118, 
    0.033215605083623, 
    0.022768865355182, 
    0.015610507580369,
    0.010704170729773, 
    0.0073406848027407,
    0.0050345207379167,
    0.0034531111260716, 
    0.0023685818000921, 
    0.0016247532026099, 
    0.0011145620697704,
    0.00076460351962958, 
    0.00052454329901607
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }
  if ((vartheta_max==360)&&(vt_min==1)){printf("Entering norm max=%d\n",vartheta_max);
  double store2[20]={
    0.40244798083753,
    0.36398332370982, 
    0.28453508994345, 
    0.20538882250424,
    0.14243347630546,
    0.097211147389497,
    0.066032733187478, 
    0.044813765175670,
    0.030414952008302, 
    0.020646582970698, 
    0.014017984549474,
    0.0095188262861504,
    0.0064644091173770,
    0.0043904832126248,
    0.0029821309043078,
    0.0020256598585520, 
    0.0013760288478192, 
    0.00093477354412463, 
    0.00063503910918492,
    0.00043142723429499
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }
  if ((vartheta_max==380)&&(vt_min==1)){printf("Entering norm max=%d\n",vartheta_max);
  double store2[20]={
    0.39921281564061,
    0.35710937928180,
    0.27654926234789, 
    0.19787557088369,
    0.13601888428298,
    0.092003640926758,
    0.061930180867685, 
    0.041647562001478,
    0.028008836334342, 
    0.018840172150314, 
    0.012675086382025, 
    0.0085285948297945,
    0.0057392045768468, 
    0.0038624602767355, 
    0.0025996044581931, 
    0.0017497499060982, 
    0.0011777846763970, 
    0.00079281834354151,
    0.00053369935252706, 
    0.00035927971220248
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }
  if ((vartheta_max==400)&&(vt_min==120)){printf("Entering norm max=%d min=%d\n",vartheta_max,vt_min);
  double store2[20]={
    8.1170735967431,
    30.806024698229,
    102.90972096846,
    340.40392260543, 
    1126.9953389013, 
    3734.6372378732, 
    12383.348364986, 
    41077.862412615,
    136303.21091981, 
    452375.68650198, 
    1.5016354942392e06,
    4.9852374980480e06,
    1.6552043829022e07, 
    5.4960824363611e07,
    1.8250894869183e08,
    6.0609314830490e08, 
    2.0128663075324e09,
    6.6850967905075e09,
    2.2203179070859e10,
    7.3745455976423e10
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }
  if ((vartheta_max==400)&&(vt_min==1)){printf("Entering norm max=%d\n",vartheta_max);
  double store2[20]={
    0.39624797808907,
    0.35081353288763,
    0.26927670488432, 
    0.19108241557446,
    0.13026484491879,
    0.087371163868121,
    0.058311690880485, 
    0.038878975480315, 
    0.025923109297229,
    0.017287916936047, 
    0.011531197662712, 
    0.0076924860382132, 
    0.0051322374468443,
    0.0034244032818082, 
    0.0022850411458931, 
    0.0015248552700360, 
    0.0010176172410441, 
    0.00067913821592976, 
    0.00045325967950598, 
    0.00030251648436782
  };
  for (j=0;j<20;j++){
    store[j]=store2[j];
  }
  }
}

void roots(int vartheta_max,int vt_min,double store[20][21])
{
  int j,k;
 
  if ((vartheta_max==500)&&(vt_min==3)){printf("Entering DES MICE 3-500 arcmin roots max=%d\n",vartheta_max);
  double store2[20][21]={{3.810201821, 4.925192939, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1.246305848, 4.020899575, 4.945380806, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.6530200054, 2.40569325, 4.214269627, 4.967069795, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.3927266699, 1.69261327, 3.107287461, 4.378416557, 4.988212852, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.2620714466, 1.209830874, 2.452585842, 3.578947419, 4.515597187, 5.007968429, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.1866329253, 0.8996113825, 1.917416156, 2.991327806, 3.908890094, 4.6276781, 5.025654343, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.1394961251, 0.6894404437, 1.525230981, 2.469951343, 3.382416766, 4.143733298, 4.716417737, 5.040744475, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.1079717803, 0.5425811172, 1.230600639, 2.053537533, 2.897066031, 3.674518842, 4.314714225, 4.785055269, 5.053084017, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.08361381652, 0.4268867232, 0.9914583471, 1.699724291, 2.466892201, 3.216110553, 3.888609812, 4.437355576, 4.835598555, 5.062491567, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }
  }
  
  
  if ((vartheta_max==299)&&(vt_min==15)){printf("Entering DES BCC 15-300 arcmin roots max=%d\n",vartheta_max);
  double store2[20][21]={{1.781530422, 2.83323302, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.5064086009, 2.086726258, 2.865096765, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.2983180693, 1.211887072, 2.30681162, 2.894000031, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.1861621416, 0.8666040441, 1.696885622, 2.474331486, 2.919378118, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.1293378671, 0.6213482345, 1.324169675, 2.021026025, 2.600157509, 2.940513687, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.09435318904, 0.4671068291, 1.03100594, 1.662361298, 2.240009704, 2.691500624, 2.957046589, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.07191008496, 0.3616240842, 0.820429869, 1.364975593, 1.913307277, 2.393582754, 2.757194699, 2.969470817, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0.05657175975, 0.2876617122, 0.6643546022, 1.131762262, 1.629556267, 2.102827282, 2.505760623, 2.805240548, 2.978729948, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.04587173888, 0.2348659359, 0.548934112, 0.9502329528, 1.395212608, 1.841023577, 2.249954072, 2.591144424, 2.841603003, 2.985772772, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }
  }
  
  
  
  
  
  if ((vartheta_max==390)&&(vt_min==5)){printf("Entering DES BCC>5arcmin roots max=%d\n",vartheta_max);
  double store2[20][21]={{3.05737, 4.16549,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.934956, 3.30016, 4.1889,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, {0.504669, 1.93793, 3.50591, 4.21288, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, {0.308251, 1.36918, 2.57158, 3.67433, 4.23558, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.208736, 0.979183, 2.02386, 2.99887, 3.81109, 4.25619, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, {0.150029, 0.730916, 1.5787, 2.49374, 3.29397, 3.91909, 4.274, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.112966, 0.56215, 1.25581, 2.05415, 2.83709, 3.50191, 4.00178, 4.28861,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.0879046, 0.443834, 1.01381, 1.70527, 2.42414, 3.09411, 3.65273, 4.06412, 4.30014,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.0664887, 0.342266, 0.804117, 1.39494, 2.04598, 2.69077, 3.2754, 3.75693, 4.10772, 4.3084,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {-0.393981, 0.0985658, 0.459438, 0.989525, 1.60696, 2.24564, 2.85277, 3.38852, 3.8229, 4.13579, 4.31381, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}; 
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }
  }
  
  
  
  if ((vartheta_max==390)&&(vt_min==2)){printf("Entering DES BCC roots max=%d\n",vartheta_max);
  double store2[20][21]={{3.93894, 5.0545,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1.30445, 4.14493, 5.07424,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, {0.680853, 2.48891, 4.33614, 5.09558, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, {0.408345, 1.75051, 3.20102, 4.4994, 5.11647, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.271821, 1.25115, 2.52776, 3.67937, 4.63643, 5.13607, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 }, {0.193259, 0.929762, 1.97704, 3.0777, 4.01462, 4.74893, 5.1537, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.144261, 0.712113, 1.57268, 2.5424, 3.47655, 4.25371, 4.83847, 5.16883,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.111533, 0.560032, 1.26866, 2.11422, 2.97887, 3.77431, 4.42792, 4.90803, 5.18127,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.0861985, 0.439971, 1.0212, 1.74911, 2.53614, 3.3035, 3.99153, 4.55223, 4.95908, 5.19074,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {-2.5213, 0.0990539, 0.495117, 1.11873, 1.86829, 2.65298, 3.40075, 4.06041, 4.593, 4.97642, 5.19405, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}}; 
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }
  }
  
  if ((vartheta_max==398)&&(vt_min==0)){printf("Entering Sato roots max=%d\n",vartheta_max);
  double store2[20][21]={{5.82911, 6.94704, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {2.33504, 5.98814, 6.96243, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1.18972, 3.81482, 6.1504, 6.97977, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.686587, 2.69621, 4.64858, 6.2989, 6.99752, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.441055, 1.92948, 3.70124, 5.20237, 6.42987, 7.01486, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.305639, 1.4236, 2.9169, 4.39631, 5.59781, 6.54318, 7.03127, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.223702, 1.08149, 2.32345, 3.65898, 4.89661, 5.88779, 6.63891, 7.04624, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.170294, 0.843122, 1.87206, 3.05451, 4.22247, 5.26863, 6.10317, 6.71761, 7.0594, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.143651, 0.71863, 1.61955, 2.68559, 3.7748, 4.78438, 5.64965, 6.3298, 6.80869, 7.07569, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {-12.37, 0.136436, 0.686446, 1.55874, 2.60428, 3.68612, 4.70099, 5.58456, 6.28513, 6.78908, 7.07206, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }
  }
  
  
  if ((vartheta_max==100)&&(vt_min==0)){printf("Entering Sato roots max=%d\n",vartheta_max);
  double store2[20][21]={{4.44859, 5.56554, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1.54935, 4.63808, 5.58373, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.798862, 2.8275, 4.82089, 5.6038, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.47394, 1.98766, 3.57833, 4.9804, 5.62374, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.312435, 1.42054, 2.83114, 4.08098, 5.11636, 5.64271, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.220669, 1.05323, 2.21842, 3.42381, 4.43559, 5.22995, 5.66008, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.163909, 0.804942, 1.76515, 2.83375, 3.85244, 4.69055, 5.32213, 5.6753, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.1266, 0.633019, 1.42556, 2.36106, 3.30823, 4.17303, 4.87773, 5.39519, 5.68812, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.103078, 0.519591, 1.18659, 2.00135, 2.86482, 3.69433, 4.43024, 5.02517, 5.45512, 5.6991, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.0993085, 0.50016, 1.14138, 1.92534, 2.75959, 3.56879, 4.29772, 4.90581, 5.36924, 5.66654, 5.8458, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }
  }
  
  if ((vartheta_max==200)&&(vt_min==0)){printf("Entering Sato roots max=%d\n",vartheta_max);
  double store2[20][21]={{5.13858, 6.25627, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1.9195, 5.31061, 6.27285, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.980306, 3.30862, 5.48266, 6.29144, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.573309, 2.32955, 4.10477, 5.63671, 6.31022, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.373056, 1.66545, 3.25713, 4.63539, 5.77045, 6.32835, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.261036, 1.23149, 2.55923, 3.90342, 5.01233, 5.88435, 6.34526, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.192499, 0.938302, 2.03717, 3.23947, 4.36947, 5.28617, 5.97885, 6.36043, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.147581, 0.734434, 1.64268, 2.70084, 3.75916, 4.71634, 5.48792, 6.05507, 6.37347, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.115389, 0.58329, 1.33592, 2.25644, 3.22936, 4.15915, 4.98019, 5.63909, 6.11529, 6.38429, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.112141, 0.568321, 1.3065, 2.21584, 3.18405, 4.11653, 4.94599, 5.61704, 6.10518, 6.38227, 11.7507, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }
  }
  
  
  
  if ((vartheta_max==195)&&(vt_min==1)){printf("Entering roots max=%d\n",vartheta_max);
  double store2[20][21]={{3.93894, 5.0545, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1.30444, 4.14492,5.07423, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.680852, 2.4889, 4.33614, 5.09557, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.408345, 1.75051, 3.20101, 4.4994,5.11646, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.271821, 1.25115, 2.52776, 3.67937, 4.63643, 5.13606, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.19326, 0.929766, 1.97705, 3.0777, 4.01462, 4.74893, 5.1537, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.144281, 0.712187, 1.57279, 2.54249, 3.47662, 4.25375, 4.83849, 5.16883, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.111699, 0.560649, 1.26953, 2.11503, 2.97948, 3.77469, 4.42813, 4.90811, 5.18129, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.0872692, 0.444047, 1.02744, 1.75607, 2.54282, 3.30921, 3.99571, 4.55482, 4.9602, 5.19095, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0.0793085, 0.407293, 0.954525, 1.65289, 2.4221, 3.1866, 3.88662, 4.47366, 4.91628, 5.18073, 5.85277, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }
  }
  
  if ((vartheta_max==388)&&(vt_min==1)){printf("Entering roots max=%d\n",vartheta_max);
  double store2[20][21]={{4.6234, 5.74062, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1.63888, 4.80797, 5.75835, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.842372, 2.94696, 4.98797, 5.77802, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.49791, 2.07198, 3.71003, 5.14612, 5.79765, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.327151, 1.48085, 2.93739, 4.22026, 5.28159, 5.81641, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.230524, 1.09714, 2.30321,  3.54409, 4.58093, 5.39537, 5.83367, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.170907, 0.837801, 1.83271, 2.93523, 3.98249, 4.84093, 5.48826, 5.84889, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.132059, 0.659406, 1.4819, 2.4488, 3.42403, 4.31184,   5.03314, 5.56269, 5.86186, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.114642, 0.574691, 1.30029, 2.16828, 3.06726, 3.91241, 4.64607, 5.23006, 5.64425, 5.87694, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {-10.5043, 0.11227, 0.563237, 1.27617, 2.13239, 3.02453, 3.86947, 4.61077, 5.20558, 5.63296, 5.87475, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }	
  }
  
  /*if ((vartheta_max==117)&&(vt_min==1)){printf("Entering roots max=%d\n",vartheta_max);
   *  double store2[20][21]={{3.43645, 4.54899, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1.08575, 3.66213,   4.57064, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.576418, 2.16946, 3.86174, 4.59341, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.34939, 1.52891, 2.83875, 4.02823, 4.61531, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.234844, 1.09306, 2.23753, 3.28957, 4.16553, 4.63551, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.168023, 0.8143, 1.74721, 2.7428, 3.60305, 4.27596, 4.6533, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.126051, 0.625156, 1.38982, 2.26194, 3.11071, 3.82498, 4.36201, 4.66821, 0, 0,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.097981, 0.493447, 1.12282, 1.88056, 2.66239, 3.38671, 3.98677, 4.42795, 4.68023, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.0800347, 0.406186, 0.936238, 1.5942, 2.30207, 2.99144, 3.60946, 4.1147, 4.48117, 4.69024, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.0783462, 0.397697, 0.917194, 1.5637, 2.2627, 2.94907, 3.57125, 4.08732, 4.
   * 46744, 4.68731, 6.57592, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
   *  for (j=0;j<20;j++){
   *    for (k=0;k<21;k++){
   *      store[j][k]=store2[j][k];
}
}	
}*/
  
  if ((vartheta_max==117)&&(vt_min==1)){printf("Entering roots max=%d\n",vartheta_max);
  double store2[20][21]={{3.43642, 4.54896, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1.08574, 3.6621, 4.57061, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.576412, 2.16944, 3.86171, 4.59338, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.349386, 1.5289, 2.83873, 4.02821, 4.61528, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.234842, 1.09305, 2.23751, 3.28955, 4.1655, 4.63548, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.168021, 0.814294, 1.7472, 2.74278, 3.60303, 4.27594, 4.65327, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.126053, 0.625161, 1.38983, 2.26194, 3.11069, 3.82497, 4.36198, 4.66818, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.0978884, 0.493026, 1.12204, 1.87958, 2.66145, 3.38599, 3.9863, 4.42774, 4.68017, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.0760503, 0.388621, 0.90441, 1.55498, 2.26391, 2.96017, 3.58775, 4.10206, 4.47583, 4.68921, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.0287829, 0.231391, 0.651411, 1.22922, 1.88959, 2.56421, 3.19855, 3.75264, 4.19808, 4.51666, 4.69706, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }	
  }
  
  if ((vartheta_max==20)&&(vt_min==1)){printf("Entering roots max=%d\n",vartheta_max);
  double store2[20][21]={{1.7595119, 2.8090662, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.49994399, 2.0655790, 2.8411155, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.29510374, 1.1998939, 2.2857259, 2.8701148, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0.18421383, 0.85817965, 1.6817682, 2.4530436, 2.8955304, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.12804079, 0.61536349, 1.3120723, 2.0037021, 2.5784760, 2.9166495, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.093431561, 0.46266409, 1.0215860, 1.6477609, 2.2211475, 2.6693726, 2.9331316, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.071222315, 0.35822860, 0.81293952, 1.3529209, 1.8969197, 2.3736598, 2.7346791, 2.9454964, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.056030409, 0.28495029, 0.65824128, 1.1216402, 1.6154049, 2.0850412, 2.4850359, 2.7824031, 2.9546988, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.045217768, 0.23168456, 0.54206637, 0.93939040, 1.3806800, 1.8233455, 2.2297480, 2.5690309, 2.8181822, 2.9616303, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.037249655, 0.19188977, 0.45313329, 0.79520810, 1.1867717, 1.5951420, 1.9895164, 2.3430826, 2.6339744, 2.8456736, 2.9669535, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.031212226, 0.16142979, 0.38385420, 0.68007907, 1.0269736, 1.3994104, 1.7726008, 2.1237914, 2.4333344, 2.6852141, 2.8672503, 2.9711235, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.026529599, 0.13762692, 0.32899378, 0.58719341, 0.89490814, 1.2327282, 1.5807927, 1.9201168, 2.2335342, 2.5062881, 2.7263392, 2.8844939, 2.9744493, 0, 0, 0, 0, 0, 0, 0, 0}, {0.022825406, 0.11868915, 0.28489993, 0.51145070, 0.78518832, 1.0910040, 1.4129845, 1.7354140, 2.0435565, 2.3242090, 2.5660509, 2.7598379, 2.8984906, 2.9771442, 0, 0, 0, 0, 0, 0, 0}, {0.019845238, 0.10338456, 0.24898013, 0.44904357, 0.69344421, 0.97031586, 1.2668673, 1.5701282, 1.8675735, 2.1476054, 2.3998931, 2.6155934, 2.7874795, 2.9100063, 2.9793580, 0, 0, 0, 0, 0, 0}, {0.017412208, 0.090845223, 0.21936279, 0.39711633, 0.61619911, 0.86720352, 1.1397960, 1.4232696, 1.7070332, 1.9810120, 2.2359558, 2.4636580, 2.6571019, 2.8105499, 2.9195941, 2.9811987, 0, 0, 0, 0, 0},
  {0.015400265, 0.080446258, 0.19467436, 0.35351037, 0.55070643, 0.77873325, 1.0291920, 1.2932245, 1.5618910, 1.8264975, 2.0788607, 2.3115115, 2.5178433, 2.6922121, 2.8300013, 2.9276610, 2.9827458, 0, 0, 0, 0}, {0.013717643, 0.071728957, 0.17389104, 0.31657875, 0.49479799, 0.70246690, 0.93271430, 1.1781851, 1.4313319, 1.6846771, 1.9310360, 2.1636974, 2.3765631, 2.5642510, 2.7221665, 2.8465508, 2.9345121, 2.9840584, 0, 0, 0}, {0.012296256, 0.064350703, 0.15623864, 0.28505245, 0.44675755, 0.63639752, 0.84831294, 1.0763702, 1.3141864, 1.5553380, 1.7935450, 2.0228275, 2.2376315, 2.4329275, 2.6042843, 2.7479217, 2.8607469, 2.9403798, 2.9851817, 0, 0},
  {0.011084733, 0.058051565, 0.14112381, 0.25794389, 0.40522016, 0.57887931, 0.77422848, 0.98612815, 1.2091658, 1.4378216, 1.6666194, 1.8902588, 2.1037266, 2.3023877, 2.4820556, 2.6390475, 2.7702233, 2.8730144, 2.9454435, 2.9861504, 0}, {0.010043728, 0.052631547, 0.12808616, 0.23447698, 0.36909403, 0.52856301, 0.70896690, 0.90597749, 1.1149901, 1.3312551, 1.5500008, 1.7665442, 1.9763873, 2.1752985, 2.3593786, 2.5251143, 2.6694181, 2.7896595, 2.8836867, 2.9498436, 2.9869915}};
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }	
  }
  
  if ((vartheta_max==40)&&(vt_min==1)){printf("Entering roots max=%d\n",vartheta_max);
  
  double store2[20][21]={{2.4079438, 3.4995002, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.70274070, 2.6832946, 3.5266715, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.39387010, 1.5592887, 2.8982608, 3.5529842, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.24344914, 1.1082313, 2.1239473, 3.0682278, 3.5771029, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.16700746, 0.79324120, 1.6659358, 2.5041276, 3.2017618, 3.5982437, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.12094960, 0.59424204, 1.2977684, 2.0715703, 2.7634513, 3.3032619, 3.6157376, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.091636741, 0.45850674, 1.0324895, 1.7034237, 2.3701626, 2.9452456, 3.3785233, 3.6295057, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.071745140, 0.36355084, 0.83505065, 1.4135505, 2.0222138, 2.5949740, 3.0776047, 3.4343599, 3.6400685, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.057676126, 0.29474441, 0.68670409, 1.1839020, 1.7306139, 2.2735838, 2.7678364, 3.1773202, 3.4765009, 3.6481559, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.047361104, 0.24350428, 0.57317452, 1.0017976, 1.4884328, 1.9917302, 2.4739288, 2.9032755, 3.2544760, 3.5090001, 3.6544193, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.039578979, 0.20439920, 0.48482179, 0.85620136, 1.2882233, 1.7487836, 2.2070800, 2.6355953, 3.0112006, 3.3154424, 3.5345783, 3.6593493, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0.033564997, 0.17392323, 0.41494489, 0.73868633, 1.1224161, 1.5411851, 1.9700483, 2.3857258, 2.7676474, 3.0985104, 3.3644561, 3.5550693, 3.6632944, 0, 0, 0, 0, 0, 0, 0, 0}, {0.028822388, 0.14973536, 0.35885889, 0.64286828, 0.98449422, 1.3642400, 1.7619683, 2.1581964, 2.5350293, 2.8767434, 3.1700977, 3.4044454, 3.5717379, 3.6664992, 0, 0, 0, 0, 0, 0, 0}, {0.025017004, 0.13023080, 0.31323390, 0.56395237, 0.86909252, 1.2133062, 1.5803199, 1.9539226, 2.3187413, 2.6607909, 2.9678199, 3.2294984, 3.4374931, 3.5854785, 3.6691377, 0, 0, 0, 0, 0, 0}, {0.021917534, 0.11428191, 0.27566570, 0.49832761, 0.77190320, 1.0842090, 1.4220490, 1.7719589, 2.1208300, 2.4563902, 2.7675423, 3.0445812, 3.2793139, 3.4651146, 3.5969384, 3.6713359, 0, 0, 0, 0, 0}, 
  {0.019359774, 0.10107893, 0.24439065, 0.44325774, 0.68950144, 0.97336591, 1.2840987, 1.6105087, 1.9414592, 2.2662724, 2.5750369, 2.8588264, 3.1098410, 3.3214902, 3.4884334, 3.6065955, 3.6731865, 0, 0, 0, 0}, {0.017224584, 0.090028889, 0.21809521, 0.39665208, 0.61917199, 0.87777478, 1.1636491, 1.4674713, 1.7797910, 2.0913636, 2.3934181, 2.6778619, 2.9374261, 3.1657633, 3.3575058, 3.5082974, 3.6148088, 3.6747590, 0, 0, 0}, {0.015423833, 0.080689867, 0.19578712, 0.35689847, 0.55875909, 0.79494959, 1.0582052, 1.3407313, 1.6345042, 1.9315397, 2.2241217, 2.5049843, 2.7674526, 3.0055428, 3.2140316, 3.3884994, 3.5253555, 3.6218522, 3.6761066, 0, 0},
  {0.013891205, 0.072727318, 0.17670671, 0.32274192, 0.50654476, 0.72284294, 0.96560929, 1.2283009, 1.5040931, 1.7860959, 2.0675441, 2.3419571, 2.6032668, 2.8459160, 3.0649307, 3.2559699, 3.4153594, 3.5401114, 3.6279375, 3.6772701, 0}, {0.012576031, 0.065884333, 0.16026525, 0.29319627, 0.46115286, 0.65977104, 0.88401922, 1.1283805, 1.3870353, 1.6540349, 1.9234582, 2.1895491, 2.4468297, 2.6901929, 2.9149729, 3.1169985, 3.2926302, 3.4387866, 3.5529607, 3.6332310, 3.6782816}};
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }	
  }
  
  if ((vartheta_max==60)&&(vt_min==1)){printf("Entering roots max=%d\n",vartheta_max);
  double store2[20][21]={{2.8008761, 3.9042400, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.83945934, 3.0561826, 3.9290136, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.45921341, 1.7857736, 3.2658501, 3.9538676, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.28186869, 1.2643863, 2.3934888, 3.4352094, 3.9771271, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.19184422, 0.90449016, 1.8814788, 2.8032977, 3.5711072, 3.9979876, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.13830571, 0.67610329, 1.4667466, 2.3265604, 3.0849411, 3.6769395, 4.0157309, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.10439654, 0.52064731, 1.1668165, 1.9150575, 2.6526915, 3.2828910, 3.7569635, 4.0300630, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.081500508, 0.41204635, 0.94313656, 1.5901396, 2.2658142, 2.8977894, 3.4268929, 3.8170022, 4.0412699, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0.065366121, 0.33349848, 0.77497875, 1.3319038, 1.9406511, 2.5417049, 3.0862070, 3.5352931, 3.8625554, 4.0499494,0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.053573181, 0.27511135, 0.64629705, 1.1268213, 1.6697171, 2.2284479, 2.7612858, 3.2338865, 3.6191798, 3.8977791, 4.0567134, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.044698851, 0.23062865, 0.54620293, 0.96271445, 1.4453074, 1.9576262, 2.4653637, 2.9383466, 3.3516251, 3.6855125, 3.9255495, 4.0620559, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.037855646, 0.19601669, 0.46709459, 0.83021742, 1.2592179, 1.7257352, 2.2018003, 2.6616783, 3.0829308, 3.4469256, 3.7388897, 3.9478280, 4.0663402, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.032469100, 0.16858565, 0.40364866, 0.72218350, 1.1043041, 1.5277915, 1.9699626, 2.4091366, 2.8256389, 3.2023775, 3.5251087, 3.7824799, 3.9659731, 4.0698261, 0, 0, 0, 0, 0, 0, 0}, {0.028153940, 0.14649451, 0.35207799, 0.63322450, 0.97462874, 1.3587725, 1.7672636, 2.1819723, 2.5858883, 2.9636963, 3.3021089, 3.5900195, 3.8185355, 3.9809476, 4.0726997, 0, 0, 0, 0, 0, 0}, {0.024644157, 0.12845166, 0.30964773, 0.55927195, 0.86539764, 1.2141064, 1.5904494, 1.9793117, 2.3661153, 2.7373393, 3.0808659, 3.3861843, 3.6444869, 3.8486957, 3.9934491, 4.0750964, 0, 0, 0, 0, 0}, {0.021751348, 0.11353097, 0.27435206, 0.49723808, 0.77278419, 1.0898408, 1.4362074, 1.7992853, 2.1666408, 2.5264513, 2.8678322, 3.1810561, 3.4576832, 3.6906274, 3.8741772, 4.0039934, 4.0771161, 0, 0, 0, 0}, {0.019339084, 0.10105519, 0.24469771, 0.44476135, 0.69374637, 0.98264615, 1.3014518, 1.6396450, 1.9866420, 2.3321663, 2.6665402, 2.9808979, 3.2673301, 3.5189718, 3.7300493, 3.8958986, 4.0129686, 4.0788338, 0, 0, 0}, 
  {0.017306628, 0.090520301, 0.21955722, 0.40001979, 0.62586446, 0.88975536, 1.1834348, 1.4980957, 1.8247302, 2.1544349, 2.4786624, 2.7894173, 3.0794005, 3.3421073, 3.5718898, 3.7639915, 3.9145638, 4.0206710, 4.0803070, 0, 0}, {0.015578298, 0.081545190, 0.19806801, 0.36159478, 0.56720763, 0.80888349, 1.0797688, 1.3724624, 1.6792885, 1.9925436, 2.3047098, 2.6086288, 2.8976372, 3.1656666, 3.4073133, 3.6178838, 3.7934212, 3.9307196, 4.0273303, 4.0815799, 0}, {0.014096351, 0.073837536, 0.17956198, 0.32837155, 0.51622788, 0.73814776, 0.98840876, 1.2607654, 1.5486629, 1.8454363, 2.1444886, 2.4394411, 2.7242573, 2.9933379, 3.2415932, 3.4644935, 3.6581031, 3.8191015, 3.9447957, 4.0331266, 4.0826872}};
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }	
  }
  
  if ((vartheta_max==80)&&(vt_min==1)){printf("Entering roots max=%d\n",vartheta_max);
  double store2[20][21]={{3.0831034, 4.1916112, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.94482232, 3.3246718, 4.2148927, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.50936156, 1.9533933, 3.5300189, 4.2387801, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.31096043, 1.3798353, 2.5895625, 3.6983304, 4.2614269, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.21046410, 0.98677699, 2.0382313, 3.0185361, 3.8351496, 4.2820158, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.15122465, 0.73648318, 1.5900197, 2.5105630, 3.3149382, 3.9433453, 4.2998246, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0.11384000, 0.56636988, 1.2648243, 2.0681781, 2.8556313, 3.5238588, 4.0262820, 4.3144572, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.088688709, 0.44763330, 1.0219701, 1.7180634, 2.4411524, 3.1146181, 3.6758449, 4.0890486, 4.3260545, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.071011849, 0.36186729, 0.83930228, 1.4391720, 2.0920630, 2.7340342, 3.3136730, 3.7901971, 4.1368816, 4.3351155, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.058120187, 0.29819649, 0.69951862, 1.2174321, 1.8005126, 2.3985047, 2.9669244, 3.4697337, 3.8786882, 4.1739488, 4.3422121, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.048436920, 0.24974767, 0.59082408, 1.0398815, 1.5586884, 2.1078036, 2.6504485, 3.1546289, 3.5942047, 3.9486917, 4.2032109, 4.3478329, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.040981538, 0.21209215, 0.50495893, 0.89649349, 1.3579673, 1.8585253, 2.3680362, 2.8590554, 3.3078767, 3.6949970, 4.0050574, 4.2267097, 4.3523478, 0, 0, 0, 0, 0, 0, 0, 0}, {0.035120924, 0.18227952, 0.43613129, 0.77957664, 1.1907761, 1.6455106, 2.1192607, 2.5887930, 3.0332046, 3.4344805, 3.7777201, 4.0511182, 4.2458655, 4.3560257, 0, 0, 0, 0, 0, 0, 0}, {0.030431376, 0.15829259, 0.38021773, 0.68331525, 1.0507775, 1.4634866, 1.9015126, 2.3453571, 2.7768639, 3.1798058, 3.5401995, 3.8464288, 4.0892412, 4.2616864, 4.3590604, 0, 0, 0, 0, 0, 0}, {0.026620894, 0.13871768, 0.33424007, 0.60330952, 0.93283277, 1.3076100, 1.7114142, 2.1279447, 2.5415850, 2.9379462, 3.3042142, 3.6293382, 3.9041067, 4.1211494, 4.2749039, 4.3615934, 0, 0, 0, 0, 0}, 
  {0.023483005, 0.12254221, 0.29601415, 0.53621619, 0.83282840, 1.1736705, 1.5454829, 1.9346503, 2.3278151, 2.7123511, 3.0767013, 3.4105933, 3.7051592, 3.9529862, 4.1481222, 4.2860592, 4.3637294, 0, 0, 0, 0}, {0.020868410, 0.10902645, 0.26391443, 0.47947648, 0.74748778, 1.0581082, 1.4004512, 1.7631319, 2.1347548, 2.5043132, 2.8614968, 3.1969078, 3.5022019, 3.7701679, 3.9947642, 4.1711263, 4.2955598, 4.3655473, 0, 0, 0}, {0.018667008, 0.097620414, 0.23671400, 0.43111543, 0.67420095, 0.95795644, 1.2733953, 1.6109742, 1.9609771, 2.3138461, 2.6604493, 2.9922844, 3.3016240, 3.5816105, 3.8263116, 4.0307482, 4.1909027, 4.3037175, 4.3671072, 0, 0}, 
  {0.016796176, 0.087908582, 0.21347462, 0.38959502, 0.61088324, 0.87076068, 1.1617665, 1.4758743, 1.8047931, 2.1402381, 2.4741585, 2.7989216, 3.1074518, 3.3933307, 3.6508639, 3.8751213, 4.0619590, 4.2080274, 4.3107737, 4.3684556, 0}, {0.015192941, 0.079572512, 0.19346991, 0.35370660, 0.55586228, 0.79449560, 1.0633763, 1.3557271, 1.6644606, 1.9823980, 2.3024612, 2.6178334, 2.9220879, 3.2092861, 3.4740482, 3.7116020, 3.9178132, 4.0892023, 4.2229536, 4.3169181, 4.3696293}};
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }	
  }
  
  if ((vartheta_max==100)&&(vt_min==1)){printf("Entering roots max=%d\n",vartheta_max);
  double store2[20][21]={{3.3033269, 4.4145895, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1.0314457, 3.5348006, 4.4368271, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.55056822, 2.0872305, 3.7366000, 4.4600058, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.33463491, 1.4721213, 2.7443510, 3.9038354, 4.4821837, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.22551074, 1.0525700, 2.1620080, 3.1872029, 4.0410278, 4.5025350, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.16160759, 0.78467372, 1.6875886, 2.6550208, 3.4944241, 4.1507031, 4.5203445, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0.12139768, 0.60279759, 1.3424121, 2.1886429, 3.0144346, 3.7115518, 4.2356434, 4.5351595, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.094422616, 0.47593099, 1.0843670, 1.8188004, 2.5785768, 3.2839144, 3.8695706, 4.3003791, 4.5470234, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.075503174, 0.38438496, 0.89018166, 1.5236795, 2.2108708, 2.8844017, 3.4909971, 3.9884249, 4.3498987, 4.5563580, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.061729377, 0.31649009, 0.74158236, 1.2888174, 1.9032115, 2.5316075, 3.1274214, 3.6533880, 4.0803945, 4.3883449, 4.5636995, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.051398438, 0.26487599, 0.62605944, 1.1006610, 1.6477430, 2.2254336, 2.7950497, 3.3232236, 3.7829512, 4.1531702, 4.4187285, 4.5695279, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.043454088, 0.22479544, 0.53483237, 0.94867484, 1.4355361, 1.9625852, 2.4980071, 3.0130515, 3.4830623, 3.8879032, 4.2117946, 4.4431474, 4.5742161, 0, 0, 0, 0, 0, 0, 0, 0},
  {0.037215599, 0.19308810, 0.46173697, 0.82474307, 1.2586954, 1.7377841, 2.2360503, 2.7290616, 3.1949832, 3.6151124, 3.9740694, 4.2597250, 4.4630666, 4.5780387, 0, 0, 0, 0, 0, 0, 0}, {0.032228163, 0.16759505, 0.40238214, 0.72271507, 1.1105769, 1.5455745, 2.0065664, 2.4729895, 2.9258090, 3.3481091, 3.7253893, 4.0456615, 4.2994144, 4.4795279, 4.5811952, 0, 0, 0, 0, 0, 0}, {0.028178803, 0.14680456, 0.35359578, 0.63793021, 0.98577555, 1.3809088, 1.8060924, 2.2440985, 2.6785078, 3.0942696, 3.4780454, 3.8183842, 4.1057790, 4.3326486, 4.4932879, 4.5838314, 0, 0, 0, 0, 0}, {0.024846492, 0.12963466, 0.31305169, 0.56684366, 0.87995433, 1.2393798, 1.6310207, 2.0404639, 2.4536332, 2.8572853, 3.2393518, 3.5891515, 3.8974988, 4.1567418, 4.3607541, 4.5049067, 4.5860557, 0, 0, 0, 0},
  {0.022071585, 0.11529560, 0.27901888, 0.50674072, 0.78965316, 1.1172497, 1.4779466, 1.8596771, 2.2504112, 2.6385783, 3.0133834, 3.3650259, 3.6848353, 3.9653439, 4.2003134, 4.3847333, 4.5148065, 4.5879495, 0, 0, 0}, {0.019736472, 0.10320059, 0.25019127, 0.45552539, 0.71211286, 1.0113972, 1.3438123, 1.6992337, 2.0673907, 2.4382167, 2.8021299, 3.1502413, 3.4744997, 3.7677831, 4.0239484, 4.2378528, 4.4053553, 4.5233103, 4.5895753, 0, 0}, {0.017752972, 0.092906615, 0.22557038, 0.41156512, 0.64512789, 0.91923582, 1.2259451, 1.5567346, 1.9028306, 2.2554947, 2.6062641, 2.9471419, 3.2707388, 3.5703724, 3.8401303, 4.0749073, 4.2704214, 4.4232182, 4.5306686, 4.5909813, 0},
  {0.016053929, 0.084074393, 0.20438348, 0.37357694, 0.58692795, 0.83862890, 1.1220453, 1.4299789, 1.7549233, 2.0892976, 2.4256502, 2.7568265, 3.0761014, 3.3772782, 3.6547587, 3.9035891, 4.1194875, 4.2988572, 4.4387924, 4.5370783, 4.5922055}};
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }	
  }
  
  if ((vartheta_max==120)&&(vt_min==1)){printf("Entering roots max=%d\n",vartheta_max);
  double store2[20][21]={{3.4838814, 4.5968115, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1.1054591, 3.7075495, 4.6182655, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.58580412, 2.1989956, 3.9063712, 4.6408915, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.35473014, 1.5493312, 2.8725367, 4.0725910, 4.6626932, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.23821348, 1.1076259, 2.2645644, 3.3261261, 4.2099032, 4.6828382, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.17033457, 0.82494896, 1.7685810, 2.7741634, 3.6417808, 4.3205829, 4.7006212, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.12772869, 0.63319996, 1.4068335, 2.2881690, 3.1450745, 3.8654091, 4.4070021, 4.7155560, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.099213223, 0.49951247, 1.1361698, 1.9020920, 2.6917785, 3.4229500, 4.0282509, 4.4732515, 4.7276151, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.079247514, 0.40312328, 0.93240217, 1.5935785, 2.3088271, 3.0080240, 3.6364500, 4.1507112, 4.5240954, 4.7371590, 0,0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.064732935, 0.33169349, 0.77646579, 1.3478645, 1.9879326, 2.6411334, 3.2591939, 3.8039036, 4.2454649, 4.5636354, 4.7446920, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.053859304, 0.27743408, 0.65526007, 1.1509271, 1.7212287, 2.3222852, 2.9138672, 3.4615103, 3.9375443, 4.3204586, 4.5949127, 4.7506846, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.045506021, 0.23532945, 0.55957229, 0.99181736, 1.4995492, 2.0482969, 2.6048661, 3.1394542, 3.6266512, 4.0458310, 4.3808909, 4.6200664, 4.7555108, 0, 0, 0, 0, 0, 0, 0, 0}, {0.038952026, 0.20204258, 0.48292831, 0.86207235, 1.3147424, 1.8138041, 2.3321128, 2.8442610, 3.3276667, 3.7630835, 4.1347602, 4.4303190, 4.6405963, 4.7594490, 0, 0, 0, 0, 0, 0, 0}, {0.033716224, 0.17529536, 0.42071381, 0.75526552, 1.1599164, 1.6132093, 2.0930004, 2.5778561, 3.0480305, 3.4860576, 3.8770291, 4.2086680, 4.4712647, 4.6575707, 4.7627029, 0, 0, 0, 0, 0, 0}, {0.029467948, 0.15349364, 0.36959493, 0.66652078, 1.0294491, 1.4413018, 1.8840034, 2.3395638, 2.7909106, 3.2224595, 3.6204546, 3.9731293, 4.2707466, 4.5055634, 4.6717657, 4.7654218, 0, 0, 0, 0, 0}, 
  {0.025973918, 0.13549727, 0.32712687, 0.59212636, 0.91882026, 1.2935158, 1.7014161, 2.1274497, 2.5569532, 2.9761825, 3.3726588, 3.7353755, 4.0548972, 4.3233852, 4.5345789, 4.6837567, 4.7677169, 0, 0, 0, 0}, {0.023065798, 0.12047453, 0.29149069, 0.52923828, 0.82441914, 1.1659691, 1.5417247, 1.9390548, 2.3454094, 2.7487580, 3.1379141, 3.5027541, 3.8343488, 4.1250282, 4.3684005, 4.5593423, 4.6939773, 4.7696718, 0, 0, 0}, {0.020619673, 0.10780783, 0.26131424, 0.47566010, 0.74336354, 1.0554141, 1.4017632, 1.7718050, 2.1548120, 2.5403022, 2.9183281, 3.2796891, 3.6160751, 3.9201541, 4.1856169, 4.4071927, 4.5806449, 4.7027595, 4.7713506, 0, 0}, 
  {0.018542706, 0.097031161, 0.23554879, 0.42968079, 0.67334814, 0.95915617, 1.2787585, 1.6232239, 1.9833803, 2.3501174, 2.7146365, 3.0686457, 3.4045037, 3.7153167, 3.9949999, 4.2383095, 4.4408556, 4.5991021, 4.7103610, 4.7728029, 0}, {0.016764241, 0.087787810, 0.21338299, 0.38995566, 0.61252175, 0.87496718, 1.1703205, 1.4910339, 1.8292547, 2.1770717, 2.5267279, 2.8707940, 3.2023038, 3.5148535, 3.8026703, 4.0606565, 4.2844135, 4.4702529, 4.6151985, 4.7169843, 4.7740677}};
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }	
  }
  
  if ((vartheta_max==140)&&(vt_min==1)){printf("Entering roots max=%d\n",vartheta_max);
  double store2[20][21]={{3.6368718, 4.7508970, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1.1703416, 3.8542750, 4.7717344, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.61673548, 2.2951474, 4.0505464, 4.7939114, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.37226718, 1.6158948, 2.9820949, 4.2158280, 4.8154007, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.24925024, 1.1551000, 2.3522747, 3.4443619, 4.3531294, 4.8353652, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.17788904, 0.85964523, 1.8379560, 2.8756705, 3.7668651, 4.4645281, 4.8531099, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.13319380, 0.65936132, 1.4620283, 2.3730878, 3.2561444, 3.9958416, 4.5520916, 4.8681270, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.10333951, 0.51977941, 1.1805505, 1.9732059, 2.7881307, 3.5409954, 4.1626853, 4.6195544, 4.8803359, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.082466814, 0.41920907, 0.96855987, 1.6532779, 2.3922670, 3.1130771, 3.7598234, 4.2881443, 4.6714801, 4.8900468, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.067311485, 0.34473069, 0.80632499, 1.3982978, 2.0601326, 2.7342767, 3.3710495, 3.9314817, 4.3852124, 4.7119219, 4.8977356, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.055969309, 0.28819244, 0.68024096, 1.1938548, 1.7838690, 2.4046916, 3.0147936, 3.5788022, 4.0685119, 4.4620492, 4.7439389, 4.9038633, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.047263524, 0.24434598, 0.58072507, 1.0286529, 1.5541193, 2.1212488, 2.6956797, 3.2467300, 3.7483683, 4.1795735, 4.5239842, 4.7697023, 4.9088036, 0, 0, 0, 0, 0, 0, 0, 0}, {0.040437930, 0.20970118, 0.50103700, 0.89393515, 1.3625198, 1.8785192, 2.4137801, 2.9420743, 3.4401983, 3.8884586, 4.2708040, 4.5746581, 4.7907396, 4.9128377, 0, 0, 0, 0, 0, 0, 0},
  {0.034988581, 0.18187671, 0.43637060, 0.78304037, 1.2019713, 1.6707908, 2.1665000, 2.6669281, 3.1517336, 3.6029944, 4.0054676, 4.3466421, 4.6166496, 4.8081407, 4.9161726, 0, 0, 0, 0, 0, 0}, {0.030569452, 0.15920716, 0.38325283, 0.69090861, 1.0666685, 1.4927182, 1.9502652, 2.4206720, 2.8863158, 3.3311669, 3.7411238,4.1041615, 4.4103563, 4.6518352, 4.8226979, 4.9189604, 0, 0, 0, 0, 0}, {0.026936647, 0.14050207, 0.33913684, 0.61368542, 0.95193624, 1.3396028, 1.7612905, 2.2013673, 2.6446727, 3.0770424, 3.4856556, 3.8592371, 4.1881463, 4.4643935, 4.6816096, 4.8349990, 4.9213144, 0, 0, 0, 0},
  {0.023914307, 0.12489346, 0.30212803, 0.54841620, 0.85403554, 1.2074411, 1.5959719, 2.0065155, 2.4260800, 2.8422476, 3.2435039, 3.6194569, 3.9609636, 4.2601883, 4.5106143, 4.7070273, 4.8454871, 4.9233202, 0, 0, 0}, {0.021373069, 0.11173717, 0.27079740, 0.49281849, 0.76997956, 1.0928788, 1.4510531, 1.8334853, 2.2290594, 2.6269419, 3.0168780, 3.3894066, 3.7360026, 4.0491609, 4.3224362, 4.5504533, 4.7288980, 4.8545017, 4.9250432, 0, 0}, {0.019216073, 0.10054738, 0.24405304, 0.44511401, 0.69737791, 0.99312948, 1.3236764, 1.6797365, 2.0517984, 2.4304362, 2.8065681, 3.1716548, 3.5178444, 3.8380682, 4.1260988, 4.3765790, 4.5850308, 4.7478517, 4.8623062, 4.9265342, 0}, 
  {0.017369645, 0.090952338, 0.22105037, 0.40390518, 0.63431033, 0.90588752, 1.2113756, 1.5429274, 1.8923949, 2.2515861, 2.6124858, 2.9674335, 3.3092611, 3.6313922, 3.9279094, 4.1935961, 4.4239583, 4.6152322, 4.7643846, 4.8691080, 4.9278329}};
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }	
  }
  
  if ((vartheta_max==160)&&(vt_min==1)){printf("Entering roots max=%d\n",vartheta_max);
  double store2[20][21]={{3.7695956, 4.8843830, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1.2282730, 3.9818268, 4.9047186, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.64439872, 2.3796420, 4.1758769, 4.9265200, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.38787591, 1.6745140, 3.0778556, 4.3402958, 4.9477434, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.25903684, 1.1969173, 2.4289944, 3.5473555, 4.4775142, 4.9675483, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.18456657, 0.89018504, 1.8987179, 2.9641689, 3.8755856, 4.5894474, 4.9852493, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.13801294, 0.68236671, 1.5103827, 2.4472181, 3.3528084, 4.1090833, 4.6779268, 5.0003249, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
  {0.10697118, 0.53758294, 1.2194297, 2.0353210, 2.8720682, 3.6436120, 4.2793361, 4.7463920, 5.0126528, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.085295838, 0.43332568, 1.0002257, 1.7054385, 2.4650040, 3.2044702, 3.8669839, 4.4073561, 4.7992265, 5.0225010, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.069574545, 0.35616141, 0.83246384, 1.4423646, 2.1230974, 2.8153603, 3.4682695, 4.0422289, 4.5064004, 4.8404323, 5.0303204, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.057819161, 0.29761731, 0.70209899, 1.2313598, 1.8385087, 2.4764598, 3.1025650, 3.6806787, 4.1821529, 4.5848104, 4.8730787, 5.0365623, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.048802926, 0.25223905, 0.59922461, 1.0608291, 1.6017230, 2.1848012, 2.7746901, 3.3399542, 3.8540355, 4.2955856, 4.6480285, 4.8993616, 5.0415996, 0, 0, 0, 0, 0, 0, 0, 0}, {0.041738415, 0.21640107, 0.51686685, 0.92176086, 1.4041969, 1.9349054, 2.4848550, 3.0271095, 3.5379349, 3.9972600, 4.3887839, 4.6997666, 4.9208319, 5.0457154, 0, 0, 0, 0, 0, 0, 0}, {0.036101404, 0.18763080, 0.45005098, 0.80728949, 1.2386533, 1.7209652, 2.2304799, 2.7443879, 3.2418357, 3.7045128, 4.1168942, 4.4662738, 4.7426521, 4.9385974, 5.0491193, 0, 0, 0, 0, 0, 0}, 
  {0.031532267, 0.16419985, 0.39518162, 0.71219467, 1.0991284, 1.5375209, 2.0079523, 2.4912223, 2.9692329, 3.4255727, 3.8458466,4.2178116, 4.5313884, 4.7785963, 4.9534643, 5.0519660, 0, 0, 0, 0, 0}, {0.027777714, 0.14487335, 0.34962223, 0.63249719, 0.98081285, 1.3797601, 1.8134204, 2.2656733, 2.7209278, 3.1646576, 3.5837498, 3.9667008, 4.3036976, 4.5866237, 4.8090202, 4.9660307, 5.0543705, 0, 0, 0, 0}, {0.024655238, 0.12875138, 0.31141167, 0.56514566, 0.87985602, 1.2435745, 1.6432035, 2.0652104, 2.4962198, 2.9234789, 3.3351919, 3.7207376, 4.0707917, 4.3773788, 4.6338777, 4.8349982, 4.9767478, 5.0564198, 0, 0, 0}, 
  {0.022030668, 0.11516632, 0.27907096, 0.50778227, 0.79317995, 1.1255174, 1.4939678, 1.8871542, 2.2936226, 2.7022349, 3.1024715, 3.4846481, 3.8400560, 4.1610435, 4.4410505, 4.6746142, 4.8573559, 4.9859614, 5.0581807, 0, 0}, {0.019803597, 0.10361489, 0.25147025, 0.45856987, 0.71832003, 1.0227231, 1.3627830, 1.7289103, 2.1112978, 2.5002459, 2.8864277, 3.2610915, 3.6162048, 3.9445500, 4.2397787, 4.4964381, 4.7099763, 4.8767353, 4.9939400, 5.0597047, 0}, {0.017897688, 0.093712156, 0.22773574, 0.41606452, 0.65329575, 0.93281842, 1.2471169, 1.5880819, 1.9473074, 2.3163574, 2.6869930, 3.0513544, 3.4020998, 3.7325053, 4.0365303, 4.3088576, 4.5449123, 4.7408677, 4.8936425, 5.0008949, 5.0610326}};
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }	
  }
  
  if ((vartheta_max==180)&&(vt_min==1)){printf("Entering roots max=%d\n",vartheta_max);
  double store2[20][21]={{3.8867912, 5.0021332, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1.2807163, 4.0946584, 5.0220496, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.66948606, 2.4550850, 4.2867466, 5.0435297, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.40197401, 1.7269641, 3.1629723, 4.4503711, 5.0645225, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.26784780, 1.2343427, 2.4972376, 3.6386412, 4.5874656, 5.0841847, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.19056169, 0.91750145, 1.9528286, 3.0426658, 3.9717670, 4.6998069, 5.1018402, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.14233053, 0.70292705, 1.5534549, 2.5130447, 3.4384165, 4.2091655, 4.7890383, 5.1169581, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.11021947, 0.55348000, 1.2540613, 2.0905065, 2.9464700, 3.7344020, 4.3823806, 4.8583506, 5.1293830, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.087822782, 0.44591979, 1.0284248, 1.7517930, 2.5295151, 3.2853852, 3.9617276, 4.5126314, 4.9119650, 5.1393469, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.071593694, 0.36635115, 0.85573263, 1.4815283, 2.1789620, 2.8871879, 3.5542738, 4.1400940, 4.6133973, 4.9538315, 5.1472779, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.059468083, 0.30601289, 0.72154907, 1.2646891, 1.8869964, 2.5400599, 3.1802498, 3.7707498, 4.2825379, 4.6931774, 4.9870241, 5.1536184, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.050174021, 0.25926556, 0.61567932, 1.0894183, 1.6439699, 2.2411350, 2.8446468, 3.4224122, 3.9474176, 4.3980371, 4.7575128, 5.0137589, 5.1587397, 0, 0, 0, 0, 0, 0, 0, 0}, {0.042895914, 0.22236194, 0.53094120, 0.94647932, 1.4411836, 1.9848944, 2.5478028, 3.1023506, 3.6243420, 4.0933800, 4.4929516, 4.8101783, 5.0356061, 5.1629266, 0, 0, 0, 0, 0, 0, 0}, 
  {0.037091282, 0.19274754, 0.46220949, 0.82882571, 1.2712045, 1.7654502, 2.2871545, 2.8129447, 3.3215189, 3.7942296, 4.2153082, 4.5718829, 4.8538432, 5.0536891, 5.1663907, 0, 0, 0, 0, 0, 0}, {0.032388259, 0.16863746, 0.40577952, 0.73109469, 1.1279298, 1.5772436, 2.0590586, 2.5536763, 3.0425811, 3.5090282, 3.9383677,4.3181690, 4.6382202, 4.8904493, 5.0688260, 5.1692886, 0, 0, 0, 0, 0}, {0.028525119, 0.14875702, 0.35893460, 0.64919620, 1.0064312, 1.4153629, 1.8596064, 2.3226078, 2.7883966, 3.2421292, 3.6704377, 4.0616209, 4.4057173, 4.6945021, 4.9214401, 5.0816239, 5.1717371, 0, 0, 0, 0}, 
  {0.025313388, 0.13217765, 0.31965410, 0.57999267, 0.90275957, 1.2756077, 1.6850505, 2.1171818, 2.5582874, 2.9953199, 3.4162370, 3.8102179, 4.1677814, 4.4808320, 4.7426590, 4.9479076, 5.0925409, 5.1738245, 0, 0, 0}, {0.022614580, 0.11821076, 0.28641445, 0.52105922, 0.81375611, 1.1544500, 1.5319896, 1.9346780, 2.3507619, 2.7688346, 3.1781440, 3.5688105, 3.9319662, 4.2598320, 4.5457492, 4.7841801, 4.9706907, 5.1019283, 5.1756184, 0, 0}, 
  {0.020325111, 0.10633743, 0.25805185, 0.47050624, 0.73689038, 1.0489539, 1.3974296, 1.7724547, 2.1639596, 2.5620029, 2.9570420, 3.3401386, 3.7031032, 4.0385882, 4.3401410, 4.6022262, 4.8202284, 4.9904423, 5.1100590, 5.1771714, 0}, {0.018366264, 0.096160905, 0.23366648, 0.42684850, 0.67012828, 0.95668652, 1.2787804, 1.6280672, 1.9959119, 2.3736625, 2.7528828, 3.1255378, 3.4841338, 3.8218182, 4.1324434, 4.4106051, 4.6516588, 4.8517235, 5.0076769, 5.1171477, 5.1785247}};
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }	
  }
  
  if ((vartheta_max==200)&&(vt_min==1)){printf("Entering roots max=%d\n",vartheta_max);
  double store2[20][21]={{3.9917095, 5.1074689, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1.3287047, 4.1958293, 5.1270284, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.69248467, 2.5232866, 4.3861649, 5.1482289, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.41485361, 1.7744766, 3.2396201, 4.5490557, 5.1690185, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.27587438, 1.2682534, 2.5587374, 3.7206455, 4.6860025, 5.1885521, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.19600960, 0.94224088, 2.0016421, 3.1132281, 4.0580289, 4.7986619, 5.2061619, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.14624670, 0.72153437, 1.5923201, 2.5722769, 3.5152669, 4.2988470, 4.8885232, 5.2213108, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.11316138, 0.56785567, 1.2853102, 2.1401868, 3.0133128, 3.8158337, 4.4746760, 4.9585642, 5.2338163, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.090108627, 0.45729995, 1.0538639, 1.7935337, 2.5875019, 3.3580041, 4.0466533, 4.6068999, 5.0128589, 5.2438792, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.073418372, 0.37555218, 0.87671753, 1.5167959, 2.2291936, 2.9516830, 3.6314049, 4.2277786, 4.7091889, 5.0553060, 5.2519074, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.060956941, 0.31358899, 0.73908379, 1.2947008, 1.9306027, 2.5971872, 3.2499509, 3.8514868, 4.3724511, 4.7901809, 5.0889802, 5.2583344, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.051411137, 0.26560257, 0.63050818, 1.1151581, 1.6819662, 2.2917471, 2.9074350, 3.4963541, 4.0310905, 4.4897790, 4.8555047, 5.1161140, 5.2635298, 0, 0, 0, 0, 0, 0, 0, 0},
  {0.043939666, 0.22773512, 0.54362032, 0.96873007, 1.4744486, 2.0298121, 2.6043142, 3.1698421, 3.7017919, 4.1794812, 4.5862133, 4.9089908, 5.1382945, 5.2677794, 0, 0, 0, 0, 0, 0, 0}, {0.037983408, 0.19735768, 0.47315891, 0.84820791, 1.3004784, 1.8054249, 2.3380426, 2.8744552, 3.3929625, 3.8746193, 4.3034443, 4.6664217, 4.9533458, 5.1566585, 5.2712967, 0, 0, 0, 0, 0, 0}, {0.033159361, 0.17263406, 0.41532043, 0.74810073, 1.1538288, 1.6129393, 2.1049518, 2.6097214, 3.1083600, 3.5838271, 4.0212484,4.4080295, 4.7338429, 4.9905382, 5.1720344, 5.2742400, 0, 0, 0, 0, 0}, {0.029198125, 0.15225344, 0.36731565, 0.66421854, 1.0294651, 1.4473554, 1.9010835, 2.3737063, 2.8489138, 3.3115800, 3.7481115, 4.1466329, 4.4970530, 4.7910529, 5.0220315, 5.1850373, 5.2767274, 0, 0, 0, 0}, {0.025905805, 0.13526124, 0.32707014, 0.59334612, 0.92334984, 1.3043910, 1.7226317, 2.1638300, 2.6139677, 3.0597345, 3.4888694, 3.8903751, 4.2546322, 4.5734403, 4.8400108, 5.0489330, 5.1961316, 5.2788485, 0, 0, 0},
  {0.023140001, 0.12094987, 0.29301995, 0.53299811, 0.83225148, 1.1804454, 1.5661353, 1.9773363, 2.4020262, 2.8285582, 3.2459730, 3.6442179, 4.0142843, 4.3482812, 4.6394636, 4.8822282, 5.0720934, 5.2056731, 5.2806717, 0, 0}, {0.020794247, 0.10878626, 0.26397061, 0.48123755, 0.75358037, 1.0725197, 1.4285433, 1.8115419, 2.2112101, 2.6173897, 3.0203459, 3.4109743, 3.7809460, 4.1227990, 4.4299889, 4.6969081, 4.9188857, 5.0921753, 5.2139388, 5.2822504, 0}, {0.018787664, 0.098362909, 0.23899871, 0.43654198, 0.68525429, 0.97812767, 1.3072139, 1.6639596, 2.0395238, 2.4250608, 2.8119580, 3.1920239, 3.5576303, 3.9018102, 4.2183219, 4.5016843, 4.7471921, 4.9509168, 5.1097005, 5.2211462, 5.2836263}};
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }	
  }
  
  if ((vartheta_max==220)&&(vt_min==1)){printf("Entering roots max=%d\n",vartheta_max);
  double store2[20][21]={{4.0866775, 5.2027599, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1.3729968, 4.2875333, 5.2220105, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.71375142, 2.5855574, 4.4762881, 5.2429643, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.42672755, 1.8179426, 3.3093643, 4.6384984, 5.2635728, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.28325546, 1.2992842, 2.6147402, 3.7951073, 4.7752840, 5.2829895, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.20100823, 0.96487034, 2.0461333, 3.1773384, 4.1362449, 4.8881953, 5.3005542, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.14983388, 0.73854387, 1.6277524, 2.6261420, 3.5850056, 4.3800995, 4.9785930, 5.3157259, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.11585254, 0.58098767, 1.3137992, 2.1853848, 3.0740132, 3.8896731, 4.5582631, 5.0492697, 5.3282994, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.092197356, 0.46768851, 1.0770519, 1.8315174, 2.6401848, 3.4238886, 4.1236191, 4.6922540, 5.1041657, 5.3384486, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.075084201, 0.38394615, 0.89584012, 1.5488908, 2.2748444, 3.0102232, 3.7013385, 4.3072132, 4.7959072, 5.1471294, 5.3465623, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.062315157, 0.32049658, 0.75505731, 1.3220112, 1.9702388, 2.6490562, 3.3131734, 3.9246563, 4.4538809, 4.8779844, 5.1812336, 5.3530662, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.052538975, 0.27137742, 0.64401230, 1.1385782, 1.7165054, 2.3377101, 2.9644044, 3.5633892, 4.1068950, 4.5728473, 4.9441935, 5.2087243, 5.3583276, 0, 0, 0, 0, 0, 0, 0, 0}, {0.044890694, 0.23262936, 0.55516301, 0.98897226, 1.5046866, 2.0706090, 2.6555997, 3.2310466, 3.7719805, 4.2574654, 4.6706440, 4.9984146, 5.2312031, 5.3626333, 0, 0, 0, 0, 0, 0, 0}, {0.038795889, 0.20155515, 0.48312386, 0.86583727, 1.3270871, 1.8417344, 2.3842319, 2.9302481, 3.4577246, 3.9474503, 4.3832555, 4.7519977, 5.0433881, 5.2498188, 5.3661983, 0, 0, 0, 0, 0, 0}, {0.033861327, 0.17627157, 0.42400098, 0.76356582, 1.1773678, 1.6453625, 2.1466114, 2.6605653, 3.1679996, 3.6516088, 4.0963187,4.4893891, 4.8203914, 5.0811063, 5.2654090, 5.3691822, 0, 0, 0, 0, 0},
  {0.029810562, 0.15543463, 0.37493885, 0.67787704, 1.0503979, 1.4764143, 1.9387364, 2.4200680, 2.9037919, 3.3745278, 3.8184803, 4.2236192, 4.5797377, 4.8784343, 5.1130503, 5.2785958, 5.3717046, 0, 0, 0, 0}, {0.026444728, 0.13806599, 0.33381389, 0.60548495, 0.94205963, 1.3305337, 1.7567486, 2.2061572, 2.6644658, 3.1181270, 3.5546830, 3.9629787, 4.3332717, 4.6572684, 4.9281115, 5.1403413, 5.2898490, 5.3738559, 0, 0, 0}, {0.023617835, 0.12344062, 0.29902523, 0.54384906, 0.84905557, 1.2040543, 1.5971332, 2.0160449, 2.4485236, 2.8827051, 3.3074436, 3.7125309, 4.0888325, 4.4283580, 4.7242856, 4.9709543, 5.1638406, 5.2995289, 5.3757054, 0, 0},
  {0.021220780, 0.11101247, 0.26935038, 0.49098925, 0.76874225, 1.0939205, 1.4567880, 1.8470109, 2.2540696, 2.6676098, 3.0777228, 3.4751549, 3.8514523, 4.1990506, 4.5113236, 4.7825997, 5.0081591, 5.1842190, 5.3079157, 5.3773071, 0}, {0.019170705, 0.10036429, 0.24384440, 0.44534910, 0.69899362, 0.99759734, 1.3330244, 1.6965296, 2.0790845, 2.4716679, 2.8655074, 3.2522710, 3.6242088, 3.9742520, 4.2960744, 4.5841266, 4.8336492, 5.0406720, 5.2020054, 5.3152299, 5.3787033}};
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }	
  }
  
  if ((vartheta_max==240)&&(vt_min==1)){printf("Entering roots max=%d\n",vartheta_max);
  double store2[20][21]={{4.1734180, 5.2897562, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1.4141676, 4.3713969, 5.3087358, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.73355619, 2.6428778, 4.5587143, 5.3294694, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.43775600, 1.8580282, 3.3733700, 4.7202910, 5.3499148, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.29009526, 1.3279094, 2.6661731, 3.8633163, 4.8569084, 5.3692245, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.20563093, 0.98573863, 2.0870280, 3.2360971, 4.2078018, 4.9700214, 5.3867453, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.15314618, 0.75422063, 1.6603281, 2.6755515, 3.6488528, 4.4543810, 5.0608818, 5.4019337, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.11833445, 0.59308310, 1.3399918, 2.2268604, 3.1296218, 3.9572277, 4.6346511, 5.1321196, 5.4145655, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.094121769, 0.47725115, 1.0983673, 1.8663806, 2.6884690, 3.4841958, 4.1939997, 4.7702402, 5.1875529, 5.4247904, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.076617724, 0.39166832, 0.91341407, 1.5783506, 2.3166953, 3.0638297, 3.7653152, 4.3798255, 4.8751276, 5.2309812, 5.4329805, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.063564638, 0.32684801, 0.76973299, 1.3470781, 2.0065813, 2.6965674, 3.3710317, 3.9915656, 4.5282978, 4.9581865, 5.2654733, 5.4395533, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.053575914, 0.27668481, 0.65641553, 1.1600720, 1.7481763, 2.3798196, 3.0165545, 3.6247087, 4.1761925, 4.6487466, 5.0251964, 5.2932864, 5.4448743, 0, 0, 0, 0, 0, 0, 0, 0}, {0.045764632, 0.23712553, 0.56576161, 1.0075469, 1.5324135, 2.1079897, 2.7025564, 3.2870469, 3.8361622, 4.3287388, 4.7477768, 5.0800828, 5.3160349, 5.4492307, 0, 0, 0, 0, 0, 0, 0}, {0.039542185, 0.20540978, 0.49227120, 0.88201162, 1.3514849, 1.8750052, 2.4265285, 2.9813072, 3.5169581, 4.0140305, 4.4561854, 4.8301680, 5.1256161, 5.3348783, 5.4528387, 0, 0, 0, 0, 0, 0}, 
  {0.034505863, 0.17961084, 0.43196720, 0.77775207, 1.1989492, 1.6750726, 2.1847634, 2.7071021, 3.2225584, 3.7135863, 4.1649314, 4.5637233, 4.8994429, 5.1638104, 5.3506625, 5.4558594, 0, 0, 0, 0, 0}, {0.030372703, 0.15835412, 0.38193300, 0.69040390, 1.0695880, 1.5030411, 1.9732208, 2.4625069, 2.9540023, 3.4320955, 3.8828084, 4.2939711, 4.6552737, 4.9582410, 5.1961629, 5.3640158, 5.4584134, 0, 0, 0, 0}, {0.026939242, 0.14063929, 0.33999974, 0.61661613, 0.95920994, 1.3544874, 1.7879951, 2.2449057, 2.7106741, 3.1715366, 3.6148569, 4.0293372, 4.4051246, 4.7338418, 5.0085704, 5.2238069, 5.3754131, 5.4605921, 0, 0, 0},
  {0.024056178, 0.12572525, 0.30453253, 0.55379761, 0.86445728, 1.2256852, 1.6255229, 2.0514824, 2.4910746, 2.9322371, 3.3636544, 3.7749775, 4.1569579, 4.5015160, 4.8017609, 5.0519803, 5.2476134, 5.3852185, 5.4624655, 0, 0}, {0.021611965, 0.11305401, 0.27428306, 0.49992855, 0.78263719, 1.1135269, 1.4826555, 1.8794831, 2.2932938, 2.7135539, 3.1301962, 3.5338317, 3.9158930, 4.2687237, 4.5856238, 4.8608644, 5.0896817, 5.2682608, 5.3937152, 5.4640881, 0}, {0.019521925, 0.10219927, 0.24828658, 0.45342128, 0.71158343, 1.0154332, 1.3566619, 1.7263478, 2.1152911, 2.5143095, 2.9144851, 3.3073577, 3.6850673, 4.0404525, 4.3671113, 4.6594327, 4.9126088, 5.1226318, 5.2862841, 5.4011260, 5.4655027}};
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }	
  }
  
  if ((vartheta_max==260)&&(vt_min==1)){printf("Entering roots max=%d\n",vartheta_max);
  double store2[20][21]={{4.2532423, 5.3697867, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1.4526637, 4.4486599, 5.3885262, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.75210839, 2.6960010, 4.6346611, 5.4090614, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.44806269, 1.8952450, 3.4325277, 4.7956456, 5.4293586, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.29647408, 1.3544933, 2.7137449, 3.9262558, 4.9320916, 5.4485698, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.20993408, 1.0051134, 2.1248817, 3.2903435, 4.2737551, 5.0453677, 5.4660480, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.15622520, 0.76876792, 1.6904883, 2.7212014, 3.7077378, 4.5228006, 5.1366316, 5.4812484, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.12063895, 0.60430064, 1.3642426, 2.2651941, 3.1809394, 4.0194930, 4.7049873, 5.2083698, 5.4939307, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.095907002, 0.48611474, 1.1180996, 1.8986093, 2.7330448, 3.5398069, 4.2588410, 4.8420339, 5.2642875, 5.5042230, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.078039267, 0.39882222, 0.92967918, 1.6055856, 2.3553420, 3.1132804, 3.8242791, 4.4467014, 4.9480475, 5.3081371, 5.5124818, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.064722147, 0.33272922, 0.78331215, 1.3702512, 2.0401460, 2.7404068, 3.4243746, 4.0532091, 4.5968193, 5.0320020, 5.3429821, 5.5191170, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.054536008, 0.28159715, 0.66788893, 1.1799400, 1.7774282, 2.4186817, 3.0646469, 3.6812191, 4.2400186, 4.7186210, 5.0997424, 5.3710890, 5.5244921, 0, 0, 0, 0, 0, 0, 0, 0}, {0.046573431, 0.24128544, 0.57556304, 1.0247142, 1.5580226, 2.1424913, 2.7458674, 3.3386673, 3.8952916, 4.3943707, 4.8187772, 5.1552358, 5.3940835, 5.5288946, 0, 0, 0, 0, 0, 0, 0}, {0.040232578, 0.20897490, 0.50072836, 0.89695832, 1.3740182, 1.9057153, 2.4655463, 3.0283816, 3.5715405, 4.0753545, 4.5233315, 4.9021159, 5.2012802, 5.4131346, 5.5325419, 0, 0, 0, 0, 0, 0}, {0.035101909, 0.18269834, 0.43933057, 0.79085949, 1.2188800, 1.7024963, 2.2199609, 2.7500131, 3.2728420, 3.7706821, 4.2281153, 4.6321532, 4.9721960, 5.2399093, 5.4290959, 5.5355962, 0, 0, 0, 0, 0}, {0.030892392, 0.16105273, 0.38839641, 0.70197627, 1.0873087, 1.5276183, 2.0050361, 2.5016431, 3.0002846, 3.4851376, 3.9420571, 4.3587466, 4.7248028, 5.0316842, 5.2726349, 5.4426013, 5.5381791, 0, 0, 0, 0},
  {0.027396286, 0.14301731, 0.34571498, 0.62689756, 0.97504557, 1.3765964, 1.8168236, 2.2806411, 2.7532719, 3.2207540, 3.6702879, 4.0904456, 4.4712738, 4.8043197, 5.0826097, 5.3006015, 5.4541302, 5.5403828, 0, 0, 0}, {0.024461205, 0.12783603, 0.30961985, 0.56298528, 0.87867688, 1.2456492, 1.6517156, 2.0841654, 2.5303039, 2.9778862, 3.4154413, 3.8324914, 4.2196845, 4.5688597, 4.8730639, 5.1265385, 5.3246888, 5.4640502, 5.5422780, 0, 0}, {0.021973339, 0.11493981, 0.27883877, 0.50818295, 0.79546431, 1.1316213, 1.5065207, 1.9094320, 2.3294580, 2.7558998, 3.1785450, 3.5878801, 3.9752346, 4.3328679, 4.6540132, 4.9328896, 5.1646940, 5.3455822, 5.4726473, 5.5439197, 0}, 
  {0.019846314, 0.10389394, 0.25238858, 0.46087397, 0.72320451, 1.0318924, 1.3784689, 1.7538490, 2.1486742, 2.5536140, 2.9596167, 3.3581043, 3.7411164, 4.1014070, 4.4325047, 4.7287430, 4.9852700, 5.1980440, 5.3638221, 5.4801467, 5.5453511}};
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }	
  }
  
  if ((vartheta_max==280)&&(vt_min==1)){printf("Entering roots max=%d\n",vartheta_max);
  double store2[20][21]={{4.3271712, 5.4438846, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1.4888396, 4.5202890, 5.4624093, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.76957406, 2.7455185, 4.7050778, 5.4827643, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.45774536, 1.9299951, 3.4875345, 4.8655068, 5.5029260, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.30245520, 1.3793222, 2.7580099, 3.9846927, 5.0017809, 5.5220460, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.21396203, 1.0232048, 2.1601291, 3.3407327, 4.3349269, 5.1151895, 5.5394831, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.15910355, 0.78234516, 1.7185780, 2.7636351, 3.7623856, 4.5862216, 5.2068088, 5.5546920, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.12279104, 0.61476466, 1.3868292, 2.3008393, 3.2285908, 4.0772446, 4.7701646, 5.2789967, 5.5674183, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.097572756, 0.49437866, 1.1364753, 1.9285834, 2.7744510, 3.5914089, 4.3189577, 4.9085498, 5.3353545, 5.5777710, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.079364748, 0.40548888, 0.94482304, 1.6309168, 2.3912494, 3.1591818, 3.8789657, 4.5086860, 5.0155982, 5.3795891, 5.5860920, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.065800804, 0.33820745, 0.79595220, 1.3918036, 2.0713358, 2.7811098, 3.4738632, 4.1103612, 4.6603154, 5.1003756, 5.4147575, 5.5927842, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.055430256, 0.28617109, 0.67856619, 1.1984169, 1.8046119, 2.4547694, 3.1092749, 3.7336263, 4.2991790, 4.7833603, 5.1687872, 5.4431340, 5.5982089, 0, 0, 0, 0, 0, 0, 0, 0}, {0.047326437, 0.24515739, 0.58468211, 1.0406776, 1.5818210, 2.1745331, 2.7860655, 3.3865501, 3.9501118, 4.4551931, 4.8845517, 5.2248383, 5.4663543, 5.6026538, 0, 0, 0, 0, 0, 0, 0}, {0.040875108, 0.21229217, 0.50859488, 0.91085489, 1.3949574, 1.9342371, 2.5017640, 3.0720549, 3.6221551, 4.1321965, 4.5855476, 4.9687618, 5.2713526, 5.4855964, 5.6063371, 0, 0, 0, 0, 0, 0}, 
  {0.035656451, 0.18557038, 0.44617816, 0.80304424, 1.2373996, 1.7279663, 2.2526349, 2.7898289, 3.3194778, 3.8236145, 4.2866709, 4.6955514, 5.0395826, 5.3103815, 5.5017204, 5.6094223, 0, 0, 0, 0, 0}, {0.031375755, 0.16356238, 0.39440583, 0.71273245, 1.1037736, 1.5504443, 2.0345717, 2.5379596, 3.0432148, 3.5343193, 3.9969748, 4.4187687, 4.7892130, 5.0997059, 5.3434504, 5.5153657, 5.6120317, 0, 0, 0, 0}, {0.027821274, 0.14522830, 0.35102775, 0.63645245, 0.98975753, 1.3971294, 1.8435870, 2.3138039, 2.7927883, 3.2663951, 3.7216740, 4.1470780, 4.5325617, 4.8696038, 5.1511803, 5.3717138, 5.5270157, 5.6142584, 0, 0, 0}, 
  {0.024837739, 0.12979814, 0.31434807, 0.57152253, 0.89188623, 1.2641892, 1.6760319, 2.1144968, 2.5666981, 3.0202224, 3.4634550, 3.8857996, 4.2778091, 4.6312485, 4.9391080, 5.1955867, 5.3960596, 5.5370412, 5.6161736, 0, 0}, {0.022309222, 0.11669244, 0.28307219, 0.51585196, 0.80737896, 1.1484240, 1.5286759, 1.9372265, 2.3630103, 2.7951755, 3.2233754, 3.6379818, 4.0302291, 4.3922998, 4.7173659, 4.9995991, 5.2341605, 5.4171795, 5.5457308, 5.6178329, 0}, {0.020147765, 0.10546867, 0.25619980, 0.46779724, 0.73399785, 1.0471758, 1.3987128, 1.7793720, 2.1796473, 2.5900710, 3.0014673, 3.4051498, 3.7930650, 4.1578895, 4.4930886, 4.7929446, 5.0525655, 5.2678790, 5.4356190, 5.5533116, 5.6192798}};
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }	
  }
  
  if ((vartheta_max==300)&&(vt_min==1)){printf("Entering roots max=%d\n",vartheta_max);
  double store2[20][21]={{4.3960153, 5.5128692, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1.5229817, 4.5870528, 5.5312000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.78608721, 2.7919037, 4.7707187, 5.5513903, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.46688260, 1.9626002, 3.5389459, 4.9306250, 5.5714273, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.30808939, 1.4026250, 2.7994101, 4.0392371, 5.0667282, 5.5904623, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.21775034, 1.0401806, 2.1931169, 3.3877860, 4.3919702, 5.1802450, 5.6078599, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.16180739, 0.79507976, 1.7448726, 2.8032855, 3.8133720, 4.6453293, 5.2721798, 5.6230745, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.12481068, 0.62457451, 1.4079729, 2.3341570, 3.2730727, 4.1310986, 4.8308915, 5.3447749, 5.6358395, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.099134769, 0.50212224, 1.1536752, 1.9566056, 2.8131163, 3.6395475, 4.3749962, 4.9705136, 5.4015352, 5.6462467, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.080606875, 0.41173296, 0.95899520, 1.6545997, 2.4247874, 3.2020162, 3.9299592, 4.5664502, 5.0785186, 5.4461238, 5.6546246, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.066811080, 0.34333633, 0.80777861, 1.4119530, 2.1004713, 2.8191018, 3.5200231, 4.1636367, 4.7194763, 5.1640566, 5.4815904, 5.6613691, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.056267424, 0.29045177, 0.68855385, 1.2156896, 1.8300065, 2.4884587, 3.1509101, 3.7824911, 4.3543135, 4.8436703, 5.2330886, 5.5102157, 5.6668394, 0, 0, 0, 0, 0, 0, 0, 0}, {0.048031100, 0.24877990, 0.59321027, 1.0555990, 1.6040531, 2.2044482, 2.8235739, 3.4312053, 4.0012127, 4.5118663, 4.9458192, 5.2896555, 5.5336448, 5.6713234, 0, 0, 0, 0, 0, 0, 0}, {0.041476179, 0.21539480, 0.51595005, 0.92384268, 1.4145178, 1.9608672, 2.5355622, 3.1127908, 3.6693445, 4.1851712, 4.6435113, 5.0308357, 5.3366042, 5.5530633, 5.6750401, 0, 0, 0, 0, 0, 0}, 
  {0.036175056, 0.18825588, 0.45257932, 0.81443063, 1.2546988, 1.7517473, 2.2831285, 2.8269711, 3.3629641, 3.8729537, 4.3412335, 4.7546097, 5.1023420, 5.3760031, 5.5693379, 5.6781537, 0, 0, 0, 0, 0}, {0.031827672, 0.16590846, 0.40002238, 0.72278250, 1.1191522, 1.5717563, 2.0621374, 2.5718406, 3.0832507, 3.5801689, 4.0481552, 4.4746905, 4.8492088, 5.1630533, 5.4093896, 5.5831126, 5.6807878, 0, 0, 0, 0}, {0.028218519, 0.14729474, 0.35599230, 0.64537889, 1.0034978, 1.4162998, 1.8685657, 2.3447444, 2.8296439, 3.3089488, 3.7695695, 4.1998488, 4.5896569, 4.9304090, 5.2150360, 5.4379278, 5.5948747, 5.6830357, 0, 0, 0}, 
  {0.025189618, 0.13163162, 0.31876565, 0.57949721, 0.90422202, 1.2814981, 1.6987267, 2.1427966, 2.6006440, 3.0596985, 3.5082121, 3.9354788, 4.3319640, 4.6893641, 5.0006175, 5.2598847, 5.4625130, 5.6049978, 5.6849695, 0, 0}, {0.022623052, 0.11832988, 0.28702687, 0.52301473, 0.81850468, 1.1641102, 1.5493533, 1.9631597, 2.3943070, 2.8318005, 3.2651688, 3.6846776, 4.0814732, 4.4476671, 4.7763750, 5.0617253, 5.2988461, 5.4838427, 5.6137730, 5.6866451, 0}, {0.020429375, 0.10693967, 0.25975957, 0.47426273, 0.74407561, 1.0614427, 1.4176059, 1.8031859, 2.2085390, 2.6240693, 3.0404858, 3.4490010, 3.8414756, 4.2105146, 4.5495244, 4.8527408, 5.1152349, 5.3329059, 5.5024670, 5.6214293, 5.6881064}};
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }	
  }
  
  if ((vartheta_max==320)&&(vt_min==1)){printf("Entering roots max=%d\n",vartheta_max);
  double store2[20][21]={{4.4604286, 5.5774007, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1.5553249, 4.6495723, 5.5955556, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.80175770, 2.8355412, 4.8321935, 5.6155943, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.47553864, 1.9933219, 3.5872120, 4.9916065, 5.6355159, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.31341807, 1.4245878, 2.8383034, 4.0903823, 5.1275409, 5.6544716, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.22132792, 1.0561775, 2.2241263, 3.4319251, 4.4454121, 5.2411460, 5.6718310, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.16435794, 0.80707523, 1.7695954, 2.8405029, 3.8611628, 4.7006766, 5.3333629, 5.6870491, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.12671411, 0.63381084, 1.4278532, 2.3654398, 3.3147870, 4.1815524, 4.8877396, 5.4063289, 5.6998484, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.10060582, 0.50940985, 1.1698456, 1.9829209, 2.8493874, 3.6846638, 4.4274787, 5.0285105, 5.4634595, 5.7103051, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.081775958, 0.41760691, 0.97231687, 1.6768410, 2.4562554, 3.2421732, 3.9777313, 4.6205352, 5.1374045, 5.5083755, 5.7187352, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.067761463, 0.34815934, 0.81889303, 1.4308755, 2.1278116, 2.8547269, 3.5632786, 4.2135319, 4.7748588, 5.2236494, 5.5441184, 5.7255281, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.057054621, 0.29447577, 0.69793823, 1.2319094, 1.8538376, 2.5200537, 3.1899336, 3.8282663, 4.4059385, 4.9001208, 5.2932580, 5.5729747, 5.7310407, 0, 0, 0, 0, 0, 0, 0, 0}, {0.048693454, 0.25218414, 0.60122161, 1.0696093, 1.6249166, 2.2325061, 2.8587348, 3.4730449, 4.0490707, 4.5649231, 5.0031599, 5.3503042, 5.5965976, 5.7355609, 0, 0, 0, 0, 0, 0, 0}, {0.042040979, 0.21830970, 0.52285808, 0.93603607, 1.4328735, 1.9858451, 2.5672485, 3.1509639, 3.7135468, 4.2347743, 4.6977691, 5.0889262, 5.3976569, 5.6161802, 5.7393085, 0, 0, 0, 0, 0, 0}, 
  {0.036662228, 0.19077825, 0.45859017, 0.82511924, 1.2709316, 1.7740532, 2.3117185, 2.8617806, 3.4037035, 3.9191600, 4.3923159, 4.8098865, 5.1610703, 5.4373997, 5.6325947, 5.7424487, 0, 0, 0, 0, 0}, {0.032252093, 0.16811154, 0.40529550, 0.73221548, 1.1335818, 1.5917461, 2.0879834, 2.6035961, 3.1207617, 3.6231128, 4.0960779, 4.5270390, 4.9053583, 5.2223288, 5.4710817, 5.6464898, 5.7451056, 0, 0, 0, 0}, {0.028591512, 0.14923484, 0.36065251, 0.65375620, 1.0163892, 1.4342804, 1.8919865, 2.3737455, 2.8641782, 3.3488100, 3.8144218, 4.2492539, 4.6430983, 4.9873123, 5.2747845, 5.4998756, 5.6583561, 5.7473733, 0, 0, 0}, 
  {0.025519950, 0.13335268, 0.32291177, 0.58698040, 0.91579481, 1.2977321, 1.7200059, 2.1693235, 2.6324538, 3.0966800, 3.5501295, 3.9819946, 4.3826592, 4.7437564, 5.0581766, 5.3200451, 5.5246836, 5.6685702, 5.7493244, 0, 0}, {0.022917612, 0.11986667, 0.29073800, 0.52973527, 0.82894137, 1.1788215, 1.5687407, 1.9874686, 2.4236355, 2.8661134, 3.3043141, 3.7284042, 4.1294485, 4.4994925, 4.8315999, 5.1198588, 5.3593673, 5.5462086, 5.6774251, 5.7510151, 0}, {0.020693650, 0.10832003, 0.26309966, 0.48032837, 0.75352845, 1.0748223, 1.4353199, 1.8255083, 2.2356147, 2.6559230, 3.0770343, 3.4900672, 3.8868023, 4.2597778, 4.6023460, 4.9086991, 5.1738744, 5.3937451, 5.5650051, 5.6851517, 5.7524898}};
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }	
  }
  
  if ((vartheta_max==340)&&(vt_min==1)){printf("Entering roots max=%d\n",vartheta_max);
  double store2[20][21]={{4.5209469, 5.6380194, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1.5860647, 4.7083565, 5.6560135, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.81667681, 2.8767484, 4.8900018, 5.6759120, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.48376667, 2.0223759, 3.6327030, 5.0489478, 5.6957263, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.31847548, 1.4453643, 2.8749842, 4.1385329, 5.1847161, 5.7146076, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.22471867, 1.0713080, 2.2533885, 3.4734960, 4.4956849, 5.2983936, 5.7319303, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.16677272, 0.81841681, 1.7929299, 2.8755749, 3.9061397, 4.7527165, 5.3908642, 5.7471501, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.12851468, 0.64254003, 1.4466177, 2.3949278, 3.3540632, 4.2290138, 4.9411768, 5.4641700, 5.7599801, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.10199643, 0.51629445, 1.1851070, 2.0077304, 2.8835486, 3.7271192, 4.4768328, 5.0830198, 5.5216428, 5.7704820, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.082880483, 0.42315385, 0.98488768, 1.6978107, 2.4858989, 3.2799725, 4.0226687, 4.6713843, 5.1927438, 5.5668629, 5.7789603, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.068658935, 0.35271223, 0.82937902, 1.4487158, 2.1535697, 2.8882670, 3.6039775, 4.2604533, 4.8269188, 5.2796485, 5.6028633, 5.7857980, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.057797692, 0.29827317, 0.70679022, 1.2472005, 1.8762908, 2.5498038, 3.2266577, 3.8713227, 4.4544766, 4.9531779, 5.3497954, 5.6319350, 5.7913499, 0, 0, 0, 0, 0, 0, 0, 0}, {0.049318462, 0.25539575, 0.60877696, 1.0828163, 1.6445737, 2.2589278, 2.8918287, 3.5124065, 4.0940757, 4.6147996, 5.0570484, 5.4072891, 5.6557388, 5.7959039, 0, 0, 0, 0, 0, 0, 0}, {0.042573772, 0.22105894, 0.52937167, 0.94752899, 1.4501673, 2.0093675, 2.5970749, 3.1868812, 3.7551208, 4.2814121, 4.7487685, 5.1435153, 5.4550190, 5.6754745, 5.7996804, 0, 0, 0, 0, 0, 0}, 
  {0.037121671, 0.19315673, 0.46425680, 0.83519264, 1.2862246, 1.7950594, 2.3386322, 2.8945363, 3.4420256, 3.9626103, 4.4403376, 4.8618386, 5.2162554, 5.4950832, 5.6920197, 5.8028454, 0, 0, 0, 0, 0}, {0.032652263, 0.17018850, 0.41026581, 0.74110447, 1.1471752, 1.6105711, 2.1123148, 2.6334803, 3.1560507, 3.6635003, 4.1411351, 4.5762453, 4.9581266, 5.2780253, 5.5290413, 5.7060273, 5.8055236, 0, 0, 0, 0}, {0.028943119, 0.15106352, 0.36504442, 0.66164948, 1.0285326, 1.4512128, 1.9140351, 2.4010391, 2.8966693, 3.3863021, 3.8565970, 4.2956991, 4.6933275, 5.0407856, 5.3309233, 5.5580743, 5.7179910, 5.8078098, 0, 0, 0}, 
  {0.025831283, 0.13497462, 0.32681861, 0.59403043, 0.92669528, 1.3130191, 1.7400385, 2.1942893, 2.6623835, 3.1314664, 3.5895490, 4.0257282, 4.4303122, 4.7948753, 5.1122634, 5.3765690, 5.5830906, 5.7282900, 5.8097770, 0, 0}, {0.023195183, 0.12131474, 0.29423448, 0.53606611, 0.83877097, 1.1926741, 1.5869920, 2.0103474, 2.4512317, 2.8983916, 3.3411294, 3.7695192, 4.1745493, 4.5482037, 4.8834980, 5.1744830, 5.4162288, 5.6047982, 5.7372194, 5.8114819, 0}, {0.020942646, 0.10962051, 0.26624617, 0.48604168, 0.76243074, 1.0874201, 1.4519955, 1.8465176, 2.2610918, 2.6858892, 3.1114097, 3.5286837, 3.9294167, 4.3060849, 4.6519900, 4.9612838, 5.2289721, 5.4509041, 5.6237557, 5.7450117, 5.8129691}};
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }	
  }
  
  if ((vartheta_max==360)&&(vt_min==1)){printf("Entering roots max=%d\n",vartheta_max);
  double store2[20][21]={{4.5780140, 5.6951727, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1.6153652, 4.7638279, 5.7130191, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.83092126, 2.9157898, 4.9445582, 5.7327875, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.49161129, 2.0499424, 3.6757266, 5.1030610, 5.7525015, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.32329027, 1.4650827, 2.9096973, 4.1840249, 5.2386666, 5.7713128, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.22794255, 1.0856659, 2.2810961, 3.5127863, 4.5431468, 5.3524030, 5.7886002, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.16906638, 0.82917561, 1.8150290, 2.9087406, 3.9486195, 4.8018245, 5.4451031, 5.8038202, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.13022356, 0.65081744, 1.4643894, 2.4228206, 3.3911749, 4.2738211, 4.9915913, 5.5187217, 5.8166777, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.10331536, 0.52282018, 1.1995597, 2.0312015, 2.9158363, 3.7672139, 4.5234131, 5.1344386, 5.5765123, 5.8272212, 0, 0,0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.083927519, 0.42840972, 0.99679054, 1.7176502, 2.5139220, 3.3156793, 4.0650918, 4.7193651, 5.2449407, 5.6220161, 5.8357440, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.069509312, 0.35702475, 0.83930607, 1.4655942, 2.1779224, 2.9199566, 3.6424084, 4.3047380, 4.8760342, 5.3324639, 5.6582571, 5.8426235, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.058501500, 0.30186900, 0.71516881, 1.2616663, 1.8975198, 2.5779160, 3.2613416, 3.9119681, 4.5002783, 5.0032280, 5.4031152, 5.6875306, 5.8482122, 0, 0, 0, 0, 0, 0, 0, 0}, {0.049910251, 0.25843606, 0.61592693, 1.0953092, 1.6631592, 2.2838969, 2.9230882, 3.5495702, 4.1365513, 4.6618574, 5.1078780, 5.4610285, 5.7115039, 5.8527978, 0, 0, 0, 0, 0, 0, 0}, {0.043078103, 0.22366090, 0.53553467, 0.95839948, 1.4665180, 2.0315976, 2.6252509, 3.2207973, 3.7943641, 4.3254212, 4.7968803, 5.1950021, 5.5091122, 5.7313831, 5.8566014, 0, 0, 0, 0, 0, 0}, 
  {0.037556461, 0.19540729, 0.46961751, 0.84471946, 1.3006829, 1.8149119, 2.3640581, 2.9254702, 3.4782038, 4.0036173, 4.4856466, 4.9108447, 5.2683014, 5.5494778, 5.7480509, 5.8597896, 0, 0, 0, 0, 0}, {0.033030877, 0.17215338, 0.41496705, 0.74951021, 1.1600259, 1.6283619, 2.1353020, 2.6617044, 3.1893687, 3.7016211, 4.1836524, 4.6226671, 5.0078989, 5.3305513, 5.5836949, 5.7621639, 5.8624878, 0, 0, 0, 0}, {0.029275721, 0.15279321, 0.36919795, 0.66911282, 1.0400118, 1.4672146, 1.9348660, 2.4268177, 2.9273481, 3.4216934, 3.8963991, 4.3395208, 4.7407101, 5.0912200, 5.3838645, 5.6129523, 5.7742190, 5.8647914, 0, 0, 0}, {0.026125737, 0.13650853, 0.33051290, 0.60069575, 0.93699875, 1.3274655, 1.7589646, 2.2178700, 2.6906453, 3.1643060, 3.6267537, 4.0669957, 4.4752692, 4.8430938, 5.1632740, 5.4298718, 5.6381641, 5.7845978, 5.8667738, 0, 0}, 
  {0.023457665, 0.12268399, 0.29754031, 0.54205087, 0.84806155, 1.2057643, 1.6042351, 2.0319572, 2.4772912, 2.9288653, 3.3758789, 3.8083189, 4.2171023, 4.5941554, 4.9324487, 5.2259984, 5.4698485, 5.6600431, 5.7935972, 5.8684919, 0}, {0.021178072, 0.11085006, 0.26922076, 0.49144215, 0.77084424, 1.0993242, 1.4677495, 1.8663616, 2.2851508, 2.7141815, 3.1438581, 3.5651282, 3.9696270, 4.3497723, 4.6988184, 5.0108796, 5.2809322, 5.5048031, 5.6791518, 5.8014512, 5.8699908}};
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }	
  }
  
  if ((vartheta_max==380)&&(vt_min==1)){printf("Entering roots max=%d\n",vartheta_max);
  double store2[20][21]={{4.6320023, 5.7492355, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1.6433658, 4.8163412, 5.7669457, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.84455620, 2.9528886, 4.9962107, 5.7865926, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.49911025, 2.0761733, 3.7165418, 5.1542920, 5.8062127, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.32788667, 1.4838509, 2.9426487, 4.2271407, 5.2897382, 5.8249579, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.23101641, 1.0993304, 2.3074111, 3.5500373, 4.5880985, 5.4035225, 5.8422112, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.17125123, 0.83941153, 1.8360214, 2.9402010, 3.9888683, 4.8483154, 5.4964310, 5.8574303, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.13185015, 0.65868968, 1.4812716, 2.4492861, 3.4263519, 4.3162582, 5.0393082, 5.5703389, 5.8703126, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.10457001, 0.52902421, 1.2132879, 2.0534749, 2.9464489, 3.8051995, 4.5675169, 5.1830998, 5.6284257, 5.8808943, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.084923017, 0.43340474, 1.0080952, 1.7364782, 2.5404962, 3.3495164, 4.1052698, 4.7647857, 5.2943338, 5.6741951, 5.8894587, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.070317488, 0.36112194, 0.84873267, 1.4816118, 2.2010183, 2.9499922, 3.6788136, 4.3466688, 4.9225218, 5.3824388, 5.7106622, 5.8963773, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.059170139, 0.30528432, 0.72312360, 1.2753935, 1.9176542, 2.6045642, 3.2942029, 3.9504606, 4.5436377, 5.0505949, 5.4535648, 5.7401251, 5.9020005, 0, 0, 0, 0, 0, 0, 0, 0}, {0.050472294, 0.26132300, 0.62271405, 1.1071634, 1.6807864, 2.3075675, 2.9527089, 3.5847709, 4.1767688, 4.7063997, 5.1559786, 5.5118728, 5.7642578, 5.9066158, 0, 0, 0, 0, 0, 0, 0}, {0.043556955, 0.22613104, 0.54138393, 0.96871316, 1.4820252, 2.0526725, 2.6519519, 3.2529260, 3.8315265, 4.3670841, 4.8424153, 5.2437212, 5.5602894, 5.7842721, 5.9104449, 0, 0, 0, 0, 0, 0}, 
  {0.037969188, 0.19754339, 0.47470450, 0.85375736, 1.3143948, 1.8337329, 2.3881545, 2.9547766, 3.5124676, 4.0424431, 4.5285346, 4.9572222, 5.3175471, 5.6009388, 5.8010552, 5.9136549, 0, 0, 0, 0, 0}, {0.033390203, 0.17401797,0.41942758, 0.75748373, 1.1722125, 1.6452283, 2.1570880, 2.6884454, 3.2209267, 3.7377182, 4.2239026, 4.6666041, 5.0549984, 5.3802491, 5.6353998, 5.8152676, 5.9163719, 0, 0, 0, 0}, {0.029591321, 0.15443435, 0.37313826, 0.67619168, 1.0508970, 1.4823846, 1.9546086, 2.4512428, 2.9564084, 3.4552088, 3.9340826, 4.3810013, 4.7855528, 5.1389433, 5.4339533, 5.6648684, 5.8274087, 5.9186919, 0, 0, 0}, 
  {0.026405092, 0.13796367, 0.33401710, 0.60701706, 0.94676849, 1.3411604, 1.7769019, 2.2402133, 2.7174175, 3.1954074, 3.6619811, 4.1060619, 4.5178203, 4.8887245, 5.2115403, 5.4803014, 5.6902643, 5.8378626, 5.9206885, 0, 0}, {0.023706650, 0.12398277, 0.30067566, 0.54772619, 0.85687025, 1.2181731, 1.6205771, 2.0524334, 2.5019780, 2.9577276, 3.4087839, 3.8450521, 4.2573815, 4.6376446, 4.9787697, 5.2747404, 5.5205767, 5.7123049, 5.8469278, 5.9224192, 0}, {0.021401363, 0.11201617, 0.27204163, 0.49656292, 0.77882080, 1.1106081, 1.4826800, 1.8851646, 2.3079430, 2.7409788, 3.1745861, 3.5996340, 4.0076916, 4.3911218, 4.7431346, 5.0578089, 5.3300935, 5.5557946, 5.7315560, 5.8548400, 5.9239291}};
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }		
  }
  
  if ((vartheta_max==400)&&(vt_min==1)){printf("Entering roots max=%d\n",vartheta_max);
  double store2[20][21]={{4.6832268, 5.8005250, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1.6701862, 4.8661968, 5.8181089, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.85763746, 2.9882349, 5.0452541, 5.8376420, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.50629586, 2.1011982, 3.7553686, 5.2029336, 5.8571738, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.33228535, 1.5017611, 2.9740133, 4.2681194, 5.3382240, 5.8758565, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.23395463, 1.1123688, 2.3324713, 3.5854541, 4.6307947, 5.4520467, 5.8930770, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.17333780, 0.84917559, 1.8560164, 2.9701266, 4.0271116, 4.8924559, 5.5451456, 5.9082941, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.13340245, 0.66619642, 1.4973522, 2.4744667, 3.4597891, 4.3565654, 5.0846028, 5.6193219, 5.9211989, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.10576666, 0.53493812, 1.2263633, 2.0746701, 2.9755547, 3.8412896, 4.6093963, 5.2292851, 5.6776858, 5.9318161, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.085872039, 0.43816463, 1.0188610, 1.7543956, 2.5657667, 3.3816725, 4.1434304, 4.8079072, 5.3412099, 5.7237049, 5.9404193, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.071087622, 0.36502509, 0.85770846, 1.4968545, 2.2229834, 2.9785406, 3.7133983, 4.3864852, 4.9666499, 5.4298641, 5.7603848, 5.9473745, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.059807086, 0.30853699, 0.73069673, 1.2884559, 1.9368036, 2.6298960, 3.3254259, 3.9870188, 4.5848036, 5.0955528, 5.5014381, 5.7900264, 5.9530302, 0, 0, 0, 0, 0, 0, 0, 0}, {0.051007539, 0.26407179, 0.62917446, 1.1184426, 1.6975514, 2.3300704, 2.9808563, 3.6182077, 4.2149578, 4.7486831, 5.2016291, 5.5601186, 5.8143095, 5.9576735, 0, 0, 0, 0, 0, 0, 0}, {0.044012859, 0.22848246, 0.54695073, 0.97852568, 1.4967736, 2.0727085, 2.6773271, 3.2834483, 3.8668193, 4.4066396, 4.8856367, 5.2899559, 5.6088495, 5.8344512, 5.9615265, 0, 0, 0, 0, 0, 0}, 
  {0.038362048, 0.19957643, 0.47954512, 0.86235527, 1.3274352, 1.8516264, 2.4110557, 2.9826202, 3.5450114, 4.0793097, 4.5692486, 5.0012398, 5.3642792, 5.6497668, 5.8513432, 5.9647572, 0, 0, 0, 0, 0}, {0.033732165, 0.17579230, 0.42367147, 0.76506833, 1.1838017, 1.6612632, 2.1777940, 2.7138534, 3.2509031, 3.7719973, 4.2621167, 4.7083098, 5.0996983, 5.4274082, 5.6844582, 5.8656493, 5.9674920, 0, 0, 0, 0}, {0.029891617, 0.15599578, 0.37688671, 0.68292459, 1.0612481, 1.4968065, 1.9733728, 2.4744513, 2.9840141, 3.4870391, 3.9698633, 4.4203793, 4.8281150, 5.1842329, 5.4814820, 5.7141263, 5.8778718, 5.9698274, 0, 0, 0}, {0.026670859, 0.13934795, 0.33735025, 0.61302887, 0.95605815, 1.3541796, 1.7939503, 2.2614442, 2.7428510, 3.2249469, 3.6954324, 4.1431513, 4.5582112, 4.9320320, 5.2573433, 5.5281521, 5.7396963, 5.8883966, 5.9718375, 0, 0}, 
  {0.023943491, 0.12521813, 0.30365763, 0.55312316, 0.86524556, 1.2299692, 1.6361091, 2.0718904, 2.5254313, 2.9851419, 3.4400320, 3.8799291, 4.2956189, 4.6789231, 5.0227300, 5.3209932, 5.5687097, 5.7618895, 5.8975241, 5.9735800, 0}, {0.021613735, 0.11312520, 0.27472423, 0.50143212, 0.78640441, 1.1213343, 1.4968701, 1.9030319, 2.3295969, 2.7664328, 3.2037683, 3.6323983, 4.0438292, 4.4303722, 4.7851954, 5.1023445, 5.3767428, 5.6041765, 5.7812753, 5.9054912, 5.9751004}};
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }
  
  }
  if ((vartheta_max==400)&&(vt_min==120)){printf("Entering roots max=%d min=%d\n",vartheta_max,vt_min);
  double store2[20][21]={{0.43058887, 1.0586370, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.14765840, 0.67742807, 1.1013503, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.096366761, 0.42994314, 0.84322366, 1.1323888, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.062691544, 0.30275755, 0.63423289, 0.94854207, 1.1533016,  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.044553356, 0.22032437, 0.48511535, 0.77139754, 1.0149905, 1.1667715, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.033175504, 0.16697336, 0.37832989, 0.62612861, 0.86597639, 1.0587776, 1.1756060, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.025652972, 0.13053095, 0.30148483, 0.51234543, 0.73279680, 0.93335028, 1.0890362, 1.1816465, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, 
  {0.020419599, 0.10468424, 0.24499304, 0.42411405, 0.62104663, 0.81398105, 0.98278234, 1.1107863, 1.1859481, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.016635362, 0.085739287, 0.20255775, 0.35536252, 0.52933963, 0.70836454, 0.87663693, 1.0200056, 1.1269288, 1.1891178, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.013811607, 0.071465007, 0.17001563, 0.30124280, 0.45448607, 0.61773520, 0.77869979, 0.92575155, 1.0486770, 1.1392302, 1.1915201, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.011649274, 0.060455380, 0.14458607, 0.25812712, 0.39325528, 0.54092507, 0.69157374, 0.83578161, 0.96483925, 1.0711996, 1.1488150, 1.1933840, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.0099570877,  0.051792244, 0.12437635, 0.22335564, 0.34288194, 0.47605711, 0.61540372, 0.75332691, 0.88253355, 0.99638811, 1.0891973,  1.1564255, 1.1948590, 0, 0, 0, 0, 0, 0, 0, 0}, {0.0086081534, 0.044856895, 0.10807158, 0.19498150, 0.30113973, 0.42121825, 0.54932716, 0.67933933, 0.80519796, 0.92119068, 1.0221813, 1.1037962, 1.1625675, 1.1960461, 0, 0, 0, 0, 0, 0, 0},
  {0.0075156278, 0.039220807, 0.094739727, 0.17157081, 0.26627943, 0.37469885, 0.49215035, 0.61367339, 0.73425118, 0.84901974,  0.95345238, 1.0435156, 1.1157953, 1.1675949, 1.1970157, 0, 0, 0, 0, 0, 0}, {0.0066184674, 0.034579920, 0.083707488, 0.15205707, 0.23693791, 0.33505270, 0.44265184, 0.55569842, 0.67003435, 0.78153912, 0.88627501, 0.98061458, 1.0613485, 1.1257735, 1.1717614, 1.1978177, 0, 0, 0, 0, 0}, {0.0058727564, 0.030713837, 0.074479816, 0.13563813, 0.21205337, 0.30108706, 0.39970760, 0.50460985, 0.61233829, 0.71940768, 0.82241666, 0.91815067, 1.0036722, 1.0763972, 1.1341579, 1.1752524, 1.1984888, 0, 0, 0, 0}, {0.0052462513, 0.027459755, 0.066686791, 0.12170368, 0.19079564, 0.27182943, 0.36233280, 0.45958232, 0.56069531, 0.66272190, 0.76273376, 0.85790658, 0.94559457, 1.0233956, 1.0892065, 1.1412691, 1.1782063, 1.1990558, 0, 0, 0},
  {0.0047148491, 0.024695358, 0.060047867, 0.10978400, 0.17251190, 0.24649010, 0.32968587, 0.41984040, 0.51453804, 0.61127692, 0.70753829, 0.80085244, 0.88885969, 0.96936512, 1.0403866, 1.1001957, 1.1473514, 1.1807276,  1.1995392, 0, 0}, {0.0042602412, 0.022327347, 0.054347375, 0.099513641, 0.15668525, 0.22442781, 0.30105688, 0.38468685,  0.47328360, 0.56471901, 0.65682525, 0.74744752, 0.83449367, 0.91597990, 0.99007171, 1.0551198, 1.1096910, 1.1525932, 1.1828968, 1.1999548, 0}, {0.0038683199, 0.020283579, 0.049417447, 0.090605265, 0.14290341, 0.20512061, 0.27585067, 0.35350961, 0.43637604, 0.52263342, 0.61041292, 0.69783552, 0.78305241, 0.86428286, 0.93984893, 1.0082066, 1.0679728, 1.1179494, 1.1571420, 1.1847763, 1.2003145}};
  
  for (j=0;j<20;j++){
    for (k=0;k<21;k++){
      store[j][k]=store2[j][k];
    }
  }
  }  

  if ((vartheta_max==488)&&(vt_min==3)){
    printf("Entering roots max=%d\n",vartheta_max);
    double store2[20][21]={{3.760152077, 4.874890308, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {1.224100018, 3.972743251, 4.895260611, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.6424044312, 2.373596954, 4.166951428, 4.917088321, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.3867529132, 1.670315933, 3.071019935, 4.331433085, 4.938330496, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.2583338656, 1.19392225, 2.423516004, 3.540013971, 4.468659536, 4.958146771, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.1840878909, 0.8879994204, 1.894377951, 2.957859394, 3.867843807, 4.580557554, 4.975851233, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.1376742987, 0.6807456982, 1.506965199, 2.441966475, 3.345948682, 4.101040486, 4.668980739, 4.990924187, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.1068349943, 0.5368141069, 1.217471495, 2.031784077, 2.866859037, 3.636849958, 4.271361563, 4.737503211, 5.003267194, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0.087214844, 0.441556818, 1.014521209, 1.722047231, 2.47968916, 3.214409772, 3.870889222, 4.405680489, 4.793015074, 5.01362487, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}, {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}};
    for (j=0;j<20;j++){
      for (k=0;k<21;k++){
	store[j][k]=store2[j][k];
      }
    }	
  }
}

  

//===========================================================
// Routines for COSEBIS =======================
//===========================================================

double integrand_Tm_log_shear_shear(double vartheta, void *p)
{
  double f,T_m,ximinus;
  
  struct cos * params = (struct cos *)p;
  int ORDER = (params->ORDER);
  double vt_max= (params->vt_max);
  double vt_min= (params->vt_min);
  int ni = (params->ni);
  int nj = (params->nj);
  T_m=T_log_minus_COSEBIS(vartheta,ORDER,vt_max,vt_min);
  ximinus=xi_pm_tomo(-1,vartheta,ni,nj);
  f=vartheta*ximinus*T_m;
  return f;
}


double integrand_Tp_log_shear_shear(double vartheta, void *p)
{
  double f,T_p,xiplus;
  
  struct cos * params = (struct cos *)p;
  int ORDER = (params->ORDER);
  double vt_max= (params->vt_max);
  double vt_min= (params->vt_min);
  int ni = (params->ni);
  int nj = (params->nj);
  
  T_p=T_log_plus_COSEBIS(vartheta,ORDER,vt_max,vt_min);
  xiplus=xi_pm_tomo(1,vartheta,ni,nj);
  f=vartheta*xiplus*T_p;
  return f;
}



// double integrand_Tp_log_shear_position(double vartheta, void *p)
// {
//   double f,T_p,xiplus;
  
//   struct cos * params = (struct cos *)p;
//   int ORDER = (params->ORDER);
//   double vt_max= (params->vt_max);
//   double vt_min= (params->vt_min);
//   int ni = (params->ni);
//   int nj = (params->nj);
  
//   T_p=T_log_plus_COSEBIS(vartheta,ORDER,vt_max,vt_min);
//   xiplus= w_gamma_t_tomo(vartheta,ni,nj);
//   f=vartheta*xiplus*T_p;
//   return f;
// }


// double integrand_Tp_log_position_position(double vartheta, void *p)
// {
//   double f,T_p,xiplus;
  
//   struct cos * params = (struct cos *)p;
//   int ORDER = (params->ORDER);
//   double vt_max= (params->vt_max);
//   double vt_min= (params->vt_min);
//   int ni = (params->ni);
//   int nj = (params->nj);
  
  
//   T_p=T_log_plus_COSEBIS(vartheta,ORDER,vt_max,vt_min);
//   xiplus= w_clustering_tomo(vartheta,ni,nj);
//   f=vartheta*xiplus*T_p;
//   return f;
// }



// double integrand_Tp_log_shear_magnification(double vartheta, void *p)
// {
//   double f,T_p,xiplus;
  
//   struct cos * params = (struct cos *)p;
//   int ORDER = (params->ORDER);
//   double vt_max= (params->vt_max);
//   double vt_min= (params->vt_min);
//   int ni = (params->ni);
//   int nj = (params->nj);
  
//   T_p=T_log_plus_COSEBIS(vartheta,ORDER,vt_max,vt_min);
//   xiplus=xi_shear_magnification_tomo(vartheta,ni,nj);
//   f=vartheta*xiplus*T_p;
//   return f;
// }


// double integrand_Tp_log_position_magnification(double vartheta, void *p)
// {
//   double f,T_p,xiplus;
  
//   struct cos * params = (struct cos *)p;
//   int ORDER = (params->ORDER);
//   double vt_max= (params->vt_max);
//   double vt_min= (params->vt_min);
//   int ni = (params->ni);
//   int nj = (params->nj);
  
//   T_p=T_log_plus_COSEBIS(vartheta,ORDER,vt_max,vt_min);
//   xiplus=xi_position_magnification_tomo(vartheta,ni,nj);
//   f=vartheta*xiplus*T_p;
//   return f;
// }


// double integrand_Tp_log_magnification_magnification(double vartheta, void *p)
// {
//   double f,T_p,xiplus;
  
//   struct cos * params = (struct cos *)p;
//   int ORDER = (params->ORDER);
//   double vt_max= (params->vt_max);
//   double vt_min= (params->vt_min);
//   int ni = (params->ni);
//   int nj = (params->nj);
  
//   T_p=T_log_plus_COSEBIS(vartheta,ORDER,vt_max,vt_min);
//   xiplus=xi_magnification_magnification_tomo(vartheta,ni,nj);
//   f=vartheta*xiplus*T_p;
//   return f;
// }


void COS_shear_shear_tomo(double *E,double *B,int ORDER,double vt_max, double vt_min,int ni, int nj)
{
  double res_tp,res_tm,error,error2;
  
  struct cos alpha;
  
  alpha.ORDER=ORDER;
  alpha.vt_min=vt_min;
  alpha.vt_max=vt_max;
  alpha.ni=ni;
  alpha.nj=nj;
  
  gsl_function integral_Tp_log;
  integral_Tp_log.function=&integrand_Tp_log_shear_shear;
  integral_Tp_log.params = &alpha;
  
  gsl_integration_cquad_workspace *w1;
  w1=gsl_integration_cquad_workspace_alloc(4000);
  gsl_integration_cquad(&integral_Tp_log,vt_min,vt_max,0,precision.low,w1,&res_tp,&error,0);
  gsl_integration_cquad_workspace_free(w1);
  
  gsl_function integral_Tm_log;
  integral_Tm_log.function=&integrand_Tm_log_shear_shear;
  integral_Tm_log.params = &alpha;
  
  gsl_integration_cquad_workspace *w2;
  w2=gsl_integration_cquad_workspace_alloc(4000);
  gsl_integration_cquad(&integral_Tm_log,vt_min,vt_max,0,precision.low,w2,&res_tm,&error2,0);
  gsl_integration_cquad_workspace_free(w2);
  
  *E=0.5*(res_tp+res_tm);
  *B=0.5*(res_tp-res_tm);
}


// void COS_pos_pos_tomo(double *E,int ORDER,double vt_max, double vt_min,int ni, int nj)
// {
//   double res_tp,error;
  
//   struct cos alpha;
  
//   alpha.ORDER=ORDER;
//   alpha.vt_min=vt_min;
//   alpha.vt_max=vt_max;
//   alpha.ni=ni;
//   alpha.nj=nj;
  
//   gsl_function integral_Tp_log;
//   integral_Tp_log.function=&integrand_Tp_log_position_position;
//   integral_Tp_log.params = &alpha;
  
//   gsl_integration_cquad_workspace *w1;
//   w1=gsl_integration_cquad_workspace_alloc(4000);
//   gsl_integration_cquad(&integral_Tp_log,vt_min,vt_max,0,precision.low,w1,&res_tp,&error,0);
//   gsl_integration_cquad_workspace_free(w1);
  
//   *E=res_tp;
// }

// void COS_shear_mag_tomo(double *E,int ORDER,double vt_max, double vt_min,int ni, int nj)
// {
//   double res_tp,error;
  
//   struct cos alpha;
  
//   alpha.ORDER=ORDER;
//   alpha.vt_min=vt_min;
//   alpha.vt_max=vt_max;
//   alpha.ni=ni;
//   alpha.nj=nj;
  
//   gsl_function integral_Tp_log;
//   integral_Tp_log.function=&integrand_Tp_log_shear_magnification;
//   integral_Tp_log.params = &alpha;
  
//   gsl_integration_cquad_workspace *w1;
//   w1=gsl_integration_cquad_workspace_alloc(4000);
//   gsl_integration_cquad(&integral_Tp_log,vt_min,vt_max,0,precision.low,w1,&res_tp,&error,0);
//   gsl_integration_cquad_workspace_free(w1);
  
//   *E=res_tp;
// }


// void COS_shear_pos_tomo(double *E,int ORDER,double vt_max, double vt_min,int ni, int nj)
// {
//   double res_tp,error;
  
//   struct cos alpha;
  
//   alpha.ORDER=ORDER;
//   alpha.vt_min=vt_min;
//   alpha.vt_max=vt_max;
//   alpha.ni=ni;
//   alpha.nj=nj;
  
//   gsl_function integral_Tp_log;
//   integral_Tp_log.function=&integrand_Tp_log_shear_position;
//   integral_Tp_log.params = &alpha;
  
//   gsl_integration_cquad_workspace *w1;
//   w1=gsl_integration_cquad_workspace_alloc(4000);
//   gsl_integration_cquad(&integral_Tp_log,vt_min,vt_max,0,precision.low,w1,&res_tp,&error,0);
//   gsl_integration_cquad_workspace_free(w1);
  
//   *E=res_tp;
// }


// void COS_mag_pos_tomo(double *E,int ORDER,double vt_max, double vt_min,int ni, int nj)
// {
//   double res_tp,error;
  
//   struct cos alpha;
  
//   alpha.ORDER=ORDER;
//   alpha.vt_min=vt_min;
//   alpha.vt_max=vt_max;
//   alpha.ni=ni;
//   alpha.nj=nj;
  
//   gsl_function integral_Tp_log;
//   integral_Tp_log.function=&integrand_Tp_log_position_magnification;
//   integral_Tp_log.params = &alpha;
  
//   gsl_integration_cquad_workspace *w1;
//   w1=gsl_integration_cquad_workspace_alloc(4000);
//   gsl_integration_cquad(&integral_Tp_log,vt_min,vt_max,0,precision.low,w1,&res_tp,&error,0);
//   gsl_integration_cquad_workspace_free(w1);
  
//   *E=res_tp;
// }

// void COS_mag_mag_tomo(double *E,int ORDER,double vt_max, double vt_min,int ni, int nj)
// {
//   double res_tp,error;
  
//   struct cos alpha;
  
//   alpha.ORDER=ORDER;
//   alpha.vt_min=vt_min;
//   alpha.vt_max=vt_max;
//   alpha.ni=ni;
//   alpha.nj=nj;
  
//   gsl_function integral_Tp_log;
//   integral_Tp_log.function=&integrand_Tp_log_magnification_magnification;
//   integral_Tp_log.params = &alpha;
  
//   gsl_integration_cquad_workspace *w1;
//   w1=gsl_integration_cquad_workspace_alloc(4000);
//   gsl_integration_cquad(&integral_Tp_log,vt_min,vt_max,0,precision.low,w1,&res_tp,&error,0);
//   gsl_integration_cquad_workspace_free(w1);
  
//   *E=res_tp;
// }




//========================================================================================
//========================================================================================
//Start of the ring statistics filter functions
//========================================================================================
//========================================================================================

/*============================================
 * Start of Y_+-functions //see SK07
 * =============================================*/
void Y_minus(double x, double eta, double *ym)
{
  double denom, etasqr, xsqr, x1;
  
  if (x>1.0+eta-precision.high || x<1.0-eta+precision.high) { *ym = 0.0; return; }
  
  assert(x>0.); assert(eta>0.);
  assert(eta<1.);
  assert(1.+eta>x); assert(x>1.-eta);
  
  etasqr = eta*eta;
  xsqr   = x*x;
  x1     = SQR(1.0-xsqr);
  
  *ym = x1 + etasqr*(-2.0*(2.0-xsqr) + etasqr*(6.0 + 2.0*xsqr + xsqr*xsqr 
  + etasqr*(-2.0*(2.0 + xsqr) + etasqr)));
  *ym = *ym/xsqr;
  
  denom  = etasqr*constants.pi*sqrt(SQR(1.0+eta) - xsqr);
  denom *= sqrt(xsqr - SQR(1.0-eta));
  
  if (denom<=0.)printf("eta=%le x=%le denom=%le\n",eta,x,denom);
  
  assert(denom>0.);
  
  *ym = *ym/denom;
}


void Y_plus(double x, double eta, double *yp)
{
  double denom, etasqr, xsqr, x1;
  
  if (x>1.0+eta-precision.high || x<1.0-eta+precision.high) { *yp = 0.0; return; }
  
  assert(x>0.); assert(eta>0.);
  assert(eta<1.);
  assert(1.+eta>x); assert(x>1.-eta);
  
  etasqr = eta*eta;
  xsqr   = x*x;
  x1     = SQR(1.0-xsqr);
  
  *yp = xsqr*(x1 + etasqr*(-2.0*xsqr + etasqr));
  
  denom  = etasqr*constants.pi*sqrt(SQR(1.0+eta) - xsqr);
  denom *= sqrt(xsqr - SQR(1.0-eta));
  
  if (denom<=0.)printf("eta=%le x=%le denom=%le\n",eta,x,denom);
  
  assert(denom>0.);
  
  *yp = *yp/denom;
}


/*========================================================================================
 * Start of Z-functions functions //
 * =========================================================================================*/


double integrand2_minus(double theta2, void *params)
{
  double res,ym,vartheta,zeta3,zeta4,theta1;
  
  double *array = (double *) params;
  zeta3 = array[2];
  zeta4 = array[3];
  vartheta=array[4];
  theta1=array[5];
  //printf("arraycheck %le %le %le %le %le %le\n",array[0], array[1],array[2], array[3],array[4], array[5]);
  //printf("eta=%le",theta1/theta2);
  Y_minus(vartheta/theta2, theta1/theta2,&ym);
  res=W2(theta2,zeta3,zeta4)*ym;
  //printf("W2%le\n",res);
  
  return res;
}


/*In this function we calculate the Integral over theta_2 MINUS*/
double integrand1_minus(double theta1, void *params)
{
  double res,f, a,b,error,zeta1,zeta2, zeta3,zeta4,vartheta,upper,lower;
  
  double *array = (double *) params;
  zeta1 = array[0];
  zeta2 = array[1];
  zeta3 = array[2];
  zeta4 = array[3];
  vartheta=array[4];
  array[5]=theta1;
  //printf("arraycheck %le %le %le %le %le %le\n",array[0], array[1],array[2], array[3],array[4], array[5]);
  a=vartheta+theta1;
  upper=FMIN(zeta4,a);
  b=vartheta-theta1;
  lower=FMAX(b,zeta3); 
  if (fabs(lower-upper)<precision.high){
    printf("integrand1_minus lower %le > upper %le \n",lower,upper);
    return 0;
  }
  //printf("theta2 integral: lower=%le upper=%le\n",lower, upper);
  gsl_function int2minus;
  int2minus.function=&integrand2_minus;
  int2minus.params = array;
  
  gsl_integration_workspace *w3;
  w3=gsl_integration_workspace_alloc(3000);
  gsl_integration_qag(&int2minus,lower,upper,1.e-21,1.e-5,3000,1,w3,&res,&error);
  gsl_integration_workspace_free(w3);
  f=res*W1(theta1,zeta1,zeta2);
  //printf("theta2 integral=%le\n",f);
  return f;
}

/*====================================================*/
/*============== Routines for RRstar ======================*/
/*====================================================*/

double integrand2_plus(double theta2, void *params)
{
  double res,yp,vartheta,zeta3,zeta4,theta1;
  
  double *array = (double *) params;
  zeta3 = array[2];
  zeta4 = array[3];
  vartheta=array[4];
  theta1=array[5];
  //printf("arraycheck %le %le %le %le %le %le\n",array[0], array[1],array[2], array[3],array[4], array[5]);
  //printf("eta=%le",theta1/theta2);
  Y_plus(vartheta/theta2, theta1/theta2,&yp);
  res=W2(theta2,zeta3,zeta4)*yp;
  //printf("W2%le\n",res);
  
  return res;
}


/*In this function we calculate the Integral over theta_2 PLUS*/
double integrand1_plus(double theta1, void *params)
{
  double res,f, a,b,error,zeta1,zeta2, zeta3,zeta4,vartheta,upper,lower;
  
  double *array = (double *) params;
  zeta1 = array[0];
  zeta2 = array[1];
  zeta3 = array[2];
  zeta4 = array[3];
  vartheta=array[4];
  array[5]=theta1;
  //printf("arraycheck %le %le %le %le %le %le\n",array[0], array[1],array[2], array[3],array[4], array[5]);
  a=vartheta+theta1;
  upper=FMIN(zeta4,a);
  b=vartheta-theta1;
  lower=FMAX(b,zeta3);
  if (fabs(lower-upper)<precision.high){
    printf("integrand1_plus lower %le > upper %le\n",lower,upper);
    return 0;
  }
  //printf("theta2 integral: lower=%le upper=%le\n",lower, upper);
  gsl_function int2plus;
  int2plus.function=&integrand2_plus;
  int2plus.params = array;
  
  gsl_integration_workspace *w3;
  w3=gsl_integration_workspace_alloc(3000);
  gsl_integration_qag(&int2plus,lower,upper,1.e-21,1.e-5,3000,1,w3,&res,&error);
  gsl_integration_workspace_free(w3);
  f=res*W1(theta1,zeta1,zeta2);
  //printf("theta2 integral=%le\n",f);
  return f;
}


double Z_minus(double vartheta,double *array)
{
  double res,a,b,zeta1,zeta2, zeta3,zeta4,upper,lower,error;
  
  zeta1 = array[0];
  zeta2 = array[1];
  zeta3 = array[2];
  zeta4 = array[3];
  array[4]=vartheta;
  
  //printf("arraycheck Z_minus%le %le %le %le %le %le\n",array[0], array[1],array[2], array[3],array[4], array[5]);
  a=vartheta-zeta4;
  b=zeta3-vartheta;
  lower=FMAX(a,b);
  lower=FMAX(zeta1,lower);
  upper=zeta2;
  if (fabs(lower-upper)<precision.high){
    //printf("Z_minus lower=%le > upper=%le\n",lower, upper);
    //printf("theta1_minus integral: lower=%le upper=%le\n",lower, upper);
    return 0;
  }
  
  gsl_function int1minus;
  int1minus.function=&integrand1_minus;
  int1minus.params = array;
  
  gsl_integration_workspace *w2;
  w2=gsl_integration_workspace_alloc(3000);
  gsl_integration_qag(&int1minus,lower,upper,1.e-21,1.e-5,3000,6,w2,&res,&error);
  gsl_integration_workspace_free(w2);
  //printf("theta1 integral=%le\n",f);
  return res;
}


double Z_plus(double vartheta, double *array)
{
  double res,a,b,zeta1,zeta2,zeta3,zeta4,upper,lower,error;
  
  zeta1 = array[0];
  zeta2 = array[1];
  zeta3 = array[2];
  zeta4 = array[3];
  array[4]=vartheta;
  
  //printf("arraycheck %le %le %le %le %le %le\n",array[0], array[1],array[2], array[3],array[4], array[5]);
  a=vartheta-zeta4;
  b=zeta3-vartheta;
  lower=FMAX(a,b);
  lower=FMAX(zeta1,lower);
  upper=zeta2;
  if (fabs(lower-upper)<precision.high){
    //printf("Z_plus lower=%le > upper=%le\n",lower, upper);
    //printf("zeta1=%le vartheta-zeta4=%le zeta3-vartheta%le",zeta1,a,b);
    //printf("theta1_plus integral: lower=%le upper=%le\n",lower, upper);
    return 0;
  }
  
  gsl_function int1plus;
  int1plus.function=&integrand1_plus;
  int1plus.params = array;
  
  gsl_integration_workspace *w2;
  w2=gsl_integration_workspace_alloc(3000);
  gsl_integration_qag(&int1plus,lower,upper,1.e-21,1.e-5,3000,6,w2,&res,&error);
  gsl_integration_workspace_free(w2);
  //printf("theta1 integral=%le\n",f);
  return res;
}


/*============================================
 * Start of W_12 functions //see SK07
 * =============================================*/


double W1(double theta_1, double zeta_1, double zeta_2)
{
  double nom,denom,res,a1,a2,a3;
  a1=theta_1-zeta_1;
  a2=zeta_2-theta_1;
  a3=zeta_2-zeta_1;
  nom=30.*SQR(a1*a2);
  denom=pow(a3,5);
  res=nom/denom;
  
  return res;
  
}


double W2(double theta_2, double zeta_3, double zeta_4)
{
  double nom,denom,res,a1,a2,a3;
  
  a1=theta_2-zeta_3;
  a2=zeta_4-theta_2;
  a3=zeta_4-zeta_3;
  nom=30.*SQR(a1*a2);
  denom=pow(a3,5);
  res=nom/denom;
  
  return res;
  
}




/*=============================================================================
 * =============================================================================
 * EB-Filter functions routines for MAP
 * ==============================================================================
 * ==============================================================================*/


double T_plus_MAP(double x)
{
  double a,b,c,d,e,f,g,res;
  res=0.;
  if (x<=2.0){
    g=x*sqrt(4.-x*x)/100./constants.pi;
    a=6.*(2.-15.*x*x)/5.*(1.-2./constants.pi*asin(x/2.));
    b=120.*g;
    c=2320.*x*x*g;
    d=754.*x*x*x*x*g;
    e=132.*pow(x,6.0)*g;
    f=9.*pow(x,8.0)*g;
    res=a+b+c-d+e-f;
  }
  return res;
}

double T_minus_MAP(double x)
{
  double res;
  res=0.;
  if (x<=2.0){
    res=192./(35.*constants.pi)*x*x*x*pow((1.-x*x/4.),3.5);
  }
  return res;
}


