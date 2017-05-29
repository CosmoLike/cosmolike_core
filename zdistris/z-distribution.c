#include <math.h>
#include <stdlib.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <stdio.h>
#include <assert.h>

#define NR_END 1
#define FREE_ARG char*

double sm2_qromb(double (*func)(double), double a, double b);
double int_for_zdistr(double z);
void set_redshift_DES_conti();
void set_redshift_DES_SV();
void set_redshift_DES_Tully_Fisher();
void set_redshift_LSST_conti();
void set_redshift_LSST_gold_conti();
void set_redshift_Euclid_conti();
void set_redshift_WFIRST_conti();


double pf(double z);
double sm2_trapzd(double (*func)(double), double a, double b, int n, double *s);
void sm2_polint(double xa[], double ya[], int n, double x, double *y, double *dy);
void sm2_error(char *s);
double *sm2_vector(long nl, long nh);
void sm2_free_vector(double *v, long nl, long nh);

double pf_LSST(double z);
void LSST_zdistricalc();
void DES_zdistricalc();


double *sm2_vector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) sm2_error("allocation failure in double vector()");
	return v-nl+NR_END;
}

void sm2_error(char *s)
{
	printf("error:%s\n ",s);
	exit(1);
}


void sm2_polint(double xa[], double ya[], int n, double x, double *y, double *dy)
{
	int i,m,ns=1;
	double den,dif,dift,ho,hp,w;
	double *c,*d;

	dif=fabs(x-xa[1]);
	c=sm2_vector(1,n);
	d=sm2_vector(1,n);
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0)
				sm2_error("Error in routine sm2_polint");
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	sm2_free_vector(d,1,n);
	sm2_free_vector(c,1,n);
}

void sm2_free_vector(double *v, long nl, long nh)
		/* free a double vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

#define FUNC(x) ((*func)(x))

double sm2_trapzd(double (*func)(double), double a, double b, int n, double *s)
{
	double x,tnm,sum,del;
	int it,j;

	if (n == 1) {
		return (*s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
	} else {
		for (it=1,j=1;j<n-1;j++) it <<= 1;
		tnm=it;
		del=(b-a)/tnm;
		x=a+0.5*del;
		for (sum=0.0,j=1;j<=it;j++,x+=del) {
			sum += FUNC(x);
		}
		*s=0.5*(*s+(b-a)*sum/tnm);
		return *s;
	}
}
#undef FUNC

typedef struct {
     double z0;		/* redshift distribution scale                    */
     double beta_p;		/* redshift distribution power index. If beta_b = 0, a distribution with a single source redshift at z0 is assumed */
     double alpha;  /* redshift distribution power index.*/
     double zdistrpar_zmin;     /*cut-offs in the redshift distribution*/
     double zdistrpar_zmax;
     char REDSHIFT_FILE[200];
     int histogram_zbins;

}redshiftclusteringpara;


redshiftclusteringpara redshiftclustering = {
  1.0,
  0.0,
  0.0,
  0.0,
  3.0,
  "",
  0
};

double int_for_zdistr(double z)
{
  double zz=z/redshiftclustering.z0;
  return pow(zz,redshiftclustering.alpha)*exp(-pow(zz,redshiftclustering.beta_p));
}

double int_for_zdistr_LSST(double z)
{
  double zz=z/redshiftclustering.z0;
  return pow(z,redshiftclustering.alpha)*exp(-pow(zz,redshiftclustering.beta_p));
}

void set_redshift_DES_conti()
{
  redshiftclustering.z0        = 0.56;
  redshiftclustering.beta_p    = 1.5;
  redshiftclustering.alpha   =   1.3;
  redshiftclustering.zdistrpar_zmin = 0.0;
  redshiftclustering.zdistrpar_zmax = 2.0;
  sprintf(redshiftclustering.REDSHIFT_FILE,"dummy");
  redshiftclustering.histogram_zbins=0;
}

void set_redshift_Euclid_conti()
{
  redshiftclustering.z0        = 0.65;
  redshiftclustering.beta_p    = 1.5;
  redshiftclustering.alpha   =   1.3;
  redshiftclustering.zdistrpar_zmin = 0.0;
  redshiftclustering.zdistrpar_zmax = 2.5;
  sprintf(redshiftclustering.REDSHIFT_FILE,"dummy");
  redshiftclustering.histogram_zbins=0;
}

void set_redshift_LSST_conti()
{
  redshiftclustering.z0        = .5;
  redshiftclustering.beta_p    = 1.02;
  redshiftclustering.alpha   =   1.27;
  redshiftclustering.zdistrpar_zmin = 0.;
  redshiftclustering.zdistrpar_zmax = 3.5;
  sprintf(redshiftclustering.REDSHIFT_FILE,"dummy");
  redshiftclustering.histogram_zbins=0;
}

void set_redshift_LSST_gold_conti()
{
  redshiftclustering.z0        = .3;
  redshiftclustering.beta_p    = 1.0;
  redshiftclustering.alpha   =   2.;
  redshiftclustering.zdistrpar_zmin = 0.;
  redshiftclustering.zdistrpar_zmax = 3.5;
  sprintf(redshiftclustering.REDSHIFT_FILE,"dummy");
  redshiftclustering.histogram_zbins=0;
}

void set_redshift_WFIRST_conti() //uses LSST routine but slightly deeper
{
  redshiftclustering.z0        = .6;
  redshiftclustering.beta_p    = 1.02;
  redshiftclustering.alpha   =   1.27;
  redshiftclustering.zdistrpar_zmin = 0.;
  redshiftclustering.zdistrpar_zmax = 4.0;
  sprintf(redshiftclustering.REDSHIFT_FILE,"dummy");
  redshiftclustering.histogram_zbins=0;
}

void set_redshift_DES_Tully_Fisher()
{
  redshiftclustering.z0        = 0.0157454;
  redshiftclustering.beta_p    = 0.63342;
  redshiftclustering.alpha   =   5.625641;
  redshiftclustering.zdistrpar_zmin = 0.;
  redshiftclustering.zdistrpar_zmax = 1.3;
  sprintf(redshiftclustering.REDSHIFT_FILE,"dummy");
  redshiftclustering.histogram_zbins=0;
}

void set_redshift_DES_SV()
{
  redshiftclustering.z0        = 0.46;
  redshiftclustering.beta_p    =1.5;
  redshiftclustering.alpha   =   1.3;
  redshiftclustering.zdistrpar_zmin = 0.;
  redshiftclustering.zdistrpar_zmax = 1.6;
  sprintf(redshiftclustering.REDSHIFT_FILE,"../../zdistris/zdistribution_DES_Y5");
  redshiftclustering.histogram_zbins=0;
  printf("Redshift set to DES_SV\n");
}

// HSC source distribution
// from Oguri & Takada 2011 (HSC whitepaper oct 2012)
void set_redshift_HSC()
{
   redshiftclustering.z0        = 0.333333333;
   redshiftclustering.beta_p    = 1.;
   redshiftclustering.alpha   =   2.;
   redshiftclustering.zdistrpar_zmin = 0.0;
   redshiftclustering.zdistrpar_zmax = 3.75; // convergence of int dz*dn/dz: 1% at z=2.8, 0.1% at z=3.75, 0.01% at z=4.65
}


/* ============================================================ *
 * qromb.c							*
 * Romberg Integration. Uses trapzd. NR p. 140			*
 * ============================================================ */

#define EPS 1.0e-6
#define JMAX 35
#define JMAXP (JMAX+1)
#define K 5

double sm2_qromb(double (*func)(double), double a, double b)
{
	double ss,dss;
	double s[JMAXP],h[JMAXP+1], strap;
	int j;

	h[1]=1.0;
	for (j=1;j<=JMAX;j++) {
		s[j]=sm2_trapzd(func,a,b,j,&strap);
		if (j >= K) {
			sm2_polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
			if (fabs(dss) <= EPS*fabs(ss)) return ss;
		}
		h[j+1]=0.25*h[j];
	}
	sm2_error("Too many steps in routine sm2_qromb");
	return 0.0;
}
#undef EPS
#undef JMAX
#undef JMAXP
#undef K

#define EPS 1.0e-7
#define JMAX 35
#define JMAXP (JMAX+1)
#define K 5




double pf(double z)
{
  static double norm = 0.;
  static double BETA_P = -42.;
  static double ALPHA  = -42.;
  static double Z0     = -42.;
  static double ZMIN   = -42.;
  static double ZMAX   = -42.;
  double x, f;
  //First, compute the normalization
  if (ALPHA != redshiftclustering.alpha || BETA_P != redshiftclustering.beta_p || Z0 != redshiftclustering.z0 || ZMIN !=redshiftclustering.zdistrpar_zmin || ZMAX !=redshiftclustering.zdistrpar_zmax)
  {
    
    if((redshiftclustering.beta_p>0.) && (redshiftclustering.histogram_zbins == 0 ))
    {
      //                      if(redshiftclustering.zdistrpar_zmin || redshiftclustering.zdistrpar_zmax)
      {       //with cuts: use numeric normalization
      norm = 1.0/(sm2_qromb(int_for_zdistr,redshiftclustering.zdistrpar_zmin,redshiftclustering.zdistrpar_zmax));
      //printf("%le\n",norm);
      }
    }
    ALPHA  = redshiftclustering.alpha;
    BETA_P = redshiftclustering.beta_p;
    Z0     = redshiftclustering.z0;
    ZMIN   = redshiftclustering.zdistrpar_zmin;
    ZMAX   = redshiftclustering.zdistrpar_zmax;
  }
  
  if((redshiftclustering.zdistrpar_zmin || redshiftclustering.zdistrpar_zmax) && (z>redshiftclustering.zdistrpar_zmax || z<redshiftclustering.zdistrpar_zmin)) return 0.0;
  if (z<0 && z<-1e-6)
  {
    printf("Negative z in %s, line%d \n",__FILE__,__LINE__);
  }
  if (z<0) x = 0;
  else x = z/redshiftclustering.z0;
  if (fabs(redshiftclustering.beta_p - 1.) < 1e-6) {
    f = x;
  } else {
    f = pow(x,redshiftclustering.beta_p);
  }
  f=exp(-f);
  return norm*pow(x,redshiftclustering.alpha)*f;
}

double pf_LSST(double z)
{
  static double norm = 0.;
  static double BETA_P = -42.;
  static double ALPHA  = -42.;
  static double Z0     = -42.;
  static double ZMIN   = -42.;
  static double ZMAX   = -42.;
  double x, f;
  //First, compute the normalization
  if (ALPHA != redshiftclustering.alpha || BETA_P != redshiftclustering.beta_p || Z0 != redshiftclustering.z0 || ZMIN !=redshiftclustering.zdistrpar_zmin || ZMAX !=redshiftclustering.zdistrpar_zmax)
  {
    
    if((redshiftclustering.beta_p>0.) && (redshiftclustering.histogram_zbins == 0 ))
    {
      //                      if(redshiftclustering.zdistrpar_zmin || redshiftclustering.zdistrpar_zmax)
      {       //with cuts: use numeric normalization
      norm = 1.0/(sm2_qromb(int_for_zdistr_LSST,redshiftclustering.zdistrpar_zmin,redshiftclustering.zdistrpar_zmax));
      //printf("%le\n",norm);
      }
    }
    ALPHA  = redshiftclustering.alpha;
    BETA_P = redshiftclustering.beta_p;
    Z0     = redshiftclustering.z0;
    ZMIN   = redshiftclustering.zdistrpar_zmin;
    ZMAX   = redshiftclustering.zdistrpar_zmax;
  }
  
  if((redshiftclustering.zdistrpar_zmin || redshiftclustering.zdistrpar_zmax) && (z>redshiftclustering.zdistrpar_zmax || z<redshiftclustering.zdistrpar_zmin)) return 0.0;
  if (z<0 && z<-1e-6)
  {
    printf("Negative z in %s, line%d \n",__FILE__,__LINE__);
  }
  if (z<0) x = 0;
  else x = z/redshiftclustering.z0;
  if (fabs(redshiftclustering.beta_p - 1.) < 1e-6) {
    f = x;
  } else {
    f = pow(x,redshiftclustering.beta_p);
  }
  f=exp(-f);
  return norm*pow(z,redshiftclustering.alpha)*f;
}



void normal_zdistricalc()
{
  int i;
  double z;
  int Nstep=300;
  double res;
  double sum[301];
  int flag[5]={0,0,0,0,0};
  set_redshift_LSST_conti();
  FILE *F;
  F=fopen("zdistribution_LSST_conti2","w");
  sum[0]=0.0;
  i=0;
  z=0.0;
  double dz=(redshiftclustering.zdistrpar_zmax-redshiftclustering.zdistrpar_zmin)/((Nstep)*1.0);
  
  for(i=0;i<Nstep;i++){
    z=0.0+(i+0.5)*dz;
    res=pf(z);
    sum[i+1]=sum[i]+pf(z);
    fprintf(F,"%le %le %le %le\n",0.0+i*dz,z,0.0+(i+1)*dz,res);   
   // printf("%le %le %le %le\n",0.0+i*dz,z,0.0+(i+1)*dz,res);
  }
  fclose(F);
  
  for(i=0;i<Nstep;i++){
    z=0.0+(i+0.5)*dz;
    if (i==0) printf("z1_min: %le\n",z); 
    if((sum[i+1]/sum[300])> 0.2 && (flag[0]==0)) {
      printf("z1_max: %le\n",z); 
      flag[0]=1;
    }
    if((sum[i+1]/sum[300]) > 0.4 && (flag[1]==0)) {
      printf("z2_max: %le\n",z); 
      flag[1]=1;
    }
    if((sum[i+1]/sum[300]) > 0.6 && (flag[2]==0)) {
      printf("z3_max: %le\n",z); 
      flag[2]=1;
    }
    if((sum[i+1]/sum[300]) > 0.8 && (flag[3]==0)) {
      printf("z4_max: %le\n",z); 
      flag[3]=1;
    }
    if((i==299) && (flag[4]==0)) {
      printf("z5_max: %le\n",z); 
      flag[4]=1;
    }
  }
}

void DES_SV_zdistricalc()
{
  int i;
  double z;
  int Nstep=300;
  double res;
  double sum[301];
  int flag[5]={0,0,0,0,0};
  set_redshift_DES_SV();
  FILE *F;
  F=fopen("zdistribution_DES_SV","w");
  sum[0]=0.0;
  i=0;
  z=0.0;
  double dz=(redshiftclustering.zdistrpar_zmax-redshiftclustering.zdistrpar_zmin)/((Nstep)*1.0);
  
  for(i=0;i<Nstep;i++){
    z=0.0+(i+0.5)*dz;
    res=pf_LSST(z);
    sum[i+1]=sum[i]+pf_LSST(z);
    fprintf(F,"%le %le %le %le\n",0.0+i*dz,z,0.0+(i+1)*dz,res);   
   // printf("%le %le %le %le\n",0.0+i*dz,z,0.0+(i+1)*dz,res);
  }
  fclose(F);
  for(i=0;i<Nstep;i++){
    z=0.0+(i+0.5)*dz;
    if (i==0) printf("z1_min: %le\n",z); 
    if((sum[i+1]/sum[300])> 0.2 && (flag[0]==0)) {
      printf("z1_max: %le\n",z); 
      flag[0]=1;
    }
    if((sum[i+1]/sum[300]) > 0.4 && (flag[1]==0)) {
      printf("z2_max: %le\n",z); 
      flag[1]=1;
    }
    if((sum[i+1]/sum[300]) > 0.6 && (flag[2]==0)) {
      printf("z3_max: %le\n",z); 
      flag[2]=1;
    }
    if((sum[i+1]/sum[300]) > 0.8 && (flag[3]==0)) {
      printf("z4_max: %le\n",z); 
      flag[3]=1;
    }
    if((i==299) && (flag[4]==0)) {
      printf("z5_max: %le\n",z); 
      flag[4]=1;
    }
  }
}

void LSST_zdistricalc()
{
  int i;
  double z;
  int Nstep=300;
  double res;
  double sum[301];
  int flag[5]={0,0,0,0,0};
  set_redshift_LSST_conti();
  FILE *F;
  F=fopen("zdistribution_LSST_conti2","w");
  sum[0]=0.0;
  i=0;
  z=0.0;
  double dz=(redshiftclustering.zdistrpar_zmax-redshiftclustering.zdistrpar_zmin)/((Nstep)*1.0);
  
  for(i=0;i<Nstep;i++){
    z=0.0+(i+0.5)*dz;
    res=pf_LSST(z);
    sum[i+1]=sum[i]+pf_LSST(z);
    fprintf(F,"%le %le %le %le\n",0.0+i*dz,z,0.0+(i+1)*dz,res);   
   // printf("%le %le %le %le\n",0.0+i*dz,z,0.0+(i+1)*dz,res);
  }
  fclose(F);
}

void LSST_gold_zdistricalc()
{
  int i;
  double z;
  int Nstep=300;
  double res;
  double sum[301];
  int flag[5]={0,0,0,0,0};
  set_redshift_LSST_gold_conti();
  FILE *F;
  F=fopen("zdistribution_LSST_gold","w");
  sum[0]=0.0;
  i=0;
  z=0.0;
  double dz=(redshiftclustering.zdistrpar_zmax-redshiftclustering.zdistrpar_zmin)/((Nstep)*1.0);
  
  for(i=0;i<Nstep;i++){
    z=0.0+(i+0.5)*dz;
    res=pf(z);
    sum[i+1]=sum[i]+pf(z);
    fprintf(F,"%le %le %le %le\n",0.0+i*dz,z,0.0+(i+1)*dz,res);   
   // printf("%le %le %le %le\n",0.0+i*dz,z,0.0+(i+1)*dz,res);
  }
  fclose(F);
}

void WFIRST_zdistricalc()
{
  int i;
  double z;
  int Nstep=300;
  double res;
  double sum[301];
  int flag[5]={0,0,0,0,0};
  set_redshift_WFIRST_conti();
  FILE *F;
  F=fopen("zdistribution_WFIRST_conti","w");
  sum[0]=0.0;
  i=0;
  z=0.0;
  double dz=(redshiftclustering.zdistrpar_zmax-redshiftclustering.zdistrpar_zmin)/((Nstep)*1.0);
  
  for(i=0;i<Nstep;i++){
    z=0.0+(i+0.5)*dz;
    res=pf_LSST(z);
    sum[i+1]=sum[i]+pf_LSST(z);
    fprintf(F,"%le %le %le %le\n",0.0+i*dz,z,0.0+(i+1)*dz,res);   
   // printf("%le %le %le %le\n",0.0+i*dz,z,0.0+(i+1)*dz,res);
  }
  fclose(F);
}

void DES_zdistricalc()
{
  int i;
  double z;
  int Nstep=300;
  double res;
  double sum[301];
  int flag[5]={0,0,0,0,0};
  set_redshift_DES_conti();
  FILE *F;
  F=fopen("zdistribution_DES_Y6","w");
  sum[0]=0.0;
  i=0;
  z=0.0;
  double dz=(redshiftclustering.zdistrpar_zmax-redshiftclustering.zdistrpar_zmin)/((Nstep)*1.0);
  
  for(i=0;i<Nstep;i++){
    z=0.0+(i+0.5)*dz;
    res=pf(z);
    sum[i+1]=sum[i]+pf(z);
    fprintf(F,"%le %le %le %le\n",0.0+i*dz,z,0.0+(i+1)*dz,res);   
   // printf("%le %le %le %le\n",0.0+i*dz,z,0.0+(i+1)*dz,res);
  }
  fclose(F);
}

void Euclid_zdistricalc()
{
  int i;
  double z;
  int Nstep=300;
  double res;
  double sum[301];
  int flag[5]={0,0,0,0,0};
  set_redshift_Euclid_conti();
  FILE *F;
  F=fopen("zdistribution_Euclid_conti","w");
  sum[0]=0.0;
  i=0;
  z=0.0;
  double dz=(redshiftclustering.zdistrpar_zmax-redshiftclustering.zdistrpar_zmin)/((Nstep)*1.0);
  
  for(i=0;i<Nstep;i++){
    z=0.0+(i+0.5)*dz;
    res=pf(z);
    sum[i+1]=sum[i]+pf(z);
    fprintf(F,"%le %le %le %le\n",0.0+i*dz,z,0.0+(i+1)*dz,res);   
   // printf("%le %le %le %le\n",0.0+i*dz,z,0.0+(i+1)*dz,res);
  }
  fclose(F);  
}


void HSC_zdistricalc()
{
   int i;
   double z;
   int Nstep=300;
   double res;
   double sum[301];
   int flag[5]={0,0,0,0,0};
   set_redshift_HSC();
   FILE *F;
   F=fopen("zdistribution_HSC","w");
   sum[0]=0.0;
   i=0;
   z=0.0;
   double dz=(redshiftclustering.zdistrpar_zmax-redshiftclustering.zdistrpar_zmin)/((Nstep)*1.0);
   
   for(i=0;i<Nstep;i++){
      z=0.0+(i+0.5)*dz;
      res=pf(z);
      sum[i+1]=sum[i]+pf(z);
      fprintf(F,"%le %le %le %le\n",0.0+i*dz,z,0.0+(i+1)*dz,res);
//      printf("%le %le %le %le\n",0.0+i*dz,z,0.0+(i+1)*dz,res);
   }
   fclose(F);
}


void multi_zdistricalc()
{
  int i,j;
  double z;
  int Nstep=300;
  double res;
  double sum[301];
  double alpha,zmin,zmax,dz;
  FILE *F;
  char filename[200];
  
  double z0[6]={0.41,0.42,0.44,0.46,0.48,0.5};
  alpha=1.27;
  double beta[6]={1.2,1.18,1.14,1.1,1.06,1.02};
  redshiftclustering.histogram_zbins=0;
  
  for(j=0;j<6;j++){  
    sprintf(filename,"../../cosmolike_light/zdistris/zdistributions_LSSTawakens/zdistri_model_%d",j);
    F=fopen(filename,"w");
    
    redshiftclustering.z0        = z0[j];
    redshiftclustering.beta_p    = beta[j];
    redshiftclustering.alpha   =   alpha;
    redshiftclustering.zdistrpar_zmin = 0.000001;
    redshiftclustering.zdistrpar_zmax = 3.5;
    dz=(redshiftclustering.zdistrpar_zmax-redshiftclustering.zdistrpar_zmin)/((Nstep)*1.0);
    
    for(i=0;i<Nstep;i++){
      z=0.0+(i+0.5)*dz;
      res=pf_LSST(z);
      fprintf(F,"%le %le %le %le\n",0.0+i*dz,z,0.0+(i+1)*dz,res);   
    }
    fclose(F);
  }
  printf("hallo\n");
}

void normalization_tomography_binning(char *filename)
{
  int i;
  double z,zmin,zmid,zmax,res[300], zmin_tot,zmax_tot,zvec[5];
  int Nstep=300;
  double sum[301];
  int flag[5]={0,0,0,0,0};
  FILE *F,*F2;
  F=fopen(filename,"r");
 
  sum[0]=0.0;
  for(i=0;i<Nstep;i++){
    fscanf(F,"%le %le %le %le\n",&zmin,&zmid,&zmax,&res[i]);    
    if(i==0)zmin_tot=zmin;
    if(i==299) zmax_tot=zmax;
    sum[i+1]=sum[i]+res[i];
   // printf("%le %le %le %le\n",0.0+i*dz,z,0.0+(i+1)*dz,res);
  }
  double dz=(zmax_tot-zmin_tot)/((Nstep)*1.0);
    
  for(i=0;i<Nstep;i++){
    z=0.0+(i+0.5)*dz;
    if (i==0) {
      printf("z1_min: %le\n",z); 
      zvec[0]=z;    
    }
    if((sum[i+1]/sum[300])> 0.2 && (flag[0]==0)) {
      zvec[1]=z;
      printf("tomo.shear_zmax[0] = %le;\n",z); 
      flag[0]=1;
    }
    if((sum[i+1]/sum[300]) > 0.4 && (flag[1]==0)) {
      zvec[2]=z;
      printf("tomo.shear_zmax[1] = %le;\n",z); 
      flag[1]=1;
    }
    if((sum[i+1]/sum[300]) > 0.6 && (flag[2]==0)) {
      zvec[3]=z;
      printf("tomo.shear_zmax[2] = %le;\n",z); 
      flag[2]=1;
    }
    if((sum[i+1]/sum[300]) > 0.8 && (flag[3]==0)) {
      zvec[4]=z;
      printf("tomo.shear_zmax[3] = %le;\n",z); 
      flag[3]=1;
    }
    if((i==299) && (flag[4]==0)) {
      printf("tomo.shear_zmax[4] = %le;\n",z); 
      flag[4]=1;
    }
  }
 printf("tomo.shear_zmin[0] = %le;\n",zvec[0]); 
 printf("tomo.shear_zmin[1] = %le;\n",zvec[1]); 
 printf("tomo.shear_zmin[2] = %le;\n",zvec[2]); 
 printf("tomo.shear_zmin[3] = %le;\n",zvec[3]); 
 printf("tomo.shear_zmin[4] = %le;\n",zvec[4]); 
  
}

int main(void)
{
//DES_SV_zdistricalc();
//BCC_zdistricalc();
//DES_zdistricalc();
//LSST_zdistricalc();
LSST_gold_zdistricalc();
//multi_zdistricalc();

//  WFIRST_zdistricalc();
//  Euclid_zdistricalc();
//  normalization_tomography_binning("zdistribution_Euclid_conti");
//  normalization_tomography_binning("zdistribution_WFIRST_conti");

  // HSC_zdistricalc();
  // normalization_tomography_binning("zdistribution_HSC");
  return 0;
}
