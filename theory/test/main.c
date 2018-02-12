#include <math.h>
#ifdef GSLLIB
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_integration.h>
#endif

typedef struct {
  double mass;
  double mass_min;
  double mass_M1;
  double intrinsic_alpha;
  double intrinsic_sigma;
  double true_lambda;
  double obs_lambda_min[10];
  double obs_lambda_max[10];
  double obs_lambda;
  double z;
  double k;
  int N_lambda;
} Par_lambda_obs;

double probability_true_richness_given_mass(const double true_lambda, 
void * params) {
  Par_lambda_obs *p=(Par_lambda_obs *)params;
  const double mass = p->mass;
  const double mass_min = p->mass_min;
  const double mass_M1 = p->mass_M1;
  const double intrinsic_alpha = p->intrinsic_alpha;
  // intrisic scatter of mass-richness relation
  const double intrinsic_sigma = p->intrinsic_sigma;

  // average true satellite richness given mass
  const double atsrgm = 
    pow((mass - mass_min)/(mass_M1 - mass_min), intrinsic_alpha);

  // Now we need skewness and variance of the skew-normal distribution
  // (which is a function of the instrinsic sigma of the mass-richness relation
  // and also the mean true satellite richness).
  static int first = 0;
  static gsl_spline2d* falpha; // skewness
  static gsl_spline2d* fsigma;
  static double* vector_intrinsic_sigma;
  static double* vector_atsrgm;
  static double* vector_alpha;
  static double* vector_sigma;

  if(first == 0) {
    first = 1;

    const int size_vector_intrinsic_sigma = 14;
    const int size_vector_atsrgm = 42;

    vector_intrinsic_sigma = 
      (double*) malloc(size_vector_intrinsic_sigma*sizeof(double));
    if(vector_intrinsic_sigma == NULL) {
      printf("malloc fail \n");
      exit(1);
    }
    vector_atsrgm = (double*) malloc(size_vector_atsrgm*sizeof(double));
    if(vector_atsrgm == NULL) {
      printf("malloc fail \n");
      exit(1);
    }
    vector_alpha = (double*) 
      malloc(size_vector_atsrgm*size_vector_intrinsic_sigma*sizeof(double));
    if(vector_alpha == NULL) {
      printf("malloc fail \n");
      exit(1);
    }
    vector_sigma = (double*) 
      malloc(size_vector_atsrgm*size_vector_intrinsic_sigma*sizeof(double));
    if(vector_sigma == NULL) {
      printf("malloc fail \n");
      exit(1);
    }

    const gsl_interp2d_type *T = gsl_interp2d_bilinear;

    falpha = gsl_spline2d_alloc(T, size_vector_intrinsic_sigma, 
      size_vector_atsrgm);

    fsigma = gsl_spline2d_alloc(T, size_vector_intrinsic_sigma, 
      size_vector_atsrgm);

    {
      FILE* file_intrinsic_sigma;
      file_intrinsic_sigma = fopen("sig_intr_grid.dat", "r");
      if(file_intrinsic_sigma == 0) {
        printf("Can't open intrinsic sigma file \n");
        exit(1);
      }
      for (int i=0; i<size_vector_intrinsic_sigma; i++) {
        fscanf(file_intrinsic_sigma, "%lf", &vector_intrinsic_sigma[i]);
      }
      fclose(file_intrinsic_sigma);
    }
    {
      FILE* file_atsrgm;
      file_atsrgm = fopen("l_sat_grid.dat", "r");
      if(file_atsrgm == 0) {
        printf("Can't open intrinsic true satellite richness file \n");
        exit(1);
      }
      for (int i=0; i<size_vector_atsrgm; i++) {
        fscanf(file_atsrgm, "%lf", &vector_atsrgm[i]);
      }
      fclose(file_atsrgm);
    }
    {
      FILE* file_alpha;
      file_alpha = fopen("skew_table.dat", "r");
      if(file_alpha == 0) {
        printf("Can't open skewness file file \n");
        exit(1);
      }
      double tmp;
      for (int i=0; i<size_vector_intrinsic_sigma; i++) {
        for (int j=0; j<size_vector_atsrgm; j++) {
          fscanf(file_alpha,"%lf",&tmp);
          gsl_spline2d_set(falpha, vector_alpha, i, j, tmp);
        }
      }
      fclose(file_alpha);
    }
    {
      FILE* file_sigma;
      file_sigma = fopen("sig_skew_table.dat", "r");
      if(file_sigma == 0) {
        printf("Can't open sigma file \n");
        exit(1);
      }
      double tmp;
      for (int i=0; i<size_vector_intrinsic_sigma; i++) {
        for (int j=0; j<size_vector_atsrgm; j++) {
          fscanf(file_sigma,"%lf",&tmp);
          gsl_spline2d_set(fsigma, vector_sigma, i, j, tmp);
        }
      }
      fclose(file_sigma);
    }

    gsl_spline2d_init(falpha, vector_intrinsic_sigma, vector_atsrgm, vector_alpha,
      size_vector_intrinsic_sigma, size_vector_atsrgm);
    if(falpha == NULL) {
      printf("malloc fail \n");
      exit(1);
    }

    gsl_spline2d_init(fsigma, vector_intrinsic_sigma, vector_atsrgm, vector_sigma,
      size_vector_intrinsic_sigma, size_vector_atsrgm);
    if(fsigma == NULL) {
      printf("malloc fail \n");
      exit(1);
    }
  }

  const double alpha=gsl_spline2d_eval(falpha,intrinsic_sigma,atsrgm, NULL,NULL);
  const double sigma=gsl_spline2d_eval(fsigma,intrinsic_sigma,atsrgm,NULL,NULL);
  const double x = 1.0/(M_SQRT2*abs(sigma));
  double result = exp(-(true_lambda-atsrgm)*(true_lambda-atsrgm)*x*x)*x/M_SQRTPI;
  result *= gsl_sf_erfc(-alpha*(true_lambda-atsrgm)*x);
  return result;
}

double
probability_observed_richness_given_true_richness(const double observed_lambda, 
void* params) {
  // The relation between observed and true richness depends on projection
  // effects and photometric noise

  // observed_lambda = observed richness
  // z = redshift
  // \mu = mean of the distribution
  // \sigma = spread of the distribution
  // tau = related to the skewness of the distribution

  // fmask = fraction of the distribution that lives in the left tail.
  // Physically: it has to do with the line of sight effects that suppress 
  // richness in small clusters near big ones

  // fprj = fraction of the distribution that lives in the right tail.
  // Physically: it has to with projection effects that may merge two halos
  // into a single one with combined richness
  Par_lambda_obs *p=(Par_lambda_obs *)params;

  const double true_lambda = p->true_lambda;
  const double z = p->z;
  static int first = 0;
  static double* vector_z;
  static double* vector_lambda;
  static double* vector_tau;
  static double* vector_mu;
  static double* vector_sigma;
  static double* vector_fmask;
  static double* vector_fprj;
  static gsl_spline2d* ftau;
  static gsl_spline2d* fmu;
  static gsl_spline2d* fsigma;
  static gsl_spline2d* ffmask;
  static gsl_spline2d* ffprj;

  if(first == 0) {
    first = 1;
    FILE* prj_params_file;
    prj_params_file = fopen("prj_params.txt", "r");
    if(prj_params_file == 0) {
      printf("Can't open prj_params file \n");
      exit(1);
    }

    const int size_vector_z = 5;
    const int size_vector_lambda = 19; 
    vector_z = (double*) malloc(size_vector_z*sizeof(double));
    if(vector_z == NULL) {
      printf("malloc fail \n");
      exit(1);
    }
    vector_lambda = (double*) malloc(size_vector_lambda*sizeof(double));
    if(vector_lambda == NULL) {
      printf("malloc fail \n");
      exit(1);
    }
    vector_tau = 
      (double*) malloc(size_vector_z*size_vector_lambda*sizeof(double));
    if(vector_tau == NULL) {
      printf("malloc fail \n");
      exit(1);
    }
    vector_mu = 
      (double*) malloc(size_vector_z*size_vector_lambda*sizeof(double));
    if(vector_mu == NULL) {
      printf("malloc fail \n");
      exit(1);
    }
    vector_sigma = 
      (double*) malloc(size_vector_z*size_vector_lambda*sizeof(double));
    if(vector_sigma == NULL) {
      printf("malloc fail \n");
      exit(1);
    }
    vector_fmask = 
      (double*) malloc(size_vector_z*size_vector_lambda*sizeof(double));
    if(vector_fmask == NULL) {
      printf("malloc fail \n");
      exit(1);
    }
    vector_fprj = 
      (double*) malloc(size_vector_z*size_vector_lambda*sizeof(double));
    if(vector_fprj == NULL) {
      printf("malloc fail \n");
      exit(1);
    }

    const gsl_interp2d_type *T = gsl_interp2d_bilinear;

    ftau = gsl_spline2d_alloc(T, size_vector_z, size_vector_lambda);
    if(ftau == NULL) {
      printf("malloc fail \n");
      exit(1);
    }
    fmu = gsl_spline2d_alloc(T, size_vector_z, size_vector_lambda);
    if(fmu == NULL) {
      printf("malloc fail \n");
      exit(1);
    }
    fsigma = gsl_spline2d_alloc(T, size_vector_z, size_vector_lambda);
    if(fsigma == NULL) {
      printf("malloc fail \n");
      exit(1);
    }
    ffmask = gsl_spline2d_alloc(T, size_vector_z, size_vector_lambda);
    if(ffmask == NULL) {
      printf("malloc fail \n");
      exit(1);
    }
    ffprj = gsl_spline2d_alloc(T, size_vector_z, size_vector_lambda);
    if(ffprj == NULL) {
      printf("malloc fail \n");
      exit(1);
    }
    //# z,mean_l,tau,mu,sig_prj,fmask_prj,fprj_prj
    double tmp, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7;
    for (int i=0; i<size_vector_z; i++) {
      for (int j=0; j<size_vector_lambda; j++) {
        fscanf(prj_params_file,"%lf %lf %lf %lf %lf %lf %lf", &tmp, &tmp2, 
          &tmp3, &tmp4, &tmp5, &tmp6, &tmp7);
        gsl_spline2d_set(ftau, vector_tau, i, j, tmp3);
        gsl_spline2d_set(fmu, vector_mu, i, j, tmp4);
        gsl_spline2d_set(fsigma, vector_sigma, i, j, tmp5);
        gsl_spline2d_set(ffmask, vector_fmask, i, j, tmp6);
        gsl_spline2d_set(ffprj, vector_fprj, i, j, tmp7);
        if(i == 0) {
          vector_lambda[j] = tmp2;
        }
        if(j == 0) {
          vector_z[i] = tmp;
        }
      }
    }

    gsl_spline2d_init(ftau, vector_z, vector_lambda, vector_tau, 
      size_vector_z, size_vector_lambda);
    if(ftau == NULL) {
      printf("malloc fail \n");
      exit(1);
    }
    gsl_spline2d_init(fmu, vector_z, vector_lambda, vector_mu, 
      size_vector_z, size_vector_lambda);
    if(fmu == NULL) {
      printf("malloc fail \n");
      exit(1);
    }
    gsl_spline2d_init(fsigma, vector_z, vector_lambda, vector_sigma, 
      size_vector_z, size_vector_lambda);
    if(fsigma == NULL) {
      printf("malloc fail \n");
      exit(1);
    }
    gsl_spline2d_init(ffmask, vector_z, vector_lambda, vector_fmask, 
      size_vector_z, size_vector_lambda);
    if(ffmask == NULL) {
      printf("malloc fail \n");
      exit(1);
    }
    gsl_spline2d_init(ffprj, vector_z, vector_lambda, vector_fprj, 
      size_vector_z, size_vector_lambda);
    if(ffprj == NULL) {
      printf("malloc fail \n");
      exit(1);
    }
  }

  const double tau = gsl_spline2d_eval(ftau, z, true_lambda, NULL, NULL);
  const double mu = gsl_spline2d_eval(fmu, z, true_lambda, NULL, NULL);
  const double sigma = gsl_spline2d_eval(fsigma, z, true_lambda, NULL, NULL);
  const double fmask = gsl_spline2d_eval(ffmask, z, true_lambda, NULL, NULL);
  const double fprj = gsl_spline2d_eval(ffprj, z, true_lambda, NULL, NULL);
  const double x = 1.0/(M_SQRT2*abs(sigma));
  const double y = 1.0/true_lambda;
  const double j = exp(0.5*tau*(2*mu+tau*sigma*sigma-2*observed_lambda));

  double r0 = exp((observed_lambda-mu)*(observed_lambda-mu)*x*x);
  r0 *= (1-fmask)*(1-fprj)*x/M_SQRTPI;

  double r1 = 0.5*((1-fmask)*fprj*tau + fmask*fprj*y)*j;
  r1 *= gsl_sf_erfc((mu+tau*sigma*sigma-observed_lambda)*x);

  double r2 = 0.5*fmask*y;
  r2 *= gsl_sf_erfc((mu-observed_lambda-true_lambda)*x) - 
  gsl_sf_erf((mu-observed_lambda)*x);

  double r3 = 0.5*fmask*fprj*y*exp(-tau*true_lambda)*j;
  r3 *= gsl_sf_erfc((mu+tau*sigma*sigma-observed_lambda-true_lambda)*x);

  return r0 + r1 + r2 + r3;
}
double int_gsl_integrate_medium_precision(double (*func)(double, void*),void *arg,double a, double b, double *error, int niter)
{
  double res, err;
  gsl_integration_cquad_workspace *w = gsl_integration_cquad_workspace_alloc(niter);
  gsl_function F;
  F.function = func;
  F.params  = arg;
 gsl_integration_cquad(&F,a,b,0,1.e-3,w,&res,&err,0);
  if(NULL!=error)
    *error=err;
  gsl_integration_cquad_workspace_free(w);
  return res;
}
double int_ltrue_probability_observed_richness_given_mass(double true_lambda, void * params){
  Par_lambda_obs *p=(Par_lambda_obs *)params;
  return probability_observed_richness_given_true_richness(p->obs_lambda, params)*probability_true_richness_given_mass(true_lambda,params);
}
double int_lobs_probability_observed_richness_given_mass(double obs_lambda, void * params){
  Par_lambda_obs *p=(Par_lambda_obs *)params;
  p->obs_lambda = obs_lambda;
  return int_gsl_integrate_medium_precision(int_ltrue_probability_observed_richness_given_mass,params,1.,500.,NULL,1000);
}
double probability_observed_richness_given_mass(int N_lambda, void * params) {
  Par_lambda_obs *p=(Par_lambda_obs *)params;
  p->N_lambda = N_lambda;
  //printf("calculating P(l_obs =[%.1f,%.1f]|M=%e)\n",p->obs_lambda_min[N_lambda],p->obs_lambda_max[N_lambda],p->mass);
  double int_l = int_gsl_integrate_medium_precision(int_lobs_probability_observed_richness_given_mass,params,p->obs_lambda_min[N_lambda],p->obs_lambda_max[N_lambda],NULL,1000);
  //printf("%.1f\n",int_l);
  return int_l;
}

int main() {
  const double true_z = 0.3;
  const double true_lambda = 6;
  const double mass_M1 = pow(10.0, 12.0593); // Mass scale with ~ one satelite
  const double mass_min = mass_M1/15.0;
  const double intrinsic_alpha = 0.7;
  const double intrinsic_sigma = 0.3;
  const double observed_lambda = 8;

  Par_lambda_obs p;
  p.mass = pow(10.0, 13);
  p.mass_min = mass_min;
  p.mass_M1 = mass_M1;
  p.intrinsic_alpha = intrinsic_alpha;
  p.intrinsic_sigma = intrinsic_sigma;
  p.z = true_z;
  p.true_lambda = true_lambda;
  p.obs_lambda_min[0] = 5;
  p.obs_lambda_max[0] = 15;
  p.obs_lambda_min[1] = 5;
  p.obs_lambda_max[1] = 20;
  printf("probability_observed_richness_given_mass(l_obs =[%.1f,%.1f]|M)\n",p.obs_lambda_min[0],p.obs_lambda_max[0]);
  for (double m = 12.0; m < 15.5; m+=0.5){
    p.mass = pow(10.,m);
    printf("%.1f %e\n",m,  probability_observed_richness_given_mass(0,&p));
  }
  printf("probability_observed_richness_given_mass(l_obs =[%.1f,%.1f]|M)\n",p.obs_lambda_min[1],p.obs_lambda_max[1]);
  for (double m = 12.0; m < 15.5; m+=0.5){
    p.mass = pow(10.,m);
    printf("%.1f %e\n",m,  probability_observed_richness_given_mass(1,&p));
  }
  {
    double r1 = probability_true_richness_given_mass(true_lambda, &p);
    printf("%lf \n",r1);

    double r2 = probability_true_richness_given_mass(0.5*true_lambda, &p);
    printf("%lf \n",r2);
  }

  {
    double r3 =   
      probability_observed_richness_given_true_richness(observed_lambda,&p);
    printf("%lf \n",r3);

    double r4 =
      probability_observed_richness_given_true_richness(observed_lambda,&p);
    printf("%lf \n",r4);
  }

  return 0;
}

// VM ends