#ifndef HEADER_FILE
#define HEADER_FILE

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

double probability_true_richness_given_mass(const double, void *);
double probability_observed_richness_given_true_richness(const double, void*);
double int_gsl_integrate_medium_precision(double (*func)(double, void*),void *,double , double , double *, int);
double int_ltrue_probability_observed_richness_given_mass(double, void *);
double int_ltrue_probability_observed_richness_given_mass(double, void *);
double int_lobs_probability_observed_richness_given_mass(double, void * );
double probability_observed_richness_given_mass(int, void * );
#endif
