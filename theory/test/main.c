#include "util.h"
int main() {
  const double true_z = 0.3;
  const double true_lambda = 6;
  const double mass_M1 = pow(10.0, 12.0593); // Mass scale with ~ one satelite
  const double mass_min = mass_M1/15.0;
  const double intrinsic_alpha = 0.7;
  const double intrinsic_sigma = 0.2;
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

