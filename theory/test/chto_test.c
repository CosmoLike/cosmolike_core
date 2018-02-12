#include "util.h"

void test_richness_true_given_mass(double M,Par_lambda_obs p){
  char result[50];
  sprintf(result, "p_richness_true_given_mass_%.3f.data", M);
  p.mass=pow(10.0,M);
  int size_of_arr=150;
  double p_richness_obs_given_true[size_of_arr][2];
  for (int i=0; i<size_of_arr; ++i){
    double lambda_obs = i+1;
    p_richness_obs_given_true[i][0] = lambda_obs;
    p_richness_obs_given_true[i][1] = probability_true_richness_given_mass(lambda_obs, &p);
  }
  char* filename = result; 
  FILE *f = fopen(filename, "wb");
  if(f == NULL){
    printf("Error");   
    exit(1);             
  }

  for(size_t n = 0; n < size_of_arr; ++n){
      for(size_t col=0; col < 2;++col){
        fprintf(f, "%.5f",p_richness_obs_given_true[n][col]);
        fprintf(f, " ");
      }
      fprintf(f, "\n");
  }
  fclose(f);
 
    
  return;

}

void test_richness_obs_given_true(double truelambda,Par_lambda_obs p){
  char result[50];
  p.true_lambda=truelambda;
  sprintf(result, "p_richness_obs_given_true_%.0f.data", truelambda);
  int size_of_arr=200;
  double p_richness_obs_given_true[size_of_arr][2];
  for (int i=0; i<size_of_arr; ++i){
    double lambda_obs = i+10;
    p_richness_obs_given_true[i][0] = lambda_obs;
    p_richness_obs_given_true[i][1] = probability_observed_richness_given_true_richness(lambda_obs, &p);
  }
  char* filename = result; 
  FILE *f = fopen(filename, "wb");
  if(f == NULL){
    printf("Error");   
    exit(1);             
  }

  for(size_t n = 0; n < size_of_arr; ++n){
      for(size_t col=0; col < 2;++col){
        fprintf(f, "%.5f",p_richness_obs_given_true[n][col]);
        fprintf(f, " ");
      }
      fprintf(f, "\n");
  }
  fclose(f);
 
    
  return;

}
int main(){
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
 
 test_richness_true_given_mass(1.8286+12.0593,p); 
 test_richness_true_given_mass(2.397+12.0593,p); 
 test_richness_true_given_mass(2.827+12.0593,p); 
 test_richness_obs_given_true(26,p);
 test_richness_obs_given_true(58,p);
 test_richness_obs_given_true(100,p);



}
