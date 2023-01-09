//#include "./covariances_binned_simple.h"
/* This new file aims to rewrite the bin-averaged real-space fullsky 
  3x2pt covariances in covariances_real_binned_fullsky.c, but also 
  provide bin-averaged mixed/real(fullsky)/fourier space 6x2pt cov 
  in a unified way
  
  Current status: 
    - binned-Fourier-space still needs testing;
    = real-space: good
    - mixed-space: not sure!
  
  -- Xiao Fang 
*/
#define __NN_SNR_FLAG__ 1
/************ ===================== ************/
/************ Function Declarations ************/
/************ ===================== ************/

// Note: Some functions are not declared separately

/******    Level 0 functions    ******/
/****** Full-sky Covariance API ******/
/****** ----------------------- ******/

/* real-space, Logarithmic bins in theta (21 blocks) */

// 3x2pt (6 blocks)
void cov_shear_shear_real_binned_fullsky(double **cov, double **covNG, 
  int z1,int z2,int z3,int z4,int pm1,int pm2,int FLAG_NG,
  double *theta,double *dtheta);
void cov_gl_shear_real_binned_fullsky(double **cov, double **covNG, 
  int z1,int z2,int z3,int z4,int pm,int FLAG_NG,double *theta,double *dtheta);
void cov_cl_shear_real_binned_fullsky(double **cov, double **covNG, 
  int z1,int z2,int z3,int z4,int pm,int FLAG_NG,double *theta,double *dtheta);
void cov_cl_gl_real_binned_fullsky(double **cov, double **covNG, 
  int z1,int z2,int z3,int z4, int FLAG_NG, double *theta, double *dtheta);
void cov_gl_gl_real_binned_fullsky(double **cov, double **covNG, 
  int z1,int z2,int z3,int z4, int FLAG_NG, double *theta, double *dtheta);
void cov_cl_cl_real_binned_fullsky(double **cov, double **covNG, 
  int z1,int z2,int z3,int z4, int FLAG_NG, double *theta, double *dtheta);

// 5x2pt (9 blocks)
void cov_ks_shear_real_binned_fullsky(double **cov, double **covNG, 
  int z1,int z3,int z4,int pm, int FLAG_NG, double *theta, double *dtheta);
void cov_gk_shear_real_binned_fullsky(double **cov, double **covNG, 
  int z1,int z3,int z4,int pm, int FLAG_NG, double *theta, double *dtheta);
void cov_gk_gl_real_binned_fullsky(double **cov, double **covNG, 
  int z1,int z3,int z4, int FLAG_NG, double *theta, double *dtheta);
void cov_ks_gl_real_binned_fullsky(double **cov, double **covNG, 
  int z1,int z3,int z4, int FLAG_NG, double *theta, double *dtheta);
void cov_gk_cl_real_binned_fullsky(double **cov, double **covNG, 
  int z1,int z3,int z4, int FLAG_NG, double *theta, double *dtheta);
void cov_ks_cl_real_binned_fullsky(double **cov, double **covNG, 
  int z1,int z3,int z4, int FLAG_NG, double *theta, double *dtheta);
void cov_gk_gk_real_binned_fullsky(double **cov, double **covNG, 
  int z1,int z2, int FLAG_NG, double *theta, double *dtheta);
void cov_ks_gk_real_binned_fullsky(double **cov, double **covNG, 
  int z1,int z2, int FLAG_NG, double *theta, double *dtheta);
void cov_ks_ks_real_binned_fullsky(double **cov, double **covNG, 
  int z1,int z2, int FLAG_NG, double *theta, double *dtheta);

// kk x {5x2pt + kk} (6 blocks)
void cov_kk_shear_real_binned_fullsky(double **cov, double **covNG, int z3,int z4,int pm, int FLAG_NG, double *theta, double *dtheta);
void cov_kk_gl_real_binned_fullsky(double **cov, double **covNG, int z3,int z4, int FLAG_NG, double *theta, double *dtheta);
void cov_kk_cl_real_binned_fullsky(double **cov, double **covNG, int z3,int z4, int FLAG_NG, double *theta, double *dtheta);
void cov_kk_gk_real_binned_fullsky(double **cov, double **covNG, int zl, int FLAG_NG, double *theta, double *dtheta);
void cov_kk_ks_real_binned_fullsky(double **cov, double **covNG, int zs, int FLAG_NG, double *theta, double *dtheta);
void cov_kk_kk_real_binned_fullsky(double **cov, double **covNG, int FLAG_NG, double *theta, double *dtheta);

/* mixed space, theta & ell both in Logarithmic bins (5 blocks) */ 

// kk x {5x2pt} (5 blocks)
void cov_kk_shear_mix_binned_fullsky(double **cov, double **covNG, 
  int z3,int z4,int pm,int FLAG_NG,double *theta, double *dtheta, double *ell);
void cov_kk_gl_mix_binned_fullsky(double **cov, double **covNG, 
  int z3,int z4, int FLAG_NG, double *theta, double *dtheta, double *ell);
void cov_kk_cl_mix_binned_fullsky(double **cov, double **covNG, 
  int z3,int z4, int FLAG_NG,double *theta, double *dtheta, double *ell);
void cov_kk_gk_mix_binned_fullsky(double **cov, double **covNG, 
  int zl, int FLAG_NG, double *theta, double *dtheta, double *ell);
void cov_kk_ks_mix_binned_fullsky(double **cov, double **covNG, 
  int zs, int FLAG_NG, double *theta, double *dtheta, double *ell);

/* mixed space, theta in Logarithmic bins, ell in band powers (5 blocks) */

// kk x {5x2pt} (5 blocks)
void cov_kk_shear_mix_banded_fullsky(double **cov, double **covNG, 
  int z3,int z4,int pm,int FLAG_NG,double *theta, double *dtheta, int **bindef);
void cov_kk_gl_mix_banded_fullsky(double **cov, double **covNG, 
  int z3,int z4, int FLAG_NG, double *theta, double *dtheta, int **bindef);
void cov_kk_cl_mix_banded_fullsky(double **cov, double **covNG, 
  int z3,int z4, int FLAG_NG, double *theta, double *dtheta, int **bindef);
void cov_kk_gk_mix_banded_fullsky(double **cov, double **covNG, 
  int zl, int FLAG_NG, double *theta, double *dtheta, int **bindef);
void cov_kk_ks_mix_banded_fullsky(double **cov, double **covNG, 
  int zs, int FLAG_NG, double *theta, double *dtheta, int **bindef);

/* Fourier space, ell in Logarithmic bins (21 blocks) */

// 3x2pt (6 blocks)
void cov_shear_shear_fourier_binned(double **cov, double **covNG, 
  int z1,int z2,int z3,int z4, int FLAG_NG, double *ell);
void cov_gl_shear_fourier_binned(double **cov, double **covNG, 
  int z1,int z2,int z3,int z4, int FLAG_NG, double *ell);
void cov_cl_shear_fourier_binned(double **cov, double **covNG, 
  int z1,int z2,int z3,int z4, int FLAG_NG, double *ell);
void cov_cl_gl_fourier_binned(double **cov, double **covNG, 
  int z1,int z2,int z3,int z4, int FLAG_NG, double *ell);
void cov_gl_gl_fourier_binned(double **cov, double **covNG, 
  int z1,int z2,int z3,int z4, int FLAG_NG, double *ell);
void cov_cl_cl_fourier_binned(double **cov, double **covNG, 
  int z1,int z2,int z3,int z4, int FLAG_NG, double *ell);

// 5x2pt (9 blocks)
void cov_ks_shear_fourier_binned(double **cov, double **covNG, 
  int z1,int z3,int z4, int FLAG_NG, double *ell);
void cov_gk_shear_fourier_binned(double **cov, double **covNG, 
  int z1,int z3,int z4, int FLAG_NG, double *ell);
void cov_gk_gl_fourier_binned(double **cov, double **covNG, 
  int z1,int z3,int z4, int FLAG_NG, double *ell);
void cov_ks_gl_fourier_binned(double **cov, double **covNG, 
  int z1,int z3,int z4, int FLAG_NG, double *ell);
void cov_gk_cl_fourier_binned(double **cov, double **covNG, 
  int z1,int z3,int z4, int FLAG_NG, double *ell);
void cov_ks_cl_fourier_binned(double **cov, double **covNG, 
  int z1,int z3,int z4, int FLAG_NG, double *ell);
void cov_gk_gk_fourier_binned(double **cov, double **covNG, 
  int z1,int z2, int FLAG_NG, double *ell);
void cov_ks_gk_fourier_binned(double **cov, double **covNG, 
  int z1,int z2, int FLAG_NG, double *ell);
void cov_ks_ks_fourier_binned(double **cov, double **covNG, 
  int z1,int z2, int FLAG_NG, double *ell);

// kk x {5x2pt + kk} (6 blocks)
void cov_kk_shear_fourier_binned(double **cov, double **covNG, 
  int z3,int z4, int FLAG_NG, double *ell);
void cov_kk_gl_fourier_binned(double **cov, double **covNG, 
  int z3,int z4, int FLAG_NG, double *ell);
void cov_kk_cl_fourier_binned(double **cov, double **covNG, 
  int z3,int z4, int FLAG_NG, double *ell);
void cov_kk_gk_fourier_binned(double **cov, double **covNG, 
  int zl, int FLAG_NG, double *ell);
void cov_kk_ks_fourier_binned(double **cov, double **covNG, 
  int zs, int FLAG_NG, double *ell);
void cov_kk_kk_fourier_binned(double **cov, double **covNG, 
  int FLAG_NG, double *ell);

/* Fourier space, ell in band powers (1 blocks)*/

void cov_kk_kk_fourier_banded(double **cov, double **covNG, 
  int FLAG_NG, double *ell, int **bindef);

/******           Level 1 functions           ******/
/****** Full-sky Covariance Template Routines ******/
/****** ------------------------------------- ******/

// Real-space template function
// theta in Logarithmic bins
void cov_real_binned_fullsky(double **cov, double **covNG, char *realcov_type, 
  int *z_ar, int FLAG_NG, double *theta, double *dtheta);

// Fourier-space template function
// ell in Logarithmic bins
void cov_fourier_binned(double **cov, double **covNG, char *cov_type, 
  int *z_ar, int FLAG_NG, double *ell);
// ell in band power and Logarithmic bins
void cov_fourier_banded(double **cov, double **covNG, char *cov_type, 
  int *z_ar, int FLAG_NG, double *ell, int **bindef);

// Mixed-space template functions
// theta & ell both in Logarithmic bins
void cov_mix_binned_fullsky(double **cov, double **covNG, char *mixcov_type, 
  int *z_ar, int FLAG_NG, double *theta, double *dtheta, double *ell);
// theta in Logarithmic bins, ell in band powers
void cov_mix_banded_fullsky(double **cov, double **covNG, char *mixcov_type, 
  int *z_ar, int FLAG_NG, double *theta, double *dtheta, int **bindef);

/******           Level 2 functions          ******/
/****** Covariance Matrix Calculation Recipe ******/
/****** ------------------------------------ ******/

/* Gaussian Covariance Matrix (21+6 functions) */

// 3x2pt without pure noise term (3 functions)
double func_for_cov_G_shear_noNN(double l, int *ar);
double func_for_cov_G_cl_noNN(double l, int *ar);
double func_for_cov_G_gl_noNN(double l, int *ar);
double func_for_cov_G_gk_noNN(double l, int *ar);
double func_for_cov_G_ks_noNN(double l, int *ar);
// 3x2pt Gaussian cov pure noise term on the diagonal, without masking effect
void pure_noise_xipm_xipm(int *z_ar, double *theta, double *dtheta, double **N);
void pure_noise_gl_gl(int *z_ar, double *theta, double *dtheta, double **N);
void pure_noise_cl_cl(int *z_ar, double *theta, double *dtheta, double **N);
void pure_noise_gk_gk(int *z_ar, double *theta, double *dtheta, double **N);
void pure_noise_ks_ks(int *z_ar, double *theta, double *dtheta, double **N);

// 3x2pt, pure noise included (6 functions)
double func_for_cov_G_shear(double l, int *ar);
double func_for_cov_G_cl(double l, int *ar);
double func_for_cov_G_gl(double l, int *ar);
double func_for_cov_G_cl_shear(double l, int *ar);
double func_for_cov_G_cl_gl(double l, int *ar);
double func_for_cov_G_gl_shear(double l, int *ar);

// 5x2pt, pure noise included (9 functions)
double func_for_cov_G_gk(double l, int *ar);
double func_for_cov_G_ks(double l, int *ar);
double func_for_cov_G_gk_shear(double l, int *ar);
double func_for_cov_G_gk_gl(double l, int *ar);
double func_for_cov_G_gk_cl(double l, int *ar);
double func_for_cov_G_ks_shear(double l, int *ar);
double func_for_cov_G_ks_gl(double l, int *ar);
double func_for_cov_G_ks_cl(double l, int *ar);
double func_for_cov_G_ks_gk(double l, int *ar);

// 6x2pt, pure noise included (6 functions)
double func_for_cov_G_kk_shear(double l, int *ar);
double func_for_cov_G_kk_gl(double l, int *ar);
double func_for_cov_G_kk_cl(double l, int *ar);
double func_for_cov_G_kk_gk(double l, int *ar);
double func_for_cov_G_kk_ks(double l, int *ar);
double func_for_cov_G_kk(double l, int *ar); // not used?

/* Non-Gaussian Covariance Matrix */

// Look up tables

// 3x2pt (6 functions)
double bin_cov_NG_shear_shear(double l1,double l2, int *z_ar);
double bin_cov_NG_gl_gl(double l1,double l2, int *z_ar);
double bin_cov_NG_cl_cl(double l1,double l2, int *z_ar);
double bin_cov_NG_cl_shear(double l1,double l2, int *z_ar);
double bin_cov_NG_cl_gl(double l1,double l2, int *z_ar);
double bin_cov_NG_gl_shear(double l1,double l2, int *z_ar);

// 5x2pt (9 functions)
double bin_cov_NG_gk_shear(double l1,double l2, int *z_ar);
double bin_cov_NG_ks_shear(double l1,double l2, int *z_ar);
double bin_cov_NG_gk_gl(double l1,double l2, int *z_ar);
double bin_cov_NG_ks_gl(double l1,double l2, int *z_ar);
double bin_cov_NG_gk_cl(double l1,double l2, int *z_ar);
double bin_cov_NG_ks_cl(double l1,double l2, int *z_ar);
double bin_cov_NG_gk_gk(double l1,double l2, int *z_ar);
double bin_cov_NG_ks_gk(double l1,double l2, int *z_ar);
double bin_cov_NG_ks_ks(double l1,double l2, int *z_ar);

// 6x2pt (6 functions)
double bin_cov_NG_kk_shear(double l1,double l2, int *z_ar);
double bin_cov_NG_kk_gl(double l1,double l2, int *z_ar);
double bin_cov_NG_kk_cl(double l1,double l2, int *z_ar);
double bin_cov_NG_kk_gk(double l1,double l2, int *z_ar);
double bin_cov_NG_kk_ks(double l1,double l2, int *z_ar);
double bin_cov_NG_kk_kk(double l1,double l2, int *z_ar);

// tabulate function
void tabulate_covNG(char *cov_type, int *z_ar, int Ntab, double **table);

// bridge functions to unpack z_ar array to call the fourier cov routines 

// 3x2pt (6 functions)
double cov_NG_shear_shear_bridge(double ll1, double ll2, int *z_ar);
double cov_NG_gl_gl_bridge(double ll1, double ll2, int *z_ar);
double cov_NG_cl_cl_bridge(double ll1, double ll2, int *z_ar);
double cov_NG_cl_gl_bridge(double ll1, double ll2, int *z_ar);
double cov_NG_cl_shear_bridge(double ll1, double ll2, int *z_ar);
double cov_NG_gl_shear_bridge(double ll1, double ll2, int *z_ar);
// 5x2pt (9 functions)
double cov_NG_gk_shear_bridge(double ll1, double ll2, int *z_ar);
double cov_NG_ks_shear_bridge(double ll1, double ll2, int *z_ar);
double cov_NG_gk_gl_bridge(double ll1, double ll2, int *z_ar);
double cov_NG_ks_gl_bridge(double ll1, double ll2, int *z_ar);
double cov_NG_gk_cl_bridge(double ll1, double ll2, int *z_ar);
double cov_NG_ks_cl_bridge(double ll1, double ll2, int *z_ar);
double cov_NG_gk_gk_bridge(double ll1, double ll2, int *z_ar);
double cov_NG_ks_gk_bridge(double ll1, double ll2, int *z_ar);
double cov_NG_ks_ks_bridge(double ll1, double ll2, int *z_ar);
// 6x2pt (5 functions)
double cov_NG_kk_shear_bridge(double ll1, double ll2, int *z_ar);
double cov_NG_kk_gl_bridge(double ll1, double ll2, int *z_ar);
double cov_NG_kk_cl_bridge(double ll1, double ll2, int *z_ar);
double cov_NG_kk_gk_bridge(double ll1, double ll2, int *z_ar);
double cov_NG_kk_ks_bridge(double ll1, double ll2, int *z_ar);
double cov_NG_kk_kk_bridge(double ll1, double ll2, int *z_ar);

/* Helper Functions */

// look-up tables for full-sky transform Legendre polys
double Pl_tab(int itheta, int ell);
double Pl2_tab(int itheta, int ell);
double Glplus_tab(int itheta, int ell);
double Glminus_tab(int itheta, int ell);
// Planck CMB lensing band power binning matrix (with correction)
double CMB_lensing_mat(int ibp, int L);
double CMB_lensing_mat_with_corr(int ibp, int L);
// CMB beam smoothing kernel
double GaussianBeam(double fwhm, int ell, double ell_left, double ell_right);
// masking effect on pure noise
double w_mask(double theta_min, int col);
// healpix pixel window function
double w_pixel(int ell);

/************ ======================== ************/
/************ Function Implementations ************/
/************ ======================== ************/

/****** Template Routines ******/

// Real-space template function
// theta in Logarithmic bins
void cov_real_binned_fullsky(double **cov, double **covNG, char *realcov_type, 
  int *z_ar, int FLAG_NG, double *theta, double *dtheta)
{
  int i,j;
  static int LMAX = 50000;
  int CMB_smooth_1 = 0, CMB_smooth_2 = 0; // flag for CMB smoothing kernel

  //double N[like.Ntheta];
  double ** N = 0;
  N = create_double_matrix(0, like.Ntheta-1, 0, like.Ntheta-1);
  for(i=0; i<like.Ntheta; i++){
    for(j=0; j<like.Ntheta; j++){
      N[i][j] = 0.;
    }
  }

  double (*func_for_cov_G)(double, int*);
  double (*func_bin_cov_NG)(double, double, int*);
  double (*func_P1)(int, int);
  double (*func_P2)(int, int);
  int Icol = 0;

  // 3x2pt
  if(strcmp(realcov_type, "xi+_xi+")==0) {
    func_for_cov_G  = &func_for_cov_G_shear_noNN;
    func_bin_cov_NG = &bin_cov_NG_shear_shear;
    func_P1 = &Glplus_tab;
    func_P2 = &Glplus_tab;
    pure_noise_xipm_xipm(z_ar, theta, dtheta, N); //C+- doesn't have the diagonal shot noise term
    Icol = 0;
  } else if(strcmp(realcov_type, "xi-_xi-")==0) {
    func_for_cov_G  = &func_for_cov_G_shear_noNN;
    func_bin_cov_NG = &bin_cov_NG_shear_shear;
    func_P1 = &Glminus_tab;
    func_P2 = &Glminus_tab;
    pure_noise_xipm_xipm(z_ar, theta, dtheta, N); //C+- doesn't have the diagonal shot noise term
    Icol = 0;
  } else if(strcmp(realcov_type, "xi+_xi-")==0) {
    func_for_cov_G  = &func_for_cov_G_shear_noNN;
    func_bin_cov_NG = &bin_cov_NG_shear_shear;
    func_P1 = &Glplus_tab;
    func_P2 = &Glminus_tab;
  } else if(strcmp(realcov_type, "xi-_xi+")==0) {
    func_for_cov_G  = &func_for_cov_G_shear_noNN;
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
    func_for_cov_G  = &func_for_cov_G_gl_noNN;
    func_bin_cov_NG = &bin_cov_NG_gl_gl;
    func_P1 = &Pl2_tab;
    func_P2 = &Pl2_tab;
    pure_noise_gl_gl(z_ar, theta, dtheta, N);
    Icol = 0;
  } else if(strcmp(realcov_type, "cl_gl")==0) {
    func_for_cov_G  = &func_for_cov_G_cl_gl;
    func_bin_cov_NG = &bin_cov_NG_cl_gl;
    func_P1 = &Pl_tab;
    func_P2 = &Pl2_tab;
  } else if(strcmp(realcov_type, "cl_cl")==0) {
    func_for_cov_G  = &func_for_cov_G_cl_noNN;
    func_bin_cov_NG = &bin_cov_NG_cl_cl;
    func_P1 = &Pl_tab;
    func_P2 = &Pl_tab;
    pure_noise_cl_cl(z_ar, theta, dtheta, N);
    Icol = 0;
  // 5x2pt
  } else if(strcmp(realcov_type, "gk_xi+")==0) {
    func_for_cov_G  = &func_for_cov_G_gk_shear;
    func_bin_cov_NG = &bin_cov_NG_gk_shear;
    func_P1 = &Pl_tab;
    func_P2 = &Glplus_tab;
    CMB_smooth_1 = 1;
    CMB_smooth_2 = 0;
  } else if(strcmp(realcov_type, "gk_xi-")==0) {
    func_for_cov_G  = &func_for_cov_G_gk_shear;
    func_bin_cov_NG = &bin_cov_NG_gk_shear;
    func_P1 = &Pl_tab;
    func_P2 = &Glminus_tab;
    CMB_smooth_1 = 1;
    CMB_smooth_2 = 0;
  } else if(strcmp(realcov_type, "ks_xi+")==0) {
    func_for_cov_G  = &func_for_cov_G_ks_shear;
    func_bin_cov_NG = &bin_cov_NG_ks_shear;
    func_P1 = &Pl2_tab;
    func_P2 = &Glplus_tab;
    CMB_smooth_1 = 1;
    CMB_smooth_2 = 0;
  } else if(strcmp(realcov_type, "ks_xi-")==0) {
    func_for_cov_G  = &func_for_cov_G_ks_shear;
    func_bin_cov_NG = &bin_cov_NG_ks_shear;
    func_P1 = &Pl2_tab;
    func_P2 = &Glminus_tab;
    CMB_smooth_1 = 1;
    CMB_smooth_2 = 0;
  } else if(strcmp(realcov_type, "gk_gl")==0) {
    func_for_cov_G  = &func_for_cov_G_gk_gl;
    func_bin_cov_NG = &bin_cov_NG_gk_gl;
    func_P1 = &Pl_tab;
    func_P2 = &Pl2_tab;
    CMB_smooth_1 = 1;
    CMB_smooth_2 = 0;
  } else if(strcmp(realcov_type, "ks_gl")==0) {
    func_for_cov_G  = &func_for_cov_G_ks_gl;
    func_bin_cov_NG = &bin_cov_NG_ks_gl;
    func_P1 = &Pl2_tab;
    func_P2 = &Pl2_tab;
    CMB_smooth_1 = 1;
    CMB_smooth_2 = 0;
  } else if(strcmp(realcov_type, "gk_cl")==0) {
    func_for_cov_G  = &func_for_cov_G_gk_cl;
    func_bin_cov_NG = &bin_cov_NG_gk_cl;
    func_P1 = &Pl_tab;
    func_P2 = &Pl_tab;
    CMB_smooth_1 = 1;
    CMB_smooth_2 = 0;
  } else if(strcmp(realcov_type, "ks_cl")==0) {
    func_for_cov_G  = &func_for_cov_G_ks_cl;
    func_bin_cov_NG = &bin_cov_NG_ks_cl;
    func_P1 = &Pl2_tab;
    func_P2 = &Pl_tab;
    CMB_smooth_1 = 1;
    CMB_smooth_2 = 0;
  } else if(strcmp(realcov_type, "gk_gk")==0) {
    func_for_cov_G  = &func_for_cov_G_gk_noNN;
    func_bin_cov_NG = &bin_cov_NG_gk_gk;
    func_P1 = &Pl_tab;
    func_P2 = &Pl_tab;
    CMB_smooth_1 = 1;
    CMB_smooth_2 = 1;
    pure_noise_gk_gk(z_ar, theta, dtheta, N);
    Icol = 2; // 0-w/o annulus; 2-w/annulus
  } else if(strcmp(realcov_type, "ks_gk")==0) {
    func_for_cov_G  = &func_for_cov_G_ks_gk;
    func_bin_cov_NG = &bin_cov_NG_ks_gk;
    func_P1 = &Pl2_tab;
    func_P2 = &Pl_tab;
    CMB_smooth_1 = 1;
    CMB_smooth_2 = 1;
  } else if(strcmp(realcov_type, "ks_ks")==0) {
    func_for_cov_G  = &func_for_cov_G_ks_noNN;
    func_bin_cov_NG = &bin_cov_NG_ks_ks;
    func_P1 = &Pl2_tab;
    func_P2 = &Pl2_tab;
    CMB_smooth_1 = 1;
    CMB_smooth_2 = 1;
    pure_noise_ks_ks(z_ar, theta, dtheta, N);
    Icol = 2; // 0-w/o annulus; 2-w/annulus
  // 6x2pt
  } else if(strcmp(realcov_type, "kk_xi+")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_shear;
    func_bin_cov_NG = &bin_cov_NG_kk_shear;
    func_P1 = &Pl_tab;
    func_P2 = &Glplus_tab;
    CMB_smooth_1 = 0;
    CMB_smooth_2 = 0;
  } else if(strcmp(realcov_type, "kk_xi-")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_shear;
    func_bin_cov_NG = &bin_cov_NG_kk_shear;
    func_P1 = &Pl_tab;
    func_P2 = &Glminus_tab;
    CMB_smooth_1 = 0;
    CMB_smooth_2 = 0;
  } else if(strcmp(realcov_type, "kk_gl")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_gl;
    func_bin_cov_NG = &bin_cov_NG_kk_gl;
    func_P1 = &Pl_tab;
    func_P2 = &Pl2_tab;
    CMB_smooth_1 = 0;
    CMB_smooth_2 = 0;
  } else if(strcmp(realcov_type, "kk_cl")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_cl;
    func_bin_cov_NG = &bin_cov_NG_kk_cl;
    func_P1 = &Pl_tab;
    func_P2 = &Pl_tab;
    CMB_smooth_1 = 0;
    CMB_smooth_2 = 0;
  } else if(strcmp(realcov_type, "kk_gk")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_gk;
    func_bin_cov_NG = &bin_cov_NG_kk_gk;
    func_P1 = &Pl_tab;
    func_P2 = &Pl_tab;
    CMB_smooth_1 = 0;
    CMB_smooth_2 = 1;
  } else if(strcmp(realcov_type, "kk_ks")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_ks;
    func_bin_cov_NG = &bin_cov_NG_kk_ks;
    func_P1 = &Pl_tab;
    func_P2 = &Pl2_tab;
    CMB_smooth_1 = 0;
    CMB_smooth_2 = 1;
  } else if(strcmp(realcov_type, "kk_kk")==0) {
    func_for_cov_G  = &func_for_cov_G_kk;
    func_bin_cov_NG = &bin_cov_NG_kk_kk;
    func_P1 = &Pl_tab;
    func_P2 = &Pl_tab;
    CMB_smooth_2 = 0;
    CMB_smooth_2 = 0;
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
  double beam1, beam2;

  for (l1 = 0; l1 < LMAX; l1++){
    l1_double = (double)l1;
    beam1 = pow(
      GaussianBeam(cmb.fwhm, l1, like.lmin_kappacmb, like.lmax_kappacmb), 
      CMB_smooth_1);

    // Gaussian Covariance Matrix
    // printf("l1,%d\n", l1);
    cov_g_l = func_for_cov_G(l1_double, z_ar);
    beam2 = pow(
      GaussianBeam(cmb.fwhm, l1, like.lmin_kappacmb, like.lmax_kappacmb),
      CMB_smooth_2);
    for(i=0; i<like.Ntheta; i++){
      covGl1 = cov_g_l * func_P1(i,l1);
      // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
      for(j=0; j<like.Ntheta ; j++){
        cov[i][j] += covGl1 * func_P2(j,l1) * beam1 * beam2;
      }
    }

    // Non-Gaussian Covariance Matrix
    if(FLAG_NG){
      for (l2 = 0; l2 < LMAX; l2++){
        tri = func_bin_cov_NG(l1_double,(double)l2,z_ar);
        beam2 = pow( 
          GaussianBeam(cmb.fwhm, l2, like.lmin_kappacmb,like.lmax_kappacmb),
          CMB_smooth_2);

        for(i=0; i<like.Ntheta ; i++){
          triP = tri * func_P1(i,l1);
          // printf("Glplus[%d][%d],%lg\n", i,l1, Glplus[i][l1]);
          for(j=0; j<like.Ntheta ; j++){
            covNG[i][j] += triP * func_P2(j,l2) * beam1 * beam2;
          }
        }
      }
    }
  }
  for(i=0; i<like.Ntheta; i++){
    for(j=0; j<like.Ntheta; j++){
      if(N[i][j]){
        cov[i][j] += N[i][j]/sqrt(w_mask(theta[i], Icol)*w_mask(theta[j], Icol));
      }
    }
  }
  free_double_matrix(N,0, like.Ntheta-1, 0, like.Ntheta-1);

}

// Fourier-space template function
// ell in Logarithmic bins
// Does not include CMB beam
void cov_fourier_binned(double **cov, double **covNG, char *cov_type, 
  int *z_ar, int FLAG_NG, double *ell)
{

  int i,j;
  static int LMAX = 50000;
  // Did not implement Gaussian smoothing since this one is a theoretical study.

  // double N[like.Ntheta];
  // for(i=0;i<like.Ntheta;i++) {N[i] = 0.;}

  double (*func_for_cov_G)(double, int*);
  double (*func_bin_cov_NG)(double, double, int*);


  if(strcmp(cov_type, "ss_ss")==0) {
    func_for_cov_G  = &func_for_cov_G_shear;
    func_bin_cov_NG = &bin_cov_NG_shear_shear;
  } else if(strcmp(cov_type, "gl_ss")==0) {
    func_for_cov_G  = &func_for_cov_G_gl_shear;
    func_bin_cov_NG = &bin_cov_NG_gl_shear;
  } else if(strcmp(cov_type, "cl_ss")==0) {
    func_for_cov_G  = &func_for_cov_G_cl_shear;
    func_bin_cov_NG = &bin_cov_NG_cl_shear;
  } else if(strcmp(cov_type, "gl_gl")==0) {
    func_for_cov_G  = &func_for_cov_G_gl;
    func_bin_cov_NG = &bin_cov_NG_gl_gl;
  } else if(strcmp(cov_type, "cl_gl")==0) {
    func_for_cov_G  = &func_for_cov_G_cl_gl;
    func_bin_cov_NG = &bin_cov_NG_cl_gl;
  } else if(strcmp(cov_type, "cl_cl")==0) {
    func_for_cov_G  = &func_for_cov_G_cl;
    func_bin_cov_NG = &bin_cov_NG_cl_cl;
  } else if(strcmp(cov_type, "gk_ss")==0) {
    func_for_cov_G  = &func_for_cov_G_gk_shear;
    func_bin_cov_NG = &bin_cov_NG_gk_shear;
  } else if(strcmp(cov_type, "ks_ss")==0) {
    func_for_cov_G  = &func_for_cov_G_ks_shear;
    func_bin_cov_NG = &bin_cov_NG_ks_shear;
  } else if(strcmp(cov_type, "gk_gl")==0) {
    func_for_cov_G  = &func_for_cov_G_gk_gl;
    func_bin_cov_NG = &bin_cov_NG_gk_gl;
  } else if(strcmp(cov_type, "ks_gl")==0) {
    func_for_cov_G  = &func_for_cov_G_ks_gl;
    func_bin_cov_NG = &bin_cov_NG_ks_gl;
  } else if(strcmp(cov_type, "gk_cl")==0) {
    func_for_cov_G  = &func_for_cov_G_gk_cl;
    func_bin_cov_NG = &bin_cov_NG_gk_cl;
  } else if(strcmp(cov_type, "ks_cl")==0) {
    func_for_cov_G  = &func_for_cov_G_ks_cl;
    func_bin_cov_NG = &bin_cov_NG_ks_cl;
  } else if(strcmp(cov_type, "gk_gk")==0) {
    func_for_cov_G  = &func_for_cov_G_gk;
    func_bin_cov_NG = &bin_cov_NG_gk_gk;
  } else if(strcmp(cov_type, "ks_gk")==0) {
    func_for_cov_G  = &func_for_cov_G_ks_gk;
    func_bin_cov_NG = &bin_cov_NG_ks_gk;
  } else if(strcmp(cov_type, "ks_ks")==0) {
    func_for_cov_G  = &func_for_cov_G_ks;
    func_bin_cov_NG = &bin_cov_NG_ks_ks;

  } else if(strcmp(cov_type, "kk_ss")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_shear;
    func_bin_cov_NG = &bin_cov_NG_kk_shear;
  } else if(strcmp(cov_type, "kk_gl")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_gl;
    func_bin_cov_NG = &bin_cov_NG_kk_gl;
  } else if(strcmp(cov_type, "kk_cl")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_cl;
    func_bin_cov_NG = &bin_cov_NG_kk_cl;
  } else if(strcmp(cov_type, "kk_gk")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_gk;
    func_bin_cov_NG = &bin_cov_NG_kk_gk;
  } else if(strcmp(cov_type, "kk_ks")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_ks;
    func_bin_cov_NG = &bin_cov_NG_kk_ks;

  } else if(strcmp(cov_type, "kk_kk")==0) {
    func_for_cov_G  = &func_for_cov_G_kk;
    func_bin_cov_NG = &bin_cov_NG_kk_kk;

  } else {
    printf("cov_fourier_binned: cov_type \"%s\" not defined!\n", cov_type); exit(1);
  }

  double l1_double,tri;
  int l1,l2;
  for(i=0; i<like.Ncl ; i++){
    for(j=0; j<like.Ncl ; j++){
      cov[i][j] = 0.;
      covNG[i][j] = 0.;
    }
  }

  double triP;
  double covGl1;

  double cov_g_l[LMAX];
    // printf("lmin,lmax, %d, %d\n", (int)floor(like.lmin),(int)ceil(like.lmax));
  for(l1=(int)ceil(like.lmin); l1<(int)ceil(like.lmax)+1; l1++){
    cov_g_l[l1] = func_for_cov_G((double)l1, z_ar);
    // printf("cov_g_l: %d, %le\n", l1,cov_g_l[l1]);
  }

  int l1_min, l1_max, l2_min, l2_max;
  int N_l1, N_l2;
  for(i=0; i<like.Ncl; i++){
    l1_min = (int)ceil(ell[i]);
    l1_max = (int)ceil(ell[i+1])-1;
    N_l1 = l1_max - l1_min + 1;
    for(l1=l1_min; l1<=l1_max; l1++){
      cov[i][i] += cov_g_l[l1];

      if(FLAG_NG){
        l1_double = (double)l1;
        for(j=0; j<like.Ncl ; j++){
          l2_min = (int)ceil(ell[j]);
          l2_max = (int)ceil(ell[j+1])-1;
          N_l2 = l2_max - l2_min + 1;
          for (l2 = l2_min; l2 <= l2_max; l2++){
            tri = func_bin_cov_NG(l1_double,(double)l2,z_ar);
            covNG[i][j] += tri/(N_l2*N_l1);
            // printf("here!!!%d, %d, %le\n", z_ar[0], z_ar[1], bin_cov_NG_shear_shear(3.366055e+01, 3.366055e+01,z_ar));
          }
          // covNG[i][j] /= (float)(N_l2);
        }
      }
    }
    cov[i][i] /= (float)(N_l1*N_l1);
  }
}

// ell in band power (5x2pt), ell in Logarithmic bins (kk)
// Does not include covariance matrix in 5x2pt, since they are Log-bin only
// Does not include CMB beam
void cov_fourier_banded(double **cov, double **covNG, char *cov_type, 
  int *z_ar, int FLAG_NG, double *ell, int **bindef)
{
  int i,j;
  static int LMAX = 50000;
  // Did not implement Gaussian smoothing since this one is a theoretical study.

  double (*func_for_cov_G)(double, int*);
  double (*func_bin_cov_NG)(double, double, int*);
  double (*func_P1)(int, int) = NULL;
  double (*func_P2)(int, int) = NULL;
  int N1 = like.Ncl;
  int N2 = like.Ncl;

  // kk x 5x2pt (5 functions, band power x Logarithmic bins)
  if(strcmp(cov_type, "kk_ss")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_shear;
    func_bin_cov_NG = &bin_cov_NG_kk_shear;
    func_P1 = &CMB_lensing_mat; N1 = like.Nbp;
  } else if(strcmp(cov_type, "kk_gl")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_gl;
    func_bin_cov_NG = &bin_cov_NG_kk_gl;
    func_P1 = &CMB_lensing_mat; N1 = like.Nbp;
  } else if(strcmp(cov_type, "kk_cl")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_cl;
    func_bin_cov_NG = &bin_cov_NG_kk_cl;
    func_P1 = &CMB_lensing_mat; N1 = like.Nbp;
  } else if(strcmp(cov_type, "kk_gk")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_gk;
    func_bin_cov_NG = &bin_cov_NG_kk_gk;
    func_P1 = &CMB_lensing_mat; N1 = like.Nbp;
  } else if(strcmp(cov_type, "kk_ks")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_ks;
    func_bin_cov_NG = &bin_cov_NG_kk_ks;
    func_P1 = &CMB_lensing_mat; N1 = like.Nbp;
  // kk x kk (band power x band power)
  } else if(strcmp(cov_type, "kk_kk")==0) {
    func_for_cov_G  = &func_for_cov_G_kk;
    func_bin_cov_NG = &bin_cov_NG_kk_kk;
    func_P1 = &CMB_lensing_mat_with_corr; N1 = like.Nbp;
    func_P2 = &CMB_lensing_mat_with_corr; N2 = like.Nbp;

  } else {
    printf("cov_fourier_banded: cov_type \"%s\" not defined!\n", cov_type); 
    exit(1);
  }

  double tri, factor1, factor2;
  int l1,l2, l1_min, l1_max, l2_min, l2_max, N_l1, N_l2;
  for(i=0; i<N1 ; i++){
    for(j=0; j<N2 ; j++){
      cov[i][j] = 0.;
      covNG[i][j] = 0.;
    }
  }

  
  // Pre-compute the Gaussian covariance matrix
  double cov_g_l[LMAX];
  printf("Pre-compute the Gaussian cov (%d <= L <= %d)\n", 
    like.lmin_bp_with_corr, like.lmax_bp_with_corr);
  for(l1=like.lmin_bp_with_corr; l1<like.lmax_bp_with_corr+1; l1++){
    cov_g_l[l1] = func_for_cov_G((double)l1, z_ar);
    // printf("cov_g_l[L=%d] = %le\n", l1, cov_g_l[l1]);
  }


  // Loop 1: Harmonic space --- band power
  for(i=0; i<N1; i++){

    // band power
    if((strcmp(cov_type, "kk_kk")!=0)){
      l1_min = bindef[i][0];
      l1_max = bindef[i][1];
    } else{
      l1_min = like.lmin_bp_with_corr;
      l1_max = like.lmax_bp_with_corr;
    }
    N_l1 = l1_max - l1_min + 1;

    // Loop 2: Harmonic space
    for(j=0; j<N2; j++){
      printf("cov index [%d, %d]\n", i+1, j+1);
      if(func_P2!=NULL){ // band power
        if((strcmp(cov_type, "kk_kk")!=0)){
          l2_min = bindef[j][0];
          l2_max = bindef[j][1];
        } else{
          l2_min = like.lmin_bp_with_corr;
          l2_max = like.lmax_bp_with_corr;
        }
      }
      else{ // Logarithmic bins
        l2_min = (int)ceil(ell[j]);
        l2_max = (int)ceil(ell[j+1])-1;
      }
      N_l2 = l2_max - l2_min + 1;

      for(l1=l1_min; l1<=l1_max; l1++){
        for(l2=l2_min; l2<=l2_max; l2++){
          
          factor1 = func_P1(i,l1);
          if(func_P2!=NULL){factor2 = func_P2(j,l2);}// band power
          else{factor2 = 1.0/N_l2;}// Log-bin average
          
          if(l1==l2){
            if(i==N1-1){
              printf("cov(bin-%d, bin-%d): L=%d f1=%e f2=%e cg=%e\n",i,j,l1,factor1,factor2, func_for_cov_G((double)l1, z_ar));
            }
            cov[i][j] += cov_g_l[l1] * factor1 * factor2;
            //cov[i][j] += func_for_cov_G((double)l1, z_ar) * factor1 * factor2;
          }

          if(FLAG_NG){
            tri = func_bin_cov_NG((double)l1, (double)l2, z_ar);
            covNG[i][j] += tri * factor1 * factor2;
          }
        }
      }
    }// End Loop 2
  }// End Loop 1
}

// Mixed-space template functions
// theta & ell both in Logarithmic bins
// Include CMB beam smoothing for real-space probes
void cov_mix_binned_fullsky(double **cov, double **covNG, char *mixcov_type, 
  int *z_ar, int FLAG_NG, double *theta, double *dtheta, double *ell)
{

  int i,j;
  static int LMAX = 50000;
  int CMB_smooth_1 = 0, CMB_smooth_2 = 0; // flag for CMB smoothing kernel

  // double N[like.Ntheta];
  // for(i=0;i<like.Ntheta;i++) {N[i] = 0.;}

  double (*func_for_cov_G)(double, int*);
  double (*func_bin_cov_NG)(double, double, int*);
  double (*func_P2)(int, int);

  if(strcmp(mixcov_type, "kk_xi+")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_shear;
    func_bin_cov_NG = &bin_cov_NG_kk_shear;
    func_P2 = &Glplus_tab;
    CMB_smooth_1 = 0;
    CMB_smooth_2 = 0;
  } else if(strcmp(mixcov_type, "kk_xi-")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_shear;
    func_bin_cov_NG = &bin_cov_NG_kk_shear;
    func_P2 = &Glminus_tab;
    CMB_smooth_1 = 0;
    CMB_smooth_2 = 0;
  } else if(strcmp(mixcov_type, "kk_gl")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_gl;
    func_bin_cov_NG = &bin_cov_NG_kk_gl;
    func_P2 = &Pl2_tab;
    CMB_smooth_1 = 0;
    CMB_smooth_2 = 0;
  } else if(strcmp(mixcov_type, "kk_cl")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_cl;
    func_bin_cov_NG = &bin_cov_NG_kk_cl;
    func_P2 = &Pl_tab;
    CMB_smooth_1 = 0;
    CMB_smooth_2 = 0;
  } else if(strcmp(mixcov_type, "kk_gk")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_gk;
    func_bin_cov_NG = &bin_cov_NG_kk_gk;
    func_P2 = &Pl_tab;
    CMB_smooth_1 = 0;
    CMB_smooth_2 = 1;
  } else if(strcmp(mixcov_type, "kk_ks")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_ks;
    func_bin_cov_NG = &bin_cov_NG_kk_ks;
    func_P2 = &Pl2_tab;
    CMB_smooth_1 = 0;
    CMB_smooth_2 = 1;
  } else {
    printf("cov_mix_binned_fullsky: mixcov_type \"%s\" not defined!\n", mixcov_type); exit(1);
  }

  double l1_double,tri;
  int l1,l2;
  double beam1, beam2;
  //double beam;
  for(i=0; i<like.Ncl ; i++){
    for(j=0; j<like.Ntheta ; j++){
      cov[i][j] = 0.;
      covNG[i][j] = 0.;
    }
  }

  double triP;
  double covGl1;

  double cov_g_l[LMAX];
  //printf("lmin,lmax, %d, %d\n", (int)floor(like.lmin),(int)ceil(like.lmax));
  for(l1=(int)ceil(like.lmin); l1<(int)ceil(like.lmax)+1; l1++){
    cov_g_l[l1] = func_for_cov_G((double)l1, z_ar);
    //printf("cov_g_l: %d, %le\n", l1,cov_g_l[l1]);
  }
  int l1_min, l1_max, l2_min, l2_max;
  int N_l1, N_l2;

  // Loop 1: ell log-bin
  for(i=0; i<like.Ncl; i++){
    l1_min = (int)ceil(ell[i]);
    l1_max = (int)ceil(ell[i+1])-1;
    N_l1 = l1_max - l1_min + 1;
    // Loop 2: real-space theta log-bin
    for(j=0; j<like.Ntheta ; j++){
      
      for(l1=l1_min; l1<=l1_max; l1++){
        
        // Gaussian Covariance Matrix
        beam1 = pow(
          GaussianBeam(cmb.fwhm, l1, like.lmin_kappacmb, like.lmax_kappacmb), 
          CMB_smooth_1 + CMB_smooth_2);
        cov[i][j] += cov_g_l[l1] * func_P2(j,l1) / N_l1 * beam1;
        
        // Non-Gaussian Covariance Matrix
        if(FLAG_NG){
          l1_double = (double)l1;
          beam1 = pow(
            GaussianBeam(cmb.fwhm,l1,like.lmin_kappacmb,like.lmax_kappacmb),
            CMB_smooth_1);
          for (l2 = 0; l2 < LMAX; l2++){
            beam2 = pow(
              GaussianBeam(cmb.fwhm,l2,like.lmin_kappacmb,like.lmax_kappacmb),
              CMB_smooth_2);
  
            tri = func_bin_cov_NG(l1_double,(double)l2,z_ar);
            covNG[i][j] += tri * func_P2(j,l2) / N_l1 * beam1 * beam2;
          }
        }
      }
    }
  }
}

// ell in band power, theta in Logarithmic bins
// Include CMB beam smoothing for real-space probes
void cov_mix_banded_fullsky(double **cov, double **covNG, char *mixcov_type, 
  int *z_ar, int FLAG_NG, double *theta, double *dtheta, int **bindef)
{
  int i,j;
  static int LMAX = 50000;
  int CMB_smooth_1 = 0, CMB_smooth_2 = 0; // flag for CMB smoothing kernel

  double (*func_for_cov_G)(double, int*);
  double (*func_bin_cov_NG)(double, double, int*);
  double (*func_P1)(int, int);
  double (*func_P2)(int, int);

  if(strcmp(mixcov_type, "kk_xi+")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_shear;
    func_bin_cov_NG = &bin_cov_NG_kk_shear;
    func_P1 = &CMB_lensing_mat;
    func_P2 = &Glplus_tab;
    CMB_smooth_1 = 0;
    CMB_smooth_2 = 0;
  } else if(strcmp(mixcov_type, "kk_xi-")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_shear;
    func_bin_cov_NG = &bin_cov_NG_kk_shear;
    func_P1 = &CMB_lensing_mat;
    func_P2 = &Glminus_tab;
    CMB_smooth_1 = 0;
    CMB_smooth_2 = 0;
  } else if(strcmp(mixcov_type, "kk_gl")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_gl;
    func_bin_cov_NG = &bin_cov_NG_kk_gl;
    func_P1 = &CMB_lensing_mat;
    func_P2 = &Pl2_tab;
    CMB_smooth_1 = 0;
    CMB_smooth_2 = 0;
  } else if(strcmp(mixcov_type, "kk_cl")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_cl;
    func_bin_cov_NG = &bin_cov_NG_kk_cl;
    func_P1 = &CMB_lensing_mat;
    func_P2 = &Pl_tab;
    CMB_smooth_1 = 0;
    CMB_smooth_2 = 0;
  } else if(strcmp(mixcov_type, "kk_gk")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_gk;
    func_bin_cov_NG = &bin_cov_NG_kk_gk;
    func_P1 = &CMB_lensing_mat;
    func_P2 = &Pl_tab;
    CMB_smooth_1 = 0;
    CMB_smooth_2 = 1;
  } else if(strcmp(mixcov_type, "kk_ks")==0) {
    func_for_cov_G  = &func_for_cov_G_kk_ks;
    func_bin_cov_NG = &bin_cov_NG_kk_ks;
    func_P1 = &CMB_lensing_mat;
    func_P2 = &Pl2_tab;
    CMB_smooth_1 = 0;
    CMB_smooth_2 = 1;
  } else {
    printf("cov_mix_banded_fullsky: mixcov_type \"%s\" not defined!\n", mixcov_type); exit(1);
  }

  double l1_double,tri;
  int l1,l2;
  double beam1, beam2;
  //double beam;
  for(i=0; i<like.Nbp ; i++){
    for(j=0; j<like.Ntheta ; j++){
      cov[i][j] = 0.;
      covNG[i][j] = 0.;
    }
  }

  double triP;
  double covGl1;

  double cov_g_l[LMAX];
  //printf("lmin,lmax, %d, %d\n", like.lmin_bp, like.lmax_bp);
  for(l1=like.lmin_bp; l1<=like.lmax_bp; l1++){
    cov_g_l[l1] = func_for_cov_G((double)l1, z_ar);
    //printf("cov_g_l: %d, %le\n", l1,cov_g_l[l1]);
  }
  int l1_min, l1_max;
  int N_l1;
  
  // Loop 1: Harmonic space
  for(i=0; i<like.Nbp; i++){
    l1_min = bindef[i][0];
    l1_max = bindef[i][1];
    N_l1 = l1_max - l1_min + 1;

    // Loop 2: Configuration space
    for(j=0; j<like.Ntheta ; j++){
      
      // sum over ell
      for(l1=l1_min; l1<=l1_max; l1++){
        
        // Gaussian
        beam1 = pow( 
          GaussianBeam(cmb.fwhm, l1, like.lmin_kappacmb, like.lmax_kappacmb),
          CMB_smooth_1 + CMB_smooth_2);

        cov[i][j] += cov_g_l[l1]*func_P1(i,l1)*func_P2(j,l1)*beam1;
        
        // Non-Gaussian
        if(FLAG_NG){
          l1_double = (double)l1;
          beam1 = pow( 
            GaussianBeam(cmb.fwhm, l1, like.lmin_kappacmb, like.lmax_kappacmb),
            CMB_smooth_1);
          for (l2 = 0; l2 < LMAX; l2++){
            beam2 = pow(
              GaussianBeam(cmb.fwhm,l2,like.lmin_kappacmb,like.lmax_kappacmb), 
              CMB_smooth_2);
            tri = func_bin_cov_NG(l1_double,(double)l2,z_ar);
            
            covNG[i][j] += tri*func_P1(i,l1)*func_P2(j,l2)*beam1*beam2;
          }
        }
      }// end sum over ell
    }// end Loop 2
  }// End Loop 1
}

/****** Full-sky Covariance API ******/

// real space
void cov_shear_shear_real_binned_fullsky(double **cov, double **covNG, 
  int z1,int z2,int z3,int z4,int pm1,int pm2, int FLAG_NG, 
  double *theta, double *dtheta)
{
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

void cov_gl_shear_real_binned_fullsky(double **cov, double **covNG, 
  int z1,int z2,int z3,int z4,int pm, int FLAG_NG, 
  double *theta, double *dtheta)
{
  static char xip[]="gl_xi+", xim[]="gl_xi-";
  char realcov_type[8];
  if(pm==1) {strcpy(realcov_type, xip);}
  else {strcpy(realcov_type, xim);}
  int z_ar[4];
  z_ar[0]=z1; z_ar[1]=z2; z_ar[2]=z3; z_ar[3]=z4;
  cov_real_binned_fullsky(cov, covNG, realcov_type, z_ar, FLAG_NG, theta, dtheta);
}

void cov_cl_shear_real_binned_fullsky(double **cov, double **covNG, 
  int z1,int z2,int z3,int z4,int pm, int FLAG_NG, 
  double *theta, double *dtheta)
{
  static char xip[]="cl_xi+", xim[]="cl_xi-";
  char realcov_type[8];
  if(pm==1) {strcpy(realcov_type, xip);}
  else {strcpy(realcov_type, xim);}
  int z_ar[4];
  z_ar[0]=z1; z_ar[1]=z2; z_ar[2]=z3; z_ar[3]=z4;
  cov_real_binned_fullsky(cov, covNG, realcov_type, z_ar, FLAG_NG, theta, dtheta);
}

void cov_cl_gl_real_binned_fullsky(double **cov, double **covNG, 
  int z1,int z2,int z3,int z4, int FLAG_NG, 
  double *theta, double *dtheta)
{
  char realcov_type[] = "cl_gl";
  int z_ar[4];
  z_ar[0]=z1; z_ar[1]=z2; z_ar[2]=z3; z_ar[3]=z4;
  cov_real_binned_fullsky(cov, covNG, realcov_type, z_ar, FLAG_NG, theta, dtheta);
}

void cov_gl_gl_real_binned_fullsky(double **cov, double **covNG, 
  int z1,int z2,int z3,int z4, int FLAG_NG, 
  double *theta, double *dtheta)
{
  char realcov_type[] = "gl_gl";
  int z_ar[4];
  z_ar[0]=z1; z_ar[1]=z2; z_ar[2]=z3; z_ar[3]=z4;
  cov_real_binned_fullsky(cov, covNG, realcov_type, z_ar, FLAG_NG, theta, dtheta);
}

void cov_cl_cl_real_binned_fullsky(double **cov, double **covNG, 
  int z1,int z2,int z3,int z4, int FLAG_NG, 
  double *theta, double *dtheta)
{
  char realcov_type[] = "cl_cl";
  int z_ar[4];
  z_ar[0]=z1; z_ar[1]=z2; z_ar[2]=z3; z_ar[3]=z4;
  cov_real_binned_fullsky(cov, covNG, realcov_type, z_ar, FLAG_NG, theta, dtheta);
}
//
void cov_ks_shear_real_binned_fullsky(double **cov, double **covNG, 
  int z1,int z3,int z4,int pm, int FLAG_NG, 
  double *theta, double *dtheta)
{
  static char xip[]="ks_xi+", xim[]="ks_xi-";
  char realcov_type[8];
  if(pm==1) {strcpy(realcov_type, xip);}
  else {strcpy(realcov_type, xim);}
  int z_ar[3];
  z_ar[0]=z1; z_ar[1]=z3; z_ar[2]=z4;
  cov_real_binned_fullsky(cov, covNG, realcov_type, z_ar, FLAG_NG, theta, dtheta);
}

void cov_gk_shear_real_binned_fullsky(double **cov, double **covNG, 
  int z1,int z3,int z4,int pm, int FLAG_NG, 
  double *theta, double *dtheta)
{
  static char xip[]="gk_xi+", xim[]="gk_xi-";
  char realcov_type[8];
  if(pm==1) {strcpy(realcov_type, xip);}
  else {strcpy(realcov_type, xim);}
  int z_ar[3];
  z_ar[0]=z1; z_ar[1]=z3; z_ar[2]=z4;
  cov_real_binned_fullsky(cov, covNG, realcov_type, z_ar, FLAG_NG, theta, dtheta);
}

void cov_gk_gl_real_binned_fullsky(double **cov, double **covNG, 
  int z1,int z3,int z4, int FLAG_NG, 
  double *theta, double *dtheta)
{
  char realcov_type[] = "gk_gl";
  int z_ar[3];
  z_ar[0]=z1; z_ar[1]=z3; z_ar[2]=z4;
  cov_real_binned_fullsky(cov, covNG, realcov_type, z_ar, FLAG_NG, theta, dtheta);
}

void cov_ks_gl_real_binned_fullsky(double **cov, double **covNG, 
  int z1,int z3,int z4, int FLAG_NG, 
  double *theta, double *dtheta)
{
  char realcov_type[] = "ks_gl";
  int z_ar[3];
  z_ar[0]=z1; z_ar[1]=z3; z_ar[2]=z4;
  cov_real_binned_fullsky(cov, covNG, realcov_type, z_ar, FLAG_NG, theta, dtheta);
}

void cov_gk_cl_real_binned_fullsky(double **cov, double **covNG, 
  int z1,int z3,int z4, int FLAG_NG, 
  double *theta, double *dtheta)
{
  char realcov_type[] = "gk_cl";
  int z_ar[3];
  z_ar[0]=z1; z_ar[1]=z3; z_ar[2]=z4;
  cov_real_binned_fullsky(cov, covNG, realcov_type, z_ar, FLAG_NG, theta, dtheta);
}

void cov_ks_cl_real_binned_fullsky(double **cov, double **covNG, 
  int z1,int z3,int z4, int FLAG_NG, 
  double *theta, double *dtheta)
{
  char realcov_type[] = "ks_cl";
  int z_ar[3];
  z_ar[0]=z1; z_ar[1]=z3; z_ar[2]=z4;
  cov_real_binned_fullsky(cov, covNG, realcov_type, z_ar, FLAG_NG, theta, dtheta);
}

void cov_gk_gk_real_binned_fullsky(double **cov, double **covNG, 
  int z1,int z2, int FLAG_NG, 
  double *theta, double *dtheta)
{
  char realcov_type[] = "gk_gk";
  int z_ar[2];
  z_ar[0]=z1; z_ar[1]=z2;
  cov_real_binned_fullsky(cov, covNG, realcov_type, z_ar, FLAG_NG, theta, dtheta);
}

void cov_ks_gk_real_binned_fullsky(double **cov, double **covNG, 
  int z1,int z2, int FLAG_NG, 
  double *theta, double *dtheta)
{
  char realcov_type[] = "ks_gk";
  int z_ar[2];
  z_ar[0]=z1; z_ar[1]=z2;
  cov_real_binned_fullsky(cov, covNG, realcov_type, z_ar, FLAG_NG, theta, dtheta);
}

void cov_ks_ks_real_binned_fullsky(double **cov, double **covNG, 
  int z1,int z2, int FLAG_NG, 
  double *theta, double *dtheta)
{
  char realcov_type[] = "ks_ks";
  int z_ar[2];
  z_ar[0]=z1; z_ar[1]=z2;
  cov_real_binned_fullsky(cov, covNG, realcov_type, z_ar, FLAG_NG, theta, dtheta);
}


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

/* mixed space, theta & ell both in Logarithmic bins (5 blocks) */ 

void cov_kk_shear_mix_binned_fullsky(double **cov, double **covNG, 
  int z3,int z4,int pm,int FLAG_NG,double *theta, double *dtheta, double *ell)
{
  static char xip[]="kk_xi+", xim[]="kk_xi-";
  char mixcov_type[8];
  if(pm==1) {strcpy(mixcov_type, xip);}
  else {strcpy(mixcov_type, xim);}
  int z_ar[2];
  z_ar[0]=z3; z_ar[1]=z4;
  cov_mix_binned_fullsky(cov,covNG,mixcov_type,z_ar,FLAG_NG,theta,dtheta,ell);
}

void cov_kk_gl_mix_binned_fullsky(double **cov, double **covNG, 
  int z3,int z4,int FLAG_NG,double *theta, double *dtheta, double *ell)
{
  char mixcov_type[] = "kk_gl";
  int z_ar[2];
  z_ar[0]=z3; z_ar[1]=z4;
  cov_mix_binned_fullsky(cov,covNG,mixcov_type,z_ar,FLAG_NG,theta,dtheta,ell);
}

void cov_kk_cl_mix_binned_fullsky(double **cov, double **covNG, 
  int z3,int z4,int FLAG_NG,double *theta, double *dtheta, double *ell)
{
  char mixcov_type[] = "kk_cl";
  int z_ar[2];
  z_ar[0]=z3; z_ar[1]=z4;
  cov_mix_binned_fullsky(cov,covNG,mixcov_type,z_ar,FLAG_NG,theta,dtheta,ell);
}

void cov_kk_gk_mix_binned_fullsky(double **cov, double **covNG, 
  int zl, int FLAG_NG, double *theta, double *dtheta, double *ell)
{
  char mixcov_type[] = "kk_gk";
  int z_ar[1];
  z_ar[0]=zl;
  cov_mix_binned_fullsky(cov,covNG,mixcov_type,z_ar,FLAG_NG,theta,dtheta,ell);
}

void cov_kk_ks_mix_binned_fullsky(double **cov, double **covNG, 
  int zs, int FLAG_NG, double *theta, double *dtheta, double *ell)
{
  char mixcov_type[] = "kk_ks";
  int z_ar[1];
  z_ar[0]=zs;
  cov_mix_binned_fullsky(cov,covNG,mixcov_type,z_ar,FLAG_NG,theta,dtheta,ell);
}

/* mixed space, theta in Logarithmic bins, ell in band powers (5 blocks) */

// kk x {5x2pt} (5 blocks)
void cov_kk_shear_mix_banded_fullsky(double **cov, double **covNG, 
  int z3,int z4,int pm,int FLAG_NG,double *theta,double *dtheta,int **bindef)
{
  static char xip[]="kk_xi+", xim[]="kk_xi-";
  char mixcov_type[8];
  if(pm==1) {strcpy(mixcov_type, xip);}
  else {strcpy(mixcov_type, xim);}
  int z_ar[2];
  z_ar[0]=z3; z_ar[1]=z4;
  printf("mixcov type = %s\n",mixcov_type);
  cov_mix_banded_fullsky(cov, covNG, mixcov_type, z_ar, FLAG_NG, 
    theta, dtheta, bindef);
}

void cov_kk_gl_mix_banded_fullsky(double **cov, double **covNG, 
  int z3,int z4, int FLAG_NG, double *theta, double *dtheta, int **bindef)
{
  char mixcov_type[] = "kk_gl";
  int z_ar[2];
  z_ar[0]=z3; z_ar[1]=z4;
  cov_mix_banded_fullsky(cov, covNG, mixcov_type, z_ar, FLAG_NG, 
    theta, dtheta, bindef);
}

void cov_kk_cl_mix_banded_fullsky(double **cov, double **covNG, 
  int z3,int z4, int FLAG_NG, double *theta, double *dtheta, int **bindef)
{
  char mixcov_type[] = "kk_cl";
  int z_ar[2];
  z_ar[0]=z3; z_ar[1]=z4;
  cov_mix_banded_fullsky(cov, covNG, mixcov_type, z_ar, FLAG_NG, 
    theta, dtheta, bindef);
}

void cov_kk_gk_mix_banded_fullsky(double **cov, double **covNG, 
  int zl, int FLAG_NG, double *theta, double *dtheta, int **bindef)
{
  char mixcov_type[] = "kk_gk";
  int z_ar[1];
  z_ar[0]=zl;
  cov_mix_banded_fullsky(cov, covNG, mixcov_type, z_ar, FLAG_NG, 
    theta, dtheta, bindef);
}

void cov_kk_ks_mix_banded_fullsky(double **cov, double **covNG, 
  int zs, int FLAG_NG, double *theta, double *dtheta, int **bindef)
{
  char mixcov_type[] = "kk_ks";
  int z_ar[1];
  z_ar[0]=zs;
  cov_mix_banded_fullsky(cov, covNG, mixcov_type, z_ar, FLAG_NG, 
    theta, dtheta, bindef);
}

/* Fourier space, ell in Logarithmic bins (21 blocks) */

void cov_shear_shear_fourier_binned(double **cov, double **covNG, 
  int z1,int z2,int z3,int z4, int FLAG_NG, double *ell)
{
  char cov_type[]= "ss_ss";
  int z_ar[4];
  z_ar[0]=z1; z_ar[1]=z2; z_ar[2]=z3; z_ar[3]=z4;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_gl_shear_fourier_binned(double **cov, double **covNG, 
  int z1,int z2,int z3,int z4, int FLAG_NG, double *ell)
{
  char cov_type[]= "gl_ss";
  int z_ar[4];
  z_ar[0]=z1; z_ar[1]=z2; z_ar[2]=z3; z_ar[3]=z4;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_cl_shear_fourier_binned(double **cov, double **covNG, 
  int z1,int z2,int z3,int z4, int FLAG_NG, double *ell)
{
  char cov_type[]= "cl_ss";
  int z_ar[4];
  z_ar[0]=z1; z_ar[1]=z2; z_ar[2]=z3; z_ar[3]=z4;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_cl_gl_fourier_binned(double **cov, double **covNG, 
  int z1,int z2,int z3,int z4, int FLAG_NG, double *ell)
{
  char cov_type[] = "cl_gl";
  int z_ar[4];
  z_ar[0]=z1; z_ar[1]=z2; z_ar[2]=z3; z_ar[3]=z4;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_gl_gl_fourier_binned(double **cov, double **covNG, 
  int z1,int z2,int z3,int z4, int FLAG_NG, double *ell)
{
  char cov_type[] = "gl_gl";
  int z_ar[4];
  z_ar[0]=z1; z_ar[1]=z2; z_ar[2]=z3; z_ar[3]=z4;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_cl_cl_fourier_binned(double **cov, double **covNG, 
  int z1,int z2,int z3,int z4, int FLAG_NG, double *ell)
{
  char cov_type[] = "cl_cl";
  int z_ar[4];
  z_ar[0]=z1; z_ar[1]=z2; z_ar[2]=z3; z_ar[3]=z4;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

// 5x2pt (9 blocks)
void cov_ks_shear_fourier_binned(double **cov, double **covNG, 
  int z1,int z3,int z4, int FLAG_NG, double *ell)
{
  char cov_type[] = "ks_ss";
  int z_ar[3];
  z_ar[0]=z1; z_ar[1]=z3; z_ar[2]=z4;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_gk_shear_fourier_binned(double **cov, double **covNG, 
  int z1,int z3,int z4, int FLAG_NG, double *ell)
{
  char cov_type[] = "gk_ss";
  int z_ar[3];
  z_ar[0]=z1; z_ar[1]=z3; z_ar[2]=z4;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_gk_gl_fourier_binned(double **cov, double **covNG, 
  int z1,int z3,int z4, int FLAG_NG, double *ell)
{
  char cov_type[] = "gk_gl";
  int z_ar[3];
  z_ar[0]=z1; z_ar[1]=z3; z_ar[2]=z4;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_ks_gl_fourier_binned(double **cov, double **covNG, 
  int z1,int z3,int z4, int FLAG_NG, double *ell)
{
  char cov_type[] = "ks_gl";
  int z_ar[3];
  z_ar[0]=z1; z_ar[1]=z3; z_ar[2]=z4;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_gk_cl_fourier_binned(double **cov, double **covNG, 
  int z1,int z3,int z4, int FLAG_NG, double *ell)
{
  char cov_type[] = "gk_cl";
  int z_ar[3];
  z_ar[0]=z1; z_ar[1]=z3; z_ar[2]=z4;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_ks_cl_fourier_binned(double **cov, double **covNG, 
  int z1,int z3,int z4, int FLAG_NG, double *ell)
{
  char cov_type[] = "ks_cl";
  int z_ar[3];
  z_ar[0]=z1; z_ar[1]=z3; z_ar[2]=z4;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_gk_gk_fourier_binned(double **cov, double **covNG, 
  int z1,int z2, int FLAG_NG, double *ell)
{
  char cov_type[] = "gk_gk";
  int z_ar[2];
  z_ar[0]=z1; z_ar[1]=z2;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_ks_gk_fourier_binned(double **cov, double **covNG, 
  int z1,int z2, int FLAG_NG, double *ell)
{
  char cov_type[] = "ks_gk";
  int z_ar[2];
  z_ar[0]=z1; z_ar[1]=z2;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_ks_ks_fourier_binned(double **cov, double **covNG, 
  int z1,int z2, int FLAG_NG, double *ell)
{
  char cov_type[] = "ks_ks";
  int z_ar[2];
  z_ar[0]=z1; z_ar[1]=z2;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

// kk x {5x2pt + kk} (6 blocks)
void cov_kk_shear_fourier_binned(double **cov, double **covNG, 
  int z3,int z4, int FLAG_NG, double *ell)
{
  char cov_type[] = "kk_ss";
  int z_ar[2];
  z_ar[0]=z3; z_ar[1]=z4;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_kk_gl_fourier_binned(double **cov, double **covNG, 
  int z3,int z4, int FLAG_NG, double *ell)
{
  char cov_type[] = "kk_gl";
  int z_ar[2];
  z_ar[0]=z3; z_ar[1]=z4;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_kk_cl_fourier_binned(double **cov, double **covNG, 
  int z3,int z4, int FLAG_NG, double *ell)
{
  char cov_type[] = "kk_cl";
  int z_ar[2];
  z_ar[0]=z3; z_ar[1]=z4;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_kk_gk_fourier_binned(double **cov, double **covNG, 
  int zl, int FLAG_NG, double *ell)
{
  char cov_type[] = "kk_gk";
  int z_ar[1];
  z_ar[0]=zl;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_kk_ks_fourier_binned(double **cov, double **covNG, 
  int zs, int FLAG_NG, double *ell)
{
  char cov_type[] = "kk_ks";
  int z_ar[1];
  z_ar[0]=zs;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

void cov_kk_kk_fourier_binned(double **cov, double **covNG, 
  int FLAG_NG, double *ell)
{
  char cov_type[] = "kk_kk";
  int z_ar[1];
  z_ar[0]=-1;
  cov_fourier_binned(cov, covNG, cov_type, z_ar, FLAG_NG, ell);
}

/* Fourier space, ell in band powers (1 blocks)*/

void cov_kk_kk_fourier_banded(double **cov, double **covNG, 
  int FLAG_NG, double *ell, int **bindef)
{
  char cov_type[] = "kk_kk";
  int z_ar[1];
  z_ar[0]=-1;
  cov_fourier_banded(cov, covNG, cov_type, z_ar, FLAG_NG, ell, bindef);
}

/****** Helper Functions ******/

// survey mask correlation function, normalized to 1 at small angular scales
double w_mask(double theta_min, int col)
{
  static int NTHETA = 0;
  static double **w_vec =0;
  static int Ncols = 0;
  int i,l;
  if (like.theta ==NULL || like.Ntheta < 1){
    printf("covariances_real_binned.c:w_mask: like.theta or like.Ntheta not initialized\nEXIT\n");
    exit(1);
  }
  if (w_vec ==0){
    NTHETA = like.Ntheta;
    FILE *F1;
    F1 = fopen(covparams.C_FOOTPRINT_FILE,"r");
    if (F1 != NULL) {
      // covparams.C_FOOTPRINT_FILE exists
      // # l [int], Cls [double *]
      // Then w_mask(theta) = Cl[Icol] * (2l+1)/4pi * P_l(cos(theta)) 
      // Cls has multiple columns showing different auto- and cross-
	  // power spectrum among several survey footprints. Users have
      // to manually set which column to use during calculation.
      fclose(F1);
      int lbins = line_count(covparams.C_FOOTPRINT_FILE);
      Ncols = column_count(covparams.C_FOOTPRINT_FILE);
      printf("MANUAL WARNING: the C_FOOTPRINT_FILE has multiple Cls columns, users have to manually specify which column to use in the source code!\n");
      printf("Reading file: %d columns, %d lines.\n", Ncols, lbins);
      double **Cl;
      w_vec = create_double_matrix(0, Ncols-2, 0, like.Ntheta-1);
      Cl = create_double_matrix(0, Ncols-2, 0, lbins-1);
      F1=fopen(covparams.C_FOOTPRINT_FILE,"r");
      for (int i = 0; i < lbins; i++){
        int tmp;
        fscanf(F1, "%d", &tmp);
        double tmp2;
        for(int j=0; j<Ncols-1; j++){
          fscanf(F1, "%le", &tmp2);
          Cl[j][i] = tmp2;
        }
      }
      fclose(F1);

      printf("\nTabulating w_mask(theta) from mask power spectrum %s\n",covparams.C_FOOTPRINT_FILE);
      for (i = 0; i < NTHETA; i++){
        printf("w_mask[%d]:", i);
        for (int j=0; j<Ncols-1; j++){
          w_vec[j][i] =0.;
          for (l = 0; l < lbins; l++){
            w_vec[j][i]+=Cl[j][l]*(2.*l+1)/(4.*M_PI)*gsl_sf_legendre_Pl(l,cos(like.theta[i]));
          }
          printf(" %le", w_vec[j][i]);
        }
        printf("\n");
      }
      free_double_matrix(Cl,0,Ncols-1,0,lbins-1);
    }
    else{
      //covparams.C_FOOTPRINT_FILE does not exit, ignore boundary effects
      Ncols = 2;
      w_vec = create_double_matrix(0,Ncols-2,0,like.Ntheta-1);
      printf("covparams.C_FOOTPRINT_FILE = %s not found\nNo boundary effect correction applied\n",covparams.C_FOOTPRINT_FILE);
      for (i = 0; i<NTHETA; i ++){
        w_vec[0][i] = 1.0;
        printf("w_mask[%d]: %e\n",i, w_vec[0][i]);
      }      
    }
  }
  i = 0;
  while(like.theta[i]< theta_min){
    i ++;
  }
  if(col < Ncols - 1){return w_vec[col][i];}
  else{
    printf("covariances_real_binned.c:w_mask: col %d/%d!\n",col,Ncols-1);
    exit(1);
  }
}

double w_pixel(int ell)
{
  static int LMAX = 50000;
  static double *cl_pixel =0;
  int i,l;

  if (cl_pixel == 0){
    cl_pixel = create_double_vector(0, LMAX -1);
    for (i = 0; i<LMAX; i ++){cl_pixel[i] = 1.0;}  
    FILE *F1;
    printf("\nTabulating Cl_pixwin(l) from pixel window function %s\n",\
        covparams.HEALPIX_WINDOW_FUNCTION_FILE);
    F1 = fopen(covparams.HEALPIX_WINDOW_FUNCTION_FILE, "r");
    if (F1 != NULL) {
      // covparams.HEALPIX_WINDOW_FUNCTION_FILE exists
      // # l [int], Cl [double]
      fclose(F1);
      int lbins = line_count(covparams.HEALPIX_WINDOW_FUNCTION_FILE);
      if(lbins <= LMAX){
        F1=fopen(covparams.HEALPIX_WINDOW_FUNCTION_FILE,"r");
        for (int i = 0; i < lbins; i++){
          int tmp;
          double tmp2;
          fscanf(F1,"%d %le\n",&tmp, &tmp2);
          cl_pixel[i] = tmp2;
        }
        fclose(F1);
        for(int i=lbins; i<LMAX; i++){cl_pixel[i] = 0.0;}
      }
      else{
        printf("LMAX = %d is too small since the file has %d lines\n",
          LMAX, lbins);
        exit(-1);
      }
    }
    else{
      //covparams.HEALPIX_WINDOW_FUNCTION_FILE does not exit, 
      // ignore healpix pixel window function
      printf("covparams.HEALPIX_WINDOW_FUNCTION_FILE = %s not found\n",
        covparams.HEALPIX_WINDOW_FUNCTION_FILE);
      printf("No healpix pixel window function correction applied\n");    
      for(int i=0; i<LMAX; i++){cl_pixel[i]=1.0;}
    }
  }
  if(ell<LMAX){return cl_pixel[ell];}
  else{
    printf("ell = %d larger than LMAX = %d!\n", ell, LMAX);
    exit(-1);
  }
}

double Glplus_tab(int itheta, int ell)
{  
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

double Glminus_tab(int itheta, int ell)
{  
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

double Pl_tab(int itheta, int ell)
{  
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

double Pl2_tab(int itheta, int ell)
{  
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

double CMB_lensing_mat(int ibp, int L)
{
  int iL;
  int NL = like.lmax_bp - like.lmin_bp + 1;
  static double **binning_matrix = NULL;
  if (binning_matrix == NULL){
    binning_matrix = create_double_matrix(0, like.Nbp-1, 
      0, NL-1);
    FILE *F1;
    F1 = fopen(covparams.BINMAT_FILE,"r");
    if (F1 != NULL) {
      fclose(F1);
      int nbp = line_count(covparams.BINMAT_FILE);
      
      F1=fopen(covparams.BINMAT_FILE, "r");
      for (int i = 0; i < like.Nbp; i++){
        for (int j = 0; j < NL; j++){
          double tmp2;
          fscanf(F1,"%le ",&tmp2);
          binning_matrix[i][j] = tmp2;
        }
      }
      fclose(F1);
    }
    else{
      printf("ERROR: Cant't open file %s\n", covparams.BINMAT_FILE);
      exit(-1);
    }  
  }
  iL = L - like.lmin_bp;
  return binning_matrix[ibp][iL];
}

double CMB_lensing_mat_with_corr(int ibp, int L)
{
  int iL;
  int NL = like.lmax_bp_with_corr - like.lmin_bp_with_corr + 1;
  static double **binning_correction_matrix = NULL;
  if (binning_correction_matrix == NULL){
    binning_correction_matrix = create_double_matrix(0, like.Nbp-1, 
      0, NL-1);
    FILE *F1;
    F1 = fopen(covparams.BINMAT_WITH_CORR_FILE,"r");
    if (F1 != NULL) {
      fclose(F1);
      int nbp = line_count(covparams.BINMAT_WITH_CORR_FILE);
      
      F1=fopen(covparams.BINMAT_WITH_CORR_FILE, "r");
      for (int i = 0; i < like.Nbp; i++){
        for (int j = 0; j < NL; j++){
          double tmp2;
          fscanf(F1,"%le ",&tmp2);
          binning_correction_matrix[i][j] = tmp2;
        }
      }
      fclose(F1);
    }
    else{
      printf("ERROR: Cant't open file %s\n", covparams.BINMAT_WITH_CORR_FILE);
      exit(-1);
    }  
  }
  iL = L - like.lmin_bp_with_corr;
  return binning_correction_matrix[ibp][iL];
}

/*  Calculate Gaussian Beam Kernel (including healpix pixwin)
    
    F(ell) = B(ell)H(ell - ell_min)H(ell_max - ell)
    where B(ell) is the Gaussian kernel
    B(ell) = exp( -ell*(ell+1)/ell_beam^2 )
    and ell_beam = sqrt(16*ln(2)) / theta_FWHM
    H(ell) is the Heaviside step function

    Inputs:
    -------
    theta_fwhm: FWHM of the CMB Gaussian beam, in units of rad
    ell:        angular mode ell where the beam kernel is evaluated
    ell_min:    large-scale cut
    ell_max:    small-scale cut

    Outputs:
    --------
    BeamKernel: The Gaussian kernel
*/
double GaussianBeam(double theta_fwhm, int ell, 
  double ell_left, double ell_right)
{
  double ell_min = ell_left - 10;
  double ell_max = ell_right + 50;
  double ell_beam = sqrt(16.0 * log(2.0)) / theta_fwhm;
  double BeamKernel, window_func, u=0;
  // Turn-On Gaussian Beam Smoothing
  if(theta_fwhm > 0.0){
  // smooth kernel
    BeamKernel = exp(-1.*ell*(ell+1.) / (ell_beam * ell_beam));
    // window function
    if((ell<ell_min) || (ell>ell_max)){
      window_func = 0.;
    }
    // tapering the edges
    else if((ell_min <= ell) && (ell<ell_left)){
      u = (ell - ell_min)/(ell_left - ell_min);
      window_func = u - sin( constants.twopi * u) / constants.twopi;
    }
    else if ((ell_right < ell) && (ell <= ell_max)){
      u = (ell_max - ell)/(ell_max - ell_right);
      window_func = u - sin( constants.twopi * u) / constants.twopi; 
    }
    else{
      window_func = 1.0;
    }
    BeamKernel *= window_func;
  }
  // Turn-Off Gaussian Beam Smoothing
  else{BeamKernel = 1.0;}
  return BeamKernel*w_pixel(ell);
}

/************ (21+6) Functions for different covariances ************/

// Pure noise term in Gaussian cov
void pure_noise_xipm_xipm(int *z_ar, double *theta, double *dtheta, double **N) {
  int z1,z2,z3,z4;
  int i;
  z1=z_ar[0]; z2=z_ar[1]; z3=z_ar[2]; z4=z_ar[3];
  if (z1 ==z3 && z2 ==z4){
    for(i=0;i<like.Ntheta;i++) {
      N[i][i] = 2.* pow(survey.sigma_e,4.0)/(M_PI*(2.*theta[i]+dtheta[i])*dtheta[i]*4.*nsource(z1)*nsource(z2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor);
    } //(number of galaxy pairs in the survey contributing to annulus of width Dtheta centered at theta1)^-1
  }
  if (z1 ==z4 && z2 ==z3){
    for(i=0;i<like.Ntheta;i++) {
      N[i][i] += 2.* pow(survey.sigma_e,4.0)/(M_PI*(2.*theta[i]+dtheta[i])*dtheta[i]*4.*nsource(z1)*nsource(z2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor);
    }
  } // factor of 2 already included
}
void pure_noise_gl_gl(int *z_ar, double *theta, double *dtheta, double **N) {
  int zl1, zs1, zl2, zs2;
  int i;
  zl1 = z_ar[0]; zs1 = z_ar[1]; zl2 = z_ar[2]; zs2 = z_ar[3];
  if (zl1 ==zl2 && zs1 ==zs2){
    for(i=0;i<like.Ntheta;i++) {
      N[i][i] = pow(survey.sigma_e,2.0)/(2.0*M_PI*(2.*theta[i]+dtheta[i])*dtheta[i]*nlens(zl1)*nsource(zs2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor);
    }
  }
}
void pure_noise_cl_cl(int *z_ar, double *theta, double *dtheta, double **N) {
  int z1,z2,z3,z4;
  int i;
  z1=z_ar[0]; z2=z_ar[1]; z3=z_ar[2]; z4=z_ar[3];
  if (z1 ==z3 && z2 ==z4){
    for(i=0;i<like.Ntheta;i++) {
      N[i][i] = 1./(M_PI*(2.*theta[i]+dtheta[i])*dtheta[i]*nlens(z1)*nlens(z2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor); //(number of galaxy pairs in the survey contributing to annulus of width Dtheta centered at theta1)^-1
    }
  }
  if (z1 ==z4 && z2 ==z3){
    for(i=0;i<like.Ntheta;i++) {
      N[i][i] += 1./(M_PI*(2.*theta[i]+dtheta[i])*dtheta[i]*nlens(z1)*nlens(z2)*pow(survey.n_gal_conversion_factor,2.0)*survey.area*survey.area_conversion_factor);
    }
  }
}
void pure_noise_gk_gk(int *z_ar, double *theta, double *dtheta, double **N){
  // N is 2d matrix
  double N13=0, N24=0, beam=0, P1=0, P2=0, ell_prefactor=0;
  double fsky = survey.area*survey.area_conversion_factor/(4.*M_PI);
  int n1,n3;
  static int LMAX = 50000;
  n1 = z_ar[0]; n3 = z_ar[1];

  N13 = 1./(nlens(n1)*survey.n_gal_conversion_factor);
  if(n1 == n3){
    
    for(int l=0; l<LMAX; l++){
      double l_d = (double)l;
      N24 = kappa_reconstruction_noise(l_d); 
      beam = GaussianBeam(cmb.fwhm, l, 
        like.lmin_kappacmb, like.lmax_kappacmb);

      for(int i=0; i<like.Ntheta; i++){
        for(int j=0; j<like.Ntheta; j++){         
          P1 = Pl_tab(i, l);
          P2 = Pl_tab(j, l);
          ell_prefactor = 1./((2.*l+1.)*fsky);
          N[i][j] += beam*P1*beam*P2*ell_prefactor*N24*N13;
        }
      }

    }

  }
}
void pure_noise_ks_ks(int *z_ar, double *theta, double *dtheta, double **N){
  // N is 2d matrix
  double N13=0, N24=0, beam=0, P1=0, P2=0, ell_prefactor=0;
  double fsky = survey.area*survey.area_conversion_factor/(4.*M_PI);
  int n1,n3;
  static int LMAX = 50000;
  n1 = z_ar[0]; n3 = z_ar[1];

  N13 = pow(survey.sigma_e,2.0)/\
    (2.0*nsource(n1)*survey.n_gal_conversion_factor);
  if(n1 == n3){
    
    for(int l=0; l<LMAX; l++){
      double l_d = (double)l;
      N24 = kappa_reconstruction_noise(l_d);
      beam = GaussianBeam(cmb.fwhm, l, 
        like.lmin_kappacmb, like.lmax_kappacmb);

      for(int i=0; i<like.Ntheta; i++){
        for(int j=0; j<like.Ntheta; j++){  
          P1 = Pl2_tab(i, l);
          P2 = Pl2_tab(j, l);
          ell_prefactor = 1./((2.*l+1.)*fsky);
          N[i][j] += beam*P1*beam*P2*ell_prefactor*N24*N13;
        }
      }

    }

  }
}

double func_for_cov_G_shear_noNN(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[0]; n2 = ar[1]; n3 = ar[2]; n4 = ar[3];
  // printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_shear_tomo(l,n1,n3) * __NN_SNR_FLAG__;
  C24 = C_shear_tomo(l,n2,n4) * __NN_SNR_FLAG__;
  C14 = C_shear_tomo(l,n1,n4) * __NN_SNR_FLAG__;
  C23 = C_shear_tomo(l,n2,n3) * __NN_SNR_FLAG__;
  
  if (n1 == n3){N13= pow(survey.sigma_e,2.0)/(2.*nsource(n1)*survey.n_gal_conversion_factor);}
  if (n1 == n4){N14= pow(survey.sigma_e,2.0)/(2.*nsource(n1)*survey.n_gal_conversion_factor);}
  if (n2 == n3){N23= pow(survey.sigma_e,2.0)/(2.*nsource(n2)*survey.n_gal_conversion_factor);}
  if (n2 == n4){N24= pow(survey.sigma_e,2.0)/(2.*nsource(n2)*survey.n_gal_conversion_factor);}
  
  return (C13*C24+C13*N24+N13*C24 + C14*C23+C14*N23+N14*C23)*4.*M_PI/(survey.area*survey.area_conversion_factor* (2.*l+1.));
}
double func_for_cov_G_shear(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[0]; n2 = ar[1]; n3 = ar[2]; n4 = ar[3];
  // printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_shear_tomo(l,n1,n3) * __NN_SNR_FLAG__;
  C24 = C_shear_tomo(l,n2,n4) * __NN_SNR_FLAG__;
  C14 = C_shear_tomo(l,n1,n4) * __NN_SNR_FLAG__;
  C23 = C_shear_tomo(l,n2,n3) * __NN_SNR_FLAG__;
  
  if (n1 == n3){N13= pow(survey.sigma_e,2.0)/(2.*nsource(n1)*survey.n_gal_conversion_factor);}
  if (n1 == n4){N14= pow(survey.sigma_e,2.0)/(2.*nsource(n1)*survey.n_gal_conversion_factor);}
  if (n2 == n3){N23= pow(survey.sigma_e,2.0)/(2.*nsource(n2)*survey.n_gal_conversion_factor);}
  if (n2 == n4){N24= pow(survey.sigma_e,2.0)/(2.*nsource(n2)*survey.n_gal_conversion_factor);}
  
  return ((C13+N13)*(C24+N24) + (C14+N14)*(C23+N23))*4.*M_PI/(survey.area*survey.area_conversion_factor* (2.*l+1.));
}
double func_for_cov_G_cl_noNN(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[0]; n2 = ar[1]; n3 = ar[2]; n4 = ar[3];
  //printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_cl_tomo(l,n1,n3) * __NN_SNR_FLAG__;
  C24 = C_cl_tomo(l,n2,n4) * __NN_SNR_FLAG__;
  C14 = C_cl_tomo(l,n1,n4) * __NN_SNR_FLAG__;
  C23 = C_cl_tomo(l,n2,n3) * __NN_SNR_FLAG__;
  
  if (n1 == n3){N13= 1./(nlens(n1)*survey.n_gal_conversion_factor);}
  if (n1 == n4){N14= 1./(nlens(n1)*survey.n_gal_conversion_factor);}
  if (n2 == n3){N23= 1./(nlens(n2)*survey.n_gal_conversion_factor);}
  if (n2 == n4){N24= 1./(nlens(n2)*survey.n_gal_conversion_factor);}

  return (C13*C24+C13*N24+N13*C24 + C14*C23+C14*N23+N14*C23)*4.*M_PI/(survey.area*survey.area_conversion_factor* (2.*l+1.));
}
double func_for_cov_G_cl(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[0]; n2 = ar[1]; n3 = ar[2]; n4 = ar[3];
  //printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_cl_tomo(l,n1,n3) * __NN_SNR_FLAG__;
  C24 = C_cl_tomo(l,n2,n4) * __NN_SNR_FLAG__;
  C14 = C_cl_tomo(l,n1,n4) * __NN_SNR_FLAG__;
  C23 = C_cl_tomo(l,n2,n3) * __NN_SNR_FLAG__;
  
  if (n1 == n3){N13= 1./(nlens(n1)*survey.n_gal_conversion_factor);}
  if (n1 == n4){N14= 1./(nlens(n1)*survey.n_gal_conversion_factor);}
  if (n2 == n3){N23= 1./(nlens(n2)*survey.n_gal_conversion_factor);}
  if (n2 == n4){N24= 1./(nlens(n2)*survey.n_gal_conversion_factor);}

  return ((C13+N13)*(C24+N24) + (C14+N14)*(C23+N23))*4.*M_PI/(survey.area*survey.area_conversion_factor* (2.*l+1.));
}
double func_for_cov_G_cl_gl(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[0]; n2 = ar[1]; n3 = ar[2]; n4 = ar[3];
  //printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_cl_tomo(l,n1,n3) * __NN_SNR_FLAG__;
  C24 = C_gl_tomo(l,n2,n4) * __NN_SNR_FLAG__;
  C14 = C_gl_tomo(l,n1,n4) * __NN_SNR_FLAG__;
  C23 = C_cl_tomo(l,n2,n3) * __NN_SNR_FLAG__;
  if (n1 == n3){N13= 1./(nlens(n1)*survey.n_gal_conversion_factor);}
  if (n2 == n3){N23= 1./(nlens(n1)*survey.n_gal_conversion_factor);}
  
  return (C13*C24+C13*N24+C24*N13+C14*C23+C14*N23+C23*N14)*4.*M_PI/(survey.area*survey.area_conversion_factor* (2.*l+1.));
}
double func_for_cov_G_cl_shear(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[0]; n2 = ar[1]; n3 = ar[2]; n4 = ar[3];
  //printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_gl_tomo(l,n1,n3) * __NN_SNR_FLAG__;
  C24 = C_gl_tomo(l,n2,n4) * __NN_SNR_FLAG__;
  C14 = C_gl_tomo(l,n1,n4) * __NN_SNR_FLAG__;
  C23 = C_gl_tomo(l,n2,n3) * __NN_SNR_FLAG__;

  return (C13*C24+ C14*C23)*4.*M_PI/(survey.area*survey.area_conversion_factor* (2.*l+1.));
}
double func_for_cov_G_gl_noNN(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[0]; n2 = ar[1]; n3 = ar[2]; n4 = ar[3];
  //printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_cl_tomo(l,n1,n3) * __NN_SNR_FLAG__;
  C24 = C_shear_tomo(l,n2,n4) * __NN_SNR_FLAG__;
  C14 = C_gl_tomo(l,n1,n4) * __NN_SNR_FLAG__;
  C23 = C_gl_tomo(l,n3,n2) * __NN_SNR_FLAG__;
  
  if (n1 == n3){N13= 1./(nlens(n1)*survey.n_gal_conversion_factor);}
  if (n2 == n4){N24= pow(survey.sigma_e,2.0)/(2.0*nsource(n2)*survey.n_gal_conversion_factor);}
  
  return (C13*C24+C13*N24+N13*C24 + C14*C23+C14*N23+N14*C23)*4.*M_PI/(survey.area*survey.area_conversion_factor* (2.*l+1.));
}
double func_for_cov_G_gl(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[0]; n2 = ar[1]; n3 = ar[2]; n4 = ar[3];
  //printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_cl_tomo(l,n1,n3) * __NN_SNR_FLAG__;
  C24 = C_shear_tomo(l,n2,n4) * __NN_SNR_FLAG__;
  C14 = C_gl_tomo(l,n1,n4) * __NN_SNR_FLAG__;
  C23 = C_gl_tomo(l,n3,n2) * __NN_SNR_FLAG__;
  
  if (n1 == n3){N13= 1./(nlens(n1)*survey.n_gal_conversion_factor);}
  if (n2 == n4){N24= pow(survey.sigma_e,2.0)/(2.0*nsource(n2)*survey.n_gal_conversion_factor);}
  
  return ((C13+N13)*(C24+N24) + (C14+N14)*(C23+N23))*4.*M_PI/(survey.area*survey.area_conversion_factor* (2.*l+1.));
}
double func_for_cov_G_gl_shear(double l, int *ar){
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  n1 = ar[0]; n2 = ar[1]; n3 = ar[2]; n4 = ar[3];
  //printf("n1=%d n2=%d n3=%d n4=%d l=%le\n",n1,n2,n3,n4,l);
  C13 = C_gl_tomo(l,n1,n3) * __NN_SNR_FLAG__;
  C24 = C_shear_tomo(l,n2,n4) * __NN_SNR_FLAG__;
  C14 = C_gl_tomo(l,n1,n4) * __NN_SNR_FLAG__;
  C23 = C_shear_tomo(l,n2,n3) * __NN_SNR_FLAG__;
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
  C13 = C_gl_tomo(l,zl, zs1) * __NN_SNR_FLAG__;
  C24 = C_ks(l, zs2) * __NN_SNR_FLAG__;
  C14 = C_gl_tomo(l, zl, zs2) * __NN_SNR_FLAG__;
  C23 = C_ks(l, zs1) * __NN_SNR_FLAG__;
  return (C13*C24+C14*C23)/((2.*l+1.)*fsky);
}
double func_for_cov_G_ks_shear(double l, int *ar){
  double fsky = survey.area*survey.area_conversion_factor/(4.*M_PI);
  double C13, C14, C23, C24;
  int z1,z2,z3;
  z1 = ar[0]; z2 = ar[1]; z3 = ar[2];
  C13 = C_ks(l,z2) * __NN_SNR_FLAG__;
  C24 = C_shear_tomo(l, z1, z3) * __NN_SNR_FLAG__;
  C14 = C_ks(l,z3) * __NN_SNR_FLAG__;
  C23 = C_shear_tomo(l, z1, z2) * __NN_SNR_FLAG__;
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
  C13 = C_cl_tomo(l,zl1, zl2) * __NN_SNR_FLAG__;
  C24 = C_ks(l, zs) * __NN_SNR_FLAG__;
  C14 = C_gl_tomo(l, zl1, zs) * __NN_SNR_FLAG__;
  C23 = C_gk(l, zl2) * __NN_SNR_FLAG__;
  double N = 0;
  if (zl1==zl2) N = 1./(nlens(zl1)*survey.n_gal_conversion_factor);
  return ((C13+N)*C24+C14*C23)/((2.*l+1.)*fsky);
}
double func_for_cov_G_ks_gl(double l, int *ar){
  double fsky = survey.area*survey.area_conversion_factor/(4.*M_PI);
  double C13, C14, C23, C24;
  int zl,zs1,zs2;
  zs1 = ar[0]; zl = ar[1]; zs2 = ar[2];
  C13 = C_gk(l, zl) * __NN_SNR_FLAG__;
  C24 = C_shear_tomo(l, zs1, zs2) * __NN_SNR_FLAG__;
  C14 = C_ks(l, zs2) * __NN_SNR_FLAG__;
  C23 = C_gl_tomo(l, zl, zs1) * __NN_SNR_FLAG__;
  double N = 0;
  if (zs1==zs2) N = pow(survey.sigma_e,2.0)/(2.0*nsource(zs1)*survey.n_gal_conversion_factor);
  return (C13*(C24+N)+C14*C23)/((2.*l+1.)*fsky);
}
double func_for_cov_G_gk_cl(double l, int *ar){
  double C13, C14, C23, C24, N13=0, N14=0, N23=0, N24=0;
  int n1,n3,n4;
  n1 = ar[0]; n3 = ar[1]; n4 = ar[2];
  double fsky = survey.area*survey.area_conversion_factor/(4.*M_PI);
  C13 = C_cl_tomo(l,n1,n3) * __NN_SNR_FLAG__;
  C24 = C_gk(l,n4) * __NN_SNR_FLAG__;
  C23 = C_gk(l,n3) * __NN_SNR_FLAG__;
  C14 = C_cl_tomo(l,n1,n4) * __NN_SNR_FLAG__;
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
  C13 = C_gk(l, n1) * __NN_SNR_FLAG__;
  C24 = C_gl_tomo(l, n2, zs) * __NN_SNR_FLAG__;
  C23 = C_gl_tomo(l, n1, zs) * __NN_SNR_FLAG__;
  C14 = C_gk(l, n2) * __NN_SNR_FLAG__;
  return (C13*C24+C14*C23)/((2.*l+1.)*fsky);
}
double func_for_cov_G_gk(double l, int *ar){
  double C13, C14, C23, C24, N13=0, N14=0, N23=0, N24=0;
  double fsky = survey.area*survey.area_conversion_factor/(4.*M_PI);
  int n1,n3;
  n1 = ar[0]; n3 = ar[1];
  C13 = C_cl_tomo(l,n1,n3) * __NN_SNR_FLAG__;
  C24 = C_kk(l) * __NN_SNR_FLAG__;
  C14 = C_gk(l,n1) * __NN_SNR_FLAG__;
  C23 = C_gk(l,n3) * __NN_SNR_FLAG__;
  if (n1 == n3){
    N13 = 1./(nlens(n1)*survey.n_gal_conversion_factor);
  }
  N24=kappa_reconstruction_noise(l);
  return ((C13+N13)*(C24+N24)+C14*C23)/((2.*l+1.)*fsky);
}
double func_for_cov_G_gk_noNN(double l, int *ar){
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
  return (C13*C24+C13*N24+N13*C24 +C14*C23)/((2.*l+1.)*fsky);
}
double func_for_cov_G_ks_gk(double l, int *ar){
  double fsky = survey.area*survey.area_conversion_factor/(4.*M_PI);
  double C13, C14, C23, C24;
  int zl,zs;
  zs = ar[0]; zl = ar[1];
  C13 = C_gk(l,zl) * __NN_SNR_FLAG__;
  C24 = C_ks(l, zs) * __NN_SNR_FLAG__;
  C23 = C_gl_tomo(l, zl, zs) * __NN_SNR_FLAG__;
  C14 = C_kk(l) * __NN_SNR_FLAG__;
  double N14 = kappa_reconstruction_noise(l);
  return (C13*C24+(C14+N14)*C23)/((2.*l+1.)*fsky);
}
double func_for_cov_G_ks(double l, int *ar){
  double C13, C14, C23, C24, N13=0, N14=0, N23=0, N24=0;
  int z1,z3;
  z1 = ar[0]; z3 = ar[1];
  double fsky = survey.area*survey.area_conversion_factor/(4.*M_PI);
  C13 = C_shear_tomo(l,z1,z3) * __NN_SNR_FLAG__;
  C24 = C_kk(l) * __NN_SNR_FLAG__;
  C14 = C_ks(l,z1) * __NN_SNR_FLAG__;
  C23 = C_ks(l,z3) * __NN_SNR_FLAG__;
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
double func_for_cov_G_ks_noNN(double l, int *ar){
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
  return (C13*C24+C13*N24+N13*C24 +C14*C23)/((2.*l+1.)*fsky);
}
double func_for_cov_G_kk_shear(double l, int *ar){
  double C13, C14, C23, C24;
  int z1,z2;
  z1 = ar[0]; z2 = ar[1];
  C13 = C_ks(l,z1) * __NN_SNR_FLAG__;
  C23 = C13 * __NN_SNR_FLAG__;
  C24 = C_ks(l, z2) * __NN_SNR_FLAG__;
  C14 = C24 * __NN_SNR_FLAG__;
  return (C13*C24+C14*C23)/((2.*l+1.)*cmb.fsky);
}
double func_for_cov_G_kk_gl(double l, int *ar){
  double C13, C14, C23, C24;
  int zl,zs;
  zl = ar[0]; zs = ar[1];
  C13 = C_gk(l, zl) * __NN_SNR_FLAG__;
  C14 = C13 * __NN_SNR_FLAG__;
  C24 = C_ks(l, zs) * __NN_SNR_FLAG__;
  C23 = C24 * __NN_SNR_FLAG__;
  return (C13*C24+C14*C23)/((2.*l+1.)*cmb.fsky);
}
double func_for_cov_G_kk_cl(double l, int *ar){
  double C13, C14, C23, C24;
  int n1, n2;
  n1=ar[0]; n2=ar[1];
  C13 = C_gk(l,n1) * __NN_SNR_FLAG__;
  C24 = C_gk(l,n2) * __NN_SNR_FLAG__;
  C14 = C13 * __NN_SNR_FLAG__;
  C23 = C24 * __NN_SNR_FLAG__;
  return (C13*C24+C14*C23)/((2.*l+1.)*cmb.fsky);
}
double func_for_cov_G_kk_gk(double l, int *ar){
  double C13, C14, C23, C24;
  int zl;
  zl=ar[0];
  C13 = C_gk(l,zl) * __NN_SNR_FLAG__;
  C14 = C13 * __NN_SNR_FLAG__;
  C24 = C_kk(l) * __NN_SNR_FLAG__;
  C23 = C24 * __NN_SNR_FLAG__;
  double N = kappa_reconstruction_noise(l);
  return (C13*(C24+N)+C14*(C23+N))/((2.*l+1.)*cmb.fsky);
}
double func_for_cov_G_kk_ks(double l, int *ar){
  double C13, C24, C14, C23;
  int zs;
  zs = ar[0];
  C13 = C_kk(l) * __NN_SNR_FLAG__;
  C23 = C13 * __NN_SNR_FLAG__;
  C14 = C_ks(l, zs) * __NN_SNR_FLAG__;
  C24 = C14 * __NN_SNR_FLAG__;
  double N = kappa_reconstruction_noise(l);
  double ans = ( (C13+N)*C24 + C14*(C23+N) )/( (2.*l+1.)*cmb.fsky );
  if(isfinite(ans)){return ans;}
  else{
    printf("cov_G_kk_ks non finite!\n");
    printf("C13=%e C23=%e C14=%e C24=%e N=%e fsky=%e\n",C13, C23, C14, C24, N, cmb.fsky);
	return 0.0;
  }
  //return ( (C13+N)*C24 + C14*(C23+N) )/( (2.*l+1.)*cmb.fsky );
}

double func_for_cov_G_kk(double l, int *ar){
   double C, N;
   C = C_kk(l) * __NN_SNR_FLAG__;
   N = kappa_reconstruction_noise(l);
   return 2.*pow((C+N), 2) / ((2.*l+1.)*cmb.fsky);
}

/****** ******/
// tabulate Fourier covNG
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
  } else if(strcmp(cov_type, "kk_kk")==0) {
    func_for_cov_NG = &cov_NG_kk_kk_bridge;
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

/**************** (21) look-up tables for covariance  *********************/
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

/************ (20) Bridge Functions  ************/
// Bridge functions to unpack z_ar array to call the Fourier cov routines
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
double cov_NG_kk_kk_bridge(double ll1, double ll2, int *z_ar){
  return cov_NG_kk_kk(ll1,ll2);
}
