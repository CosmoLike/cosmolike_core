#ifndef __COV_BINNED_SIMPLE__
#define __COV_BINNED_SIMPLE__
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

#endif