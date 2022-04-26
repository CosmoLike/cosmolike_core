/* run_covariance_real_fullsky_6x2pt.c

  This file contains the wrapper functions to calculate 6x2pt - 3x2pt 
  covariance matrix in configuration / Fourier / mix space. 
  Assume curved sky and bin-averaged

*/

/************ Function Declarations ************/

// configuration space: {gk, ks} x {shear, ggl, clustering}
void run_cov_gk_shear_real_bin(char *OUTFILE, char *PATH, 
  double *theta, double *dtheta, int Ntheta, 
  int n1, int n2, int pm, int start);
void run_cov_ks_shear_real_bin(char *OUTFILE, char *PATH, 
  double *theta, double *dtheta, int Ntheta, 
  int n1, int n2, int pm, int start);
void run_cov_gk_ggl_real_bin(char *OUTFILE, char *PATH, 
  double *theta, double *dtheta,int Ntheta, 
  int n1, int n2, int start);
void run_cov_ks_ggl_real_bin(char *OUTFILE, char *PATH, 
  double *theta, double *dtheta,int Ntheta, 
  int n1, int n2, int start);
void run_cov_gk_clustering_real_bin(char *OUTFILE, char *PATH, 
  double *theta, double *dtheta, int Ntheta, 
  int n1, int n2, int start);
void run_cov_ks_clustering_real_bin(char *OUTFILE, char *PATH, 
  double *theta, double *dtheta, int Ntheta, 
  int n1, int n2, int start);
void run_cov_gk_gk_real_bin(char *OUTFILE, char *PATH, 
  double *theta, double *dtheta, int Ntheta, 
  int n1, int n2, int start);
void run_cov_ks_gk_real_bin(char *OUTFILE, char *PATH, 
  double *theta, double *dtheta, int Ntheta, 
  int n1, int n2, int start);
void run_cov_ks_ks_real_bin(char *OUTFILE, char *PATH, 
  double *theta, double *dtheta, int Ntheta, 
  int n1, int n2, int start);

// mix space: {kk} x {shear, ggl, clustering, gk, ks}
// kk ell in Logarithmic bins
void run_cov_kk_shear_mix_bin(char *OUTFILE, char *PATH, 
  double *theta, double *dtheta, int Ntheta, 
  double *ell, int Ncl, int n2, int pm, int start);
void run_cov_kk_ggl_mix_bin(char *OUTFILE, char *PATH, 
  double *theta, double *dtheta, int Ntheta, 
  double *ell, int Ncl, int n2, int start);
void run_cov_kk_clustering_mix_bin(char *OUTFILE, char *PATH, 
  double *theta, double *dtheta, int Ntheta, 
  double *ell, int Ncl, int n2, int start);
void run_cov_kk_gk_mix_bin(char *OUTFILE, char *PATH, 
  double *theta, double *dtheta, int Ntheta, 
  double *ell, int Ncl, int n2, int start);
void run_cov_kk_ks_mix_bin(char *OUTFILE, char *PATH, 
  double *theta, double *dtheta, int Ntheta, 
  double *ell, int Ncl, int n2, int start);

// kkk ell in band power
void run_cov_kk_shear_mix_band(char *OUTFILE, char *PATH, 
  double *theta, double *dtheta, int Ntheta,
  int **bindef, int Nbp, int n2, int pm, int start);
void run_cov_kk_ggl_mix_band(char *OUTFILE, char *PATH, 
  double *theta, double *dtheta, int Ntheta, 
  int **bindef, int Nbp, int n2, int start);
void run_cov_kk_clustering_mix_band(char *OUTFILE, char *PATH, 
  double *theta, double *dtheta, int Ntheta, 
  int **bindef, int Nbp, int n2, int start);
void run_cov_kk_gk_mix_band(char *OUTFILE, char *PATH, 
  double *theta, double *dtheta, int Ntheta, 
  int **bindef, int Nbp, int n2, int start);
void run_cov_kk_ks_mix_band(char *OUTFILE, char *PATH, 
  double *theta, double *dtheta, int Ntheta, 
  int **bindef, int Nbp, int n2, int start);

// Fourier space: kk x kk
void run_cov_kk_kk(char *OUTFILE, char *PATH, 
  double *ell, double *dell, int start);
void run_cov_kk_kk_fourier_band(char *OUTFILE, char *PATH, double *ell,
  int **bindef, int Nbp, int start);

/************ Function Implementations ************/

/****** =================== ******/
/****** configuration space ******/
/****** =================== ******/

void run_cov_gk_shear_real_bin(char *OUTFILE, char *PATH, double *theta, double *dtheta, int Ntheta, int n1, int n2, int pm, int start)
{
  int z1,z3,z4,nl1,nl2;
  int zk=-1;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1;
  printf("\nN_gk = %d \n", n1);
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear = %d (%d, %d)\n", n2,z3,z4);

  double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
  like.vtmin = theta[0];
  like.vtmax = theta[Ntheta];
  cov_fullsky_G = create_double_matrix(0, like.Ntheta-1, 0, like.Ntheta-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ntheta-1, 0, like.Ntheta-1);
  cov_gk_shear_real_binned_fullsky(cov_fullsky_G,cov_fullsky_NG,z1,z3,z4,pm, NG, theta, dtheta);

  for (nl1 = 0; nl1 < Ntheta; nl1 ++){
    double t1 = 2./3.*(pow(theta[nl1+1],3.)-pow(theta[nl1],3.))/(pow(theta[nl1+1],2.)-pow(theta[nl1],2.));
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
    double t2 = 2./3.*(pow(theta[nl2+1],3.)-pow(theta[nl2],3.))/(pow(theta[nl2+1],2.)-pow(theta[nl2],2.));
      c_ng = 0.; c_g = 0.;
      // if (test_zoverlap_cov(z1,z3)*test_zoverlap_cov(z1,z4) && NG){ c_ng = cov_NG_cl_shear_real_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,z3,z4,pm);}
      if (test_zoverlap(z1, z3)*test_zoverlap(z1, z4) && covparams.ng){ c_ng = cov_fullsky_NG[nl1][nl2];}
      c_g = cov_fullsky_G[nl1][nl2];
      if(pm==1)fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", 
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + tomo.clustering_Nbin + n1) + nl1,
        Ntheta*(n2) + nl2,
        t1,t2,z1,zk,z3,z4,c_g,c_ng);
      if(pm==0)fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", 
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + tomo.clustering_Nbin + n1) + nl1,
        Ntheta*(tomo.shear_Npowerspectra + n2) + nl2,
        t1,t2,z1,zk,z3,z4,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ntheta-1, 0, like.Ntheta-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ntheta-1, 0, like.Ntheta-1);
  fclose(F1);
}

void run_cov_ks_shear_real_bin(char *OUTFILE, char *PATH, double *theta, double *dtheta, int Ntheta, int n1, int n2, int pm, int start)
{
  int z1,z2,z3,z4,nl1,nl2;
  int zk=-1;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1;
  printf("\nN_ks = %d \n", n1);
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear = %d (%d, %d)\n", n2,z3,z4);

  double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
  like.vtmin = theta[0];
  like.vtmax = theta[Ntheta];
  cov_fullsky_G = create_double_matrix(0, like.Ntheta-1, 0, like.Ntheta-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ntheta-1, 0, like.Ntheta-1);
  cov_ks_shear_real_binned_fullsky(cov_fullsky_G,cov_fullsky_NG,z1,z3,z4,pm, NG, theta, dtheta);

  for (nl1 = 0; nl1 < Ntheta; nl1 ++){
    double t1 = 2./3.*(pow(theta[nl1+1],3.)-pow(theta[nl1],3.))/(pow(theta[nl1+1],2.)-pow(theta[nl1],2.));
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
    double t2 = 2./3.*(pow(theta[nl2+1],3.)-pow(theta[nl2],3.))/(pow(theta[nl2+1],2.)-pow(theta[nl2],2.));
      c_ng = 0.; c_g = 0.;
      // if (test_zoverlap_cov(z1,z3)*test_zoverlap_cov(z1,z4) && NG){ c_ng = cov_NG_cl_shear_real_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,z3,z4,pm);}
      if (covparams.ng){ c_ng = cov_fullsky_NG[nl1][nl2];}
      c_g = cov_fullsky_G[nl1][nl2];
      if(pm==1)fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", 
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + tomo.clustering_Nbin + tomo.clustering_Nbin + n1) + nl1,
        Ntheta*(n2)+nl2,
        t1,t2,z1,zk,z3,z4,c_g,c_ng);
      if(pm==0)fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", 
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + tomo.clustering_Nbin + tomo.clustering_Nbin + n1) + nl1,
        Ntheta*(tomo.shear_Npowerspectra + n2) + nl2,
        t1,t2,z1,zk,z3,z4,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ntheta-1, 0, like.Ntheta-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ntheta-1, 0, like.Ntheta-1);
  fclose(F1);
}

void run_cov_gk_ggl_real_bin(char *OUTFILE, char *PATH, double *theta, double *dtheta, int Ntheta, int n1, int n2, int start)
{
  int z1,z2,zl,zs,nl1,nl2;
  int zk=-1;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1;
  printf("\nN_gk = %d \n", n1);
  zl = ZL(n2); zs = ZS(n2);
  printf("N_tomo_2 = %d (%d, %d)\n", n2,zl,zs);

  double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
  like.vtmin = theta[0];
  like.vtmax = theta[Ntheta];
  cov_fullsky_G = create_double_matrix(0, like.Ntheta-1, 0, like.Ntheta-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ntheta-1, 0, like.Ntheta-1);
  cov_gk_gl_real_binned_fullsky(cov_fullsky_G,cov_fullsky_NG,z1,zl,zs, NG, theta, dtheta);

  for (nl1 = 0; nl1 < Ntheta; nl1 ++){
    double t1 = 2./3.*(pow(theta[nl1+1],3.)-pow(theta[nl1],3.))/(pow(theta[nl1+1],2.)-pow(theta[nl1],2.));
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
    double t2 = 2./3.*(pow(theta[nl2+1],3.)-pow(theta[nl2],3.))/(pow(theta[nl2+1],2.)-pow(theta[nl2],2.));
      c_ng = 0.; 
      c_g = cov_fullsky_G[nl1][nl2];
      // if (z1 == zl && NG && test_zoverlap_cov(z1,zs)*test_zoverlap_cov(zl,zs)){c_ng = cov_NG_cl_gl_real_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,zl,zs);}
      if (z1 == zl && covparams.ng){c_ng = cov_fullsky_NG[nl1][nl2];}
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", 
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + tomo.clustering_Nbin + n1) + nl1,
        Ntheta*(2*tomo.shear_Npowerspectra + n2) + nl2,
        t1,t2,z1,zk,zl,zs,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ntheta-1, 0, like.Ntheta-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ntheta-1, 0, like.Ntheta-1);
  fclose(F1);
}


void run_cov_ks_ggl_real_bin(char *OUTFILE, char *PATH, double *theta, double *dtheta, int Ntheta, int n1, int n2, int start)
{
  int z1,z2,zl,zs,nl1,nl2;
  int zk=-1;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1;
  printf("\nN_ks = %d \n", n1);
  zl = ZL(n2); zs = ZS(n2);
  printf("N_tomo_2 = %d (%d, %d)\n", n2,zl,zs);

  double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
  like.vtmin = theta[0];
  like.vtmax = theta[Ntheta];
  cov_fullsky_G = create_double_matrix(0, like.Ntheta-1, 0, like.Ntheta-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ntheta-1, 0, like.Ntheta-1);
  cov_ks_gl_real_binned_fullsky(cov_fullsky_G,cov_fullsky_NG,z1,zl,zs, NG, theta, dtheta);

  for (nl1 = 0; nl1 < Ntheta; nl1 ++){
    double t1 = 2./3.*(pow(theta[nl1+1],3.)-pow(theta[nl1],3.))/(pow(theta[nl1+1],2.)-pow(theta[nl1],2.));
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
    double t2 = 2./3.*(pow(theta[nl2+1],3.)-pow(theta[nl2],3.))/(pow(theta[nl2+1],2.)-pow(theta[nl2],2.));
      c_ng = 0.; 
      c_g = cov_fullsky_G[nl1][nl2];
      // if (z1 == zl && NG && test_zoverlap_cov(z1,zs)*test_zoverlap_cov(zl,zs)){c_ng = cov_NG_cl_gl_real_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,zl,zs);}
      if (covparams.ng && test_zoverlap_cov(z1,zl)){c_ng = cov_fullsky_NG[nl1][nl2];}
      // fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,Ntheta*(2*tomo.shear_Npowerspectra+n2)+nl2,t1,t2,z1,z2,zl,zs,c_g,c_ng);
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", 
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + tomo.clustering_Nbin + tomo.clustering_Nbin + n1) + nl1,
        Ntheta*(2*tomo.shear_Npowerspectra + n2) + nl2,
        t1,t2,z1,zk,zl,zs,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ntheta-1, 0, like.Ntheta-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ntheta-1, 0, like.Ntheta-1);
  fclose(F1);
}

void run_cov_gk_clustering_real_bin(char *OUTFILE, char *PATH, double *theta, double *dtheta, int Ntheta, int n1, int n2, int start)
{
  int z1,z2,z3,nl1,nl2;
  int zk=-1;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1;
  printf("\nN_gk = %d \n", n1);
  z2 = n2; z3 = n2;
  printf("N_tomo_2 = %d (%d, %d)\n", n2,z2,z3);

  double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
  like.vtmin = theta[0];
  like.vtmax = theta[Ntheta];
  cov_fullsky_G = create_double_matrix(0, like.Ntheta-1, 0, like.Ntheta-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ntheta-1, 0, like.Ntheta-1);
  cov_gk_cl_real_binned_fullsky(cov_fullsky_G,cov_fullsky_NG,z1,z2,z3, NG, theta, dtheta);

  for (nl1 = 0; nl1 < Ntheta; nl1 ++){
    double t1 = 2./3.*(pow(theta[nl1+1],3.)-pow(theta[nl1],3.))/(pow(theta[nl1+1],2.)-pow(theta[nl1],2.));
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
    double t2 = 2./3.*(pow(theta[nl2+1],3.)-pow(theta[nl2],3.))/(pow(theta[nl2+1],2.)-pow(theta[nl2],2.));
      c_ng = 0.; 
      c_g = cov_fullsky_G[nl1][nl2];
      // if (z1 == zl && NG && test_zoverlap_cov(z1,zs)*test_zoverlap_cov(zl,zs)){c_ng = cov_NG_cl_gl_real_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,zl,zs);}
      if (z1 == z2 && covparams.ng){c_ng = cov_fullsky_NG[nl1][nl2];}
      // fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,Ntheta*(2*tomo.shear_Npowerspectra+n2)+nl2,t1,t2,z1,z2,zl,zs,c_g,c_ng);
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", 
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + tomo.clustering_Nbin + n1) + nl1,
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + n2) + nl2,
        t1,t2,z1,zk,z2,z3,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ntheta-1, 0, like.Ntheta-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ntheta-1, 0, like.Ntheta-1);
  fclose(F1);
}

void run_cov_ks_clustering_real_bin(char *OUTFILE, char *PATH, double *theta, double *dtheta, int Ntheta, int n1, int n2, int start)
{
  int z1,z2,zs,nl1,nl2;
  int zk=-1;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  zs = n1;
  printf("\nN_ks = %d \n", n1);
  z1 = n2; z2 = n2;
  printf("N_cl_2 = %d (%d, %d)\n", n2,z1,z2);

  double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
  like.vtmin = theta[0];
  like.vtmax = theta[Ntheta];
  cov_fullsky_G = create_double_matrix(0, like.Ntheta-1, 0, like.Ntheta-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ntheta-1, 0, like.Ntheta-1);
  cov_ks_cl_real_binned_fullsky(cov_fullsky_G,cov_fullsky_NG,zs,z1,z2, NG, theta, dtheta);

  for (nl1 = 0; nl1 < Ntheta; nl1 ++){
    double t1 = 2./3.*(pow(theta[nl1+1],3.)-pow(theta[nl1],3.))/(pow(theta[nl1+1],2.)-pow(theta[nl1],2.));
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
    double t2 = 2./3.*(pow(theta[nl2+1],3.)-pow(theta[nl2],3.))/(pow(theta[nl2+1],2.)-pow(theta[nl2],2.));
      c_ng = 0.; 
      c_g = cov_fullsky_G[nl1][nl2];
      // if (z1 == zl && NG && test_zoverlap_cov(z1,zs)*test_zoverlap_cov(zl,zs)){c_ng = cov_NG_cl_gl_real_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,zl,zs);}
      if (covparams.ng && test_zoverlap_cov(z1,zs)){c_ng = cov_fullsky_NG[nl1][nl2];}
      // fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,Ntheta*(2*tomo.shear_Npowerspectra+n2)+nl2,t1,t2,z1,z2,zl,zs,c_g,c_ng);
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", 
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + tomo.clustering_Nbin + tomo.clustering_Nbin + n1) + nl1,
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + n2) + nl2,
        t1,t2,zs,zk,z1,z2,c_g,c_ng);

    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ntheta-1, 0, like.Ntheta-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ntheta-1, 0, like.Ntheta-1);
  fclose(F1);
}

void run_cov_gk_gk_real_bin(char *OUTFILE, char *PATH, double *theta, double *dtheta, int Ntheta, int n1, int n2, int start)
{
  int z1,z2,nl1,nl2;
  int zk=-1;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1; z2 = n2;
  printf("\nN_gk_1 = %d \n", n1);
  printf("N_gk_2 = %d\n", n2);

  double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
  like.vtmin = theta[0];
  like.vtmax = theta[Ntheta];
  cov_fullsky_G = create_double_matrix(0, like.Ntheta-1, 0, like.Ntheta-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ntheta-1, 0, like.Ntheta-1);
  cov_gk_gk_real_binned_fullsky(cov_fullsky_G,cov_fullsky_NG,z1,z2, NG, theta, dtheta);

  for (nl1 = 0; nl1 < Ntheta; nl1 ++){
    double t1 = 2./3.*(pow(theta[nl1+1],3.)-pow(theta[nl1],3.))/(pow(theta[nl1+1],2.)-pow(theta[nl1],2.));
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
    double t2 = 2./3.*(pow(theta[nl2+1],3.)-pow(theta[nl2],3.))/(pow(theta[nl2+1],2.)-pow(theta[nl2],2.));
      c_ng = 0.; 
      c_g = cov_fullsky_G[nl1][nl2];
      // if (z1 == zl && NG && test_zoverlap_cov(z1,zs)*test_zoverlap_cov(zl,zs)){c_ng = cov_NG_cl_gl_real_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,zl,zs);}
      if (z1 == z2 && covparams.ng){c_ng = cov_fullsky_NG[nl1][nl2];}
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", 
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + tomo.clustering_Nbin + n1) + nl1,
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + tomo.clustering_Nbin + n2) + nl2,
        t1,t2,z1,zk,z2,zk,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ntheta-1, 0, like.Ntheta-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ntheta-1, 0, like.Ntheta-1);
  fclose(F1);
}

void run_cov_ks_gk_real_bin(char *OUTFILE, char *PATH, double *theta, double *dtheta, int Ntheta, int n1, int n2, int start)
{
  int z1,z2,zl,zs,nl1,nl2;
  int zk=-1;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  zs = n1;
  printf("\nN_ks = %d \n", n1);
  zl = n2;
  printf("N_gk = %d\n", n2);

  double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
  like.vtmin = theta[0];
  like.vtmax = theta[Ntheta];
  cov_fullsky_G = create_double_matrix(0, like.Ntheta-1, 0, like.Ntheta-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ntheta-1, 0, like.Ntheta-1);
  cov_ks_gk_real_binned_fullsky(cov_fullsky_G,cov_fullsky_NG,zs,zl, NG, theta, dtheta);

  for (nl1 = 0; nl1 < Ntheta; nl1 ++){
    double t1 = 2./3.*(pow(theta[nl1+1],3.)-pow(theta[nl1],3.))/(pow(theta[nl1+1],2.)-pow(theta[nl1],2.));
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
    double t2 = 2./3.*(pow(theta[nl2+1],3.)-pow(theta[nl2],3.))/(pow(theta[nl2+1],2.)-pow(theta[nl2],2.));
      c_ng = 0.; 
      c_g = cov_fullsky_G[nl1][nl2];
      // if (z1 == zl && NG && test_zoverlap_cov(z1,zs)*test_zoverlap_cov(zl,zs)){c_ng = cov_NG_cl_gl_real_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,zl,zs);}
      if (covparams.ng && test_zoverlap_cov(zl,zs)){c_ng = cov_fullsky_NG[nl1][nl2];}
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", 
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + tomo.clustering_Nbin + tomo.clustering_Nbin + n1) + nl1,
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + tomo.clustering_Nbin + n2) + nl2,
        t1,t2,zs,zk,zl,zk,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ntheta-1, 0, like.Ntheta-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ntheta-1, 0, like.Ntheta-1);
  fclose(F1);
}


void run_cov_ks_ks_real_bin(char *OUTFILE, char *PATH, double *theta, double *dtheta, int Ntheta, int n1, int n2, int start)
{
  int z1,z2,nl1,nl2;
  int zk=-1;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1;
  printf("\nN_ks_1 = %d \n", n1);
  z2 = n2;
  printf("N_ks_2 = %d \n", n2);

  double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
  like.vtmin = theta[0];
  like.vtmax = theta[Ntheta];
  cov_fullsky_G = create_double_matrix(0, like.Ntheta-1, 0, like.Ntheta-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ntheta-1, 0, like.Ntheta-1);
  cov_ks_ks_real_binned_fullsky(cov_fullsky_G,cov_fullsky_NG,z1,z2, NG, theta, dtheta);

  for (nl1 = 0; nl1 < Ntheta; nl1 ++){
    double t1 = 2./3.*(pow(theta[nl1+1],3.)-pow(theta[nl1],3.))/(pow(theta[nl1+1],2.)-pow(theta[nl1],2.));
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
    double t2 = 2./3.*(pow(theta[nl2+1],3.)-pow(theta[nl2],3.))/(pow(theta[nl2+1],2.)-pow(theta[nl2],2.));
      c_ng = 0.; 
      c_g = cov_fullsky_G[nl1][nl2];
      // if (z1 == zl && NG && test_zoverlap_cov(z1,zs)*test_zoverlap_cov(zl,zs)){c_ng = cov_NG_cl_gl_real_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,zl,zs);}
      if (covparams.ng){c_ng = cov_fullsky_NG[nl1][nl2];}
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", 
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + tomo.clustering_Nbin + tomo.clustering_Nbin + n1) + nl1,
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + tomo.clustering_Nbin + tomo.clustering_Nbin + n2) + nl2,
        t1,t2,z1,zk,z2,zk,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ntheta-1, 0, like.Ntheta-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ntheta-1, 0, like.Ntheta-1);
  fclose(F1);
}

/****** ========= ******/
/****** mix space ******/
/****** ========= ******/

// ell in Logarithmic bins
void run_cov_kk_shear_mix_bin(char *OUTFILE, char *PATH, 
  double *theta, double *dtheta, int Ntheta, 
  double *ell, int Ncl, int n2, int pm, int start)
{
  int z1,z2,z3,z4,nl1,nl2;
  int zk=-1;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear = %d (%d, %d)\n", n2,z3,z4);

  double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
  like.vtmin = theta[0];
  like.vtmax = theta[Ntheta];
  cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ntheta-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ntheta-1);
  cov_kk_shear_mix_binned_fullsky(cov_fullsky_G,cov_fullsky_NG,z3,z4,pm, NG, theta, dtheta, ell);

  for (nl1 = 0; nl1 < Ncl; nl1 ++){
    double t1 = 2./3.*(pow(ell[nl1+1],3.)-pow(ell[nl1],3.))/(pow(ell[nl1+1],2.)-pow(ell[nl1],2.));
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
      double t2 = 2./3.*(pow(theta[nl2+1],3.)-pow(theta[nl2],3.))/(pow(theta[nl2+1],2.)-pow(theta[nl2],2.));
      c_ng = 0.; c_g = 0.;
      // if (test_zoverlap_cov(z1,z3)*test_zoverlap_cov(z1,z4) && NG){ c_ng = cov_NG_cl_shear_real_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,z3,z4,pm);}
      if (covparams.ng){ c_ng = cov_fullsky_NG[nl1][nl2];}
      c_g = cov_fullsky_G[nl1][nl2];
      if(pm==1)fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", 
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + tomo.clustering_Nbin + tomo.clustering_Nbin + tomo.shear_Nbin) + nl1,
        Ntheta*(n2)+nl2,
        t1,t2,zk,zk,z3,z4,c_g,c_ng);
      if(pm==0)fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", 
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + tomo.clustering_Nbin + tomo.clustering_Nbin + tomo.shear_Nbin) + nl1,
        Ntheta*(tomo.shear_Npowerspectra+n2)+nl2, 
        t1,t2,zk,zk,z3,z4,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ntheta-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ntheta-1);
  fclose(F1);
}

void run_cov_kk_ggl_mix_bin(char *OUTFILE, char *PATH, 
  double *theta, double *dtheta, int Ntheta, 
  double *ell, int Ncl, int n2, int start)
{
  int z1,z2,z3,z4,nl1,nl2;
  int zk=-1;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z3 = ZL(n2); z4 = ZS(n2);
  printf("N_ggl = %d (%d, %d)\n", n2,z3,z4);

  double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
  like.vtmin = theta[0];
  like.vtmax = theta[Ntheta];
  cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ntheta-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ntheta-1);
  cov_kk_gl_mix_binned_fullsky(cov_fullsky_G,cov_fullsky_NG,z3,z4, NG, theta, dtheta, ell);

  for (nl1 = 0; nl1 < Ncl; nl1 ++){
    double t1 = 2./3.*(pow(ell[nl1+1],3.)-pow(ell[nl1],3.))/(pow(ell[nl1+1],2.)-pow(ell[nl1],2.));
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
      double t2 = 2./3.*(pow(theta[nl2+1],3.)-pow(theta[nl2],3.))/(pow(theta[nl2+1],2.)-pow(theta[nl2],2.));
      c_ng = 0.; c_g = 0.;
      // if (test_zoverlap_cov(z1,z3)*test_zoverlap_cov(z1,z4) && NG){ c_ng = cov_NG_cl_shear_real_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,z3,z4,pm);}
      if (covparams.ng){ c_ng = cov_fullsky_NG[nl1][nl2];}
      c_g = cov_fullsky_G[nl1][nl2];
      // fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+n1)+nl1,Ntheta*(n2)+nl2,t1,t2,z1,zk,z3,z4,c_g,c_ng);
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", 
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + tomo.clustering_Nbin + tomo.clustering_Nbin + tomo.shear_Nbin) + nl1,
        Ntheta*(2*tomo.shear_Npowerspectra + n2)+nl2,
        t1,t2,zk,zk,z3,z4,c_g,c_ng);

    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ntheta-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ntheta-1);
  fclose(F1);
}

void run_cov_kk_clustering_mix_bin(char *OUTFILE, char *PATH, 
  double *theta, double *dtheta, int Ntheta, 
  double *ell, int Ncl, int n2, int start)
{
  int z1,z2,z3,z4,nl1,nl2;
  int zk=-1;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z3 = n2; z4 = n2;
  printf("N_cl = %d (%d, %d)\n", n2,z3,z4);

  double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
  like.vtmin = theta[0];
  like.vtmax = theta[Ntheta];
  cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ntheta-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ntheta-1);
  cov_kk_cl_mix_binned_fullsky(cov_fullsky_G,cov_fullsky_NG,z3,z4, NG, theta, dtheta, ell);

  for (nl1 = 0; nl1 < Ncl; nl1 ++){
    double t1 = 2./3.*(pow(ell[nl1+1],3.)-pow(ell[nl1],3.))/(pow(ell[nl1+1],2.)-pow(ell[nl1],2.));
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
      double t2 = 2./3.*(pow(theta[nl2+1],3.)-pow(theta[nl2],3.))/(pow(theta[nl2+1],2.)-pow(theta[nl2],2.));
      c_ng = 0.; c_g = 0.;
      // if (test_zoverlap_cov(z1,z3)*test_zoverlap_cov(z1,z4) && NG){ c_ng = cov_NG_cl_shear_real_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,z3,z4,pm);}
      if (covparams.ng){ c_ng = cov_fullsky_NG[nl1][nl2];}
      c_g = cov_fullsky_G[nl1][nl2];
      // fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+n1)+nl1,Ntheta*(n2)+nl2,t1,t2,z1,zk,z3,z4,c_g,c_ng);
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", 
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + tomo.clustering_Nbin + tomo.clustering_Nbin + tomo.shear_Nbin) + nl1,
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + n2) + nl2,
        t1,t2,zk,zk,z3,z4,c_g,c_ng);

    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ntheta-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ntheta-1);
  fclose(F1);
}

void run_cov_kk_gk_mix_bin(char *OUTFILE, char *PATH, 
  double *theta, double *dtheta, int Ntheta, 
  double *ell, int Ncl, int n2, int start)
{
  int z1,z2,z3,z4,nl1,nl2;
  int zk=-1;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z3 = n2;
  printf("N_gk_2 = %d\n", n2);

  double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
  like.vtmin = theta[0];
  like.vtmax = theta[Ntheta];
  cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ntheta-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ntheta-1);
  cov_kk_gk_mix_binned_fullsky(cov_fullsky_G,cov_fullsky_NG,z3, NG, theta, dtheta, ell);

  for (nl1 = 0; nl1 < Ncl; nl1 ++){
    double t1 = 2./3.*(pow(ell[nl1+1],3.)-pow(ell[nl1],3.))/(pow(ell[nl1+1],2.)-pow(ell[nl1],2.));
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
      double t2 = 2./3.*(pow(theta[nl2+1],3.)-pow(theta[nl2],3.))/(pow(theta[nl2+1],2.)-pow(theta[nl2],2.));
      c_ng = 0.; c_g = 0.;
      // if (test_zoverlap_cov(z1,z3)*test_zoverlap_cov(z1,z4) && NG){ c_ng = cov_NG_cl_shear_real_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,z3,z4,pm);}
      if (covparams.ng){ c_ng = cov_fullsky_NG[nl1][nl2];}
      c_g = cov_fullsky_G[nl1][nl2];
      // fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+n1)+nl1,Ntheta*(n2)+nl2,t1,t2,z1,zk,z3,z4,c_g,c_ng);
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", 
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + tomo.clustering_Nbin + tomo.clustering_Nbin + tomo.shear_Nbin) + nl1,
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + tomo.clustering_Nbin + n2) + nl2,
        t1,t2,zk,zk,z3,zk,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ntheta-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ntheta-1);
  fclose(F1);
}

void run_cov_kk_ks_mix_bin(char *OUTFILE, char *PATH, 
  double *theta, double *dtheta, int Ntheta, 
  double *ell, int Ncl, int n2, int start)
{
  int z1,z2,z3,z4,nl1,nl2;
  int zk=-1;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z3 = n2;
  printf("N_ks_2 = %d\n", n2);

  double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
  like.vtmin = theta[0];
  like.vtmax = theta[Ntheta];
  cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ntheta-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ntheta-1);
  cov_kk_ks_mix_binned_fullsky(cov_fullsky_G,cov_fullsky_NG,z3, NG, theta, dtheta, ell);

  for (nl1 = 0; nl1 < Ncl; nl1 ++){
    double t1 = 2./3.*(pow(ell[nl1+1],3.)-pow(ell[nl1],3.))/(pow(ell[nl1+1],2.)-pow(ell[nl1],2.));
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
      double t2 = 2./3.*(pow(theta[nl2+1],3.)-pow(theta[nl2],3.))/(pow(theta[nl2+1],2.)-pow(theta[nl2],2.));
      c_ng = 0.; c_g = 0.;
      // if (test_zoverlap_cov(z1,z3)*test_zoverlap_cov(z1,z4) && NG){ c_ng = cov_NG_cl_shear_real_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,z3,z4,pm);}
      if (covparams.ng){ c_ng = cov_fullsky_NG[nl1][nl2];}
      c_g = cov_fullsky_G[nl1][nl2];
      // fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+n1)+nl1,Ntheta*(n2)+nl2,t1,t2,z1,zk,z3,z4,c_g,c_ng);
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", 
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + tomo.clustering_Nbin + tomo.clustering_Nbin + tomo.shear_Nbin) + nl1,
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + tomo.clustering_Nbin + tomo.clustering_Nbin + n2) + nl2,
        t1,t2,zk,zk,z3,zk,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ntheta-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ntheta-1);
  fclose(F1);
}

// ell in band powers
void run_cov_kk_shear_mix_band(char *OUTFILE, char *PATH, 
  double *theta, double *dtheta, int Ntheta, 
  int **bindef, int Nbp, int n2, int pm, int start)
{
  int z1,z2,z3,z4,nl1,nl2;
  int zk=-1;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear = %d (%d, %d)\n", n2,z3,z4);

  double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
  like.vtmin = theta[0];
  like.vtmax = theta[Ntheta];
  cov_fullsky_G = create_double_matrix(0, like.Nbp-1, 0, like.Ntheta-1);
  cov_fullsky_NG = create_double_matrix(0, like.Nbp-1, 0, like.Ntheta-1);
  cov_kk_shear_mix_banded_fullsky(cov_fullsky_G,cov_fullsky_NG,z3,z4,pm, NG, theta, dtheta, bindef);

  for (nl1 = 0; nl1 < Nbp; nl1 ++){
    double t1 = 2./3.*(pow(bindef[nl1][1],3.)-pow(bindef[nl1][0],3.))/(pow(bindef[nl1][1],2.)-pow(bindef[nl1][1],2.));
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
      double t2 = 2./3.*(pow(theta[nl2+1],3.)-pow(theta[nl2],3.))/(pow(theta[nl2+1],2.)-pow(theta[nl2],2.));
      c_ng = 0.; c_g = 0.;
      // if (test_zoverlap_cov(z1,z3)*test_zoverlap_cov(z1,z4) && NG){ c_ng = cov_NG_cl_shear_real_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,z3,z4,pm);}
      if (covparams.ng){ c_ng = cov_fullsky_NG[nl1][nl2];}
      c_g = cov_fullsky_G[nl1][nl2];
      if(pm==1)fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", 
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + tomo.clustering_Nbin + tomo.clustering_Nbin + tomo.shear_Nbin) + nl1,
        Ntheta*(n2)+nl2,
        t1,t2,zk,zk,z3,z4,c_g,c_ng);
      if(pm==0)fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", 
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + tomo.clustering_Nbin + tomo.clustering_Nbin + tomo.shear_Nbin) + nl1,
        Ntheta*(tomo.shear_Npowerspectra+n2)+nl2, 
        t1,t2,zk,zk,z3,z4,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Nbp-1, 0, like.Ntheta-1);
  free_double_matrix(cov_fullsky_NG,0, like.Nbp-1, 0, like.Ntheta-1);
  fclose(F1);
}

void run_cov_kk_ggl_mix_band(char *OUTFILE, char *PATH, 
  double *theta, double *dtheta, int Ntheta, 
  int **bindef, int Nbp, int n2, int start)
{
  int z1,z2,z3,z4,nl1,nl2;
  int zk=-1;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z3 = ZL(n2); z4 = ZS(n2);
  printf("N_ggl = %d (%d, %d)\n", n2,z3,z4);

  double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
  like.vtmin = theta[0];
  like.vtmax = theta[Ntheta];
  cov_fullsky_G = create_double_matrix(0, like.Nbp-1, 0, like.Ntheta-1);
  cov_fullsky_NG = create_double_matrix(0, like.Nbp-1, 0, like.Ntheta-1);
  cov_kk_gl_mix_banded_fullsky(cov_fullsky_G,cov_fullsky_NG,z3,z4, NG, theta, dtheta, bindef);

  for (nl1 = 0; nl1 < Nbp; nl1 ++){
    double t1 = 2./3.*(pow(bindef[nl1][1],3.)-pow(bindef[nl1][0],3.))/(pow(bindef[nl1][1],2.)-pow(bindef[nl1][0],2.));
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
      double t2 = 2./3.*(pow(theta[nl2+1],3.)-pow(theta[nl2],3.))/(pow(theta[nl2+1],2.)-pow(theta[nl2],2.));
      c_ng = 0.; c_g = 0.;
      // if (test_zoverlap_cov(z1,z3)*test_zoverlap_cov(z1,z4) && NG){ c_ng = cov_NG_cl_shear_real_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,z3,z4,pm);}
      if (covparams.ng){ c_ng = cov_fullsky_NG[nl1][nl2];}
      c_g = cov_fullsky_G[nl1][nl2];
      // fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+n1)+nl1,Ntheta*(n2)+nl2,t1,t2,z1,zk,z3,z4,c_g,c_ng);
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", 
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + tomo.clustering_Nbin + tomo.clustering_Nbin + tomo.shear_Nbin) + nl1,
        Ntheta*(2*tomo.shear_Npowerspectra + n2)+nl2,
        t1,t2,zk,zk,z3,z4,c_g,c_ng);

    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Nbp-1, 0, like.Ntheta-1);
  free_double_matrix(cov_fullsky_NG,0, like.Nbp-1, 0, like.Ntheta-1);
  fclose(F1);
}

void run_cov_kk_clustering_mix_band(char *OUTFILE, char *PATH, 
  double *theta, double *dtheta, int Ntheta, 
  int **bindef, int Nbp, int n2, int start)
{
  int z1,z2,z3,z4,nl1,nl2;
  int zk=-1;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z3 = n2; z4 = n2;
  printf("N_cl = %d (%d, %d)\n", n2,z3,z4);

  double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
  like.vtmin = theta[0];
  like.vtmax = theta[Ntheta];
  cov_fullsky_G = create_double_matrix(0, like.Nbp-1, 0, like.Ntheta-1);
  cov_fullsky_NG = create_double_matrix(0, like.Nbp-1, 0, like.Ntheta-1);
  cov_kk_cl_mix_banded_fullsky(cov_fullsky_G,cov_fullsky_NG,z3,z4, NG, theta, dtheta, bindef);

  for (nl1 = 0; nl1 < Nbp; nl1 ++){
    double t1 = 2./3.*(pow(bindef[nl1][1],3.)-pow(bindef[nl1][0],3.))/(pow(bindef[nl1][1],2.)-pow(bindef[nl1][0],2.));
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
      double t2 = 2./3.*(pow(theta[nl2+1],3.)-pow(theta[nl2],3.))/(pow(theta[nl2+1],2.)-pow(theta[nl2],2.));
      c_ng = 0.; c_g = 0.;
      // if (test_zoverlap_cov(z1,z3)*test_zoverlap_cov(z1,z4) && NG){ c_ng = cov_NG_cl_shear_real_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,z3,z4,pm);}
      if (covparams.ng){ c_ng = cov_fullsky_NG[nl1][nl2];}
      c_g = cov_fullsky_G[nl1][nl2];
      // fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+n1)+nl1,Ntheta*(n2)+nl2,t1,t2,z1,zk,z3,z4,c_g,c_ng);
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", 
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + tomo.clustering_Nbin + tomo.clustering_Nbin + tomo.shear_Nbin) + nl1,
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + n2) + nl2,
        t1,t2,zk,zk,z3,z4,c_g,c_ng);

    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ntheta-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ntheta-1);
  fclose(F1);
}

void run_cov_kk_gk_mix_band(char *OUTFILE, char *PATH, 
  double *theta, double *dtheta, int Ntheta, 
  int **bindef, int Nbp, int n2, int start)
{
  int z1,z2,z3,z4,nl1,nl2;
  int zk=-1;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z3 = n2;
  printf("N_gk_2 = %d\n", n2);

  double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
  like.vtmin = theta[0];
  like.vtmax = theta[Ntheta];
  cov_fullsky_G = create_double_matrix(0, like.Nbp-1, 0, like.Ntheta-1);
  cov_fullsky_NG = create_double_matrix(0, like.Nbp-1, 0, like.Ntheta-1);
  cov_kk_gk_mix_banded_fullsky(cov_fullsky_G,cov_fullsky_NG,z3, NG, theta, dtheta, bindef);

  for (nl1 = 0; nl1 < Nbp; nl1 ++){
    double t1 = 2./3.*(pow(bindef[nl1][1],3.)-pow(bindef[nl1][0],3.))/(pow(bindef[nl1][1],2.)-pow(bindef[nl1][0],2.));
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
      double t2 = 2./3.*(pow(theta[nl2+1],3.)-pow(theta[nl2],3.))/(pow(theta[nl2+1],2.)-pow(theta[nl2],2.));
      c_ng = 0.; c_g = 0.;
      // if (test_zoverlap_cov(z1,z3)*test_zoverlap_cov(z1,z4) && NG){ c_ng = cov_NG_cl_shear_real_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,z3,z4,pm);}
      if (covparams.ng){ c_ng = cov_fullsky_NG[nl1][nl2];}
      c_g = cov_fullsky_G[nl1][nl2];
      // fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+n1)+nl1,Ntheta*(n2)+nl2,t1,t2,z1,zk,z3,z4,c_g,c_ng);
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", 
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + tomo.clustering_Nbin + tomo.clustering_Nbin + tomo.shear_Nbin) + nl1,
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + tomo.clustering_Nbin + n2) + nl2,
        t1,t2,zk,zk,z3,zk,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Nbp-1, 0, like.Ntheta-1);
  free_double_matrix(cov_fullsky_NG,0, like.Nbp-1, 0, like.Ntheta-1);
  fclose(F1);
}

void run_cov_kk_ks_mix_band(char *OUTFILE, char *PATH, 
  double *theta, double *dtheta, int Ntheta, 
  int **bindef, int Nbp, int n2, int start)
{
  int z1,z2,z3,z4,nl1,nl2;
  int zk=-1;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z3 = n2;
  printf("N_ks_2 = %d\n", n2);

  double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
  like.vtmin = theta[0];
  like.vtmax = theta[Ntheta];
  cov_fullsky_G = create_double_matrix(0, like.Nbp-1, 0, like.Ntheta-1);
  cov_fullsky_NG = create_double_matrix(0, like.Nbp-1, 0, like.Ntheta-1);
  cov_kk_ks_mix_banded_fullsky(cov_fullsky_G,cov_fullsky_NG,z3, NG, theta, dtheta, bindef);

  for (nl1 = 0; nl1 < Nbp; nl1 ++){
    double t1 = 2./3.*(pow(bindef[nl1][1],3.)-pow(bindef[nl1][0],3.))/(pow(bindef[nl1][1],2.)-pow(bindef[nl1][0],2.));
    for (nl2 = 0; nl2 < Ntheta; nl2 ++){
      double t2 = 2./3.*(pow(theta[nl2+1],3.)-pow(theta[nl2],3.))/(pow(theta[nl2+1],2.)-pow(theta[nl2],2.));
      c_ng = 0.; c_g = 0.;
      // if (test_zoverlap_cov(z1,z3)*test_zoverlap_cov(z1,z4) && NG){ c_ng = cov_NG_cl_shear_real_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,z3,z4,pm);}
      if (covparams.ng){ c_ng = cov_fullsky_NG[nl1][nl2];}
      c_g = cov_fullsky_G[nl1][nl2];
      // fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+n1)+nl1,Ntheta*(n2)+nl2,t1,t2,z1,zk,z3,z4,c_g,c_ng);
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", 
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + tomo.clustering_Nbin + tomo.clustering_Nbin + tomo.shear_Nbin) + nl1,
        Ntheta*(2*tomo.shear_Npowerspectra + tomo.ggl_Npowerspectra + tomo.clustering_Nbin + tomo.clustering_Nbin + n2) + nl2,
        t1,t2,zk,zk,z3,zk,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Nbp-1, 0, like.Ntheta-1);
  free_double_matrix(cov_fullsky_NG,0, like.Nbp-1, 0, like.Ntheta-1);
  fclose(F1);
}


/****** ============= ******/
/****** Fourier space ******/
/****** ============= ******/

// Note that this one does not use template function
/*
void run_cov_kk_kk(char *OUTFILE, char *PATH, 
  double *ell, double *dell, int start)
{
  int nl1,nl2,i,j,weight;
  int zk=-1;
  double c_ng, c_g;
  FILE *F1;
  char filename[300];
  // sprintf(filename,"%scov/%s_%s_cov_kkkk_Nell%d_Ns%d_Ng%d_%d",PATH,survey.name, cmb.name,like.Ncl, tomo.shear_Nbin, tomo.clustering_Nbin, start);
  // printf("Saving to: %s\n",filename);
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");

  for (nl1 = 0; nl1 < like.Nbp; nl1 ++){
    int l1_min = (int)ceil(bindef[nl1][0]);
    int l1_max = (int)ceil(bindef[nl1][1]);
    int N_l1 = l1_max - l1_min + 1;
    double t1 = 2./3.*(pow(l1_max,3.)-pow(l1_min,3.))/(pow(l1_max,2.)-pow(l1_min,2.));

    for (nl2 = 0; nl2 < like.Nbp; nl2 ++){
      int l2_min = (int)ceil(bindef[nl2][0]);
      int l2_max = (int)ceil(bindef[nl2][1]);
      int N_l2 = l2_max - l2_min + 1;
      double t2 = 2./3.*(pow(l2_max,3.)-pow(l2_min,3.))/(pow(l2_max,2.)-pow(l2_min,2.));

      // This routine skip the template function
      // Need to apply binning matrix
      c_ng = 0.; c_g = 0.;
      for (int L = l1_min; L <= l1_max; L++){
        // Diagonal blocks
        if (nl1==nl2){
          //c_g = cov_G_kk_kk(ell[nl1], dell[nl1]);
          c_g += cov_G_kk_kk((double)L, 1.0) * ;
        }
      }
      
      if (covparams.ng){
        c_ng = cov_NG_kk_kk(ell[nl1],ell[nl2]);
      }
      
      i=like.Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin)+nl1;
      j=like.Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin)+nl2;
      fprintf(F1, "%d %d %e %e %d %d %d %d %e %e\n",i,j,t1,t2,zk,zk,zk,zk,c_g,c_ng);
      //printf("l1, l2, Cg, Cng = %le, %le, %le, %le\n", ell[nl1], ell[nl2], c_g, c_ng);
      }
   }
   fclose(F1);
}*/
void run_cov_kk_kk(char *OUTFILE, char *PATH, 
  double *ell, double *dell,int start)
{
  int nl1,nl2,i,j;
  int zk=-1;
  double c_ng, c_g;
  FILE *F1;
  char filename[300];
  // sprintf(filename,"%scov/%s_%s_cov_kkkk_Nell%d_Ns%d_Ng%d_%d",PATH,survey.name, cmb.name,like.Ncl, tomo.shear_Nbin, tomo.clustering_Nbin, start);
  // printf("Saving to: %s\n",filename);
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    double t1 = 2./3.*(pow(ell[nl1+1],3.)-pow(ell[nl1],3.))/(pow(ell[nl1+1],2.)-pow(ell[nl1],2.));
    for (nl2 = 0; nl2 < like.Ncl; nl2 ++){
      double t2 = 2./3.*(pow(ell[nl2+1],3.)-pow(ell[nl2],3.))/(pow(ell[nl2+1],2.)-pow(ell[nl2],2.));
      c_ng = 0.; c_g = 0.;
      if (covparams.ng){
        c_ng = cov_NG_kk_kk(ell[nl1],ell[nl2]);
      }
      if (nl1==nl2){
        c_g = cov_G_kk_kk(ell[nl1], dell[nl1]);
      }
      i=like.Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin)+nl1;
      j=like.Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin)+nl2;
      fprintf(F1, "%d %d %e %e %d %d %d %d %e %e\n",i,j,t1,t2,zk,zk,zk,zk,c_g,c_ng);
      //printf("l1, l2, Cg, Cng = %le, %le, %le, %le\n", ell[nl1], ell[nl2], c_g, c_ng);
      }
   }
   fclose(F1);
}


void run_cov_kk_kk_fourier_band(char *OUTFILE, char *PATH, double *ell,
  int **bindef, int Nbp, int start)
{
  int nl1,nl2,i,j;
  int zk=-1;
  double c_ng, c_g;
  //double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");

  double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;

  cov_fullsky_G = create_double_matrix(0, like.Nbp-1, 0, like.Nbp-1);
  cov_fullsky_NG = create_double_matrix(0, like.Nbp-1, 0, like.Nbp-1);
  cov_kk_kk_fourier_banded(cov_fullsky_G, cov_fullsky_NG, NG, ell, bindef);

  for (nl1 = 0; nl1 < Nbp; nl1 ++){
    double t1 = 2./3.*(pow(bindef[nl1][1],3.)-pow(bindef[nl1][0],3.))/(pow(bindef[nl1][1],2.)-pow(bindef[nl1][1],2.));
    for (nl2 = 0; nl2 < Nbp; nl2 ++){
      double t2 = 2./3.*(pow(bindef[nl2][1],3.)-pow(bindef[nl2][0],3.))/(pow(bindef[nl2][1],2.)-pow(bindef[nl2][1],2.));

      c_ng = 0.; c_g = 0.;

      if (covparams.ng){ c_ng = cov_fullsky_NG[nl1][nl2];}
      c_g = cov_fullsky_G[nl1][nl2];

      i=like.Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin)+nl1;
      j=like.Ntheta*(2*tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin)+nl2;
      fprintf(F1, "%d %d %e %e %d %d %d %d %e %e\n",i,j,t1,t2,zk,zk,zk,zk,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Nbp-1, 0, like.Nbp-1);
  free_double_matrix(cov_fullsky_NG,0, like.Nbp-1, 0, like.Nbp-1);
  fclose(F1);
}