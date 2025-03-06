void run_cov_shear_shear_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n1, int n2, int start);
void run_cov_ggl_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n1, int n2, int start);
void run_cov_clustering_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n1, int n2, int start);
void run_cov_ggl_shear_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n1, int n2, int start);
void run_cov_clustering_shear_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n1, int n2, int start);
void run_cov_clustering_ggl_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n1, int n2, int start);

void run_cov_gk_shear_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n1, int n2, int start);
void run_cov_ks_shear_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n1, int n2, int start);
void run_cov_gk_ggl_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n1, int n2, int start);
void run_cov_ks_ggl_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n1, int n2, int start);
void run_cov_gk_clustering_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n1, int n2, int start);
void run_cov_ks_clustering_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n1, int n2, int start);
void run_cov_gk_gk_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n1, int n2, int start);
void run_cov_ks_gk_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n1, int n2, int start);
void run_cov_ks_ks_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n1, int n2, int start);

void run_cov_kk_shear_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n2, int start);
void run_cov_kk_ggl_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n2, int start);
void run_cov_kk_clustering_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n2, int start);
void run_cov_kk_gk_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n2, int start);
void run_cov_kk_ks_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n2, int start);

void run_cov_kk_kk_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int start);

/////////////////

void run_cov_clustering_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n1, int n2, int start)
{
  int NG = covparams.ng;
  int z1,z2,z3,z4,nl1,nl2;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1; z2 = n1;
  //  printf("\nN_cl_1 = %d \n", n1);
  z3 = n2; z4 = n2;
  printf("N_cl_1 = %d, N_cl_2 = %d\n", n1,n2);

  double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
  like.lmin = ell[0];
  like.lmax = ell[Ncl];
  cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_cl_cl_fourier_binned(cov_fullsky_G,cov_fullsky_NG,z1,z2,z3,z4, NG, ell);

  for (nl1 = 0; nl1 < Ncl; nl1 ++){
    double t1 = (floor(ell[nl1+1]) + ceil(ell[nl1]))/2.; // average ell in band
    for (nl2 = 0; nl2 < Ncl; nl2 ++){
    double t2 = (floor(ell[nl2+1]) + ceil(ell[nl2]))/2.; // average ell in band
      c_ng = 0.; 
      c_g = cov_fullsky_G[nl1][nl2];
      if (z1 == z3 && NG){c_ng = cov_fullsky_NG[nl1][nl2];}
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra + n2)+nl2,t1,t2,z1,z2,z3,z4,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ncl-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ncl-1);
  fclose(F1);
}

void run_cov_ggl_fourier_bin(char *OUTFILE, char *PATH, double *ell,int Ncl, int n1, int n2, int start)
{
  int NG = covparams.ng;
  int zl1,zl2,zs1,zs2,nl1,nl2;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");

  zl1 = ZL(n1); zs1 = ZS(n1);
  printf("\nN_tomo_1 = %d (%d, %d)\n", n1,zl1,zs1);
  zl2 = ZL(n2); zs2 = ZS(n2);
  printf("N_tomo_2 = %d (%d, %d)\n", n2,zl2,zs2);

  double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
  like.lmin = ell[0];
  like.lmax = ell[Ncl];
  cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_gl_gl_fourier_binned(cov_fullsky_G,cov_fullsky_NG,zl1,zs1,zl2,zs2, NG, ell);

  for (nl1 = 0; nl1 < Ncl; nl1 ++){
    double t1 = (floor(ell[nl1+1]) + ceil(ell[nl1]))/2.; // average ell in band
    for (nl2 = 0; nl2 < Ncl; nl2 ++){
    double t2 = (floor(ell[nl2+1]) + ceil(ell[nl2]))/2.; // average ell in band
      c_ng = 0.;
      // if(NG && zl1 == zl2 && test_zoverlap_cov(zl1,zs1)*test_zoverlap_cov(zl2,zs2)){c_ng = cov_NG_gl_gl_fourier_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],zl1,zs1,zl2,zs2);}
      if(NG && zl1 == zl2 && test_zoverlap_cov(zl1,zs1)*test_zoverlap_cov(zl2,zs2)){c_ng = cov_fullsky_NG[nl1][nl2];}
      c_g = cov_fullsky_G[nl1][nl2];
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+n1)+nl1,Ncl*(tomo.shear_Npowerspectra+n2)+nl2,t1,t2,zl1,zs1,zl2,zs2,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ncl-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ncl-1);
  fclose(F1);
}
void run_cov_shear_shear_fourier_bin(char *OUTFILE, char *PATH,double *ell, int Ncl, int n1, int n2, int start)
{
  int NG = covparams.ng;
  int z1,z2,z3,z4,nl1,nl2;
  double c_ng, c_g,sn;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  z1 = Z1(n1); z2 = Z2(n1);
  printf("N_shear1 = %d (%d,%d)\n", n1,z1,z2);
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear2 = %d (%d, %d)\n",n2,z3,z4);
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");


  double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
  like.lmin = ell[0];
  like.lmax = ell[Ncl];
  cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_shear_shear_fourier_binned(cov_fullsky_G,cov_fullsky_NG,z1,z2,z3,z4, NG, ell);

  for (nl1 = 0; nl1 < Ncl; nl1 ++){
    double t1 = (floor(ell[nl1+1]) + ceil(ell[nl1]))/2.; // average ell in band
    for (nl2 = 0; nl2 < Ncl; nl2 ++){
    double t2 = (floor(ell[nl2+1]) + ceil(ell[nl2]))/2.; // average ell in band
      c_ng = 0.; c_g = 0.;sn = 0.;
      if(NG){ c_ng = cov_fullsky_NG[nl1][nl2];}
      c_g = cov_fullsky_G[nl1][nl2];
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*n1+nl1,Ncl*(n2)+nl2, t1,t2,z1,z2,z3,z4,c_g,c_ng); 

    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ncl-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ncl-1);
  fclose(F1);
}

void run_cov_ggl_shear_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n1, int n2, int start)
{
  int NG = covparams.ng;
  int zl,zs,z3,z4,nl1,nl2;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w"); 
  zl = ZL(n1); zs = ZS(n1);
  printf("\nN_ggl = %d (%d, %d)\n", n1,zl,zs);
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear = %d (%d, %d)\n", n2,z3,z4);

  double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
  like.lmin = ell[0];
  like.lmax = ell[Ncl];
  cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_gl_shear_fourier_binned(cov_fullsky_G,cov_fullsky_NG,zl,zs,z3,z4, NG, ell);

  for (nl1 = 0; nl1 < Ncl; nl1 ++){
    double t1 = (floor(ell[nl1+1]) + ceil(ell[nl1]))/2.; // average ell in band
    for (nl2 = 0; nl2 < Ncl; nl2 ++){
    double t2 = (floor(ell[nl2+1]) + ceil(ell[nl2]))/2.; // average ell in band
      c_ng = 0.; c_g = 0.;
      // if (test_zoverlap_cov(zl,zs)*test_zoverlap_cov(zl,z3)*test_zoverlap_cov(zl,z4) && NG){ c_ng = cov_NG_gl_shear_fourier_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],zl,zs,z3,z4,pm);}
      if (test_zoverlap_cov(zl,zs)*test_zoverlap_cov(zl,z3)*test_zoverlap_cov(zl,z4) && NG){ c_ng = cov_fullsky_NG[nl1][nl2];}
      c_g = cov_fullsky_G[nl1][nl2];
      fprintf(F1, "%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+n1)+nl1,Ncl*(n2)+nl2,t1,t2,zl,zs,z3,z4,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ncl-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ncl-1);
  fclose(F1);
}
void run_cov_clustering_shear_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n1, int n2, int start)
{
  int NG = covparams.ng;
  int z1,z2,z3,z4,nl1,nl2;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1; z2 = n1;
  printf("\nN_cl = %d \n", n1);
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear = %d (%d, %d)\n", n2,z3,z4);

  double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
  like.lmin = ell[0];
  like.lmax = ell[Ncl];
  cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_cl_shear_fourier_binned(cov_fullsky_G,cov_fullsky_NG,z1,z2,z3,z4, NG, ell);

  for (nl1 = 0; nl1 < Ncl; nl1 ++){
    double t1 = (floor(ell[nl1+1]) + ceil(ell[nl1]))/2.; // average ell in band
    for (nl2 = 0; nl2 < Ncl; nl2 ++){
    double t2 = (floor(ell[nl2+1]) + ceil(ell[nl2]))/2.; // average ell in band
      c_ng = 0.; c_g = 0.;
      // if (test_zoverlap_cov(z1,z3)*test_zoverlap_cov(z1,z4) && NG){ c_ng = cov_NG_cl_shear_fourier_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,z3,z4,pm);}
      if (test_zoverlap_cov(z1,z3)*test_zoverlap_cov(z1,z4) && NG){ c_ng = cov_fullsky_NG[nl1][nl2];}
      c_g = cov_fullsky_G[nl1][nl2];
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,Ncl*(n2)+nl2,t1,t2,z1,z2,z3,z4,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ncl-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ncl-1);
fclose(F1);
}

void run_cov_clustering_ggl_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n1, int n2, int start)
{
  int NG = covparams.ng;
  int z1,z2,zl,zs,nl1,nl2;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1; z2 = n1;
  printf("\nN_cl_1 = %d \n", n1);
  zl = ZL(n2); zs = ZS(n2);
  printf("N_tomo_2 = %d (%d, %d)\n", n2,zl,zs);

  double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
  like.lmin = ell[0];
  like.lmax = ell[Ncl];
  cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_cl_gl_fourier_binned(cov_fullsky_G,cov_fullsky_NG,z1,z2,zl,zs, NG, ell);

  for (nl1 = 0; nl1 < Ncl; nl1 ++){
    double t1 = (floor(ell[nl1+1]) + ceil(ell[nl1]))/2.; // average ell in band
    for (nl2 = 0; nl2 < Ncl; nl2 ++){
    double t2 = (floor(ell[nl2+1]) + ceil(ell[nl2]))/2.; // average ell in band
      c_ng = 0.; 
      c_g = cov_fullsky_G[nl1][nl2];
      // if (z1 == zl && NG && test_zoverlap_cov(z1,zs)*test_zoverlap_cov(zl,zs)){c_ng = cov_NG_cl_gl_fourier_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,zl,zs);}
      if (z1 == zl && NG && test_zoverlap_cov(z1,zs)*test_zoverlap_cov(zl,zs)){c_ng = cov_fullsky_NG[nl1][nl2];}
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,Ncl*(tomo.shear_Npowerspectra+n2)+nl2,t1,t2,z1,z2,zl,zs,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ncl-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ncl-1);
  fclose(F1);
}

//////////

void run_cov_gk_shear_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n1, int n2, int start)
{
  int NG = covparams.ng;
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
  like.lmin = ell[0];
  like.lmax = ell[Ncl];
  cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_gk_shear_fourier_binned(cov_fullsky_G,cov_fullsky_NG,z1,z3,z4, NG, ell);

  for (nl1 = 0; nl1 < Ncl; nl1 ++){
    double t1 = (floor(ell[nl1+1]) + ceil(ell[nl1]))/2.; // average ell in band
    for (nl2 = 0; nl2 < Ncl; nl2 ++){
      double t2 = (floor(ell[nl2+1]) + ceil(ell[nl2]))/2.; // average ell in band
      c_ng = 0.; c_g = 0.;
      // if (test_zoverlap_cov(z1,z3)*test_zoverlap_cov(z1,z4) && NG){ c_ng = cov_NG_cl_shear_fourier_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,z3,z4,pm);}
      if (test_zoverlap(z1, z3)*test_zoverlap(z1, z4) && covparams.ng){ c_ng = cov_fullsky_NG[nl1][nl2];}
      c_g = cov_fullsky_G[nl1][nl2];
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+n1)+nl1,Ncl*(n2)+nl2,t1,t2,z1,zk,z3,z4,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ncl-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ncl-1);
fclose(F1);
}

void run_cov_ks_shear_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n1, int n2, int start)
{
  int NG = covparams.ng;
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
  like.lmin = ell[0];
  like.lmax = ell[Ncl];
  cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_ks_shear_fourier_binned(cov_fullsky_G,cov_fullsky_NG,z1,z3,z4, NG, ell);

  for (nl1 = 0; nl1 < Ncl; nl1 ++){
    double t1 = (floor(ell[nl1+1]) + ceil(ell[nl1]))/2.; // average ell in band
    for (nl2 = 0; nl2 < Ncl; nl2 ++){
    double t2 = (floor(ell[nl2+1]) + ceil(ell[nl2]))/2.; // average ell in band
      c_ng = 0.; c_g = 0.;
      // if (test_zoverlap_cov(z1,z3)*test_zoverlap_cov(z1,z4) && NG){ c_ng = cov_NG_cl_shear_fourier_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,z3,z4,pm);}
      if (covparams.ng){ c_ng = cov_fullsky_NG[nl1][nl2];}
      c_g = cov_fullsky_G[nl1][nl2];
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+n1)+nl1,Ncl*(n2)+nl2,t1,t2,z1,zk,z3,z4,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ncl-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ncl-1);
fclose(F1);
}

void run_cov_gk_ggl_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n1, int n2, int start)
{
  int NG = covparams.ng;
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
  like.lmin = ell[0];
  like.lmax = ell[Ncl];
  cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_gk_gl_fourier_binned(cov_fullsky_G,cov_fullsky_NG,z1,zl,zs, NG, ell);

  for (nl1 = 0; nl1 < Ncl; nl1 ++){
    double t1 = (floor(ell[nl1+1]) + ceil(ell[nl1]))/2.; // average ell in band
    for (nl2 = 0; nl2 < Ncl; nl2 ++){
      double t2 = (floor(ell[nl2+1]) + ceil(ell[nl2]))/2.; // average ell in band
      c_ng = 0.; 
      c_g = cov_fullsky_G[nl1][nl2];
      // if (z1 == zl && NG && test_zoverlap_cov(z1,zs)*test_zoverlap_cov(zl,zs)){c_ng = cov_NG_cl_gl_fourier_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,zl,zs);}
      if (z1 == zl && covparams.ng){c_ng = cov_fullsky_NG[nl1][nl2];}
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+n1)+nl1,Ncl*(tomo.shear_Npowerspectra+n2)+nl2,t1,t2,z1,zk,zl,zs,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ncl-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ncl-1);
  fclose(F1);
}


void run_cov_ks_ggl_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n1, int n2, int start)
{
  int NG = covparams.ng;
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
  like.lmin = ell[0];
  like.lmax = ell[Ncl];
  cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_ks_gl_fourier_binned(cov_fullsky_G,cov_fullsky_NG,z1,zl,zs, NG, ell);

  for (nl1 = 0; nl1 < Ncl; nl1 ++){
    double t1 = (floor(ell[nl1+1]) + ceil(ell[nl1]))/2.; // average ell in band
    for (nl2 = 0; nl2 < Ncl; nl2 ++){
    double t2 = (floor(ell[nl2+1]) + ceil(ell[nl2]))/2.; // average ell in band
      c_ng = 0.; 
      c_g = cov_fullsky_G[nl1][nl2];
      // if (z1 == zl && NG && test_zoverlap_cov(z1,zs)*test_zoverlap_cov(zl,zs)){c_ng = cov_NG_cl_gl_fourier_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,zl,zs);}
      if (covparams.ng && test_zoverlap_cov(z1,zl)){c_ng = cov_fullsky_NG[nl1][nl2];}
      // fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,Ncl*(tomo.shear_Npowerspectra+n2)+nl2,t1,t2,z1,z2,zl,zs,c_g,c_ng);
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+n1)+nl1,Ncl*(tomo.shear_Npowerspectra+n2)+nl2,t1,t2,z1,zk,zl,zs,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ncl-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ncl-1);
  fclose(F1);
}

void run_cov_gk_clustering_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n1, int n2, int start)
{
  int NG = covparams.ng;
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
  like.lmin = ell[0];
  like.lmax = ell[Ncl];
  cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_gk_cl_fourier_binned(cov_fullsky_G,cov_fullsky_NG,z1,z2,z3, NG, ell);

  for (nl1 = 0; nl1 < Ncl; nl1 ++){
    double t1 = (floor(ell[nl1+1]) + ceil(ell[nl1]))/2.; // average ell in band
    for (nl2 = 0; nl2 < Ncl; nl2 ++){
      double t2 = (floor(ell[nl2+1]) + ceil(ell[nl2]))/2.; // average ell in band
      c_ng = 0.; 
      c_g = cov_fullsky_G[nl1][nl2];
      // if (z1 == zl && NG && test_zoverlap_cov(z1,zs)*test_zoverlap_cov(zl,zs)){c_ng = cov_NG_cl_gl_fourier_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,zl,zs);}
      if (z1 == z2 && covparams.ng){c_ng = cov_fullsky_NG[nl1][nl2];}
      // fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,Ncl*(tomo.shear_Npowerspectra+n2)+nl2,t1,t2,z1,z2,zl,zs,c_g,c_ng);
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+n1)+nl1,Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n2)+nl2,t1,t2,z1,zk,z2,z3,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ncl-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ncl-1);
  fclose(F1);
}

void run_cov_ks_clustering_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n1, int n2, int start)
{
  int NG = covparams.ng;
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
  like.lmin = ell[0];
  like.lmax = ell[Ncl];
  cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_ks_cl_fourier_binned(cov_fullsky_G,cov_fullsky_NG,zs,z1,z2, NG, ell);

  for (nl1 = 0; nl1 < Ncl; nl1 ++){
    double t1 = (floor(ell[nl1+1]) + ceil(ell[nl1]))/2.; // average ell in band
    for (nl2 = 0; nl2 < Ncl; nl2 ++){
    double t2 = (floor(ell[nl2+1]) + ceil(ell[nl2]))/2.; // average ell in band
      c_ng = 0.; 
      c_g = cov_fullsky_G[nl1][nl2];
      // if (z1 == zl && NG && test_zoverlap_cov(z1,zs)*test_zoverlap_cov(zl,zs)){c_ng = cov_NG_cl_gl_fourier_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,zl,zs);}
      if (covparams.ng && test_zoverlap_cov(z1,zs)){c_ng = cov_fullsky_NG[nl1][nl2];}
      // fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,Ncl*(tomo.shear_Npowerspectra+n2)+nl2,t1,t2,z1,z2,zl,zs,c_g,c_ng);
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+n1)+nl1,Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n2)+nl2,t1,t2,zs,zk,z1,z2,c_g,c_ng);

    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ncl-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ncl-1);
  fclose(F1);
}

void run_cov_gk_gk_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n1, int n2, int start)
{
  int NG = covparams.ng;
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
  like.lmin = ell[0];
  like.lmax = ell[Ncl];
  cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_gk_gk_fourier_binned(cov_fullsky_G,cov_fullsky_NG,z1,z2, NG, ell);

  for (nl1 = 0; nl1 < Ncl; nl1 ++){
    double t1 = (floor(ell[nl1+1]) + ceil(ell[nl1]))/2.; // average ell in band
    for (nl2 = 0; nl2 < Ncl; nl2 ++){
    double t2 = (floor(ell[nl2+1]) + ceil(ell[nl2]))/2.; // average ell in band
      c_ng = 0.; 
      c_g = cov_fullsky_G[nl1][nl2];
      // if (z1 == zl && NG && test_zoverlap_cov(z1,zs)*test_zoverlap_cov(zl,zs)){c_ng = cov_NG_cl_gl_fourier_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,zl,zs);}
      if (z1 == z2 && covparams.ng){c_ng = cov_fullsky_NG[nl1][nl2];}
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+n1)+nl1,Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+n2)+nl2,t1,t2,z1,zk,z2,zk,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ncl-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ncl-1);
  fclose(F1);
}

void run_cov_ks_gk_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n1, int n2, int start)
{
  int NG = covparams.ng;
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
  like.lmin = ell[0];
  like.lmax = ell[Ncl];
  cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_ks_gk_fourier_binned(cov_fullsky_G,cov_fullsky_NG,zs,zl, NG, ell);

  for (nl1 = 0; nl1 < Ncl; nl1 ++){
    double t1 = (floor(ell[nl1+1]) + ceil(ell[nl1]))/2.; // average ell in band
    for (nl2 = 0; nl2 < Ncl; nl2 ++){
    double t2 = (floor(ell[nl2+1]) + ceil(ell[nl2]))/2.; // average ell in band
      c_ng = 0.; 
      c_g = cov_fullsky_G[nl1][nl2];
      // if (z1 == zl && NG && test_zoverlap_cov(z1,zs)*test_zoverlap_cov(zl,zs)){c_ng = cov_NG_cl_gl_fourier_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,zl,zs);}
      if (covparams.ng && test_zoverlap_cov(zl,zs)){c_ng = cov_fullsky_NG[nl1][nl2];}
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+n1)+nl1,Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+n2)+nl2,t1,t2,zs,zk,zl,zk,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ncl-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ncl-1);
  fclose(F1);
}


void run_cov_ks_ks_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n1, int n2, int start)
{
  int NG = covparams.ng;
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
  like.lmin = ell[0];
  like.lmax = ell[Ncl];
  cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_ks_ks_fourier_binned(cov_fullsky_G,cov_fullsky_NG,z1,z2, NG, ell);

  for (nl1 = 0; nl1 < Ncl; nl1 ++){
    double t1 = (floor(ell[nl1+1]) + ceil(ell[nl1]))/2.; // average ell in band
    for (nl2 = 0; nl2 < Ncl; nl2 ++){
    double t2 = (floor(ell[nl2+1]) + ceil(ell[nl2]))/2.; // average ell in band
      c_ng = 0.; 
      c_g = cov_fullsky_G[nl1][nl2];
      // if (z1 == zl && NG && test_zoverlap_cov(z1,zs)*test_zoverlap_cov(zl,zs)){c_ng = cov_NG_cl_gl_fourier_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,zl,zs);}
      if (covparams.ng){c_ng = cov_fullsky_NG[nl1][nl2];}
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+n1)+nl1,Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+n2)+nl2,t1,t2,z1,zk,z2,zk,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ncl-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ncl-1);
  fclose(F1);
}

// Mix covs
void run_cov_kk_shear_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n2, int start)
{
  int NG = covparams.ng;
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
  like.lmin = ell[0];
  like.lmax = ell[Ncl];
  cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_kk_shear_fourier_binned(cov_fullsky_G,cov_fullsky_NG,z3,z4, NG, ell);

  for (nl1 = 0; nl1 < Ncl; nl1 ++){
    double t1 = (floor(ell[nl1+1]) + ceil(ell[nl1]))/2.; // average ell in band
    for (nl2 = 0; nl2 < Ncl; nl2 ++){
      double t2 = (floor(ell[nl2+1]) + ceil(ell[nl2]))/2.; // average ell in band
      c_ng = 0.; c_g = 0.;
      // if (test_zoverlap_cov(z1,z3)*test_zoverlap_cov(z1,z4) && NG){ c_ng = cov_NG_cl_shear_fourier_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,z3,z4,pm);}
      if (covparams.ng){ c_ng = cov_fullsky_NG[nl1][nl2];}
      c_g = cov_fullsky_G[nl1][nl2];
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin)+nl1,Ncl*(n2)+nl2,t1,t2,zk,zk,z3,z4,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ncl-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ncl-1);
  fclose(F1);
}

void run_cov_kk_ggl_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n2, int start)
{
  int NG = covparams.ng;
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
  like.lmin = ell[0];
  like.lmax = ell[Ncl];
  cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_kk_gl_fourier_binned(cov_fullsky_G,cov_fullsky_NG,z3,z4, NG, ell);

  for (nl1 = 0; nl1 < Ncl; nl1 ++){
    double t1 = (floor(ell[nl1+1]) + ceil(ell[nl1]))/2.; // average ell in band
    for (nl2 = 0; nl2 < Ncl; nl2 ++){
      double t2 = (floor(ell[nl2+1]) + ceil(ell[nl2]))/2.; // average ell in band
      c_ng = 0.; c_g = 0.;
      // if (test_zoverlap_cov(z1,z3)*test_zoverlap_cov(z1,z4) && NG){ c_ng = cov_NG_cl_shear_fourier_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,z3,z4,pm);}
      if (covparams.ng){ c_ng = cov_fullsky_NG[nl1][nl2];}
      c_g = cov_fullsky_G[nl1][nl2];
      // fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+n1)+nl1,Ncl*(n2)+nl2,t1,t2,z1,zk,z3,z4,c_g,c_ng);
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin)+nl1,Ncl*(tomo.shear_Npowerspectra+n2)+nl2,t1,t2,zk,zk,z3,z4,c_g,c_ng);

    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ncl-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ncl-1);
  fclose(F1);
}

void run_cov_kk_clustering_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n2, int start)
{
  int NG = covparams.ng;
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
  like.lmin = ell[0];
  like.lmax = ell[Ncl];
  cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_kk_cl_fourier_binned(cov_fullsky_G,cov_fullsky_NG,z3,z4, NG, ell);

  for (nl1 = 0; nl1 < Ncl; nl1 ++){
    double t1 = (floor(ell[nl1+1]) + ceil(ell[nl1]))/2.; // average ell in band
    for (nl2 = 0; nl2 < Ncl; nl2 ++){
      double t2 = (floor(ell[nl2+1]) + ceil(ell[nl2]))/2.; // average ell in band
      c_ng = 0.; c_g = 0.;
      // if (test_zoverlap_cov(z1,z3)*test_zoverlap_cov(z1,z4) && NG){ c_ng = cov_NG_cl_shear_fourier_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,z3,z4,pm);}
      if (covparams.ng){ c_ng = cov_fullsky_NG[nl1][nl2];}
      c_g = cov_fullsky_G[nl1][nl2];
      // fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+n1)+nl1,Ncl*(n2)+nl2,t1,t2,z1,zk,z3,z4,c_g,c_ng);
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin)+nl1,Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n2)+nl2,t1,t2,zk,zk,z3,z4,c_g,c_ng);

    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ncl-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ncl-1);
  fclose(F1);
}

void run_cov_kk_gk_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n2, int start)
{
  int NG = covparams.ng;
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
  like.lmin = ell[0];
  like.lmax = ell[Ncl];
  cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_kk_gk_fourier_binned(cov_fullsky_G,cov_fullsky_NG,z3, NG, ell);

  for (nl1 = 0; nl1 < Ncl; nl1 ++){
    double t1 = (floor(ell[nl1+1]) + ceil(ell[nl1]))/2.; // average ell in band
    for (nl2 = 0; nl2 < Ncl; nl2 ++){
      double t2 = (floor(ell[nl2+1]) + ceil(ell[nl2]))/2.; // average ell in band
      c_ng = 0.; c_g = 0.;
      // if (test_zoverlap_cov(z1,z3)*test_zoverlap_cov(z1,z4) && NG){ c_ng = cov_NG_cl_shear_fourier_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,z3,z4,pm);}
      if (covparams.ng){ c_ng = cov_fullsky_NG[nl1][nl2];}
      c_g = cov_fullsky_G[nl1][nl2];
      // fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+n1)+nl1,Ncl*(n2)+nl2,t1,t2,z1,zk,z3,z4,c_g,c_ng);
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin)+nl1,Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+n2)+nl2,t1,t2,zk,zk,z3,zk,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ncl-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ncl-1);
  fclose(F1);
}

void run_cov_kk_ks_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n2, int start)
{
  int NG = covparams.ng;
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
  like.lmin = ell[0];
  like.lmax = ell[Ncl];
  cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_kk_ks_fourier_binned(cov_fullsky_G,cov_fullsky_NG,z3, NG, ell);

  for (nl1 = 0; nl1 < Ncl; nl1 ++){
    double t1 = (floor(ell[nl1+1]) + ceil(ell[nl1]))/2.; // average ell in band
    for (nl2 = 0; nl2 < Ncl; nl2 ++){
      double t2 = (floor(ell[nl2+1]) + ceil(ell[nl2]))/2.; // average ell in band
      c_ng = 0.; c_g = 0.;
      // if (test_zoverlap_cov(z1,z3)*test_zoverlap_cov(z1,z4) && NG){ c_ng = cov_NG_cl_shear_fourier_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,z3,z4,pm);}
      if (covparams.ng){ c_ng = cov_fullsky_NG[nl1][nl2];}
      c_g = cov_fullsky_G[nl1][nl2];
      // fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+n1)+nl1,Ncl*(n2)+nl2,t1,t2,z1,zk,z3,z4,c_g,c_ng);
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin)+nl1,Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+n2)+nl2,t1,t2,zk,zk,z3,zk,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ncl-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ncl-1);
  fclose(F1);
}


// kk_kk Fourier
void run_cov_kk_kk_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int start)
{
  int NG = covparams.ng;
  int nl1,nl2;
  int zk=-1;
  double c_ng, c_g;
  double fsky = survey.area/41253.0;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");

  double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
  like.lmin = ell[0];
  like.lmax = ell[Ncl];
  cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
  cov_kk_kk_fourier_binned(cov_fullsky_G,cov_fullsky_NG, NG, ell);

  int i,j;
  for (nl1 = 0; nl1 < Ncl; nl1 ++){
    double t1 = (floor(ell[nl1+1]) + ceil(ell[nl1]))/2.; // average ell in band
    for (nl2 = 0; nl2 < Ncl; nl2 ++){
    double t2 = (floor(ell[nl2+1]) + ceil(ell[nl2]))/2.; // average ell in band
      c_ng = 0.; c_g = 0.;
      // if (test_zoverlap_cov(z1,z3)*test_zoverlap_cov(z1,z4) && NG){ c_ng = cov_NG_cl_shear_fourier_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,z3,z4,pm);}
      if (covparams.ng){ c_ng = cov_fullsky_NG[nl1][nl2];}
      c_g = cov_fullsky_G[nl1][nl2];
      i = Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin)+nl1;
      j = Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin)+nl2;
      // fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+n1)+nl1,Ncl*(n2)+nl2,t1,t2,z1,zk,z3,z4,c_g,c_ng);
      fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", i,j,t1,t2,zk,zk,zk,zk,c_g,c_ng);
    }
  }
  free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ncl-1);
  free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ncl-1);
  fclose(F1);
}

// ///////////////
// ///////////////


// void run_cov_kk_shear_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n2, int start)
// {
//   int z1,z2,z3,z4,nl1,nl2;
//   int zk=-1;
//   double c_ng, c_g;
//   double fsky = survey.area/41253.0;
//   FILE *F1;
//   char filename[300];
//   sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
//   F1 =fopen(filename,"w");
//   z3 = Z1(n2); z4 = Z2(n2);
//   printf("N_shear = %d (%d, %d)\n", n2,z3,z4);

//   double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
//   like.lmin = ell[0];
//   like.lmax = ell[Ncl];
//   cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
//   cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
//   cov_kk_shear_fourier_binned(cov_fullsky_G,cov_fullsky_NG,z3,z4, NG, ell);

//   for (nl1 = 0; nl1 < Ncl; nl1 ++){
//     double t1 = (floor(ell[nl1+1]) + ceil(ell[nl1]))/2.; // average ell in band
//     for (nl2 = 0; nl2 < Ncl; nl2 ++){
//     double t2 = (floor(ell[nl2+1]) + ceil(ell[nl2]))/2.; // average ell in band
//       c_ng = 0.; c_g = 0.;
//       // if (test_zoverlap_cov(z1,z3)*test_zoverlap_cov(z1,z4) && NG){ c_ng = cov_NG_cl_shear_fourier_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,z3,z4,pm);}
//       if (covparams.ng){ c_ng = cov_fullsky_NG[nl1][nl2];}
//       c_g = cov_fullsky_G[nl1][nl2];
//       fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin)+nl1,Ncl*(n2)+nl2,t1,t2,zk,zk,z3,z4,c_g,c_ng);
//     }
//   }
//   free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ncl-1);
//   free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ncl-1);
//   fclose(F1);
// }

// void run_cov_kk_ggl_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n2, int start)
// {
//   int z1,z2,z3,z4,nl1,nl2;
//   int zk=-1;
//   double c_ng, c_g;
//   double fsky = survey.area/41253.0;
//   FILE *F1;
//   char filename[300];
//   sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
//   F1 =fopen(filename,"w");
//   z3 = ZL(n2); z4 = ZS(n2);
//   printf("N_ggl = %d (%d, %d)\n", n2,z3,z4);

//   double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
//   like.lmin = ell[0];
//   like.lmax = ell[Ncl];
//   cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
//   cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
//   cov_kk_gl_fourier_binned(cov_fullsky_G,cov_fullsky_NG,z3,z4, NG, ell);

//   for (nl1 = 0; nl1 < Ncl; nl1 ++){
//     double t1 = (floor(ell[nl1+1]) + ceil(ell[nl1]))/2.; // average ell in band
//     for (nl2 = 0; nl2 < Ncl; nl2 ++){
//     double t2 = (floor(ell[nl2+1]) + ceil(ell[nl2]))/2.; // average ell in band
//       c_ng = 0.; c_g = 0.;
//       // if (test_zoverlap_cov(z1,z3)*test_zoverlap_cov(z1,z4) && NG){ c_ng = cov_NG_cl_shear_fourier_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,z3,z4,pm);}
//       if (covparams.ng){ c_ng = cov_fullsky_NG[nl1][nl2];}
//       c_g = cov_fullsky_G[nl1][nl2];
//       // fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+n1)+nl1,Ncl*(n2)+nl2,t1,t2,z1,zk,z3,z4,c_g,c_ng);
//       fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin)+nl1,Ncl*(tomo.shear_Npowerspectra+n2)+nl2,t1,t2,zk,zk,z3,z4,c_g,c_ng);

//     }
//   }
//   free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ncl-1);
//   free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ncl-1);
//   fclose(F1);
// }

// void run_cov_kk_clustering_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n2, int start)
// {
//   int z1,z2,z3,z4,nl1,nl2;
//   int zk=-1;
//   double c_ng, c_g;
//   double fsky = survey.area/41253.0;
//   FILE *F1;
//   char filename[300];
//   sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
//   F1 =fopen(filename,"w");
//   z3 = n2; z4 = n2;
//   printf("N_cl = %d (%d, %d)\n", n2,z3,z4);

//   double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
//   like.lmin = ell[0];
//   like.lmax = ell[Ncl];
//   cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
//   cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
//   cov_kk_cl_fourier_binned(cov_fullsky_G,cov_fullsky_NG,z3,z4, NG, ell);

//   for (nl1 = 0; nl1 < Ncl; nl1 ++){
//     double t1 = (floor(ell[nl1+1]) + ceil(ell[nl1]))/2.; // average ell in band
//     for (nl2 = 0; nl2 < Ncl; nl2 ++){
//     double t2 = (floor(ell[nl2+1]) + ceil(ell[nl2]))/2.; // average ell in band
//       c_ng = 0.; c_g = 0.;
//       // if (test_zoverlap_cov(z1,z3)*test_zoverlap_cov(z1,z4) && NG){ c_ng = cov_NG_cl_shear_fourier_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,z3,z4,pm);}
//       if (covparams.ng){ c_ng = cov_fullsky_NG[nl1][nl2];}
//       c_g = cov_fullsky_G[nl1][nl2];
//       // fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+n1)+nl1,Ncl*(n2)+nl2,t1,t2,z1,zk,z3,z4,c_g,c_ng);
//       fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin)+nl1,Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n2)+nl2,t1,t2,zk,zk,z3,z4,c_g,c_ng);

//     }
//   }
//   free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ncl-1);
//   free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ncl-1);
//   fclose(F1);
// }

// void run_cov_kk_gk_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n2, int start)
// {
//   int z1,z2,z3,z4,nl1,nl2;
//   int zk=-1;
//   double c_ng, c_g;
//   double fsky = survey.area/41253.0;
//   FILE *F1;
//   char filename[300];
//   sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
//   F1 =fopen(filename,"w");
//   z3 = n2;
//   printf("N_gk_2 = %d\n", n2);

//   double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
//   like.lmin = ell[0];
//   like.lmax = ell[Ncl];
//   cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
//   cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
//   cov_kk_gk_fourier_binned(cov_fullsky_G,cov_fullsky_NG,z3, NG, ell);

//   for (nl1 = 0; nl1 < Ncl; nl1 ++){
//     double t1 = (floor(ell[nl1+1]) + ceil(ell[nl1]))/2.; // average ell in band
//     for (nl2 = 0; nl2 < Ncl; nl2 ++){
//     double t2 = (floor(ell[nl2+1]) + ceil(ell[nl2]))/2.; // average ell in band
//       c_ng = 0.; c_g = 0.;
//       // if (test_zoverlap_cov(z1,z3)*test_zoverlap_cov(z1,z4) && NG){ c_ng = cov_NG_cl_shear_fourier_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,z3,z4,pm);}
//       if (covparams.ng){ c_ng = cov_fullsky_NG[nl1][nl2];}
//       c_g = cov_fullsky_G[nl1][nl2];
//       // fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+n1)+nl1,Ncl*(n2)+nl2,t1,t2,z1,zk,z3,z4,c_g,c_ng);
//       fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin)+nl1,Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+n2)+nl2,t1,t2,zk,zk,z3,zk,c_g,c_ng);
//     }
//   }
//   free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ncl-1);
//   free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ncl-1);
//   fclose(F1);
// }

// void run_cov_kk_ks_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int n2, int start)
// {
//   int z1,z2,z3,z4,nl1,nl2;
//   int zk=-1;
//   double c_ng, c_g;
//   double fsky = survey.area/41253.0;
//   FILE *F1;
//   char filename[300];
//   sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
//   F1 =fopen(filename,"w");
//   z3 = n2;
//   printf("N_ks_2 = %d\n", n2);

//   double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
//   like.lmin = ell[0];
//   like.lmax = ell[Ncl];
//   cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
//   cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
//   cov_kk_ks_fourier_binned(cov_fullsky_G,cov_fullsky_NG,z3, NG, ell);

//   for (nl1 = 0; nl1 < Ncl; nl1 ++){
//     double t1 = (floor(ell[nl1+1]) + ceil(ell[nl1]))/2.; // average ell in band
//     for (nl2 = 0; nl2 < Ncl; nl2 ++){
//     double t2 = (floor(ell[nl2+1]) + ceil(ell[nl2]))/2.; // average ell in band
//       c_ng = 0.; c_g = 0.;
//       // if (test_zoverlap_cov(z1,z3)*test_zoverlap_cov(z1,z4) && NG){ c_ng = cov_NG_cl_shear_fourier_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,z3,z4,pm);}
//       if (covparams.ng){ c_ng = cov_fullsky_NG[nl1][nl2];}
//       c_g = cov_fullsky_G[nl1][nl2];
//       // fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+n1)+nl1,Ncl*(n2)+nl2,t1,t2,z1,zk,z3,z4,c_g,c_ng);
//       fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin)+nl1,Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+n2)+nl2,t1,t2,zk,zk,z3,zk,c_g,c_ng);
//     }
//   }
//   free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ncl-1);
//   free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ncl-1);
//   fclose(F1);
// }


// // kk_kk Fourier
// void run_cov_kk_kk_fourier_bin(char *OUTFILE, char *PATH, double *ell, int Ncl, int start)
// {
//   int nl1,nl2;
//   int zk=-1;
//   double c_ng, c_g;
//   double fsky = survey.area/41253.0;
//   FILE *F1;
//   char filename[300];
//   sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
//   F1 =fopen(filename,"w");

//   double **cov_fullsky_G = 0, **cov_fullsky_NG = 0;
//   like.lmin = ell[0];
//   like.lmax = ell[Ncl];
//   cov_fullsky_G = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
//   cov_fullsky_NG = create_double_matrix(0, like.Ncl-1, 0, like.Ncl-1);
//   cov_kk_kk_fourier_binned(cov_fullsky_G,cov_fullsky_NG, NG, ell);

//   int i,j;
//   for (nl1 = 0; nl1 < Ncl; nl1 ++){
//     double t1 = (floor(ell[nl1+1]) + ceil(ell[nl1]))/2.; // average ell in band
//     for (nl2 = 0; nl2 < Ncl; nl2 ++){
//     double t2 = (floor(ell[nl2+1]) + ceil(ell[nl2]))/2.; // average ell in band
//       c_ng = 0.; c_g = 0.;
//       // if (test_zoverlap_cov(z1,z3)*test_zoverlap_cov(z1,z4) && NG){ c_ng = cov_NG_cl_shear_fourier_binned(theta[nl1],theta[nl1+1],theta[nl2],theta[nl2+1],z1,z2,z3,z4,pm);}
//       if (covparams.ng){ c_ng = cov_fullsky_NG[nl1][nl2];}
//       c_g = cov_fullsky_G[nl1][nl2];
//       i = Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin)+nl1;
//       j = Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+tomo.shear_Nbin)+nl2;
//       // fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Nbin+tomo.clustering_Nbin+n1)+nl1,Ncl*(n2)+nl2,t1,t2,z1,zk,z3,z4,c_g,c_ng);
//       fprintf(F1,"%d %d %e %e %d %d %d %d %e %e\n", i,j,t1,t2,zk,zk,zk,zk,c_g,c_ng);
//     }
//   }
//   free_double_matrix(cov_fullsky_G,0, like.Ncl-1, 0, like.Ncl-1);
//   free_double_matrix(cov_fullsky_NG,0, like.Ncl-1, 0, like.Ncl-1);
//   fclose(F1);
// }






