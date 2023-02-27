#include <math.h>
#include <stdlib.h>
#if !defined(__APPLE__)
#include <malloc.h>
#endif
#include <stdio.h>
#include <assert.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_sf_erf.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf_expint.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include <fftw3.h>

#include "../cfftlog/cfftlog.h"
#include "../cfftlog/utils.h"
#include "../cfastpt/cfastpt.h"
#include "../cfastpt/utils.h"


#include "basics.c"
#include "structs.c"
#include "parameters.c"
#include "../emu17/P_cb/emu.c"
#include "recompute.c"
#include "cosmo3D.c"
#include "redshift_spline.c"
#include "halo.c"
#include "HOD.c"
#include "pt_cfastpt.c"
#include "cosmo2D_fourier.c"
#include "IA.c"
#include "cluster.c"
#include "BAO.c"
//#include "external_prior.c"
#include "covariances_3D.c"
#include "covariances_fourier.c"
#include "covariances_cluster.c"
#include "../../des_y3/init_y3_fourier.c"
//#include "init.c"
#include "../cfftlog/utils_complex.h"
#include "../cfastpt/utils_complex.h"

void run_cov_N_N (char *OUTFILE, char *PATH, int nzc1, int nzc2,int start);
void run_cov_cgl_N (char *OUTFILE, char *PATH, double *ell_Cluster, double *dell_Cluster,int N1, int nzc2, int start);
void run_cov_cgl_cgl (char *OUTFILE, char *PATH, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int start);
void run_cov_cgl_cgl_all (char *OUTFILE, char *PATH, double *ell_Cluster, double *dell_Cluster);
void run_cov_shear_N (char *OUTFILE, char *PATH, double *ell, double *dell, int N1, int nzc2, int start);
void run_cov_shear_cgl (char *OUTFILE, char *PATH, double *ell, double *dell, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int nl1, int start);
void run_cov_ggl_N (char *OUTFILE, char *PATH, double *ell, double *dell, int N1, int nzc2, int start);
void run_cov_ggl_cgl (char *OUTFILE, char *PATH, double *ell, double *dell, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int nl1, int start);
void run_cov_cl_N (char *OUTFILE, char *PATH, double *ell, double *dell, int N1, int nzc2, int start);
void run_cov_cl_cgl (char *OUTFILE, char *PATH, double *ell, double *dell, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int nl1, int start);

void run_cov_ggl_shear(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);
void run_cov_clustering_shear(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);
void run_cov_clustering_ggl(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);
void run_cov_clustering(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);
void run_cov_ggl(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);
void run_cov_shear_shear(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start);


void run_cov_N_N (char *OUTFILE, char *PATH, int nzc1, int nzc2,int start)
{
  int nN1, nN2,i,j;
  double cov;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  for (nN1 = 0; nN1 < Cluster.N200_Nbin; nN1 ++){
    for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
      i = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
      i += Cluster.N200_Nbin*nzc1+nN1;
      j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
      j += Cluster.N200_Nbin*nzc2+nN2;

      cov =cov_N_N(nzc1,nN1, nzc2, nN2);
      fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j,0.0,0.0, nzc1, nN1, nzc2, nN2,cov,0.0);
    }
  }
  fclose(F1);
}

void run_cov_cgl_N (char *OUTFILE, char *PATH, double *ell_Cluster, double *dell_Cluster,int N1, int nzc2, int start)
{
  int nN1, nN2, nl1, nzc1, nzs1,i,j;
  double cov;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  nzc1 = ZC(N1);
  nzs1 = ZSC(N1);
  for (nN1 = 0; nN1 < Cluster.N200_Nbin; nN1 ++){
    for( nl1 = 0; nl1 < Cluster.lbin; nl1 ++){
     for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
       i = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
       i += (N1*Cluster.N200_Nbin+nN1)*Cluster.lbin +nl1;
       j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
       j += Cluster.N200_Nbin*nzc2+nN2;

       cov =cov_cgl_N(ell_Cluster[nl1],nzc1,nN1, nzs1, nzc2, nN2);
       fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j, ell_Cluster[nl1], 0., nzc1, nzs1, nzc2, nN2,cov,0.);
     }
    }
  }
  fclose(F1);
}

void run_cov_cgl_cgl (char *OUTFILE, char *PATH, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int start)
{
  int nN1, nN2, nl1, nzc1, nzs1, nl2, nzc2, nzs2,i,j;
  double c_g, c_ng;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  nzc1 = ZC(N1);
  nzs1 = ZSC(N1);
  nzc2 = ZC(N2);
  nzs2 = ZSC(N2);
  for (nN1 = 0; nN1 < Cluster.N200_Nbin; nN1 ++){
    for( nl1 = 0; nl1 < Cluster.lbin; nl1 ++){
      for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
        for( nl2 = 0; nl2 < Cluster.lbin; nl2 ++){
          i = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
          i += (N1*Cluster.N200_Nbin+nN1)*Cluster.lbin +nl1;
          j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
          j += (N2*Cluster.N200_Nbin+nN2)*Cluster.lbin +nl2;

          c_g = 0;
          c_ng = cov_NG_cgl_cgl(ell_Cluster[nl1],ell_Cluster[nl2],nzc1,nN1, nzs1, nzc2, nN2,nzs2);
          if (nl2 == nl1){c_g =cov_G_cgl_cgl(ell_Cluster[nl1],dell_Cluster[nl1],nzc1,nN1, nzs1, nzc2, nN2,nzs2);}
          fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j,ell_Cluster[nl1],ell_Cluster[nl2], nzc1, nzs1, nzc2, nzs2,c_g, c_ng);
        }
      }
    }
  }
  fclose(F1);
}

void run_cov_cgl_cgl_all (char *OUTFILE, char *PATH, double *ell_Cluster, double *dell_Cluster)
{
  int nN1, nN2, nl1, nzc1, nzs1, nl2, nzc2, nzs2,i,j, N1,N2;
  double c_g, c_ng;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s",PATH,OUTFILE);
  F1 =fopen(filename,"w");
  for (N1 = 0; N1 < tomo.cgl_Npowerspectra; N1 ++){
    for (N2 = 0; N2 < tomo.cgl_Npowerspectra; N2 ++){
      printf("%d %d\n",N1, N2);
      nzc1 = ZC(N1);
      nzs1 = ZSC(N1);
      nzc2 = ZC(N2);
      nzs2 = ZSC(N2);
      for (nN1 = 0; nN1 < Cluster.N200_Nbin; nN1 ++){
        for( nl1 = 0; nl1 < Cluster.lbin; nl1 ++){
          for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
            for( nl2 = 0; nl2 < Cluster.lbin; nl2 ++){
              i = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
              i += (N1*Cluster.N200_Nbin+nN1)*Cluster.lbin +nl1;
              j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
              j += (N2*Cluster.N200_Nbin+nN2)*Cluster.lbin +nl2;
          
              c_g = 0;
              c_ng = cov_NG_cgl_cgl(ell_Cluster[nl1],ell_Cluster[nl2],nzc1,nN1, nzs1, nzc2, nN2,nzs2);
              if (nl2 == nl1){c_g =cov_G_cgl_cgl(ell_Cluster[nl1],dell_Cluster[nl1],nzc1,nN1, nzs1, nzc2, nN2,nzs2);}
              fprintf(F1,"%d %d %e %e %d %d %d  %d %d %d  %e %e\n",i,j,ell_Cluster[nl1],ell_Cluster[nl2], nzc1, nN1,nzs1, nzc2, nN2,nzs2,c_g, c_ng);
            }
          }
        }
      }
    }
  }
  fclose(F1);
}

void run_cov_shear_N (char *OUTFILE, char *PATH, double *ell, double *dell, int N1, int nzc2, int start)
{
  int nz1,nz2, nN2, nl1, nzc1, i,j;
  double cov;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  nz1 = Z1(N1);
  nz2 = Z2(N1);
  for( nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
      cov = 0.;
      i = like.Ncl*N1+nl1;
      j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
      j += Cluster.N200_Nbin*nzc2+nN2;
        
      if (ell[nl1] < like.lmax_shear){cov =cov_shear_N(ell[nl1],nz1,nz2, nzc2, nN2);}
      fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j, ell[nl1], 0., nz1, nz2, nzc2, nN2,cov,0.);
    }
  }
  fclose(F1);
}

void run_cov_shear_cgl (char *OUTFILE, char *PATH, double *ell, double *dell, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int nl1, int start)
{
  int nN1, nN2, nzs1, nzs2,nl2, nzc2, nzs3,i,j;
  double c_g, c_ng;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  nzs1 = Z1(N1);
  nzs2 = Z2(N1);
  nzc2 = ZC(N2);
  nzs3 = ZSC(N2);
  for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
    for( nl2 = 0; nl2 < Cluster.lbin; nl2 ++){
      i =  like.Ncl*N1+nl1;
      j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
      j += (N2*Cluster.N200_Nbin+nN2)*Cluster.lbin +nl2;          
      c_g = 0.;
      c_ng = 0.;
      if (ell[nl1] < like.lmax_shear){
        c_ng = cov_NG_shear_cgl(ell[nl1],ell_Cluster[nl2],nzs1, nzs2, nzc2, nN2,nzs3);
        if (fabs(ell[nl1]/ell_Cluster[nl2] -1.) < 0.1){ c_g =cov_G_shear_cgl(ell[nl1],dell_Cluster[nl2],nzs1,nzs2, nzc2, nN2,nzs3);}
      }
      fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j,ell[nl1],ell_Cluster[nl2], nzs1, nzs2, nzc2, nzs3,c_g, c_ng);
    }
  }
  fclose(F1);
}

void run_cov_ggl_N (char *OUTFILE, char *PATH, double *ell, double *dell, int N1, int nzc2, int start)
{
  int zl,zs, nN2, nl1, nzc1, i,j;
  double cov,weight;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  zl = ZL(N1);
  zs = ZS(N1);
  for( nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
      i = like.Ncl*(tomo.shear_Npowerspectra+N1)+nl1;
      j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
      j += Cluster.N200_Nbin*nzc2+nN2;
      cov = 0.;
      weight = test_kmax(ell[nl1],zl);
      if (weight){
        cov =cov_ggl_N(ell[nl1],zl,zs, nzc2, nN2);
      }
      fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j, ell[nl1], 0., zl, zs, nzc2, nN2,cov,0.);
    }
  }
  fclose(F1);
}

void run_cov_ggl_cgl (char *OUTFILE, char *PATH, double *ell, double *dell, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int nl1, int start)
{
  int nN2, zl, zs, nzs1, nl2, nzc2, nzs3,i,j;
  double c_g, c_ng,weight;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  zl = ZL(N1);
  zs = ZS(N1);
  nzc2 = ZC(N2);
  nzs3 = ZSC(N2);
  for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
    for( nl2 = 0; nl2 < Cluster.lbin; nl2 ++){
      i =   like.Ncl*(tomo.shear_Npowerspectra+N1)+nl1;
      j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
      j += (N2*Cluster.N200_Nbin+nN2)*Cluster.lbin +nl2;
        
      c_g = 0; c_ng = 0.;
      weight = test_kmax(ell[nl1],zl);
      if (weight){
        c_ng = cov_NG_ggl_cgl(ell[nl1],ell_Cluster[nl2],zl,zs, nzc2, nN2,nzs3);
        if (fabs(ell[nl1]/ell_Cluster[nl2] -1.) < 0.1){
          c_g =cov_G_ggl_cgl(ell[nl1],dell_Cluster[nl2],zl,zs, nzc2, nN2,nzs3);
        }
      }
      fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j,ell[nl1],ell_Cluster[nl2], zl, zs, nzc2, nzs3,c_g, c_ng);
    }
  }
  fclose(F1);
}

void run_cov_cl_N (char *OUTFILE, char *PATH, double *ell, double *dell,int N1, int nzc2, int start)
{
  int zl,zs, nN2, nl1, nzc1, i,j;
  double cov,weight;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  for( nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
      i = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+N1)+nl1;
      j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra);
      j += Cluster.N200_Nbin*nzc2+nN2;
      cov = 0.;
      weight = test_kmax(ell[nl1],N1);
      if (weight){
        cov =cov_cl_N(ell[nl1],N1,N1, nzc2, nN2);
      }
      fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j, ell[nl1], 0., N1, N1, nzc2, nN2,cov,0.);
    }
  }
  fclose(F1);
}

void run_cov_cl_cgl (char *OUTFILE, char *PATH, double *ell, double *dell, double *ell_Cluster, double *dell_Cluster,int N1, int N2, int nl1, int start)
{
  int nN2,nzc2, nzs3,i,j,nl2;
  double c_g, c_ng,weight;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  nzc2 = ZC(N2);
  nzs3 = ZSC(N2);
    for (nN2 = 0; nN2 < Cluster.N200_Nbin; nN2 ++){
      for( nl2 = 0; nl2 < Cluster.lbin; nl2 ++){
        i = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+N1)+nl1;
        j = like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+tomo.clustering_Npowerspectra)+Cluster.N200_Nbin*tomo.cluster_Nbin;
        j += (N2*Cluster.N200_Nbin+nN2)*Cluster.lbin +nl2;
        
        c_g = 0; c_ng = 0.;
        weight = test_kmax(ell[nl1],N1);
        if (weight){
          c_ng = cov_NG_cl_cgl(ell[nl1],ell_Cluster[nl2],N1,N1, nzc2, nN2,nzs3);
          if (fabs(ell[nl1]/ell_Cluster[nl2] -1.) < 0.1){c_g =cov_G_cl_cgl(ell[nl1],dell_Cluster[nl2],N1,N1, nzc2, nN2,nzs3);}
        }
        fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n",i,j,ell[nl1],ell_Cluster[nl2], N1,N1, nzc2, nzs3,c_g, c_ng);
        //printf("%d %d %e %e %d %d %d %d  %e %e\n",i,j,ell[nl1],ell_Cluster[nl2], N1,N1, nzc2, nzs3,c_g, c_ng);
      }
    }
  fclose(F1);
}

void run_cov_ggl_shear(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
  int zl,zs,z3,z4,nl1,nl2,weight;
  double c_ng, c_g;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w"); 
  zl = ZL(n1); zs = ZS(n1);
  printf("\nN_ggl = %d (%d, %d)\n", n1,zl,zs);
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear = %d (%d, %d)\n", n2,z3,z4);
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 <like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      weight = test_kmax(ell[nl1],zl);
      if (weight && ell[nl2] < like.lmax_shear){
        if (test_zoverlap(zl,z3)*test_zoverlap(zl,z4)){c_ng = cov_NG_gl_shear_tomo(ell[nl1],ell[nl2],zl,zs,z3,z4);
        }
        if (nl1 == nl2){
          c_g =  cov_G_gl_shear_tomo(ell[nl1],dell[nl1],zl,zs,z3,z4);
        }
      }
      fprintf(F1, "%d %d %e %e %d %d %d %d  %e %e\n", like.Ncl*(tomo.shear_Npowerspectra+n1)+nl1,like.Ncl*(n2)+nl2, ell[nl1],ell[nl2],zl,zs,z3,z4,c_g,c_ng);
    }
  }
  fclose(F1);
}

void run_cov_clustering_shear(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
  int z1,z2,z3,z4,nl1,nl2,weight;
  double c_ng, c_g;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1; z2 = n1;
  printf("\nN_cl = %d \n", n1);
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear = %d (%d, %d)\n", n2,z3,z4);
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 <like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      weight = test_kmax(ell[nl1],z1);
      if (weight && ell[nl2] < like.lmax_shear){
        if (test_zoverlap(z1,z3)*test_zoverlap(z1,z4)){
        c_ng = cov_NG_cl_shear_tomo(ell[nl1],ell[nl2],z1,z2,z3,z4);
      }
      if (nl1 == nl2){
        c_g =  cov_G_cl_shear_tomo(ell[nl1],dell[nl1],z1,z2,z3,z4);
      }
    }
    fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n", like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,like.Ncl*(n2)+nl2, ell[nl1],ell[nl2],z1,z2,z3,z4,c_g,c_ng);
    }
  }
  fclose(F1);
}

void run_cov_clustering_ggl(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
  int z1,z2,zl,zs,nl1,nl2,weight;
  double c_ng, c_g;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1; z2 = n1;
  printf("\nN_cl_1 = %d \n", n1);
  zl = ZL(n2); zs = ZS(n2);
  printf("N_tomo_2 = %d (%d, %d)\n", n2,zl,zs);
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 < like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      if (z1 == zl){
        weight = test_kmax(ell[nl1],z1)*test_kmax(ell[nl2],zl);
        if (weight){
          c_ng = cov_NG_cl_gl_tomo(ell[nl1],ell[nl2],z1,z2,zl,zs);
          if (nl1 == nl2){
            c_g =  cov_G_cl_gl_tomo(ell[nl1],dell[nl1],z1,z2,zl,zs);
          }
        }
      }
      fprintf(F1, "%d %d %e %e %d %d %d %d  %e %e\n", like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,like.Ncl*(tomo.shear_Npowerspectra+n2)+nl2, ell[nl1],ell[nl2],z1,z2,zl,zs,c_g,c_ng);
    }
  }
  fclose(F1);
}

void run_cov_clustering(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
  int z1,z2,z3,z4,nl1,nl2,weight;
  double c_ng, c_g;
  FILE *F1;
  char filename[300];
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  F1 =fopen(filename,"w");
  z1 = n1; z2 = n1;
  printf("\nN_cl_1 = %d \n", n1);
  z3 = n2; z4 = n2;
  printf("N_cl_2 = %d\n", n2);
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 < like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      if (z1 == z3){
        weight = test_kmax(ell[nl1],z1)*test_kmax(ell[nl2],z3);
        if (weight) {
          c_ng = cov_NG_cl_cl_tomo(ell[nl1],ell[nl2],z1,z2,z3,z4);
        }
        if (nl1 == nl2){
          c_g =  cov_G_cl_cl_tomo(ell[nl1],dell[nl1],z1,z2,z3,z4);
        }
      }
      fprintf(F1, "%d %d %e %e %d %d %d %d %e %e\n", like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra+n1)+nl1,like.Ncl*(tomo.shear_Npowerspectra+tomo.ggl_Npowerspectra + n2)+nl2, ell[nl1],ell[nl2],z1,z2,z3,z4,c_g,c_ng);
    }
  }
  fclose(F1);
}

void run_cov_ggl(char *OUTFILE, char *PATH, double *ell, double *dell, int n1, int n2,int start)
{
  int zl1,zl2,zs1,zs2,nl1,nl2, weight;
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
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 < like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      weight = test_kmax(ell[nl1],zl1)*test_kmax(ell[nl2],zl2);
      if (weight && zl1 == zl2) {
        c_ng = cov_NG_gl_gl_tomo(ell[nl1],ell[nl2],zl1,zs1,zl2,zs2);
      }
      if (nl1 == nl2){
        c_g =  cov_G_gl_gl_tomo(ell[nl1],dell[nl1],zl1,zs1,zl2,zs2);
      }
      if (weight ==0 && n2 != n1){
        c_g = 0;
      }
      fprintf(F1,"%d %d %e %e %d %d %d %d  %e %e\n", like.Ncl*(tomo.shear_Npowerspectra+n1)+nl1,like.Ncl*(tomo.shear_Npowerspectra+n2)+nl2, ell[nl1],ell[nl2],zl1,zs1,zl2,zs2,c_g,c_ng);   
    }
  }
  fclose(F1);
}


void run_cov_shear_shear(char *OUTFILE, char *PATH, double *ell, double *dell,int n1, int n2,int start)
{
  int z1,z2,z3,z4,nl1,nl2,weight;
  double c_ng, c_g;
  FILE *F1;
  char filename[300];
  z1 = Z1(n1); z2 = Z2(n1);
  printf("N_shear = %d\n", n1);
  z3 = Z1(n2); z4 = Z2(n2);
  printf("N_shear = %d (%d, %d)\n",n2,z3,z4);
  sprintf(filename,"%s%s_%d",PATH,OUTFILE,start);
  printf("%s\n",filename);
  F1 =fopen(filename,"w");
  for (nl1 = 0; nl1 < like.Ncl; nl1 ++){
    for (nl2 = 0; nl2 < like.Ncl; nl2 ++){
      c_ng = 0.; c_g = 0.;
      if (ell[nl1] < like.lmax_shear && ell[nl2] < like.lmax_shear){c_ng = cov_NG_shear_shear_tomo(ell[nl1],ell[nl2],z1,z2,z3,z4);}
      if (nl1 == nl2){
        c_g =  cov_G_shear_shear_tomo(ell[nl1],dell[nl1],z1,z2,z3,z4);
        if (ell[nl1] > like.lmax_shear && n1!=n2){c_g = 0.;} 
      }         
      fprintf(F1, "%d %d %e %e %d %d %d %d %e %e\n", like.Ncl*n1+nl1,like.Ncl*(n2)+nl2, ell[nl1],ell[nl2],z1,z2,z3,z4,c_g,c_ng);
      //printf("%d %d %e %e %d %d %d %d %e %e\n", like.Ncl*n1+nl1,like.Ncl*(n2)+nl2, ell[nl1],ell[nl2],z1,z2,z3,z4,c_g,c_ng);
    }
  }
  fclose(F1);
}


int main(int argc, char** argv)
{
  int i,l,m,n,o,s,p,nl1;
  int hit=1;//atoi(argv[1]);
  char OUTFILE[400],PATH[400];
    nuisance.A_ia = 0.7;
  //nuisance.A2_ia = -1.36;
  //nuisance.b_ta_z[0] = 1.0;

  nuisance.eta_ia = -1.7;
  double b1[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}, b2[10] ={0.0,0.0,0.0,0.0,0.,0.,0.,0.,0.,0.},b_mag[10] ={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}, stretch[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  b1[0] = 1.7;
  b1[1] = 1.7;
  b1[2] = 1.7;
  b1[3] = 2.0;
  b1[4] = 2.0;
  survey.area   = 12300.0;
  survey.n_gal   = 10.7;
  survey.sigma_e   = 0.26;  
  survey.m_lim=24.1;
  survey.n_lens=13.1;
  //RUN MODE setup
  double mean_m[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double sigma_m[10]={0.005,0.005,0.005,0.005,0.005,0.0,0.0,0.0,0.0,0.0};
  double bias_photoz_s[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double sigma_b_photoz_s[10]={0.05,0.05,0.05,0.05,0.05,0.0,0.0,0.0,0.0,0.0};
  double bias_photoz_l[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
  double sigma_b_photoz_l[10]={0.03,0.03,0.03,0.03,0.03,0.0,0.0,0.0,0.0,0.0};
  double mean_stretch[10]={1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
  double sigma_stretch[10]={0.03,0.03,0.03,0.03,0.003,0.0,0.0,0.0,0.0,0.0};
  
  double pm[10]={0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
//set_stretch_priors_mpp(mean_stretch, sigma_stretch);

  double A_1[10]={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}, A_2[10] ={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.},b_ta[10] ={0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  A_1[0]= nuisance.A_ia; A_1[1] = nuisance.eta_ia;
  A_2[0]= nuisance.A2_ia; A_2[1] = nuisance.eta_ia_tt;
  b_ta[0] = nuisance.b_ta_z[0];

  init_IA_mpp(6);
  init_cosmo_runmode("CLASS");
  init_survey_mpp("LSST Y1", 12300.0, 0.26);
  
  init_source_sample_mpp("cov_LSSTY1/source_LSSTY1.nz",5);
  init_lens_sample_mpp("cov_LSSTY1/lens_LSSTY1.nz",5,mean_stretch, bias_photoz_l, b1,b2,10.0, 0, 1);
  //init_source_sample_mpp("../DESC_SRD/zdistris/zdistri_model_z0=1.300000e-01_beta=7.800000e-01_Y1_source",5);
  //init_lens_sample_mpp("../DESC_SRD/zdistris/zdistri_model_z0=2.600000e-01_beta=9.400000e-01_Y1_lens",5,stretch, b_mag, b1,b2,10.0, 0.0, 1.0);

  init_binning_fourier_mpp(20,20,15000,15000, 0.0001);
  init_probes_fourier_mpp("3x2pt");
  tomo.n_source[0] = survey.n_gal;
  tomo.n_source[1] = survey.n_gal;
  tomo.n_source[2] = survey.n_gal;
  tomo.n_source[3] = survey.n_gal;
  tomo.n_source[4] = survey.n_gal;

  tomo.n_lens[0] = survey.n_lens;
  tomo.n_lens[1] = survey.n_lens;
  tomo.n_lens[2] = survey.n_lens;
  tomo.n_lens[3] = survey.n_lens;
  tomo.n_lens[4] = survey.n_lens;

  //init_binning_fourier_mpp(20,20,15000,15000, 0.00001);
  //init_galaxies("../zdistris/zdistribution_Euclid","../zdistris/zdistribution_const_comoving", "none", "gaussian", "redmagic");
  //init_clusters();

  //init_probes("all_2pt");
  //sprintf(PATH,"/home/teifler/Dropbox/cosmolike/top-level/WFIRST_SIT/covparallel/");  
  

  sprintf(PATH,"/home/paul/des_y3/cov_LSSTY1/nz/");
  //printf("HIHIHIHIHIHIHIH\n");
  int k;

  //set l-bins for shear, ggl, clustering
  double logdl=(log(like.lmax)-log(like.lmin))/like.Ncl;
  double *ell, *dell;
  ell=create_double_vector(0,like.Ncl-1);
  dell=create_double_vector(0,like.Ncl-1);
  for(i=0;i<like.Ncl;i++){
    ell[i]=exp(log(like.lmin)+(i+0.5)*logdl);
    dell[i]=exp(log(like.lmin)+(i+1)*logdl) - exp(log(like.lmin)+(i*logdl));
  } 
  //set l-bins for cluster lensing
  double *ell_Cluster, *dell_Cluster;
  logdl=(log(Cluster.l_max)-log(Cluster.l_min))/Cluster.lbin;
  ell_Cluster=create_double_vector(0,Cluster.lbin-1);
  dell_Cluster=create_double_vector(0,Cluster.lbin-1);
  for(i=0;i<Cluster.lbin;i++){
    ell_Cluster[i]=exp(log(Cluster.l_min)+(i+0.5)*logdl);
    dell_Cluster[i]=exp(log(Cluster.l_min)+(i+1)*logdl) - exp(log(Cluster.l_min)+(i*logdl));
  }

  k=1;
/*
  sprintf(OUTFILE,"%s_ssss_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.shear_Npowerspectra; l++){
    for (m=l;m<tomo.shear_Npowerspectra; m++){
      //if(k==hit) 
      run_cov_shear_shear(OUTFILE,PATH,ell,dell,l,m,k);
      k=k+1;
      printf("%d\n",k);
    }
  }

//k=1;
  sprintf(OUTFILE,"%s_lsls_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.ggl_Npowerspectra; l++){
    for (m=l;m<tomo.ggl_Npowerspectra; m++){
      //if(k==hit) 
      run_cov_ggl(OUTFILE,PATH,ell,dell,l,m,k);
     //printf("%d\n",k);
      k=k+1;
    }
  }
    
  //k=1;
  sprintf(OUTFILE,"%s_llll_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.clustering_Npowerspectra; l++){ //auto bins only for now!
    for (m=l;m<tomo.clustering_Npowerspectra; m++){
      //if(k==hit) 
      run_cov_clustering(OUTFILE,PATH,ell,dell,l,m,k);
      k=k+1;
      printf("%d %d %d\n",l,m,k);
    }
  }

  sprintf(OUTFILE,"%s_llss_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.clustering_Npowerspectra; l++){
    for (m=0;m<tomo.shear_Npowerspectra; m++){
      //if(k==hit) 
      run_cov_clustering_shear(OUTFILE,PATH,ell,dell,l,m,k);
      k=k+1;
      //printf("%d\n",k);
    }
  }
  sprintf(OUTFILE,"%s_llls_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.clustering_Npowerspectra; l++){
    for (m=0;m<tomo.ggl_Npowerspectra; m++){
      //if(k==hit) 
      run_cov_clustering_ggl(OUTFILE,PATH,ell,dell,l,m,k);
      k=k+1;
      //printf("%d\n",k);
    }
  }
  sprintf(OUTFILE,"%s_lsss_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
  for (l=0;l<tomo.ggl_Npowerspectra; l++){
    for (m=0;m<tomo.shear_Npowerspectra; m++){
      //if(k==hit) 
      run_cov_ggl_shear(OUTFILE,PATH,ell,dell,l,m,k);
      k=k+1;
      //printf("%d\n",k);
    }
  }


*/
// /********** cluster covariance ************/
//   sprintf(OUTFILE,"%s_nn_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
//   for (l=0;l<tomo.cluster_Nbin; l++){
//     for (m=0;m<tomo.cluster_Nbin; m++){
//       if(k==hit) run_cov_N_N (OUTFILE,PATH,l,m,k);
//       k=k+1;
//       //printf("%d\n",k);
//     }
//   }
//   sprintf(OUTFILE,"%s_cscs_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
//   for (l=0;l<tomo.cgl_Npowerspectra; l++){
//     for (m=0;m<tomo.cgl_Npowerspectra; m++){
//       if(k==hit) run_cov_cgl_cgl (OUTFILE,PATH,ell_Cluster,dell_Cluster,l,m,k);
//       k=k+1;
//     }
//   }
//   sprintf(OUTFILE,"%s_csn_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
//   for (l=0;l<tomo.cgl_Npowerspectra; l++){
//     for (m=0;m<tomo.cluster_Nbin; m++){
//       if(k==hit) run_cov_cgl_N (OUTFILE,PATH,ell_Cluster,dell_Cluster,l,m,k);
//       k=k+1;
//     }
//   }
// // shear X cluster
//   sprintf(OUTFILE,"%s_ssn_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
//   for (l=0;l<tomo.shear_Npowerspectra; l++){
//     for (m=0;m<tomo.cluster_Nbin; m++){
//       if(k==hit) run_cov_shear_N (OUTFILE,PATH,ell,dell,l,m,k);
//       k=k+1;
//     }
//   }  
//   sprintf(OUTFILE,"%s_sscs_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
//   for (l=0;l<tomo.shear_Npowerspectra; l++){
//     for (m=0;m<tomo.cgl_Npowerspectra; m++){
//       for(nl1 = 0; nl1 < like.Ncl; nl1 ++){
//         if(k==hit) run_cov_shear_cgl (OUTFILE,PATH,ell,dell,ell_Cluster,dell_Cluster,l,m,nl1,k);
//         k=k+1;
//       }
//     }
//   }
//   // ggl X cluster
//   sprintf(OUTFILE,"%s_lsn_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
//   for (l=0;l<tomo.ggl_Npowerspectra; l++){
//     for (m=0;m<tomo.cluster_Nbin; m++){
//       if(k==hit) run_cov_ggl_N (OUTFILE,PATH,ell,dell,l,m,k);
//       k=k+1;
//     }
//   }
//   sprintf(OUTFILE,"%s_lscs_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
//   for (l=0;l<tomo.ggl_Npowerspectra; l++){
//     for (m=0;m<tomo.cgl_Npowerspectra; m++){
//       for(nl1 = 0; nl1 < like.Ncl; nl1 ++){
//         if(k==hit) run_cov_ggl_cgl (OUTFILE,PATH,ell,dell,ell_Cluster,dell_Cluster,l,m,nl1,k);
//         k=k+1;
//       }
//     }
//   }
//   // clustering X cluster
//   sprintf(OUTFILE,"%s_lln_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
//   for (l=0;l<tomo.clustering_Npowerspectra; l++){
//     for (m=0;m<tomo.cluster_Nbin; m++){
//       if(k==hit) run_cov_cl_N (OUTFILE,PATH,ell,dell,l,m,k);
// //      printf("%d\n",k);
//       k=k+1;
//     }
//   }
//   sprintf(OUTFILE,"%s_llcs_cov_Ncl%d_Ntomo%d",survey.name,like.Ncl,tomo.shear_Nbin);
//   for (l=0;l<tomo.clustering_Npowerspectra; l++){
//     for (m=0;m<tomo.cgl_Npowerspectra; m++){
//       for(nl1 = 0; nl1 < like.Ncl; nl1 ++){
//         if(k==hit) run_cov_cl_cgl (OUTFILE,PATH,ell,dell,ell_Cluster,dell_Cluster,l,m,nl1,k);
//         k=k+1;
//       }
//     }
//   }
  printf("number of cov blocks for parallelization: %d\n",k-1); 
  
  printf("-----------------\n");
  printf("PROGRAM EXECUTED\n");
  printf("-----------------\n");
  return 0;   
}

