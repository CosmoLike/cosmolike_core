//================================================================================================
// Covariances of observables involving CMB lensing
// no look up tables for now

// Naming convention:
// g = galaxy positions ("g" as in "galaxy")
// k = kappa CMB ("k" as in "kappa")
// s = kappa from source galaxies ("s" as in "shear")
// And alphabetical order

double cov_NG_gg_ks(double l1,double l2, int n1, int n2, int zs);
double cov_NG_gs_ks(double l1,double l2, int zl, int zs1, int zs2);
double cov_G_ks_ks(double l, double delta_l, int z1, int z3);
double cov_NG_ks_ss(double l1,double l2, int z1, int z2, int z3);

double cov_NG_gk_ks(double l1,double l2, int zl, int zs);
double cov_NG_gg_gk(double l1, double l2, int n1, int n2, int n3);
double cov_NG_gk_gk(double l1,double l2, int n1, int n2);
double cov_NG_gk_gs(double l1,double l2, int zl1, int zl2, int zs);
double cov_NG_gk_ss(double l1,double l2, int zl, int zs1, int zs2);
//================================================================================================
// utility routines

// Reads in the noise N_mv for mv quadratic estimator for d,
// without reduction due to average over ell-bin.
// Returns value at nearest upper ell

// Assuming planck fsky is much larger and covers the galaxy survey area
// This is going to affect all covariances of C_kk and C_xy
double fsky_planck = 0.673;
double area_planck = 27763.269;

double kappa_reconstruction_noise(double l){
   
   static double *ell;
   static double *noise;
   static int nEll;

   printf("run %s\n", cmb.name);
   if (noise==0){
      // count lines
      nEll = line_count(cmb.pathLensRecNoise)-1;
      printf("Reading CMB lensing noise: %s\n", cmb.pathLensRecNoise);

      // allocate ell and Nlkk
      ell = create_double_vector(0, nEll-1);
      noise = create_double_vector(0, nEll-1);
      
      // read each line
      FILE *file = fopen(cmb.pathLensRecNoise, "r");
      int iEll;
      for (iEll=0; iEll<nEll; iEll++) {
         fscanf(file, "%le %le", &ell[iEll], &noise[iEll]);
      }
      fclose(file);
   }
   // if l is in the range
   if ((l>=ell[0]) &&(l<=ell[nEll-1])){
      // find value of ell just above l
      int iEll = 0;
      while (ell[iEll] < l) {
         iEll ++;
      }
      // evaluate at that ell
      // C_ell^kk = l*(l+1)/4 * C_ell^dd
      return noise[iEll];
   }
   return 0.;
}


//================================================================================================
// Covariances

/********** gg gk ***************/

double cov_G_gg_gk(double l, double delta_l, int n1, int n2, int n3){
   double C13, C14, C23, C24, N13=0, N14=0, N23=0, N24=0;
   double fsky = survey.area/41253.0;
   C13 = C_cl_tomo_nointerp(l,n1,n3);
   C24 = C_gk_nointerp(l,n2);
   C14 = C_gk_nointerp(l,n1);
   C23 = C_cl_tomo_nointerp(l,n2,n3);
   if (n1==n3){
      N13 = 1./(nlens(n1)*survey.n_gal_conversion_factor);
   }
   if (n2==n3){
      N23 = 1./(nlens(n2)*survey.n_gal_conversion_factor);
   }
   return ((C13+N13)*C24+C14*(C23+N23))/((2.*l+1.)*delta_l*fsky);
}

double inner_project_tri_cov_gg_gk(double a,void *params)
{
   double k1,k2,fK,weights,res = 0.;
   double *ar = (double *) params;
   fK = f_K(chi(a));
   k1 = (ar[0]+0.5)/fK;
   k2 = (ar[1]+0.5)/fK;
   weights = W_gal(a,ar[2])*W_gal(a,ar[3])*W_gal(a,ar[4])*W_k(a,fK)*dchi_da(a);
   if (weights >0.){
      res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);
      res += delP_SSC(k1,a)*delP_SSC(k2,a)*survey_variance(a,survey.area/41253.0)*pow(fK,-4.); //super-sample covariance
   }
   res *= weights;
   return res;
}

double cov_NG_gg_gk(double l1, double l2, int n1, int n2, int n3){
   double a1,a2,array[5];
   int zmax, zmin;
   zmax = (n2 > n1 ? n2 : n1);
   a1 = amin_lens(zmax);
   zmin = (n2 < n1 ? n2 : n1);
   a2 = amax_lens(zmin);
   array[0] = l1;
   array[1] = l2;
   array[2] = (double) n1;
   array[3] = (double) n2;
   array[4] = (double) n3;
   return int_gsl_integrate_low_precision(inner_project_tri_cov_gg_gk,(void*)array,a1,a2,NULL,1000);
}

/********** gg kk ***************/

double cov_G_gg_kk(double l, double delta_l, int n1, int n2){
   double C13, C14, C23, C24;
   // double fsky = survey.area/41253.0;
   C13 = C_gk_nointerp(l,n1);
   C24 = C_gk_nointerp(l,n2);
   C14 = C13;
   C23 = C24;
   return (C13*C24+C14*C23)/((2.*l+1.)*delta_l*fsky_planck);
}

double inner_project_tri_cov_gg_kk(double a,void *params)
{
   double k1,k2,fK,weights,res = 0.;
   double *ar = (double *) params;
   fK = f_K(chi(a));
   k1 = (ar[0]+0.5)/fK;
   k2 = (ar[1]+0.5)/fK;
   weights = W_gal(a,ar[2])*W_gal(a,ar[3])*W_k(a,fK)*W_k(a,fK)*dchi_da(a);
   if (weights >0.){
      res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(area_planck*survey.area_conversion_factor);
      res += delP_SSC(k1,a)*delP_SSC(k2,a)*cross_survey_variance(a,survey.area/41253.0,fsky_planck)*pow(fK,-4.); //super-sample covariance
   }
   res *= weights;
   return res;
}

double cov_NG_gg_kk(double l1,double l2, int n1, int n2){
   double a1,a2,array[4];
   int zmax, zmin;
   zmax = (n2 > n1 ? n2 : n1);
   a1 = amin_lens(zmax);
   zmin = (n2 < n1 ? n2 : n1);
   a2 = amax_lens(zmin);
   array[0] = l1;
   array[1] = l2;
   array[2] = (double) n1;
   array[3] = (double) n2;
   return int_gsl_integrate_low_precision(inner_project_tri_cov_gg_kk,(void*)array,a1,a2,NULL,1000);
}

/********** gg ks ***************/

double cov_G_gg_ks(double l, double delta_l, int n1, int n2, int zs){
   double fsky = survey.area/41253.0;
   double C13, C14, C23, C24;
   C13 = C_gk_nointerp(l, n1);
   C24 = C_gl_tomo_nointerp(l, n2, zs);
   C14 = C_gl_tomo_nointerp(l, n1, zs);
   C23 = C_gk_nointerp(l, n2);
   return (C13*C24+C14*C23)/((2.*l+1.)*delta_l*fsky);
}

double inner_project_tri_cov_gg_ks(double a,void *params)
{
   double k1,k2,fK,weights,res = 0.;
   double *ar = (double *) params;
   fK = f_K(chi(a));
   k1 = (ar[0]+0.5)/fK;
   k2 = (ar[1]+0.5)/fK;
   weights = W_gal(a,ar[2])*W_gal(a,ar[3])*W_k(a,fK)*W_kappa(a,fK,ar[4])*dchi_da(a);
   if (weights >0.){
      res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);
      res += delP_SSC(k1,a)*delP_SSC(k2,a)*survey_variance(a,survey.area/41253.0)*pow(fK,-4.); //super-sample covariance
   }
   res *= weights;
   return res;
}

double cov_NG_gg_ks(double l1,double l2, int n1, int n2, int zs){
   double a1,a2,array[5];
   int zmin, zmax;
   zmax = (n1 > n2 ? n1 : n2);
   zmin = (n1 < n2 ? n1 : n2);
   a1 = amin_lens(zmax);
   a2 = amax_lens(zmin);
   array[0] = l1;
   array[1] = l2;
   array[2] = (double) n1;
   array[3] = (double) n2;
   array[4] = (double) zs;
   return int_gsl_integrate_low_precision(inner_project_tri_cov_gg_ks,(void*)array,a1,a2,NULL,1000);
}

/********** gk gk ***************/

double cov_G_gk_gk(double l, double delta_l, int n1, int n3){
   double C13, C14, C23, C24, N13=0, N14=0, N23=0, N24=0;
   double fsky = survey.area/41253.0;
   C13 = C_cl_tomo_nointerp(l,n1,n3);
   C24 = C_kk_nointerp(l);
   C14 = C_gk_nointerp(l,n1);
   C23 = C_gk_nointerp(l,n3);
   if (n1 == n3){
      N13 = 1./(nlens(n1)*survey.n_gal_conversion_factor);
   }
   N24=kappa_reconstruction_noise(l);
   return ((C13+N13)*(C24+N24)+C14*C23)/((2.*l+1.)*delta_l*fsky);
}

double inner_project_tri_cov_gk_gk(double a,void *params)
{
   double k1,k2,fK,weights,res = 0.;
   double *ar = (double *) params;
   fK = f_K(chi(a));
   k1 = (ar[0]+0.5)/fK;
   k2 = (ar[1]+0.5)/fK;
   weights = W_gal(a,ar[2])*W_k(a,fK)*W_gal(a,ar[3])*W_k(a,fK)*dchi_da(a);
   if (weights >0.){
      res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);
      res += delP_SSC(k1,a)*delP_SSC(k2,a)*survey_variance(a,survey.area/41253.0)*pow(fK,-4.); //super-sample covariance
   }
   res *= weights;
   return res;
}

double cov_NG_gk_gk(double l1,double l2, int n1, int n2){
   double a1,a2,array[4];
   int zmax, zmin;
   zmax = (n2 > n1 ? n2 : n1);
   a1 = amin_lens(zmax);
   zmin = (n2 < n1 ? n2 : n1);
   a2 = amax_lens(zmin);
   array[0] = l1;
   array[1] = l2;
   array[2] = (double) n1;
   array[3] = (double) n2;
   return int_gsl_integrate_low_precision(inner_project_tri_cov_gk_gk,(void*)array,a1,a2,NULL,1000);
}

/********** gk gs ***************/

double cov_G_gk_gs(double l, double delta_l, int zl1, int zl2, int zs){
   double fsky = survey.area/41253.0;
   double C13, C14, C23, C24;
   C13 = C_cl_tomo_nointerp(l,zl1, zl2);
   C24 = C_ks_nointerp(l, zs);
   C14 = C_gl_tomo_nointerp(l, zl1, zs);
   C23 = C_gk_nointerp(l, zl2);
   double N = 0;
   if (zl1==zl2) N = 1./(nlens(zl1)*survey.n_gal_conversion_factor);
   return ((C13+N)*C24+C14*C23)/((2.*l+1.)*delta_l*fsky);
}

double inner_project_tri_cov_gk_gs(double a, void *params)
{
   double k1,k2,fK,weights,res = 0.;
   double *ar = (double *) params;
   fK = f_K(chi(a));
   k1 = (ar[0]+0.5)/fK;
   k2 = (ar[1]+0.5)/fK;
   weights = W_gal(a,ar[2])*W_k(a,fK)*W_gal(a,ar[3])*W_kappa(a, fK, ar[4])*dchi_da(a);
   if (weights >0.){
      res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);
      res += delP_SSC(k1,a)*delP_SSC(k2,a)*survey_variance(a,survey.area/41253.0)*pow(fK,-4.); //super-sample covariance
   }
   res *= weights;
   return res;
}

double cov_NG_gk_gs(double l1,double l2, int zl1, int zl2, int zs){
   double a1,a2,array[5];
   a1 = amin_lens(zl1);
   a2 = amax_lens(zl1);
   array[0] = l1;
   array[1] = l2;
   array[2] = (double) zl1;
   array[3] = (double) zl2;
   array[4] = (double) zs;
   return int_gsl_integrate_low_precision(inner_project_tri_cov_gk_gs,(void*)array,a1,a2,NULL,1000);
}

/********** gk kk ***************/

double cov_G_gk_kk(double l, double delta_l, int zl){
   // double fsky = survey.area/41253.0;
   double C13, C14, C23, C24;
   C13 = C_gk_nointerp(l,zl);
   C14 = C13;
   C24 = C_kk_nointerp(l);
   C23 = C24;
   double N = kappa_reconstruction_noise(l);
   return (C13*(C24+N)+C14*(C23+N))/((2.*l+1.)*delta_l*fsky_planck);
}

double inner_project_tri_cov_gk_kk(double a, void *params)
{
   double k1,k2,fK,weights,res = 0.;
   double *ar = (double *) params;
   fK = f_K(chi(a));
   k1 = (ar[0]+0.5)/fK;
   k2 = (ar[1]+0.5)/fK;
   weights = W_gal(a,ar[2])*pow(W_k(a,fK), 3)*dchi_da(a);
   if (weights >0.){
      res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(area_planck*survey.area_conversion_factor);
      res += delP_SSC(k1,a)*delP_SSC(k2,a)*cross_survey_variance(a,survey.area/41253.0,fsky_planck)*pow(fK,-4.); //super-sample covariance
   }
   res *= weights;
   return res;
}

double cov_NG_gk_kk(double l1,double l2, int zl){
   double a1,a2,array[3];
   a1 = amin_lens(zl);
   a2 = amax_lens(zl);
   array[0] = l1;
   array[1] = l2;
   array[2] = (double) zl;
   return int_gsl_integrate_low_precision(inner_project_tri_cov_gk_kk,(void*)array,a1,a2,NULL,1000);
}

/********** gk ks ***************/

double cov_G_gk_ks(double l, double delta_l, int zl, int zs){
   double fsky = survey.area/41253.0;
   double C13, C14, C23, C24;
   C13 = C_gk_nointerp(l,zl);
   C24 = C_ks_nointerp(l, zs);
   C14 = C_gl_tomo_nointerp(l, zl, zs);
   C23 = C_kk_nointerp(l);
   
   double N23 = kappa_reconstruction_noise(l);
   return (C13*C24+C14*(C23+N23))/((2.*l+1.)*delta_l*fsky);
}

double inner_project_tri_cov_gk_ks(double a, void *params)
{
   double k1,k2,fK,weights,res = 0.;
   double *ar = (double *) params;
   fK = f_K(chi(a));
   k1 = (ar[0]+0.5)/fK;
   k2 = (ar[1]+0.5)/fK;
   weights = W_gal(a,ar[2])*pow(W_k(a,fK), 2)*W_kappa(a,fK,ar[3])*dchi_da(a);
   if (weights >0.){
      res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);
      res += delP_SSC(k1,a)*delP_SSC(k2,a)*survey_variance(a,survey.area/41253.0)*pow(fK,-4.); //super-sample covariance
   }
   res *= weights;
   return res;
}

double cov_NG_gk_ks(double l1,double l2, int zl, int zs){
   double a1,a2,array[4];
   a1 = amin_lens(zl);
   a2 = amax_lens(zl);
   array[0] = l1;
   array[1] = l2;
   array[2] = (double) zl;
   array[3] = (double) zs;
   return int_gsl_integrate_low_precision(inner_project_tri_cov_gk_ks,(void*)array,a1,a2,NULL,1000);
}


/********** gk ss ***************/

double cov_G_gk_ss(double l, double delta_l, int zl, int zs1, int zs2){
   double fsky = survey.area/41253.0;
   double C13, C14, C23, C24;
   C13 = C_gl_tomo_nointerp(l,zl, zs1);
   C24 = C_ks_nointerp(l, zs2);
   C14 = C_gl_tomo_nointerp(l, zl, zs2);
   C23 = C_ks_nointerp(l, zs1);
   return (C13*C24+C14*C23)/((2.*l+1.)*delta_l*fsky);
}

double inner_project_tri_cov_gk_ss(double a, void *params)
{
   double k1,k2,fK,weights,res = 0.;
   double *ar = (double *) params;
   fK = f_K(chi(a));
   k1 = (ar[0]+0.5)/fK;
   k2 = (ar[1]+0.5)/fK;
   weights = W_gal(a,ar[2])*W_k(a,fK)*W_kappa(a, fK, ar[3])*W_kappa(a, fK, ar[4])*dchi_da(a);
   if (weights >0.){
      res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);
      res += delP_SSC(k1,a)*delP_SSC(k2,a)*survey_variance(a,survey.area/41253.0)*pow(fK,-4.); //super-sample covariance
   }
   res *= weights;
   return res;
}

double cov_NG_gk_ss(double l1,double l2, int zl, int zs1, int zs2){
   double a1,a2,array[5];
   a1 = amin_lens(zl);
   a2 = amax_lens(zl);
   array[0] = l1;
   array[1] = l2;
   array[2] = (double) zl;
   array[3] = (double) zs1;
   array[4] = (double) zs2;
   return int_gsl_integrate_low_precision(inner_project_tri_cov_gk_ss,(void*)array,a1,a2,NULL,1000);
}

/********** gs kk ***************/

double cov_G_gs_kk(double l, double delta_l, int zl, int zs){
   // double fsky = survey.area/41253.0;
   double C13, C14, C23, C24;
   C13 = C_gk_nointerp(l, zl);
   C14 = C13;
   C24 = C_ks_nointerp(l, zs);
   C23 = C24;
   return (C13*C24+C14*C23)/((2.*l+1.)*delta_l*fsky_planck);
}

double inner_project_tri_cov_gs_kk(double a,void *params)
{
   double fsky = survey.area/41253.0;
   double k1,k2,fK,weights,res = 0.;
   double *ar = (double *) params;
   fK = f_K(chi(a));
   k1 = (ar[0]+0.5)/fK;
   k2 = (ar[1]+0.5)/fK;
   weights = W_gal(a,ar[2])*W_kappa(a,fK,ar[3])*pow(W_k(a,fK), 2)*dchi_da(a);
   if (weights >0.){
      res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(area_planck*survey.area_conversion_factor);
      res += delP_SSC(k1,a)*delP_SSC(k2,a)*cross_survey_variance(a,fsky,fsky_planck)*pow(fK,-4.); //super-sample covariance
   }
   res *= weights;
   return res;
}

double cov_NG_gs_kk(double l1,double l2, int zl, int zs){
   double a1,a2,array[4];
   a1 = amin_source(zl);
   a2 = amax_source(zl);
   array[0] = l1;
   array[1] = l2;
   array[2] = (double) zl;
   array[3] = (double) zs;
   return int_gsl_integrate_low_precision(inner_project_tri_cov_gs_kk,(void*)array,a1,a2,NULL,1000);
}

/********** gs ks ***************/

double cov_G_gs_ks(double l, double delta_l, int zl, int zs1, int zs2){
   
   
   double fsky = survey.area/41253.0;
   double C13, C14, C23, C24;
   C13 = C_gk_nointerp(l, zl);
   C24 = C_shear_tomo_nointerp(l, zs1, zs2);
   C14 = C_gl_tomo_nointerp(l, zl, zs2);
   C23 = C_ks_nointerp(l, zs1);
   double N = 0;
   if (zs1==zs2) N = pow(survey.sigma_e,2.0)/(2.0*nsource(zs1)*survey.n_gal_conversion_factor);
   return (C13*(C24+N)+C14*C23)/((2.*l+1.)*delta_l*fsky);
}

double inner_project_tri_cov_gs_ks(double a,void *params)
{
   double k1,k2,fK,weights,res = 0.;
   double *ar = (double *) params;
   fK = f_K(chi(a));
   k1 = (ar[0]+0.5)/fK;
   k2 = (ar[1]+0.5)/fK;
   weights = W_gal(a,ar[2])*W_kappa(a,fK,ar[3])*W_k(a,fK)*W_kappa(a,fK,ar[4])*dchi_da(a);
   if (weights >0.){
      res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);
      res += delP_SSC(k1,a)*delP_SSC(k2,a)*survey_variance(a,survey.area/41253.0)*pow(fK,-4.); //super-sample covariance
   }
   res *= weights;
   return res;
}

double cov_NG_gs_ks(double l1,double l2, int zl, int zs1, int zs2){
   double a1,a2,array[5];
   a1 = amin_source(zl);
   a2 = amax_source(zl);
   array[0] = l1;
   array[1] = l2;
   array[2] = (double) zl;
   array[3] = (double) zs1;
   array[4] = (double) zs2;
   return int_gsl_integrate_low_precision(inner_project_tri_cov_gs_ks,(void*)array,a1,a2,NULL,1000);
}

/********** kk kk ***************/

double cov_G_kk_kk(double l, double delta_l){
   double C, N;
   // double fsky = survey.area/41253.0;
   C = C_kk_nointerp(l);
   N = kappa_reconstruction_noise(l);
//   printf("l, C, N = %le, %le, %le\n", l, C, N);
   return 2.*pow((C+N), 2) / ((2.*l+1.)*delta_l*fsky_planck);
}

double inner_project_tri_cov_kk_kk(double a,void *params)
{
   double fsky = survey.area/41253.0;
   double k1,k2,fK,weights,res = 0.;
   double *ar = (double *) params;
   fK = f_K(chi(a));
   k1 = (ar[0]+0.5)/fK;
   k2 = (ar[1]+0.5)/fK;
   weights = pow(W_k(a,fK), 4)*dchi_da(a);
   if (weights>0.){
// MANUWARNING: put an if to switch between low-z and high-z cov?
//      res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);
//      res += delP_SSC(k1,a)*delP_SSC(k2,a)*survey_variance(a,fsky)*pow(fK,-4.); //SSC
      res = tri_1h_cov(k1, k2, a)*pow(fK,-6.)/(area_planck*survey.area_conversion_factor); // T1h only
      res += delPlin_SSC(k1, a)*delPlin_SSC(k2, a)*cross_survey_variance(a,fsky,fsky_planck)*pow(fK,-4.); //SSC for Plin
   }
   res *= weights;
   return res;
}

double cov_NG_kk_kk(double l1,double l2){
   double a1, a2, array[2];
   array[0] = l1;
   array[1] = l2;
   a1 = limits.a_min*(1.+1.e-5);
   a2 = 1.-1.e-5;
   return int_gsl_integrate_low_precision(inner_project_tri_cov_kk_kk,(void*)array,a1,a2,NULL,1000);
}

/********** kk ks ***************/

double cov_G_kk_ks(double l, double delta_l, int zs){
   // double fsky = survey.area/41253.0;
   double C13, C24, C14, C23;
   C13 = C_kk_nointerp(l);
   C23 = C13;
   C14 = C_ks_nointerp(l, zs);
   C24 = C14;
   double N = kappa_reconstruction_noise(l);
   return ((C13+N)*C24+C14*(C23+N))/((2.*l+1.)*delta_l*fsky_planck);
}

double inner_project_tri_cov_kk_ks(double a,void *params)
{
   double fsky = survey.area/41253.0;
   double k1,k2,fK,weights,res = 0.;
   double *ar = (double *) params;
   fK = f_K(chi(a));
   k1 = (ar[0]+0.5)/fK;
   k2 = (ar[1]+0.5)/fK;
   weights = pow(W_k(a,fK), 3) * W_kappa(a, fK, ar[2]) *dchi_da(a);
   if (weights >0.){
      res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(area_planck*survey.area_conversion_factor);
      res += delP_SSC(k1,a)*delP_SSC(k2,a)*cross_survey_variance(a,fsky,fsky_planck)*pow(fK,-4.); //super-sample covariance
   }
   res *= weights;
   return res;
}

double cov_NG_kk_ks(double l1,double l2, int zs){
   double a1, a2, array[3];
   array[0] = l1;
   array[1] = l2;
   array[2] = (double) zs;
   a1 = amin_source(zs);
   a2 = 0.9999;
   return int_gsl_integrate_low_precision(inner_project_tri_cov_kk_ks,(void*)array,a1,a2,NULL,1000);
}

/********** kk ss ***************/

double cov_G_kk_ss(double l, double delta_l, int z1, int z2){
   // double fsky = survey.area/41253.0;
   double C13, C14, C23, C24;
   C13 = C_ks_nointerp(l,z1);
   C23 = C13;
   C24 = C_ks_nointerp(l, z2);
   C14 = C24;
   return (C13*C24+C14*C23)/((2.*l+1.)*delta_l*fsky_planck);
}

double inner_project_tri_cov_kk_ss(double a,void *params)
{
   double k1,k2,fK,weights,res = 0.;
   double *ar = (double *) params;
   fK = f_K(chi(a));
   k1 = (ar[0]+0.5)/fK;
   k2 = (ar[1]+0.5)/fK;
   weights = pow(W_k(a,fK), 2)*W_kappa(a,fK,ar[2])*W_kappa(a,fK,ar[3])*dchi_da(a);
   if (weights >0.){
      res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(area_planck*survey.area_conversion_factor);
      res += delP_SSC(k1,a)*delP_SSC(k2,a)*cross_survey_variance(a,survey.area/41253.0,fsky_planck)*pow(fK,-4.); //super-sample covariance
   }
   res *= weights;
   return res;
}

double cov_NG_kk_ss(double l1,double l2, int z1, int z2){
   double a1,a2,array[5];
   int zmin;
   zmin = (z2 < z1 ? z2 : z1);
   a1 = amin_source(zmin);
   a2 = amax_source(zmin);
   array[0] = l1;
   array[1] = l2;
   array[2] = (double) z1;
   array[3] = (double) z2;
   return int_gsl_integrate_low_precision(inner_project_tri_cov_kk_ss,(void*)array,a1,a2,NULL,1000);
}

/********** ks ks ***************/

double cov_G_ks_ks(double l, double delta_l, int z1, int z3){
   double C13, C14, C23, C24, N13=0, N14=0, N23=0, N24=0;
   double fsky = survey.area/41253.0;
   C13 = C_shear_tomo_nointerp(l,z1,z3);
   C24 = C_kk_nointerp(l);
   C14 = C_ks_nointerp(l,z1);
   C23 = C_ks_nointerp(l,z3);
   if (z1 == z3){
      N13= pow(survey.sigma_e,2.0)/(2.0*nsource(z1)*survey.n_gal_conversion_factor);
   }
   N24=kappa_reconstruction_noise(l);
   return ((C13+N13)*(C24+N24)+C14*C23)/((2.*l+1.)*delta_l*fsky);
}

double inner_project_tri_cov_ks_ks(double a,void *params)
{
   double k1,k2,fK,weights,res = 0.;
   double *ar = (double *) params;
   fK = f_K(chi(a));
   k1 = (ar[0]+0.5)/fK;
   k2 = (ar[1]+0.5)/fK;
   weights = W_kappa(a,fK,ar[2])*W_k(a,fK)*W_kappa(a,fK,ar[3])*W_k(a,fK)*dchi_da(a);
   if (weights >0.){
      res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);
      res += delP_SSC(k1,a)*delP_SSC(k2,a)*survey_variance(a,survey.area/41253.0)*pow(fK,-4.); //super-sample covariance
   }
   res *= weights;
   return res;
}

double cov_NG_ks_ks(double l1,double l2, int z1, int z2){
   double a1,a2,array[4];
   int zmin;
   zmin = (z2 < z1 ? z2 : z1);
   a1 = amin_source(zmin);
   a2 = amax_source(zmin);
   array[0] = l1;
   array[1] = l2;
   array[2] = (double) z1;
   array[3] = (double) z2;
   return int_gsl_integrate_low_precision(inner_project_tri_cov_ks_ks,(void*)array,a1,a2,NULL,1000);
}

/********** ks ss ***************/

double cov_G_ks_ss(double l, double delta_l, int z1, int z2, int z3){
   double fsky = survey.area/41253.0;
   double C13, C14, C23, C24;
   C13 = C_ks_nointerp(l,z2);
   C24 = C_shear_tomo_nointerp(l, z1, z3);
   C14 = C_ks_nointerp(l,z3);
   C23 = C_shear_tomo_nointerp(l, z1, z2);
   double N24 = 0, N23 = 0;
   if (z1 == z3){
      N24= pow(survey.sigma_e,2.0)/(2.0*nsource(z1)*survey.n_gal_conversion_factor);
   }
   if (z1 == z2){
      N23 = pow(survey.sigma_e,2.0)/(2.0*nsource(z1)*survey.n_gal_conversion_factor);
   }
   return (C13*(C24+N24)+C14*(C23+N23))/((2.*l+1.)*delta_l*fsky);
}

double inner_project_tri_cov_ks_ss(double a,void *params)
{
   double k1,k2,fK,weights,res = 0.;
   double *ar = (double *) params;
   fK = f_K(chi(a));
   k1 = (ar[0]+0.5)/fK;
   k2 = (ar[1]+0.5)/fK;
   weights = W_k(a,fK)*W_kappa(a,fK,ar[2])*W_kappa(a,fK,ar[3])*W_kappa(a,fK,ar[4])*dchi_da(a);
   if (weights >0.){
      res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);
      res += delP_SSC(k1,a)*delP_SSC(k2,a)*survey_variance(a,survey.area/41253.0)*pow(fK,-4.); //super-sample covariance
   }
   res *= weights;
   return res;
}

double cov_NG_ks_ss(double l1,double l2, int z1, int z2, int z3){
   double a1,a2,array[5];
   int zmin;
   zmin = (z2 < z1 ? z2 : z1);
   zmin = (z3 < zmin ? z3 : zmin);
   a1 = amin_source(zmin);
   a2 = amax_source(zmin);
   array[0] = l1;
   array[1] = l2;
   array[2] = (double) z1;
   array[3] = (double) z2;
   array[4] = (double) z3;
   return int_gsl_integrate_low_precision(inner_project_tri_cov_ks_ss,(void*)array,a1,a2,NULL,1000);
}
















