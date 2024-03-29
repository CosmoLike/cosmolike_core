/************************* covariance routines for angular power spectra *********************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
//     No T(l1,l2) look-up tables; most useful for small number of l-bins
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/************************* routines for angular trispectrum terms ************************/
/********************           AUTO-COVARIANCE BLOCKS              **********************/
double cov_NG_shear_shear_tomo(double l1,double l2, int z1, int z2, int z3, int z4);
double cov_G_shear_shear_tomo(double l, double delta_l, int z1, int z2, int z3, int z4);

double cov_NG_gl_gl_tomo(double l1,double l2, int z1l, int z1s, int z2l, int z2s);
double cov_G_gl_gl_tomo(double l, double delta_l, int z1l, int z1s, int z2l, int z2s);

double cov_NG_cl_cl_tomo(double l1,double l2, int z1, int z2, int z3, int z4);
double cov_G_cl_cl_tomo(double l, double delta_l, int z1, int z2, int z3, int z4);

/********************          CROSS-COVARIANCE BLOCKS              **********************/
double cov_NG_gl_shear_tomo(double l1,double l2, int zl, int zs, int z1, int z2);//zl,zs g-g lensing bins; z3,z4 shear bins
double cov_G_gl_shear_tomo(double l, double delta_l, int zl, int zs, int z1, int z2);

double cov_NG_cl_shear_tomo(double l1,double l2, int z1, int z2, int z3, int z4); //z1,z2 clustering bins; z3,z4 shear bins
double cov_G_cl_shear_tomo(double l, double delta_l, int z1, int z2, int z3, int z4);

double cov_NG_cl_gl_tomo(double l1,double l2, int z1, int z2, int zl, int zs); //z1,z2 clustering bins; zl,zs g-g lensing bins
double cov_G_cl_gl_tomo(double l, double delta_l, int z1, int z2, int zl, int zs);


/************** shear x shear routines ***************/
double inner_project_tri_cov_shear_shear_tomo(double a,void *params)
{
  double k1,k2,fK,weights,res = 0.;
  double *ar = (double *) params;
  fK = f_K(chi(a));
  k1 = (ar[0]+0.5)/fK;
  k2 = (ar[1]+0.5)/fK;
  weights = W_kappa(a,fK,ar[2])*W_kappa(a,fK,ar[3])*W_kappa(a,fK,ar[4])*W_kappa(a,fK,ar[5])*dchi_da(a);
  if (weights >0.){
    //res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);
    res = delP_SSC(k1,a)*delP_SSC(k2,a)*survey_variance(a,ar[6])*pow(fK,-4.); //SSC
  }
  res *= weights;
  return res;
}

double cov_NG_shear_shear_tomo(double l1,double l2, int z1, int z2, int z3, int z4){
  double a1,a2,array[7];
  int zmin;
  zmin = (z2 < z1 ? z2 : z1);
  zmin = (z3 < zmin ? z3 : zmin);
  zmin = (z4 < zmin ? z4 : zmin);
  a1 = amin_source(zmin);
  a2 = amax_source(zmin);
  array[0] = l1;
  array[1] = l2;
  array[2] = (double) z1;
  array[3] = (double) z2;
  array[4] = (double) z3;
  array[5] = (double) z4;
  array[6] = survey.area/41253.0;
  return int_gsl_integrate_low_precision(inner_project_tri_cov_shear_shear_tomo,(void*)array,a1,a2,NULL,1000);
}

double cov_G_shear_shear_tomo(double l, double delta_l, int z1, int z2, int z3, int z4){
  double C13, C14, C23, C24, N13 =0, N14=0, N23=0, N24=0;
   double fsky = survey.area/41253.0;
  C13 = C_shear_tomo_nointerp(l,z1,z3);C24 = C_shear_tomo_nointerp(l,z2,z4);
  C14 = C_shear_tomo_nointerp(l,z1,z4);C23 = C_shear_tomo_nointerp(l,z2,z3);
  if (z1 == z3){N13= pow(survey.sigma_e,2.0)/(2.0*nsource(z1)*survey.n_gal_conversion_factor);}
  if (z1 == z4){N14= pow(survey.sigma_e,2.0)/(2.0*nsource(z1)*survey.n_gal_conversion_factor);}
  if (z2 == z3){N23= pow(survey.sigma_e,2.0)/(2.0*nsource(z2)*survey.n_gal_conversion_factor);}
  if (z2 == z4){N24=pow(survey.sigma_e,2.0)/(2.0*nsource(z2)*survey.n_gal_conversion_factor);}
  
  return (C13*C24+ C13*N24+N13*C24 + C14*C23+C14*N23+N14*C23+N13*N24+N14*N23)/((2.*l+1.)*delta_l*fsky);
}

/***** gg-lensing x gg-lensing routines ****/
double bgal_a(double a, double nz){
  return gbias.b1_function(1./a-1.,(int)nz);
}
double inner_project_tri_cov_gl_gl_tomo(double a,void *params)
{
  double k1,k2,fK,weights,res = 0.;
  double *ar = (double *) params;
  fK = f_K(chi(a));
  k1 = (ar[0]+0.5)/fK;
  k2 = (ar[1]+0.5)/fK;
  weights = W_gal(a,ar[2])*W_kappa(a,fK,ar[3])*W_gal(a,ar[4])*W_kappa(a,fK,ar[5])*dchi_da(a);
  if (weights >0.){
    //res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);
    res =(delP_SSC(k1,a)-bgal_a(a,ar[2])*Pdelta(k1,a))*(delP_SSC(k2,a)-bgal_a(a,ar[4])*Pdelta(k2,a))*survey_variance(a,ar[6])*pow(fK,-4.); //SSC
  }
  res *= weights;
  return res;
}

double cov_NG_gl_gl_tomo(double l1,double l2, int z1l, int z1s, int z2l, int z2s){
  double a1,a2,array[7];
  int zmax, zmin;
  zmax = (z2l > z1l ? z2l : z1l);
  a1 = amin_lens(zmax);
  zmin = (z2l < z1l ? z2l : z1l);
  a2 = amax_lens(zmin);
//  printf("%le %le\n",l1,l2);
  array[0] = l1;
  array[1] = l2;
  array[2] = (double) z1l;
  array[3] = (double) z1s;
  array[4] = (double) z2l;
  array[5] = (double) z2s;
  array[6] = survey.area/41253.0;

  return int_gsl_integrate_low_precision(inner_project_tri_cov_gl_gl_tomo,(void*)array,a1,a2,NULL,1000);
}


double cov_G_gl_gl_tomo(double l, double delta_l, int z1, int z2, int z3, int z4){
  double C13 = 0., C14, C23, C24, N13 =0, N24=0;
   double fsky = survey.area/41253.0;
  if (z1 ==z3){C13 = C_cl_tomo_nointerp(l,z1,z3);}
  C24 = C_shear_tomo_nointerp(l,z2,z4);
  C14 = C_gl_tomo_nointerp(l,z1,z4);
  C23 = C_gl_tomo_nointerp(l,z3,z2);
  if (z1 == z3){N13= 1./(nlens(z1)*survey.n_gal_conversion_factor);}
  if (z2 == z4){N24= pow(survey.sigma_e,2.0)/(2.0*nsource(z2)*survey.n_gal_conversion_factor);}
  
  return (C13*C24+ C13*N24+N13*C24 + N13*N24+C14*C23)/((2.*l+1.)*delta_l*fsky);
}


/***** clustering x clustering routines ****/

double inner_project_tri_cov_cl_cl_tomo(double a,void *params)
{
  double k1,k2,fK,weights,res = 0.;
  double *ar = (double *) params;
  fK = f_K(chi(a));
  k1 = (ar[0]+0.5)/fK;
  k2 = (ar[1]+0.5)/fK;
  weights = W_gal(a,ar[2])*W_gal(a,ar[3])*W_gal(a,ar[4])*W_gal(a, ar[5])*dchi_da(a);
  if (weights >0.){
    //res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);
    res = (delP_SSC(k1,a)-2.*bgal_a(a,ar[2])*Pdelta(k1,a))*(delP_SSC(k2,a)-2.*bgal_a(a,ar[4])*Pdelta(k2,a))*survey_variance(a,ar[6])*pow(fK,-4.); //SSC
  }
  res *= weights;
  return res;
}


double cov_NG_cl_cl_tomo(double l1,double l2, int z1, int z2, int z3, int z4){
  double a1,a2,array[7];
  int zmin, zmax;
  zmax = (z2 > z1 ? z2 : z1);
  zmax = (z3 > zmax ? z3 : zmax);
  zmax = (z4 > zmax ? z4 : zmax);
  a1 = amin_lens(zmax);
  zmin = (z2 < z1 ? z2 : z1);
  zmin = (z3 < zmin ? z3 : zmin);
  zmin = (z4 < zmin ? z4 : zmin);
  a2 = amax_lens(zmin);
  array[0] = l1;
  array[1] = l2;
  array[2] = (double) z1;
  array[3] = (double) z2;
  array[4] = (double) z3;
  array[5] = (double) z4;
  array[6] = survey.area/41253.0;

  return int_gsl_integrate_low_precision(inner_project_tri_cov_cl_cl_tomo,(void*)array,a1,a2,NULL,1000);
}


double cov_G_cl_cl_tomo(double l, double delta_l, int z1, int z2, int z3, int z4){
  double C13, C14, C23, C24, N13 =0, N14=0, N23=0, N24=0;
   double fsky = survey.area/41253.0;
  C13 = C_cl_tomo_nointerp(l,z1,z3);C24 = C_cl_tomo_nointerp(l,z2,z4);
  C14 = C_cl_tomo_nointerp(l,z1,z4);C23 = C_cl_tomo_nointerp(l,z2,z3);
  if (z1 == z3){N13= 1./(nlens(z1)*survey.n_gal_conversion_factor);}
  if (z1 == z4){N14= 1./(nlens(z1)*survey.n_gal_conversion_factor);}
  if (z2 == z3){N23= 1./(nlens(z2)*survey.n_gal_conversion_factor);}
  if (z2 == z4){N24= 1./(nlens(z2)*survey.n_gal_conversion_factor);}
  return ((C13+N13)*(C24+N24)+(C14+N14)*(C23+N23))/((2.*l+1.)*delta_l*fsky);
}


/************** g-g lensing x shear routines ***************/

double inner_project_tri_cov_gl_shear_tomo(double a,void *params)
{
  double k1,k2,fK,weights,res = 0.;
  double *ar = (double *) params;
  fK = f_K(chi(a));
  k1 = (ar[0]+0.5)/fK;
  k2 = (ar[1]+0.5)/fK;
  weights = W_gal(a, ar[2])*W_kappa(a,fK,ar[3])*W_kappa(a,fK,ar[4])*W_kappa(a,fK,ar[5])*dchi_da(a);
  if (weights >0.){
    //res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);
    res = (delP_SSC(k1,a)-bgal_a(a,ar[2])*Pdelta(k1,a))*(delP_SSC(k2,a))*survey_variance(a,ar[6])*pow(fK,-4.); //SSC
  }
  res *= weights;
  return res;
}

double cov_NG_gl_shear_tomo(double l1,double l2, int zl, int zs, int z3, int z4){ //zl,zs g-g lensing bins; z3,z4 shear bins
  double a1,a2,array[7];
  a1 = amin_lens(zl);
  a2 = amax_lens(zl);
  array[0] = l1;
  array[1] = l2;
  array[2] = (double) zl;
  array[3] = (double) zs;
  array[4] = (double) z3;
  array[5] = (double) z4;
  array[6] = survey.area/41253.0;
  return int_gsl_integrate_low_precision(inner_project_tri_cov_gl_shear_tomo,(void*)array,a1,a2,NULL,1000);
}

double cov_G_gl_shear_tomo(double l, double delta_l, int zl, int zs, int z3, int z4){
  double C13, C14, C23, C24, N13 =0, N14=0, N23=0, N24=0;
   double fsky = survey.area/41253.0;
  C13 = C_gl_tomo_nointerp(l,zl,z3);C24 = C_shear_tomo_nointerp(l,zs,z4);
  C14 = C_gl_tomo_nointerp(l,zl,z4);C23 = C_shear_tomo_nointerp(l,zs,z3);
  if (zs == z3){N23= pow(survey.sigma_e,2.0)/(2.0*nsource(zs)*survey.n_gal_conversion_factor);}
  if (zs == z4){N24=pow(survey.sigma_e,2.0)/(2.0*nsource(zs)*survey.n_gal_conversion_factor);}
  
  return (C13*(C24+N24) + C14*(C23+N23))/((2.*l+1.)*delta_l*fsky);
    
}

/************** clustering x shear routines ***************/

double inner_project_tri_cov_cl_shear_tomo(double a,void *params)
{
  double k1,k2,fK,weights,res = 0.;
  double *ar = (double *) params;
  fK = f_K(chi(a));
  k1 = (ar[0]+0.5)/fK;
  k2 = (ar[1]+0.5)/fK;
  weights = W_gal(a,ar[2])*W_gal(a, ar[3])*W_kappa(a,fK,ar[4])*W_kappa(a,fK,ar[5])*dchi_da(a);
  if (weights >0.){
    //res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);
    res = (delP_SSC(k1,a)-2.*bgal_a(a,ar[2])*Pdelta(k1,a))*(delP_SSC(k2,a))*survey_variance(a,ar[6])*pow(fK,-4.); //SSC
  }
  res *= weights;
  return res;
}

double cov_NG_cl_shear_tomo(double l1,double l2, int z1, int z2, int z3, int z4){ //z1,z2 clustering bins; z3,z4 shear bins
  double a1,a2,array[7];
  int zmin, zmax;
  zmax = (z1 > z2 ? z1 : z2);
  zmin = (z1 < z2 ? z1 : z2);
  a1 = amin_lens(zmax);
  a2 = amax_lens(zmin);
  array[0] = l1;
  array[1] = l2;
  array[2] = (double) z1;
  array[3] = (double) z2;
  array[4] = (double) z3;
  array[5] = (double) z4;
  array[6] = survey.area/41253.0;
  return int_gsl_integrate_low_precision(inner_project_tri_cov_cl_shear_tomo,(void*)array,a1,a2,NULL,1000);
}

double cov_G_cl_shear_tomo(double l, double delta_l, int z1, int z2, int z3, int z4){
  double C13, C14, C23, C24, N13 =0, N14=0, N23=0, N24=0;
   double fsky = survey.area/41253.0;
  C13 = C_gl_tomo_nointerp(l,z1,z3);C24 = C_gl_tomo_nointerp(l,z2,z4);
  C14 = C_gl_tomo_nointerp(l,z1,z4);C23 = C_gl_tomo_nointerp(l,z2,z3);
  
  return (C13*C24+  C14*C23)/((2.*l+1.)*delta_l*fsky);
  
}

/************** clustering x shear routines ***************/

double inner_project_tri_cov_cl_gl_tomo(double a,void *params)
{
  double k1,k2,fK,weights,res = 0.;
  double *ar = (double *) params;
  fK = f_K(chi(a));
  k1 = (ar[0]+0.5)/fK;
  k2 = (ar[1]+0.5)/fK;
  weights = W_gal(a,ar[2])*W_gal(a,ar[3])*W_gal(a,ar[4])*W_kappa(a,fK,ar[5])*dchi_da(a);
  if (weights >0.){
    //res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);
    res = (delP_SSC(k1,a)-2.*bgal_a(a,ar[2])*Pdelta(k1,a))*(delP_SSC(k2,a)-bgal_a(a,ar[4])*Pdelta(k2,a))*survey_variance(a,ar[6])*pow(fK,-4.); //SSC
  }
  res *= weights;
  return res;
}

double cov_NG_cl_gl_tomo(double l1,double l2, int z1, int z2, int zl, int zs){ //z1,z2 clustering bins; zl,zs g-g lensing bins
  double a1,a2,array[7];
  a1 = amin_lens(z1);
  a2 = amax_lens(z2);
  array[0] = l1;
  array[1] = l2;
  array[2] = (double) z1;
  array[3] = (double) z2;
  array[4] = (double) zl;
  array[5] = (double) zs;
  array[6] = survey.area/41253.0;
  return int_gsl_integrate_low_precision(inner_project_tri_cov_cl_gl_tomo,(void*)array,a1,a2,NULL,1000);
}

double cov_G_cl_gl_tomo(double l, double delta_l, int z1, int z2, int zl, int zs){
  double C13, C14, C23, C24, N13 =0, N14=0, N23=0, N24=0;
   double fsky = survey.area/41253.0;
  C13 = C_cl_tomo_nointerp(l,z1,zl);C24 = C_gl_tomo_nointerp(l,z2,zs);
  C14 = C_gl_tomo_nointerp(l,z1,zs);C23 = C_cl_tomo_nointerp(l,z2,zl);
  if (z1 == zl){N13= 1./(nlens(z1)*survey.n_gal_conversion_factor);}
  if (z2 == zl){N23= 1./(nlens(z1)*survey.n_gal_conversion_factor);}
  
  return ((C13+N13)*C24 + C14*(C23 +N23))/((2.*l+1.)*delta_l*fsky);
  
}

