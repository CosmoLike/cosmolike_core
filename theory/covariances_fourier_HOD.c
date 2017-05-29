/************************* covariance routines for angular power spectra *********************/
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
//     No T(l1,l2) look-up tables; most useful for small number of l-bins
/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!*/
/************************* routines for angular trispectrum terms ************************/
/********************           AUTO-COVARIANCE BLOCKS              **********************/
double cov_NG_gl_gl_tomo_HOD(double l1,double l2, int z1l, int z1s, int z2l, int z2s);
double cov_G_gl_gl_tomo_HOD(double l, double delta_l, int z1l, int z1s, int z2l, int z2s);

double cov_NG_cl_cl_tomo(double l1,double l2, int z1, int z2, int z3, int z4);
double cov_G_cl_cl_tomo(double l, double delta_l, int z1, int z2, int z3, int z4);

/********************          CROSS-COVARIANCE BLOCKS              **********************/
double cov_NG_gl_shear_tomo_HOD(double l1,double l2, int zl, int zs, int z1, int z2);//zl,zs g-g lensing bins; z3,z4 shear bins
double cov_G_gl_shear_tomo_HOD(double l, double delta_l, int zl, int zs, int z1, int z2);

double cov_NG_cl_shear_tomo_HOD(double l1,double l2, int z1, int z2, int z3, int z4); //z1,z2 clustering bins; z3,z4 shear bins
double cov_G_cl_shear_tomo_HOD(double l, double delta_l, int z1, int z2, int z3, int z4);

double cov_NG_cl_gl_tomo_HOD(double l1,double l2, int z1, int z2, int zl, int zs); //z1,z2 clustering bins; zl,zs g-g lensing bins
double cov_G_cl_gl_tomo_HOD(double l, double delta_l, int z1, int z2, int zl, int zs);

/********************      CROSS-COVARIANCE BLOCKS - CLUSTERS       **********************/
double cov_G_cl_cgl_HOD (double l, double delta_l,int nzl1,int nzl2, int nzc, int nN, int nzs);
double cov_NG_cl_cgl_HOD (double l1,double l2,int nzl1,int nzl2, int nzc, int nN, int nzs);
double cov_G_ggl_cgl_HOD (double l, double delta_l,int nzl,int nzs1, int nzc, int nN, int nzs2);
double cov_NG_ggl_cgl_HOD (double l1, double l2, int nzl,int nzs1, int nzc, int nN, int nzs2);

double cov_cl_N_HOD (double l, double delta_l,int nzl1,int nzl2, int nzc, int nN);
double cov_ggl_N_HOD (double l, double delta_l,int nzl,int nzs, int nzc, int nN);


/***** gg-lensing x gg-lensing routines ****/
double W_rm(double a, double nz){
  return bgal_rm(a)*pf_photoz(1./a-1,(int)nz)*hoverh0(a);
}
double inner_SSC_gg (double logm, void *para){
	double *array = (double *) para;
	double m = exp(logm);
	double u,ns,k = array[0],a = array[1];
	u = u_g_rm(k,m,a);
	ns = n_s_rm(m);
	return massfunc(m,a)*m*B1(m,a)*(u*u*ns*ns*1.0+2.0*u*ns*n_c_rm(m));
}

double inner_SSC_gm (double logm, void *para){
	double *array = (double *) para;
	double m = exp(logm);
	double k = array[0],a = array[1];
	return massfunc(m,a)*m*m*B1(m,a)/(cosmology.rho_crit*cosmology.Omega_m)*u_nfw_c(conc(m,a),k,m,a)*(u_g_rm(k,m,a)*n_s_rm(m)+n_c_rm(m));
}

double I12_SSC_gm (double k,double a){//one-halo term contribution to super sample covariance
  static cosmopara C;
  static double amin = 0, amax = 0,logkmin = 0., logkmax = 0., dk = 0., da = 0.,res = 0.;
  
  static double **table_I1=0;
  
  double aa,klog;
  int i,j;
  
  if (recompute_cosmo3D(C)){ //extend this by halo model parameters if these become part of the model
    update_cosmopara(&C);
    if (table_I1==0) table_I1 = create_double_matrix(0, Ntable.N_a-1, 0, Ntable.N_k_nlin-1);
    double result,error,array[2];
    
    amin= 1./(1+redshift.shear_zdistrpar_zmax);
    amax = 1.;
    da = (amax - amin)/(Ntable.N_a);
    aa = amin;
    logkmin = log(limits.k_min_cH0);
    logkmax = log(limits.k_max_cH0);
    dk = (logkmax - logkmin)/(Ntable.N_k_nlin);
    for (i=0; i<Ntable.N_a; i++, aa +=da) {
      klog  = logkmin;
      for (j=0; j<Ntable.N_k_nlin; j++, klog += dk) {
        array[0]= exp(klog);array[1]= aa;
        table_I1[i][j] = log(int_gsl_integrate_medium_precision(inner_SSC_gm,(void*)array,log(limits.M_min),log(limits.M_max),NULL, 2000)/ngal_rm(aa));
      }
    }
  }
  if (log(k) <logkmin){return 1.0;};
  if (log(k) >logkmax){return 0.0;};
  aa = fmin(a,amax-1.1*da);
  return exp(interpol2d(table_I1, Ntable.N_a, amin, amax, da, aa, Ntable.N_k_nlin, logkmin, logkmax, dk, log(k), 1.0, 1.0));
}
double I12_SSC_gg (double k,double a){//one-halo term contribution to super sample covariance
  static cosmopara C;
  static double amin = 0, amax = 0,logkmin = 0., logkmax = 0., dk = 0., da = 0.,res = 0.;
  
  static double **table_I1=0;
  
  double aa,klog;
  int i,j;
  
  if (recompute_cosmo3D(C)){ //extend this by halo model parameters if these become part of the model
    update_cosmopara(&C);
    if (table_I1==0) table_I1 = create_double_matrix(0, Ntable.N_a-1, 0, Ntable.N_k_nlin-1);
    double result,error,array[2];
    
    amin= 1./(1+redshift.shear_zdistrpar_zmax);
    amax = 1.;
    da = (amax - amin)/(Ntable.N_a);
    aa = amin;
    logkmin = log(limits.k_min_cH0);
    logkmax = log(limits.k_max_cH0);
    dk = (logkmax - logkmin)/(Ntable.N_k_nlin);
    for (i=0; i<Ntable.N_a; i++, aa +=da) {
      klog  = logkmin;
      for (j=0; j<Ntable.N_k_nlin; j++, klog += dk) {
        array[0]= exp(klog);array[1]= aa;
        table_I1[i][j] = log(int_gsl_integrate_medium_precision(inner_SSC_gg,(void*)array,log(limits.M_min),log(limits.M_max),NULL, 2000)/pow(ngal_rm(aa),2.));
      }
    }
  }
  if (log(k) <logkmin){return 1.0;};
  if (log(k) >logkmax){return 0.0;};
  aa = fmin(a,amax-1.1*da);
  return exp(interpol2d(table_I1, Ntable.N_a, amin, amax, da, aa, Ntable.N_k_nlin, logkmin, logkmax, dk, log(k), 1.0, 1.0));
}
double delPgg_SSC(double k, double a){
  double b = bgal_rm(a);
  return delP_SSC(k,a)-2.*b*Pdelta(k,a);
//  return (68./21.+LD_term(k)-b*b)*p_2h(k,a)+I12_SSC_gg(k,a)/(b*b);
}

double delPgm_SSC(double k, double a){
  double b = bgal_rm(a);
  return delP_SSC(k,a)-b*Pdelta(k,a);
//  return (68./21.+LD_term(k)-b)*p_2h(k,a)+I12_SSC_gm(k,a)/b;
}

double inner_project_tri_cov_gl_gl_tomo_HOD(double a,void *params)
{
  double k1,k2,fK,weights,res = 0.;
  double *ar = (double *) params;
  fK = f_K(chi(a));
  k1 = (ar[0]+0.5)/fK;
  k2 = (ar[1]+0.5)/fK;
  weights = W_rm(a,ar[2])*W_kappa(a,fK,ar[3])*W_rm(a,ar[4])*W_kappa(a,fK,ar[5])*dchi_da(a);
  if (weights >0.){
    res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);
    res += delPgm_SSC(k1,a)*delPgm_SSC(k2,a)*survey_variance(a,ar[6])*pow(fK,-4.); //SSC
  }
  res *= weights;
  return res;
}

double cov_NG_gl_gl_tomo_HOD(double l1,double l2, int z1l, int z1s, int z2l, int z2s){
  double a1,a2,array[7];
  int zmax, zmin;
  if (z1l != z2l){return 0.;}
  a1 = amin_lens(z1l);
  a2 = amax_lens(z2l);
  array[0] = l1;
  array[1] = l2;
  array[2] = (double) z1l;
  array[3] = (double) z1s;
  array[4] = (double) z2l;
  array[5] = (double) z2s;
  array[6] = survey.area/41253.0;

  return int_gsl_integrate_low_precision(inner_project_tri_cov_gl_gl_tomo_HOD,(void*)array,a1,a2,NULL,1000);
}


double cov_G_gl_gl_tomo_HOD(double l, double delta_l, int z1, int z2, int z3, int z4){
  double C13 = 0., C14, C23, C24, N13 =0, N24=0;
   double fsky = survey.area/41253.0;
  if (z1 ==z3){C13 = C_cl_HOD_rm_tomo(l,z1);}
  C24 = C_shear_tomo_nointerp(l,z2,z4);
  C14 = C_gl_HOD_rm_tomo(l,z1,z4);
  C23 = C_gl_HOD_rm_tomo(l,z3,z2);
  if (z1 == z3){N13= 1./(nlens(z1)*survey.n_gal_conversion_factor);}
  if (z2 == z4){N24= pow(survey.sigma_e,2.0)/(2.0*nsource(z2)*survey.n_gal_conversion_factor);}
  
  return (C13*C24+ C13*N24+N13*C24 + N13*N24+C14*C23)/((2.*l+1.)*delta_l*fsky);
}


/***** clustering x clustering routines ****/

double inner_project_tri_cov_cl_cl_tomo_HOD(double a,void *params)
{
  double k1,k2,fK,weights,res =0.;
  double *ar = (double *) params;
  fK = f_K(chi(a));
  k1 = (ar[0]+0.5)/fK;
  k2 = (ar[1]+0.5)/fK;
  weights = W_rm(a,ar[2])*W_rm(a,ar[3])*W_rm(a,ar[4])*W_rm(a, ar[5])*dchi_da(a);
  if (weights >0.){
    res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);
    res += delPgg_SSC(k1,a)*delPgg_SSC(k2,a)*survey_variance(a,ar[6])*pow(fK,-4.); //SSC
  }
  res *= weights;
  return res;
}


double cov_NG_cl_cl_tomo_HOD(double l1,double l2, int z1, int z2, int z3, int z4){
  double a1,a2,array[7];
  if (z1 != z2){return 0.;}
  if (z3 != z4){return 0.;}
  if (z1 != z3){return 0.;}
  a1 = amin_lens(z1);
  a2 = amax_lens(z1);
  array[0] = l1;
  array[1] = l2;
  array[2] = (double) z1;
  array[3] = (double) z2;
  array[4] = (double) z3;
  array[5] = (double) z4;
  array[6] = survey.area/41253.0;

  return int_gsl_integrate_low_precision(inner_project_tri_cov_cl_cl_tomo_HOD,(void*)array,a1,a2,NULL,1000);
}


double cov_G_cl_cl_tomo_HOD(double l, double delta_l, int z1, int z2, int z3, int z4){
  double C13 = 0., C14 = 0., C23 = 0., C24 = 0., N13 =0, N14=0, N23=0, N24=0;
   double fsky = survey.area/41253.0;
  if (z1 == z3){
    N13= 1./(nlens(z1)*survey.n_gal_conversion_factor);
    C13 = C_cl_HOD_rm_tomo(l,z1);
  }
  if (z1 == z4){
    N14= 1./(nlens(z1)*survey.n_gal_conversion_factor);
    C14 = C_cl_HOD_rm_tomo(l,z1);
  }
  if (z2 == z3){
    N23= 1./(nlens(z2)*survey.n_gal_conversion_factor);
    C23 = C_cl_HOD_rm_tomo(l,z2);
  }
  if (z2 == z4){
    N24= 1./(nlens(z2)*survey.n_gal_conversion_factor);
    C24 = C_cl_HOD_rm_tomo(l,z2);
  }
  
  return (C13*C24+ C13*N24+N13*C24 + C14*C23+C14*N23+N14*C23+N13*N24+N14*N23)/((2.*l+1.)*delta_l*fsky);
}


/************** g-g lensing x shear routines ***************/

double inner_project_tri_cov_gl_shear_tomo_HOD(double a,void *params)
{
  double k1,k2,fK,weights,res =0.;
  double *ar = (double *) params;
  fK = f_K(chi(a));
  k1 = (ar[0]+0.5)/fK;
  k2 = (ar[1]+0.5)/fK;
  weights = W_rm(a, ar[2])*W_kappa(a,fK,ar[3])*W_kappa(a,fK,ar[4])*W_kappa(a,fK,ar[5])*dchi_da(a);
  if (weights >0.){
    res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);
    res += delPgm_SSC(k1,a)*delP_SSC(k2,a)*survey_variance(a,ar[6])*pow(fK,-4.); //SSC
  }
  res *= weights;
  return res;
}

double cov_NG_gl_shear_tomo_HOD(double l1,double l2, int zl, int zs, int z3, int z4){ //zl,zs g-g lensing bins; z3,z4 shear bins
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
  return int_gsl_integrate_low_precision(inner_project_tri_cov_gl_shear_tomo_HOD,(void*)array,a1,a2,NULL,1000);
}

double cov_G_gl_shear_tomo_HOD(double l, double delta_l, int zl, int zs, int z3, int z4){
  double C13, C14, C23, C24, N13 =0, N14=0, N23=0, N24=0;
   double fsky = survey.area/41253.0;
  C13 = C_gl_HOD_rm_tomo(l,zl,z3);C24 = C_shear_tomo_nointerp(l,zs,z4);
  C14 = C_gl_HOD_rm_tomo(l,zl,z4);C23 = C_shear_tomo_nointerp(l,zs,z3);
  if (zs == z3){N23= pow(survey.sigma_e,2.0)/(2.0*nsource(zs)*survey.n_gal_conversion_factor);}
  if (zs == z4){N24=pow(survey.sigma_e,2.0)/(2.0*nsource(zs)*survey.n_gal_conversion_factor);}
  
  return (C13*(C24+N24) + C14*(C23+N23))/((2.*l+1.)*delta_l*fsky);
    
}

/************** clustering x shear routines ***************/

double inner_project_tri_cov_cl_shear_tomo_HOD(double a,void *params)
{
  double k1,k2,fK,weights,res =0.;
  double *ar = (double *) params;
  fK = f_K(chi(a));
  k1 = (ar[0]+0.5)/fK;
  k2 = (ar[1]+0.5)/fK;
  weights = W_rm(a,ar[2])*W_rm(a, ar[3])*W_kappa(a,fK,ar[4])*W_kappa(a,fK,ar[5])*dchi_da(a);
  if (weights >0.){
    res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);
    res += delPgg_SSC(k2,a)*delP_SSC(k2,a)*survey_variance(a,ar[6])*pow(fK,-4.); //SSC
  }
  res *= weights;
  return res;
}

double cov_NG_cl_shear_tomo_HOD(double l1,double l2, int z1, int z2, int z3, int z4){ //z1,z2 clustering bins; z3,z4 shear bins
  double a1,a2,array[7];
  int zmin, zmax;
  if (z1 != z2){return 0.;}
  a1 = amin_lens(z1);
  a2 = amax_lens(z1);
  array[0] = l1;
  array[1] = l2;
  array[2] = (double) z1;
  array[3] = (double) z2;
  array[4] = (double) z3;
  array[5] = (double) z4;
  array[6] = survey.area/41253.0;
  return int_gsl_integrate_low_precision(inner_project_tri_cov_cl_shear_tomo_HOD,(void*)array,a1,a2,NULL,1000);
}

double cov_G_cl_shear_tomo_HOD(double l, double delta_l, int z1, int z2, int z3, int z4){
  double C13, C14, C23, C24, N13 =0, N14=0, N23=0, N24=0;
   double fsky = survey.area/41253.0;
  C13 = C_gl_HOD_rm_tomo(l,z1,z3);C24 = C_gl_HOD_rm_tomo(l,z2,z4);
  C14 = C_gl_HOD_rm_tomo(l,z1,z4);C23 = C_gl_HOD_rm_tomo(l,z2,z3);
  
  return (C13*C24+  C14*C23)/((2.*l+1.)*delta_l*fsky);
  
}

/************** clustering x shear routines ***************/

double inner_project_tri_cov_cl_gl_tomo_HOD(double a,void *params)
{
  double k1,k2,fK,weights,res =0.;
  double *ar = (double *) params;
  fK = f_K(chi(a));
  k1 = (ar[0]+0.5)/fK;
  k2 = (ar[1]+0.5)/fK;
  weights = W_rm(a,ar[2])*W_rm(a,ar[3])*W_rm(a,ar[4])*W_kappa(a,fK,ar[5])*dchi_da(a);
  if (weights >0.){
    res = tri_matter_cov(k1,k2,a)*pow(fK,-6.)/(survey.area*survey.area_conversion_factor);
    res += delPgg_SSC(k1,a)*delPgm_SSC(k2,a)*survey_variance(a,ar[6])*pow(fK,-4.); //SSC
  }
  res *= weights;
  return res;
}

double cov_NG_cl_gl_tomo_HOD(double l1,double l2, int z1, int z2, int zl, int zs){ //z1,z2 clustering bins; zl,zs g-g lensing bins
  double a1,a2,array[7];
  if (z1 != z2){return 0.;}
  if (z1 != zl){return 0.;}
  a1 = amin_lens(z1);
  a2 = amax_lens(z1);
  array[0] = l1;
  array[1] = l2;
  array[2] = (double) z1;
  array[3] = (double) z2;
  array[4] = (double) zl;
  array[5] = (double) zs;
  array[6] = survey.area/41253.0;
  return int_gsl_integrate_low_precision(inner_project_tri_cov_cl_gl_tomo_HOD,(void*)array,a1,a2,NULL,1000);
}

double cov_G_cl_gl_tomo_HOD(double l, double delta_l, int z1, int z2, int zl, int zs){
  double C13 =0., C14, C23 = 0., C24, N13 =0, N14=0, N23=0, N24=0;
   double fsky = survey.area/41253.0;
  C24 = C_gl_HOD_rm_tomo(l,z2,zs);
  C14 = C_gl_HOD_rm_tomo(l,z1,zs);
  if (z1 == zl){
    N13= 1./(nlens(z1)*survey.n_gal_conversion_factor);
    C13 = C_cl_HOD_rm_tomo(l,z1);
  }
  if (z2 == zl){
    N23= 1./(nlens(z1)*survey.n_gal_conversion_factor);
    C23 =C_cl_HOD_rm_tomo(l,z2);
  }
  
  return ((C13+N13)*C24 + C14*(C23 +N23))/((2.*l+1.)*delta_l*fsky);
  
}

/************ HOD x cluster routines ****************/
double project_cov_cl_N_HOD(double a,void *params)
{
  double k,res,weights,wa;
  double *ar = (double *) params;
  wa = f_K(chi(a));
  k = ar[0]/wa;
  weights =W_rm(a,(int)ar[1]); //normalized redshift distribution of galaxy sample
	weights *= W_rm(a,(int)ar[1]);
  weights *= dN200_dchi(a,(int)ar[3]);//dN/dV of number counts cluster bin
  weights *= dchi_da(a);
  if (weights > 0){
    res= delPgg_SSC(k,a)*b_cluster((int)ar[2],(int)ar[3])*survey_variance(a,survey.area/41253.0);
  }
  return res*weights;
}
double project_cov_ggl_N_HOD(double a,void *params)
{
  double k,res,weights,wa;
  double *ar = (double *) params;
  wa = f_K(chi(a));
  k = ar[0]/wa;
  weights =W_rm(a,(int)ar[1]); //normalized redshift distribution of galaxy sample
	weights *= W_kappa(a,wa,ar[2]);
  weights *= dN200_dchi(a,(int)ar[4]);//dN/dV of number counts cluster bin
  weights *= dchi_da(a);
  if (weights > 0){
    res= delPgm_SSC(k,a)*b_cluster((int)ar[3],(int)ar[4])*survey_variance(a,survey.area/41253.0);
  }
  return res*weights;
}

double cov_cl_N_HOD(double l,int nl1,int nl2, int nzc, int nN){
  if (nl1 != nl2){return 0.;}
  if ((tomo.clustering_zmin[nl1]+tomo.clustering_zmax[nl1])/2. > tomo.cluster_zmax[nzc]){return 0.;}
  if ((tomo.clustering_zmin[nl1]+tomo.clustering_zmax[nl1])/2. < tomo.cluster_zmin[nzc]){return 0.;}
  double array[4] = {l,(double)nl1, (double) nzc, (double) nN};
  return int_gsl_integrate_medium_precision(project_cov_cl_N_HOD,(void*) array, 1./(1+tomo.cluster_zmax[nzc]),1./(1+tomo.cluster_zmin[nzc]), NULL,1000);
}
double cov_ggl_N_HOD(double l,int nzl,int nzs, int nzc, int nN){
  if ((tomo.clustering_zmin[nzl]+tomo.clustering_zmax[nzl])/2. > tomo.cluster_zmax[nzc]){return 0.;}
  if ((tomo.clustering_zmin[nzl]+tomo.clustering_zmax[nzl])/2. < tomo.cluster_zmin[nzc]){return 0.;}
  double array[5] = {l,(double)nzl,(double)nzs, (double) nzc, (double) nN};
  return int_gsl_integrate_medium_precision(project_cov_ggl_N_HOD,(void*) array, 1./(1+tomo.cluster_zmax[nzc]),1./(1+tomo.cluster_zmin[nzc]), NULL,1000);
}

double project_tri_cl_cgl_HOD(double a,void *params)
{
  double k1,k2,res =0.,wa,weights;
  double *ar = (double *) params;
  wa = f_K(chi(a));
  k1 = ar[0]/wa;
  k2 = ar[1]/wa;
	weights = W_rm(a,ar[2]);
  weights *= W_rm(a,ar[2]);
  weights *= W_cluster(a,(int)ar[3],(int)ar[4]); //normalized redshift distribution of lensing cluster bin
  weights *= W_kappa(a,wa,ar[5]);
  weights *= dchi_da(a);
  if (weights > 0){
    //res = (b_cluster((int)ar[3],(int)ar[4])*tri_multih_cov(k1,k2,a)+tri_1h_mmcm(k1,k2,a,(int)ar[3],(int)ar[4]))*pow(wa,-6.)/(survey.area*survey.area_conversion_factor);
    res = (b_cluster((int)ar[3],(int)ar[4])*tri_matter_cov(k1,k2,a))*pow(wa,-6.)/(survey.area*survey.area_conversion_factor);
    res += delPgg_SSC(k1, a)*delP_SSC_cm(k2, a, ar[3], ar[4])*survey_variance(a,survey.area/41253.0)*pow(wa,-4.);
  }
  return weights*res;
}


double cov_NG_cl_cgl_HOD (double l1,double l2, int nzl1,int nzl2, int nzc, int nN,int nzs){
  if (nzl1 != nzl2){return 0.;}
  if ((tomo.clustering_zmin[nzl1]+tomo.clustering_zmax[nzl1])/2. > tomo.cluster_zmax[nzc]){return 0.;}
  if ((tomo.clustering_zmin[nzl1]+tomo.clustering_zmax[nzl1])/2. < tomo.cluster_zmin[nzc]){return 0.;}
  double array[6] = {l1,l2,(double)nzl1, (double) nzc, (double) nN, (double) nzs};
  return int_gsl_integrate_medium_precision(project_tri_cl_cgl_HOD,(void*) array, 1./(1+tomo.cluster_zmax[nzc]),1./(1+tomo.cluster_zmin[nzc]), NULL,1000);
}
double int_for_C_cl_g_rm_tomo(double a, void *params){
  double *ar = (double *) params;
  double res,ell, fK, k;
  
  ell       = ar[3]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  res= b_cluster((int)(ar[0]),(int)(ar[1]))*W_cluster(a,ar[0],ar[1])*W_rm(a,ar[2])*dchi_da(a)/fK/fK;
  res= res*p_lin(k,a);
  return res;
}

double C_cl_g_rm_tomo(double l, int nzc1, int nN1, int nzl1){
  double array[4] ={(double)nzc1, (double)nN1,  (double)nzl1, l};
  return int_gsl_integrate_low_precision(int_for_C_cl_g_rm_tomo, (void*) array, 1./(1+fmin(tomo.cluster_zmax[nzc1],tomo.clustering_zmax[nzl1])),1./(1+fmax(tomo.cluster_zmin[nzc1],tomo.clustering_zmin[nzl1])), NULL,1000);
}

double cov_G_cl_cgl_HOD (double l, double delta_l,int nzl1,int nzl2, int nzc, int nN,int nzs){
  double C13, C14, C23, C24;
   double fsky = survey.area/41253.0;
  C13 = C_cl_g_rm_tomo(l,nzc, nN, nzl1);C24 = C_gl_HOD_rm_tomo(l,nzl2,nzs);
  C23 = C_cl_g_rm_tomo(l,nzc, nN, nzl2);C14 = C_gl_HOD_rm_tomo(l,nzl1,nzs);
  
  return (C13*C24+  C14*C23)/((2.*l+1.)*delta_l*fsky);
}

double cov_G_ggl_cgl_HOD (double l, double delta_l,int nzl,int nzs1, int nzc, int nN,int nzs2){
  double C13, C14, C23, C24, N13 =0, N14=0, N23=0, N24=0;
   double fsky = survey.area/41253.0;
  C13 = C_cl_g_tomo(l,nzc,nN,nzl);C24 = C_shear_tomo_nointerp(l,nzs1,nzs2);
  C14 = C_gl_HOD_rm_tomo(l,nzl,nzs2);C23 = C_cgl_tomo_nointerp(l,nzc,nN,nzs1);
  if (nzs1 == nzs2){N24=pow(survey.sigma_e,2.0)/(2.0*nsource(nzs1)*survey.n_gal_conversion_factor);}
  
  return (C13*C24+ C13*N24+N13*C24 + C14*C23+C14*N23+N14*C23+N13*N24+N14*N23)/((2.*l+1.)*delta_l*fsky);
  
}
double project_tri_ggl_cgl_HOD(double a,void *params)
{
  double k1,k2,res = 0.,wa,weights;
  double *ar = (double *) params;
  wa = f_K(chi(a));
  k1 = ar[0]/wa;
  k2 = ar[1]/wa;
	weights = W_rm(a,ar[2]);
  weights *= W_kappa(a,wa,ar[3]);
  weights *= W_cluster(a,(int)ar[4],(int)ar[5]); //normalized redshift distribution of lensing cluster bin
  weights *= W_kappa(a,wa,ar[6]);
  weights *= dchi_da(a);
  if (weights > 0){
    //res = (b_cluster((int)ar[4],(int)ar[5])*tri_multih_cov(k1,k2,a)+tri_1h_mmcm(k1,k2,a,(int)ar[4],(int)ar[5]))*pow(wa,-6.)/(survey.area*survey.area_conversion_factor);
    res = (b_cluster((int)ar[4],(int)ar[5])*tri_matter_cov(k1,k2,a))*pow(wa,-6.)/(survey.area*survey.area_conversion_factor);
    res += delPgm_SSC(k1, a)*delP_SSC_cm(k2, a, ar[4], ar[5])*survey_variance(a,survey.area/41253.0)*pow(wa,-4.);
  }
  return weights*res;
}


double cov_NG_ggl_cgl_HOD (double l1,double l2,int nzl,int nzs1, int nzc, int nN,int nzs2){
  if ((tomo.clustering_zmin[nzl]+tomo.clustering_zmax[nzl])/2. > tomo.cluster_zmax[nzc]){return 0.;}
  if ((tomo.clustering_zmin[nzl]+tomo.clustering_zmax[nzl])/2. < tomo.cluster_zmin[nzc]){return 0.;}
  double array[7] = {l1,l2,(double)nzl, (double) nzs1, (double) nzc, (double) nN, (double) nzs2};
  return int_gsl_integrate_medium_precision(project_tri_ggl_cgl_HOD,(void*) array, 1./(1+tomo.cluster_zmax[nzc]),1./(1+tomo.cluster_zmin[nzc]), NULL,1000);
}
