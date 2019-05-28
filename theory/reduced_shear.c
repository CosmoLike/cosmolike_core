double Delta_C_shear_reduced_shear_nointerp(double l, int ni, int nj);
double Delta_C_reduced_shear_tomo(double l, int ni, int nj);
double Delta_C_shear_gl_shear_nointerp(double l, int nl, int ns);
double Delta_C_gl_reduced_shear_tomo(double l, int ni, int nj);
double Delta_xi_pm_reduced_shear_tomo(int pm, double theta, int ni, int nj);
double Delta_w_gamma_t_reduced_shear_tomo(double theta,int ni, int nj);
//typedef double (*C_tomo_pointer)(double l, int n1, int n2);

double fs_2(double k1x, double k1y, double k2x, double k2y)
{
  double k1,k2,res;
  k1 = sqrt(k1x*k1x+k1y*k1y);
  k2 = sqrt(k2x*k2x+k2y*k2y);
  res = 5./7.+2./7.*pow((k1x*k2x+k1y*k2y)/(k1*k2),2.)+.5*((k1x*k2x+k1y*k2y)/(k1*k2))*(k1/k2+k2/k1);
  if(isnan(res)) {res=0.0;}
  //  printf("fs_2 %le\n",res);
  return res;
}

double bi_tree_nl (double k1x, double k1y, double k2x, double k2y, double a)
{
	double k1,k2,k3,k3x,k3y;
	k3x = -k1x-k2x;
	k3y = -k1y-k2y;
	k1 = sqrt(k1x*k1x+k1y*k1y);
	k2 = sqrt(k2x*k2x+k2y*k2y);
	k3 = sqrt(k3x*k3x+k3y*k3y);
	double res = fs_2(k1x,k1y,k2x,k2y)*Pdelta(k1,a)*Pdelta(k2,a);
	//if (k3 > limits.k_min_cH0 && k3 < limits.k_max_cH0){
		res += fs_2(k2x,k2y,k3x,k3y)*Pdelta(k2,a)*Pdelta(k3,a)+fs_2(k3x,k3y,k1x,k1y)*Pdelta(k3,a)*Pdelta(k1,a);
	//}
	return 2.*res;
}


double bi_eff_mu(double k1,double k2, double mu, double a){
	double bi_eff = bi_tree_nl(k1,0.,k2*mu,k2*sqrt(1.-mu*mu),a);
 // bi_eff += bi_1h(k1,k2,sqrt(k1*k1 +2*k1*k2*mu*mu + k2*k2),a);
	return bi_eff;
}

double int_bi_reduced_shear_phi(double phi, void *params){
	double *ar = (double *) params;
	double mu = cos(phi);
	return cos(2.*phi)*bi_eff_mu(ar[0],ar[1],cos(phi),ar[2]);
}
double bi_3D_reduced_shear_angle_weighted(double k1, double k2, double a){ //linear bispectrum, averaged over angle between k1 and k2
 	double array[3] = {k1,k2,a};
 	double eps = 0.0;
	return int_gsl_integrate_low_precision(int_bi_reduced_shear_phi,(void*)array,0.+eps,M_PI-eps,NULL,1000)/M_PI;
}

double int_Delta_P_reduced_shear_dk2(double logk2, void *params){
	double *ar = (double *) params;
	double k2 = exp(logk2);
	return k2*k2*bi_3D_reduced_shear_angle_weighted(ar[0],k2,ar[1]);
}
// perform integral over d^2l' before Limber projection
double Delta_P_reduced_shear(double k1, double a){
	double array[2] = {k1,a};
	return 2.*pow(f_K(chi(a)),2.)*int_gsl_integrate_low_precision(int_Delta_P_reduced_shear_dk2, (void*)array,log(limits.k_min_cH0),log(limits.k_max_cH0), NULL, 1000)/(2.*M_PI);
}

double Delta_P_reduced_shear_tab(double k, double a){
	static cosmopara C;
 	static double amin = 0, amax = 0,logkmin = 0., logkmax = 0., dk = 0., da = 0.;
  static int NA = 20, NK =50;
  	static double **table_P=0;
  
	double aa,klog;
	int i,j;
  
	if (recompute_cosmo3D(C)){
    	update_cosmopara(&C);
		if (table_P==0) table_P = create_double_matrix(0, NA-1, 0, NK-1);

		amin = 1./(1+tomo.shear_zmax[tomo.shear_Nbin -1]);
		amax = 1.;
		da = (amax - amin)/(NA-1);
		aa = amin;
		logkmin = log(limits.k_min_cH0);
		logkmax = log(limits.k_max_cH0);
		dk = (logkmax - logkmin)/(NK-1);
		for (i=0; i<NA; i++, aa +=da) {
			klog  = logkmin;
			for (j=0; j<NK; j++, klog += dk) {
				table_P[i][j] = log(Delta_P_reduced_shear(exp(klog),aa));
			}
		}
 	}
	aa = fmin(a,amax-1.1*da);//to avoid interpolation errors near z=0
	return exp(interpol2d(table_P, NA, amin, amax, da, aa, NK, logkmin, logkmax, dk, log(k), 0.0, 0.0));
}

double int_for_Delta_C_shear_reduced_shear(double a, void *params){
  double *ar = (double *) params;
  double res,ell, fK, k,w1,w2;
  if (a >= 1.0) error("a>=1 in int_for_C_shear_tomo");
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  w1 = W_kappa(a,fK,ar[0]);
  w2 = W_kappa(a,fK,ar[1]);
  res= 0.5*w1*w2*(w1+w2)*dchi_da(a)*pow(fK,-4.);
  res *= Delta_P_reduced_shear_tab(k,a); 
  return res;
}

double Delta_C_shear_reduced_shear_nointerp(double l, int ni, int nj){
  double array[3] = {(double) ni, (double) nj,l};
  int j,k;
  if (ni <= nj){j =nj; k = ni;}
  else{j = ni; k = nj;}
  return int_gsl_integrate_low_precision(int_for_Delta_C_shear_reduced_shear,(void*)array,amin_source(j),amax_source(k),NULL,1000);
}

double Delta_C_reduced_shear_tomo(double l, int ni, int nj)  //shear power spectrum of source galaxies in bins ni, nj
{
  static cosmopara C;
  static nuisancepara N;
  
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  if (ni < 0 || ni >= tomo.shear_Nbin ||nj < 0 || nj >= tomo.shear_Nbin){
    printf("C_shear_tomo(l,%d,%d) outside tomo.shear_Nbin range\nEXIT\n",ni,nj); exit(1);
  }
  
  if (recompute_shear(C,N)){
    if (table==0) {
      table   = create_double_matrix(0, tomo.shear_Npowerspectra-1, 0, Ntable.N_ell-1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell-1.);
    }
    
    double llog;
    int i,k;
    
    for (k=0; k<tomo.shear_Npowerspectra; k++) {
      llog = logsmin;
      for (i=0; i<Ntable.N_ell; i++, llog+=ds) {
        table[k][i]= Delta_C_shear_reduced_shear_nointerp(exp(llog),Z1(k),Z2(k));
      }
    }
    update_cosmopara(&C); update_nuisance(&N); 
  }
  double f1 = interpol_fitslope(table[N_shear(ni,nj)], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1.);
  if (isnan(f1)){f1 = 0.;}
  return f1;
}


double int_for_Delta_C_gl_reduced_shear(double a, void *params){
  double *ar = (double *) params;
  double res,ell, fK, k,w1,w2;
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  w1 = W_gal(a,ar[0]);
  w2 = W_kappa(a,fK,ar[1]);
  res= w1*w2*w2*dchi_da(a)*pow(fK,-4.);
  res *= Delta_P_reduced_shear_tab(k,a); 
  return res;
}

double Delta_C_shear_gl_shear_nointerp(double l, int nl, int ns){
  double array[3] = {(double) nl, (double) ns,l};
  return int_gsl_integrate_low_precision(int_for_Delta_C_gl_reduced_shear,(void*)array,amin_lens(nl),amax_lens(nl),NULL,1000);
}

double Delta_C_gl_reduced_shear_tomo(double l, int ni, int nj)  //G-G lensing power spectrum, lens bin ni, source bin nj
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  double f1 = 0.;
  
  if (ni < 0 || ni >= tomo.clustering_Nbin ||nj < 0 || nj >= tomo.shear_Nbin){
    printf("C_gl_tomo(l,%d,%d) outside tomo.X_Nbin range\nEXIT\n",ni,nj); exit(1);
  }
  
  if (recompute_ggl(C,G,N,ni)){
    if (table==0){
      table   = create_double_matrix(0, tomo.ggl_Npowerspectra-1, 0, Ntable.N_ell-1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell-1.);
    }
    int i,k;
    double llog;
    
    for (k=0; k<tomo.ggl_Npowerspectra; k++) {
      llog = logsmin;
      for (i=0; i<Ntable.N_ell; i++, llog+=ds) {
          table[k][i]= Delta_C_shear_gl_shear_nointerp(exp(llog),ZL(k),ZS(k));
      }
    }
    
    update_cosmopara(&C); update_nuisance(&N); update_galpara(&G);
    
  }
  if(test_zoverlap(ni,nj)){f1 = interpol_fitslope(table[N_ggl(ni,nj)], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1.);}
  if (isnan(f1)){f1 = 0;}
  return f1;
}


double Delta_xi_pm_reduced_shear_tomo(int pm, double theta, int ni, int nj) //shear tomography correlation functions
{
  static cosmopara C;
  static nuisancepara N; 
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;
  if (recompute_shear(C,N)){
	C_tomo_pointer C_pointer = &C_reduced_shear_tomo;
    update_cosmopara(&C); update_nuisance(&N);
    double **tab;
    int i,k;
    tab   = create_double_matrix(0, 1, 0, Ntable.N_thetaH-1);
    if (table==0) table   = create_double_matrix(0, 2*tomo.shear_Npowerspectra-1, 0, Ntable.N_thetaH-1);
    for (i = 0; i < tomo.shear_Npowerspectra; i++){
    	xipm_via_hankel(tab, &logthetamin, &logthetamax,C_pointer, Z1(i),Z2(i));
      for (k = 0; k < Ntable.N_thetaH; k++){
        table[2*i][k] = tab[0][k];
          table[2*i+1][k] = tab[1][k];
      }
    }
    dlogtheta = (logthetamax-logthetamin)/((double)Ntable.N_thetaH);
    free_double_matrix(tab,0, 1, 0, Ntable.N_thetaH-1);
  }
  return interpol(table[2*N_shear(ni,nj)+(1-pm)/2], Ntable.N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 0.0, 0.0);
}

double Delta_w_gamma_t_reduced_shear_tomo(double theta,int ni, int nj) //G-G lensing, lens bin ni, source bin nj
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;
  double res =0.;
  if (recompute_ggl(C,G,N,ni)){
    C_tomo_pointer C_gl_pointer = &C_gl_reduced_shear_tomo;
    double **tab;
    int i, k;
    tab   = create_double_matrix(0, 1, 0, Ntable.N_thetaH-1);
    if (table==0) table  = create_double_matrix(0, tomo.ggl_Npowerspectra, 0, Ntable.N_thetaH-1);
    for (i = 0; i <tomo.ggl_Npowerspectra; i++){
      twopoint_via_hankel(tab, &logthetamin, &logthetamax,C_gl_pointer, ZL(i),ZS(i),2);
      for (k = 0; k < Ntable.N_thetaH; k++){table[i][k] = tab[0][k];}
    }
    free_double_matrix(tab,0, 1, 0, Ntable.N_thetaH);
    dlogtheta = (logthetamax-logthetamin)/((double)Ntable.N_thetaH);
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }
  if(test_zoverlap(ni,nj)) {
    res = interpol(table[N_ggl(ni,nj)], Ntable.N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 0.0, 0.0);}
  return res;
}

