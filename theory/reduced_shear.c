double Delta_P_reduced_shear_tab(double k, double a);
double Delta_CC_reduced_shear_tomo(double l, int ni, int nj);

//double Delta_C_shear_reduced_shear_nointerp(double l, int ni, int nj);
//double Delta_C_reduced_shear_tomo(double l, int ni, int nj);
//double Delta_C_shear_gl_shear_nointerp(double l, int nl, int ns);
//double Delta_C_gl_reduced_shear_tomo(double l, int ni, int nj);
//double Delta_xi_pm_reduced_shear_tomo(int pm, double theta, int ni, int nj);
//double Delta_w_gamma_t_reduced_shear_tomo(double theta,int ni, int nj);
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
	if (k3 > limits.k_min_cH0 && k3 < limits.k_max_cH0){
		res += fs_2(k2x,k2y,k3x,k3y)*Pdelta(k2,a)*Pdelta(k3,a)+fs_2(k3x,k3y,k1x,k1y)*Pdelta(k3,a)*Pdelta(k1,a);
	}
	return 2.*res;
}


double bi_eff_mu(double k1,double k2, double mu, double a){
	double bi_eff = bi_tree_nl(k1,0.,k2*mu,k2*sqrt(1.-mu*mu),a);
  //bi_eff += bi_1h(k1,k2,sqrt(k1*k1 +2*k1*k2*mu*mu + k2*k2),a);
	return bi_eff;
}

double int_bi_reduced_shear_phi(double phi, void *params){
	double *ar = (double *) params;
	double mu = cos(phi);
	return cos(2.*phi)*bi_eff_mu(ar[0],ar[1],cos(phi),ar[2]);
}
double bi_3D_reduced_shear_angle_weighted(double k1, double k2, double a){ //linear bispectrum, averaged over angle between k1 and k2
 	double array[3] = {k1,k2,a};
 	double eps = 0.00;
	return int_gsl_integrate_low_precision(int_bi_reduced_shear_phi,(void*)array,0.+eps,M_PI-eps,NULL,1000)/(M_PI);
}

double int_Delta_P_reduced_shear_dk2(double logl2, void *params){
	double *ar = (double *) params;
	double l2 = exp(logl2);
	return l2*l2*bi_3D_reduced_shear_angle_weighted(ar[0],(l2+0.5)/chi(ar[1]),ar[1]);
}
// perform integral over d^2l' before Limber projection
double Delta_P_reduced_shear(double k1, double a){
	double array[2] = {k1,a};
  double fk = chi(a);
	return int_gsl_integrate_low_precision(int_Delta_P_reduced_shear_dk2, (void*)array,log(limits.k_min_cH0*fk),log(limits.k_max_cH0/10.*fk), NULL, 1000)/(2.*M_PI);
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
		amax = .99;
		da = (amax - amin)/(NA-1);
		aa = amin;
		logkmin = log(limits.k_min_cH0);
		logkmax = log(limits.k_max_cH0/10.);
		dk = (logkmax - logkmin)/(NK-1);
		for (i=0; i<NA; i++, aa +=da) {
			klog  = logkmin;
			for (j=0; j<NK; j++, klog += dk) {
				table_P[i][j] = log(Delta_P_reduced_shear(exp(klog),aa));
//        printf("%e %e %e  %e\n",exp(klog)/cosmology.coverH0, aa, table_P[i][j],table_P[i][j]/pow(chi(aa),2.)/Pdelta(exp(klog),aa));
			}
		}
 	}
	aa = fmin(a,amax-1.1*da);//to avoid interpolation errors near z=0
	return exp(interpol2d(table_P, NA, amin, amax, da, aa, NK, logkmin, logkmax, dk, log(k), 0.0, 0.0));
}

double int_for_C_source_clustering(double a, void *params){
  double *ar = (double *) params;
  double res,ell, fK, k,w2,ws1;
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  w2 = W_source(a,ar[1]);
  ws1 =  W_source(a,ar[0]);
  res= ws1*w2*Pdelta(k,a)*dchi_da(a)*pow(fK,-2.);
  return res;
}

double C_source_clustering_no_interp(double l, int ni, int nj){
  double array[3] = {(double) ni, (double) nj,l};
  return int_gsl_integrate_low_precision(int_for_C_source_clustering,(void*)array,amin_source(ni),amax_source(ni),NULL,1000);
}

double int_for_C_source_lensing(double a, void *params){
  double *ar = (double *) params;
  double res,ell, fK, k,w2,ws1;
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  w2 = W_kappa(a,fK,ar[1]);
  ws1 =  W_source(a,ar[0]);
  res= ws1*w2*Pdelta(k,a)*dchi_da(a)*pow(fK,-2.);
  return res;
}

double C_source_lensing_no_interp(double l, int ni, int nj){
  double array[3] = {(double) ni, (double) nj,l};
  return int_gsl_integrate_low_precision(int_for_C_source_lensing,(void*)array,amin_source(ni),amax_source(ni),NULL,1000);
}
double C_source_lensing(double l, int ni, int nj)  //shear power spectrum of source galaxies in bins ni, nj
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
      table   = create_double_matrix(0, tomo.shear_Nbin*tomo.shear_Nbin, 0, Ntable.N_ell-1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell-1.);
    }

    double llog;
    int i,k;

    for (k=0; k<tomo.shear_Nbin; k++) {
      for(int j = 0; j < tomo.shear_Nbin; j++){
        llog = logsmin;
        for (i=0; i<Ntable.N_ell; i++, llog+=ds) {
          table[tomo.shear_Nbin*k+j][i]= C_source_lensing_no_interp(exp(llog),k,j);
        }
      }
    }
    update_cosmopara(&C); update_nuisance(&N);
  }
  double f1 = interpol_fitslope(table[tomo.shear_Nbin*ni+nj], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1.);
  if (isnan(f1)){f1 = 0.;}
  return f1;
}

double int_phi1_Delta_CC_shear(double phi1,void *params){
  double *ar = (double *) params;
  double l1 = ar[0];
  double l = ar[1];
  int ni = (int) ar[2];
  int nj = (int) ar[3];
  double cosphi1_l1_l = (cos(phi1)*l1-l)/sqrt(pow(cos(phi1)*l1-l,2.)+pow(sin(phi1)*l1,2.));
  double cos2phi_l1_l = 2.*cosphi1_l1_l*cosphi1_l1_l-1.;
  return cos2phi_l1_l*cos(2.*phi1)*C_source_lensing(sqrt(pow(cos(phi1)*l1-l,2.)+pow(sin(phi1)*l1,2.)),ni,nj);
}
double int_l1_Delta_CC_shear(double ll1,void *params){
  double l1 = exp(ll1);
  double *ar = (double *) params;
  double l = ar[1];
  int ni = (int) ar[2];
  int nj = (int) ar[3];
  ar[0] = l1;
  return l1*l1*C_source_lensing(l1,nj,ni)*int_gsl_integrate_low_precision(int_phi1_Delta_CC_shear,(void*)ar,0.,2.*M_PI,NULL,1000)/(2.*M_PI);
}

double Delta_CC_shear_reduced_shear_nointerp(double l,int n1,int n2){
  double ar[4] = {0.,l,1.*n1,1.*n2};
  return int_gsl_integrate_low_precision(int_l1_Delta_CC_shear,(void*)ar,log(limits.P_2_s_min),log(limits.P_2_s_max),NULL,1000)/(2.*M_PI);
}
double Delta_CC_reduced_shear_tomo(double l, int ni, int nj)  //shear power spectrum of source galaxies in bins ni, nj
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
      table   = create_double_matrix(0, tomo.shear_Npowerspectra, 0, Ntable.N_ell-1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell-1.);
    }

    double llog;
    int i,k;

    for (k=0; k<tomo.shear_Npowerspectra; k++) {
      printf("tabulating %d\n",k);
      llog = logsmin;
        for (i=0; i<Ntable.N_ell; i++, llog+=ds) {
          double Ckk= C_shear_tomo_nointerp(exp(llog),Z1(k),Z2(k));
          table[k][i]= Delta_CC_shear_reduced_shear_nointerp(exp(llog),Z1(k),Z2(k));
          printf("%e %d %d   %e %e    %e %e %e\n",exp(llog),Z1(k),Z2(k),table[k][i],
           Delta_CC_shear_reduced_shear_nointerp(exp(llog),Z2(k),Z1(k)),Ckk,C_source_lensing(exp(llog),Z1(k),Z2(k)),C_source_lensing(exp(llog),Z2(k),Z1(k)));
      //  }
      }
    }
    update_cosmopara(&C); update_nuisance(&N);
  }

  double f1 = interpol_fitslope(table[N_shear(ni,nj)], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1.);
  if (isnan(f1)){f1 = 0.;}
  return f1;
}
