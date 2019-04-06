double Delta_C_kappa_reduced_shear(double l, int ni, int nj);

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
	//  bi_eff += bi_1h(k1,k2,sqrt(k1*k1 +2*k1*k2*mu*mu + k2*k2),a);
	return bi_eff;
}

double int_bi_reduced_shear_phi(double phi, void *params){
	double *ar = (double *) params;
	double mu = cos(phi);
	return cos(2.*phi)*bi_eff_mu(ar[0],ar[1],cos(phi),ar[2]);
}
double bi_3D_reduced_shear_angle_weighted(double k1, double k2, double a){ //linear bispectrum, averaged over angle between k1 and k2
 	double array[3] = {k1,k2,a};
	return int_gsl_integrate_medium_precision(int_bi_reduced_shear_phi,(void*)array,0.,M_PI,NULL,1000)/M_PI;
}

double int_Delta_P_reduced_shear_dk2(double logk2, void *params){
	double *ar = (double *) params;
	double k2 = exp(logk2);
	return k2*bi_3D_reduced_shear_angle_weighted(ar[0],k2,ar[1]);
}
double Delta_P_reduced_shear(double k1, double a){
	double array[2] = {k1,a};
	return 2.*int_gsl_integrate_medium_precision(int_Delta_P_reduced_shear_dk2, (void*)array,log(limits.k_min_cH0),log(limits.k_max_cH0), NULL, 1000)(2.*M_PI);
}

double int_for_Delta_C_kappa_reduced_shear(double a, void *params){
  double *ar = (double *) params;
  double res,ell, fK, k,w1,w2;
  if (a >= 1.0) error("a>=1 in int_for_C_shear_tomo");
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  w1 = W_kappa(a,fK,ar[0]);
  w2 = W_kappa(a,fK,ar[1]);
  res= 0.5*w1*w2*(w1+w2)*dchi_da(a)*pow(fK,-4.);
  res *= Delta_P_reduced_shear(k,a); 
  return res;
}

double Delta_C_kappa_reduced_shear(double l, int ni, int nj){
  double array[3] = {(double) ni, (double) nj,l};
  int j,k;
  if (ni <= nj){j =nj; k = ni;}
  else{j = ni; k = nj;}
  return int_gsl_integrate_medium_precision(int_for_Delta_C_kappa_reduced_shear,(void*)array,amin_source(j),amax_source(k),NULL,1000);
}