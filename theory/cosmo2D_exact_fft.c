double C_cl_non_Limber(int l, int ni, int nj); //includes RSD
double C_cl_RSD(int l, int ni, int nj); //C_cl_Limber + non-Limber RSD terms only
double w_tomo_nonLimber(int nt, int ni, int nj); //w(theta) including non-Limber+RSD
void C_cl_mixed(int L, int LMAX, int ni, int nj, double *Cl, double dev, double tolerance);
void f_chi_for_Psi_cl(double* chi_ar, int Nchi, double* f_chi_ar, int ni);

double G_taper(double k){
	double s_bao = 5.5/cosmology.coverH0;
	return exp(-k*k*s_bao*s_bao);
}
double f_growth(double z){
	double aa = 1./(1+z);
	double gamma = 0.55;
	return pow(cosmology.Omega_m /(cosmology.Omega_m +omv_vareos(aa) *aa*aa*aa),gamma);
}

double Psi_RSD_z(double z,void *params){
	double *ar = (double *) params;
	int l = (int) ar[0];
	//if (l>50){return 0.;}
	double x = ar[1]*f_K(chi(1./(1.+z)));
	double WRSD =0.;
	if (x < 0.1*sqrt(2.*l) && (2*l+5) < GSL_SF_DOUBLEFACT_NMAX){ //small-x limit for j_l(x) in Eq. 27 of https://arxiv.org/pdf/astro-ph/0605302v2.pdf
		WRSD =(2.*l*l+2.*l-1.)/((2.*l+3.)*(2.*l-1))*pow(x,l)/gsl_sf_doublefact((unsigned int)(2*l+1))*(1-0.5*x*x/(2.*l+3));
		WRSD -= (l+1.)*(l+2.)/((2.*l+1.)*(2.*l+3.))*pow(x,l+2.)/gsl_sf_doublefact((unsigned int)(2*l+5))*(1-0.5*x*x/(2.*l+7));
		if (l>1){WRSD -= l*(l-1.)/((2.*l-1.)*(2.*l+1.))*pow(x,l-2.)/gsl_sf_doublefact((unsigned int)(2*l-3))*(1-0.5*x*x/(2.*l-1));}
	}
	else{ // Eq. 27 of https://arxiv.org/pdf/astro-ph/0605302v2.pdf
		WRSD =(2.*l*l+2.*l-1.)/((2.*l+3.)*(2.*l-1))*gsl_sf_bessel_jl(l,x);
		WRSD -= (l+1.)*(l+2.)/((2.*l+1.)*(2.*l+3.))*gsl_sf_bessel_jl(l+2,x);
		if (l>1){WRSD -= l*(l-1.)/((2.*l-1.)*(2.*l+1.))*gsl_sf_bessel_jl(l-2,x);}
	}
	return WRSD*f_growth(z);
}
double int_for_Psi_RSD(double z, void *params){
	static double g0 = 0.;
	if (g0 ==0){g0 =1./growfac(1.);}
	double *ar = (double *) params;
	return Psi_RSD_z(z,params)*pf_photoz(z,(int)ar[2])*growfac(1./(1+z))*g0;
}
double int_for_Psi_cl (double z, void *params){
	static double g0 = 0.;
	if (g0 ==0){g0 =1./growfac(1.);}
	double *ar = (double *) params;
	int l = (int) ar[0];
	double res, x = ar[1]*f_K(chi(1./(1.+z)));
	if (x < 0.1*sqrt(2.*l) && (2*l+1) < GSL_SF_DOUBLEFACT_NMAX){ // small-x limit for j_l(x), to avoid gsl underflow error
		res =pow(x,l)/gsl_sf_doublefact((unsigned int)(2*l+1))*(1-0.5*x*x/(2.*l+3));
	}
	else if (x > 100.*l*l){// large-x limit for j_l(x), to avoid gsl error
		res = sin(x-l*M_PI*0.5)/x;
	}
	else{
		 res = gsl_sf_bessel_jl(l,x);
	}
	return pf_photoz(z,(int)ar[2])*growfac(1./(1+z))*g0*(res*gbias.b1_function(z,(int)ar[2])+Psi_RSD_z(z,params));
}


double Psi_cl(double k, int l,int ni){
	double ar[3] = {(double)l,k,(double) ni};
	return int_gsl_integrate_low_precision(int_for_Psi_cl,(void*)ar,tomo.clustering_zmin[ni],tomo.clustering_zmax[ni],NULL,100);
}
double Psi_RSD(double k, int l,int ni){
	double ar[3] = {(double)l,k,(double) ni};
	return int_gsl_integrate_low_precision(int_for_Psi_RSD,(void*)ar,tomo.clustering_zmin[ni],tomo.clustering_zmax[ni],NULL,100);
}

double int_for_C_cl_nonLimber (double lk, void *params){
	double k = exp(lk);
	int *ar = (int *) params;
	return k*k*k*pow(Psi_cl(k,ar[0],ar[1]),2.)*p_lin(k,1.0)*G_taper(k);
}
double int_for_C_cl_RSD (double klog, void *params){
	int *ar = (int *) params;
	double k = exp(klog);
	return k*k*k*pow(Psi_RSD(k,ar[0],ar[1]),2.)*p_lin(k,1.0)*G_taper(k);
}
double int_for_C_cl_nonLimber_tomo (double klog, void *params){
	int *ar = (int *) params;
	double k = exp(klog);
	return k*k*k*Psi_cl(k,ar[0],ar[1])*Psi_cl(k,ar[0],ar[2])*p_lin(k,1.0)*G_taper(k);
}

double int_for_C_cl_lin(double a, void *params)
{
	double res,ell, fK, k;
	double *ar = (double *) params;
	ell       = ar[2]+0.5;
	fK     = f_K(chi(a));
	k      = ell/fK;
	
	res=W_gal(a,ar[0])*W_gal(a,ar[1])*dchi_da(a)/fK/fK;
	res= res*p_lin(k,a)*G_taper(k);
	return res;
}


double C_cl_lin_nointerp(double l, int ni, int nj)  //galaxy clustering power spectrum of galaxy bins ni, nj
{
	double array[3] = {1.0*ni,1.0*nj,l};
	return int_gsl_integrate_medium_precision(int_for_C_cl_lin,(void*)array,fmax(amin_lens(ni),amin_lens(nj)),fmin(amax_lens(ni),amax_lens(nj)),NULL,1000);
}


double C_cl_non_Limber(int l, int ni, int nj){ //includes RSD too!
	int ar[3] ={l,ni,nj};
	double res;
	if (2*l+5 > GSL_SF_DOUBLEFACT_NMAX){
		printf("cosmo2D_exact.c:C_cl_non_Limber: l = %d too large for stable evaluation of Bessel functions\nSwitching to Limber approximation\n", l);
		return C_cl_tomo_nointerp(1.*l,ni,nj);
	}
	//gsl_set_error_handler_off ();
	if (ni == nj){
		//checked that these boundaries give better than 1% accuracy for l > 2
		double kmin = fmax(limits.k_min_cH0, 0.25*l/f_K(chi(1./(1.+tomo.clustering_zmax[ni]))));
		double kmax = fmin(limits.k_max_cH0, 20.*l/f_K(chi(1./(1.+tomo.clustering_zmin[ni]))));
		res = 2.*int_gsl_integrate_low_precision(int_for_C_cl_nonLimber,(void*)ar,log(kmin),log(kmax),NULL,100)/M_PI;
	}
	else res = 2.*int_gsl_integrate_low_precision(int_for_C_cl_nonLimber_tomo,(void*)ar,log(limits.k_min_cH0),log(limits.k_max_cH0),NULL,1000)/M_PI;
	//gsl_set_error_handler (NULL);

	// non add non-linear power spectrum contributions in Limber approximation
	res = res + C_cl_tomo_nointerp(1.*l,ni,nj) - C_cl_lin_nointerp(1.*l,ni,nj);
	return res;
}


double C_cl_RSD(int l, int ni, int nj){
	int ar[3] ={l,ni,nj};
	double res;
	if (2*l+5 > GSL_SF_DOUBLEFACT_NMAX){
		printf("cosmo2D_exact.c:C_cl_RSD: l = %d too large for stable evaluation of Bessel functions\nSwitching to Limber approximation\n", l);
		return C_cl_tomo_nointerp(1.*l,ni,nj);
	}
	gsl_set_error_handler_off ();
	if (ni == nj){
		double kmin = fmax(limits.k_min_cH0, 0.25*l/f_K(chi(1./(1.+tomo.clustering_zmax[ni]))));
		double kmax = fmin(limits.k_max_cH0, 20.*l/f_K(chi(1./(1.+tomo.clustering_zmin[ni]))));
		res = 2.*int_gsl_integrate_low_precision(int_for_C_cl_RSD,(void*)ar,log(kmin),log(kmax),NULL,1000)/M_PI;
	}
	else res = 0.;
	gsl_set_error_handler (NULL);
	//add regular galaxy power spectrum contributions in Limber approximation
	res = res + C_cl_tomo_nointerp(1.*l,ni,nj) - C_cl_lin_nointerp(1.*l,ni,nj);
	return res;
}



///////
void f_chi_for_Psi_cl(double* chi_ar, int Nchi, double* f_chi_ar, int ni){
	double g0 =1./growfac(1.);
	double a, z;
	int i;
	double real_coverH0 = cosmology.coverH0 / cosmology.h0;
	double pf;
	for(i=0;i<Nchi;i++) {
		a = a_chi(chi_ar[i] / real_coverH0) ; // first convert unit of chi from Mpc to c/H0
		z = 1./a - 1.;
		pf = (pf_photoz(z,ni)<0.)? 0:pf_photoz(z,ni); // get rid of unphysical negatives
		// pf =pf_photoz(z,ni);
		f_chi_ar[i] = chi_ar[i] * pf*growfac(a)*g0*gbias.b1_function(z,ni)*hoverh0(a)/real_coverH0;
		// if(f_chi_ar[i]<1.) {f_chi_ar[i]=0.;}
		// printf("%lf, %lf, %lf, %lf, %lf \n", f_chi_ar[i],pf,growfac(a),g0,gbias.b1_function(z,ni));
	}
}

///////
void f_chi_for_Psi_cl_RSD(double* chi_ar, int Nchi, double* f_chi_RSD_ar, int ni){
	double g0 =1./growfac(1.);
	double a, z;
	int i;
	double real_coverH0 = cosmology.coverH0 / cosmology.h0;
	double pf;
	for(i=0;i<Nchi;i++) {
		a = a_chi(chi_ar[i] / real_coverH0) ; // first convert unit of chi from Mpc to c/H0
		z = 1./a - 1.;
		pf = (pf_photoz(z,ni)<0.)? 0:pf_photoz(z,ni); // get rid of unphysical negatives
		// pf =pf_photoz(z,ni);
		f_chi_RSD_ar[i] = -chi_ar[i] * pf*growfac(a)*g0*f_growth(z)*hoverh0(a)/real_coverH0;
		// if(f_chi_ar[i]<1.) {f_chi_ar[i]=0.;}
		// printf("%lf, %lf, %lf, %lf, %lf \n", f_chi_ar[i],pf,growfac(a),g0,gbias.b1_function(z,ni));
	}
}

// Mixture of non-Limber and Limber
void C_cl_mixed(int L, int LMAX, int ni, int nj, double *Cl, double dev, double tolerance) {
	ni = 4;
	nj = 4;
	int i,j,i_block;
	long l;
	// run 100 ells at a time, and see if switching to Limber is needed.
	// Save runtime for Limber, and save re-creation time of fftw_plan.
	int Nell_block = 100, Nchi = 1000;
	int ell_ar[Nell_block];
	double **k1_ar, **k2_ar, **Fk1_ar, **Fk2_ar;
	double **Fk1_RSD_ar, **Fk2_RSD_ar;
	// ell_ar = malloc(Nell_block * sizeof(int));
	k1_ar = malloc(Nell_block * sizeof(double *));
	k2_ar = malloc(Nell_block * sizeof(double *));
	Fk1_ar = malloc(Nell_block * sizeof(double *));
	Fk2_ar = malloc(Nell_block * sizeof(double *));
	Fk1_RSD_ar = malloc(Nell_block * sizeof(double *));
	Fk2_RSD_ar = malloc(Nell_block * sizeof(double *));
	for(i=0;i<Nell_block;i++) {
		k1_ar[i] = malloc(Nchi * sizeof(double));
		k2_ar[i] = malloc(Nchi * sizeof(double));
		Fk1_ar[i] = malloc(Nchi * sizeof(double));
		Fk2_ar[i] = malloc(Nchi * sizeof(double));
		Fk1_RSD_ar[i] = malloc(Nchi * sizeof(double));
		Fk2_RSD_ar[i] = malloc(Nchi * sizeof(double));
	}

	double chi_ar[Nchi], f1_chi_ar[Nchi], f2_chi_ar[Nchi];
	double f1_chi_RSD_ar[Nchi], f2_chi_RSD_ar[Nchi];
	double chi_min = 60., chi_max = 6000.;
	double dlnchi = log(chi_max/chi_min) / (Nchi - 1.);
	double dlnk = dlnchi;

	for(i=0; i<Nchi; i++) {
		chi_ar[i] = chi_min * exp(dlnchi*i);
	}

	f_chi_for_Psi_cl(chi_ar, Nchi, f1_chi_ar, ni);
	if(ni != nj) {f_chi_for_Psi_cl(chi_ar, Nchi, f2_chi_ar, nj);}

	f_chi_for_Psi_cl_RSD(chi_ar, Nchi, f1_chi_RSD_ar, ni);
	if(ni != nj) {f_chi_for_Psi_cl_RSD(chi_ar, Nchi, f2_chi_RSD_ar, nj);}

	// char outfilename[] = "f1_chi4.txt";
	// char outfilename[] = "f1_chi4_rsd.txt";
	// FILE *OUT = fopen(outfilename, "w");
	
	// for(i=0; i<Nchi; i++) {
	// 	// fprintf(OUT, "%lg %lg", chi_ar[i], f1_chi_ar[i]);
	// 	fprintf(OUT, "%lg %lg", chi_ar[i], f1_chi_RSD_ar[i]);
	// 	fprintf(OUT, "\n");
	// }
	// for(i=0; i<Nchi; i++) {
	// 	printf("f_chi_ar: %d, %lf\n", i, f1_chi_ar[i]);
	// }
	// exit(0);
	char outfilename[] = "c_cl4_rsd.txt";
	FILE *OUT = fopen(outfilename, "w");


	i_block = 0;
	double cl_temp;

	config my_config, my_config_RSD;
	my_config.nu = 1.;
	my_config.c_window_width = 0.25;
	my_config.derivative = 0;
	my_config_RSD.nu = 1.1;
	my_config_RSD.c_window_width = 0.25;
	my_config_RSD.derivative = 2;

	double real_coverH0 = cosmology.coverH0 / cosmology.h0;
	double k1_cH0;

	while (fabs(dev) > tolerance){
		//Cl[L] = C_cl_RSD(L,nz,nz);
		for(i=0;i<Nell_block;i++) {ell_ar[i]=i+i_block*Nell_block;}

		cfftlog_ells(chi_ar, f1_chi_ar, Nchi, &my_config, ell_ar, Nell_block, k1_ar, Fk1_ar);
		if(ni != nj) {cfftlog_ells(chi_ar, f2_chi_ar, Nchi, &my_config, ell_ar, Nell_block, k2_ar, Fk2_ar);}

		cfftlog_ells(chi_ar, f1_chi_RSD_ar, Nchi, &my_config_RSD, ell_ar, Nell_block, k1_ar, Fk1_RSD_ar);
		if(ni != nj) {cfftlog_ells(chi_ar, f2_chi_RSD_ar, Nchi, &my_config_RSD, ell_ar, Nell_block, k2_ar, Fk2_RSD_ar);}

		for(i=0;i<Nell_block;i++) {
			cl_temp = 0.;
			for(j=0;j<Nchi;j++) {
				// printf("k,Fk: %d,%d, %lf,%lf\n", i,j, k1_ar[i][j], Fk1_ar[i][j]);
				k1_cH0 = k1_ar[i][j] * real_coverH0;
				if(ni == nj) {
					// cl_temp += (Fk1_ar[i][j])*(Fk1_ar[i][j]) *k1_cH0*k1_cH0*k1_cH0 *p_lin(k1_cH0,1.0)*G_taper(k1_cH0);
					cl_temp += (Fk1_ar[i][j]+Fk1_RSD_ar[i][j]) * (Fk1_ar[i][j]+Fk1_RSD_ar[i][j]) *k1_cH0*k1_cH0*k1_cH0 *p_lin(k1_cH0,1.0)*G_taper(k1_cH0);
					// printf("plin,%lg, %lg\n", k1_ar[i][j],p_lin(k1_cH0,1.0));
				}
				else {
					// cl_temp += (Fk1_ar[i][j])*(Fk2_ar[i][j]) *k1_cH0*k1_cH0*k1_cH0 *p_lin(k1_cH0,1.0)*G_taper(k1_cH0);
					cl_temp += (Fk1_ar[i][j]+Fk1_RSD_ar[i][j]) * (Fk2_ar[i][j]+Fk2_RSD_ar[i][j]) *k1_cH0*k1_cH0*k1_cH0 *p_lin(k1_cH0,1.0)*G_taper(k1_cH0);
				}
			}
			Cl[ell_ar[i]] = cl_temp * dlnk * 2./M_PI + C_cl_tomo_nointerp(1.*ell_ar[i],ni,nj) - C_cl_lin_nointerp(1.*ell_ar[i],ni,nj);
			// printf("cl_t/emp: %d, %lg\n", i, cl_temp);
			fprintf(OUT, "%d %lg %lg %lg", ell_ar[i], Cl[ell_ar[i]], C_cl_tomo_nointerp(1.*ell_ar[i],ni,nj), C_cl_lin_nointerp(1.*ell_ar[i],ni,nj));
			fprintf(OUT, "\n");
		}
		// if(i_block==9) {exit(0);}

		i_block++;
		L = i_block*Nell_block -1 ;
		dev = Cl[L]/C_cl_tomo_nointerp((double)L,ni,nj)-1.;
	   // printf("ni,L,Cl[L],dev=%d %d %e %e\n",ni,L,Cl[L],dev);
		// printf("i_block: %d\n", i_block);

	}
	L++;
	printf("switching to Limber calculation at l = %d\n",L);
	for (l = L; l < LMAX; l++){
		Cl[l]=C_cl_tomo((double)l,ni,nj);
	}
	fclose(OUT);
	exit(0);
}

double w_tomo_nonLimber(int nt, int ni, int nj){

	static int LMAX = 100000;
	static int NTHETA = 0;
	static double ** Pl =0;
	static double *Cl =0;
	static double *w_vec =0;
	static cosmopara C;
	static nuisancepara N;
	static galpara G;
	int i,l,nz;
	if (like.Ntheta ==0){
		printf("cosmo2D_real.c:w_tomo_exact: like.Ntheta not initialized\nEXIT\n"); exit(1);
	}
	if (ni != nj){
		printf("cosmo2D_real.c:w_tomo_exact: ni != nj tomography not supported\nEXIT\n"); exit(1);    
	}
	if (Pl ==0){
		Pl =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
		Cl = create_double_vector(0,LMAX-1);
		w_vec = create_double_vector(0,tomo.clustering_Nbin*like.Ntheta-1);
		NTHETA = like.Ntheta;
		double *xmin, *xmax, *Pmin, *Pmax;
		xmin= create_double_vector(0, like.Ntheta-1);
		xmax= create_double_vector(0, like.Ntheta-1);
		double logdt=(log(like.vtmax)-log(like.vtmin))/like.Ntheta;
		Pmin= create_double_vector(0, LMAX+1);
		Pmax= create_double_vector(0, LMAX+1);
		for(i=0; i<like.Ntheta ; i++){
			xmin[i]=cos(exp(log(like.vtmin)+(i+0.0)*logdt));
			xmax[i]=cos(exp(log(like.vtmin)+(i+1.0)*logdt));
		}

		for (i = 0; i<NTHETA; i ++){
			printf("Tabulating Legendre coefficients %d/%d\n",i+1, NTHETA);
			gsl_sf_legendre_Pl_array(LMAX, xmin[i],Pmin);
			gsl_sf_legendre_Pl_array(LMAX, xmax[i],Pmax);
			for (int l = 1; l < LMAX; l ++){
				Pl[i][l] = (2*l+1.)/(4.*M_PI)*Pmin[l];
				//Pl[i][l] = 1./(4.*M_PI)*(Pmin[l+1]-Pmax[l+1]-Pmin[l-1]+Pmax[l-1])/(xmin[i]-xmax[i]);
			}
		}
		free_double_vector(xmin,0,like.Ntheta-1);
		free_double_vector(xmax,0,like.Ntheta-1);
		free_double_vector(Pmin,0,LMAX+1);
		free_double_vector(Pmax,0,LMAX+1);
		}
		if (recompute_clustering(C,G,N,ni,nj)){
		//required fractional accuracy in C(l)
		double tolerance= 0.01;
		//dev will be the actual difference between exact and Limber calcuation
		double dev;

		for (nz = 0; nz <tomo.clustering_Nbin; nz ++){
			int L = 1;
			// initialize to large value in order to start while loop
			dev=10.*tolerance;
		//   while (fabs(dev) > tolerance){
		//     //Cl[L] = C_cl_RSD(L,nz,nz);
		//     Cl[L] = C_cl_non_Limber(L,nz,nz);
		//     dev = Cl[L]/C_cl_tomo_nointerp((double)L,nz,nz)-1.;
		// //    printf("nL, nz=%d %d %e %e\n",nz,L,Cl[L],dev);
		//     L = L+1;
		//   }
		//   printf("switching to Limber calculation at l = %d\n",L);
		//   for (l = L; l < LMAX; l++){
		//     Cl[l]=C_cl_tomo((double)l,nz,nz);
		//   }
			C_cl_mixed(L, LMAX, nz,nz, Cl, dev, tolerance);
			//////////////
			for (i = 0; i < NTHETA; i++){
				w_vec[nz*like.Ntheta+i] =0;
				for (l = 1; l < LMAX; l++){
					w_vec[nz*like.Ntheta+i]+=Pl[i][l]*Cl[l];
				}
			}
		}
		update_cosmopara(&C);
		update_galpara(&G);
		update_nuisance(&N);
	}
	return w_vec[ni*like.Ntheta+nt];  
}

