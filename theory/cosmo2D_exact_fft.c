void f_chi_for_Psi_cl(double* chi_ar, int Nchi, double* f_chi_ar, int ni);
void f_chi_for_Psi_cl_RSD(double* chi_ar, int Nchi, double* f_chi_RSD_ar, int ni);
void f_chi_for_Psi_cl_Mag(double* chi_ar, int Nchi, double* f_chi_Mag_ar, int ni);
void C_cl_mixed(int L, int LMAX, int ni, int nj, double *Cl, double dev, double tolerance);
double w_tomo_nonLimber(int nt, int ni, int nj); //w(theta) including non-Limber+RSD


void f_chi_for_Psi_sh(double* chi_ar, int Nchi, double* f_chi_ar, int nz);
void f_chi_for_Psi_sh_IA(double* chi_ar, int Nchi, double* f_chi_IA_ar, int nz);
void C_gl_mixed(int L, int LMAX, int ni, int nj, double *Cl, double dev, double tolerance);
double w_gamma_t_nonLimber(int nt, int ni, int nj);


/////////////////////
double G_taper(double k){
	double s_bao = 5.5/cosmology.coverH0;
	return exp(-k*k*s_bao*s_bao);
}
double f_growth(double z){
	double aa = 1./(1+z);
	double gamma = 0.55;
	return pow(cosmology.Omega_m /(cosmology.Omega_m +omv_vareos(aa) *aa*aa*aa),gamma);
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

/////// Integrand for galaxy density
void f_chi_for_Psi_cl(double* chi_ar, int Nchi, double* f_chi_ar, int ni){
	double g0 =1./growfac(1.);
	double a, z;
	int i;
	double real_coverH0 = cosmology.coverH0 / cosmology.h0; // unit Mpc
	double pf;
	for(i=0;i<Nchi;i++) {
		a = a_chi(chi_ar[i] / real_coverH0) ; // first convert unit of chi from Mpc to c/H0
		z = 1./a - 1.;
		pf = (pf_photoz(z,ni)<0.)? 0:pf_photoz(z,ni); // get rid of unphysical negatives
		f_chi_ar[i] = chi_ar[i] * pf*growfac(a)*g0*gbias.b1_function(z,ni)*hoverh0(a)/real_coverH0;
	}
}

// Integrand for galaxy density RSD
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
		f_chi_RSD_ar[i] = -chi_ar[i] * pf*growfac(a)*g0*f_growth(z)*hoverh0(a)/real_coverH0;
	}
}

// Integrand for lensing magnification of galaxy density
void f_chi_for_Psi_cl_Mag(double* chi_ar, int Nchi, double* f_chi_Mag_ar, int ni){
	double g0 =1./growfac(1.);
	double a, z, fK;
	int i;
	double real_coverH0 = cosmology.coverH0 / cosmology.h0;
	double window_M, wmag;
	for(i=0;i<Nchi;i++) {
		a = a_chi(chi_ar[i] / real_coverH0) ; // first convert unit of chi from Mpc to c/H0
		z = 1./a - 1.;
		fK = f_K(chi_ar[i]/real_coverH0);
		// printf("Here! a, fK, ni: %lg,%lg,%d\n", a, fK, ni);
		wmag = W_mag(a, fK, (double)ni);
		window_M = (-gbias.b_mag[ni]) * wmag/ fK / (real_coverH0*real_coverH0);
		// printf("bmag, wkappa, f_K, real_coverH0, %lg %lg %lg %lg\n", gbias.b_mag[ni], wkappa, fK,real_coverH0);
		// pf = (pf_photoz(z,ni)<0.)? 0:pf_photoz(z,ni); // get rid of unphysical negatives
		// f_chi_Mag_ar[i] = chi_ar[i]/a * window_M*growfac(a)*g0;
		f_chi_Mag_ar[i] = window_M*growfac(a)*g0; // unit [Mpc^-2]
	}
}

// Mixture of non-Limber and Limber of C_cl (galaxy clustering)
void C_cl_mixed(int L, int LMAX, int ni, int nj, double *Cl, double dev, double tolerance) {
	// ni = 4;
	// nj = 4;
	int i,j,i_block;
	long l;
	// run 100 ells at a time, and see if switching to Limber is needed.
	// Save runtime for Limber, and save re-creation time of fftw_plan.
	int Nell_block = 100, Nchi = 1000;
	int ell_ar[Nell_block];
	double **k1_ar, **k2_ar, **Fk1_ar, **Fk2_ar;
	double **Fk1_Mag_ar, **Fk2_Mag_ar;

	k1_ar = malloc(Nell_block * sizeof(double *));
	k2_ar = malloc(Nell_block * sizeof(double *));
	Fk1_ar = malloc(Nell_block * sizeof(double *));
	Fk2_ar = malloc(Nell_block * sizeof(double *));

	Fk1_Mag_ar = malloc(Nell_block * sizeof(double *));
	Fk2_Mag_ar = malloc(Nell_block * sizeof(double *));
	for(i=0;i<Nell_block;i++) {
		k1_ar[i] = malloc(Nchi * sizeof(double));
		k2_ar[i] = malloc(Nchi * sizeof(double));
		Fk1_ar[i] = malloc(Nchi * sizeof(double));
		Fk2_ar[i] = malloc(Nchi * sizeof(double));
		Fk1_Mag_ar[i] = malloc(Nchi * sizeof(double));
		Fk2_Mag_ar[i] = malloc(Nchi * sizeof(double));
		for(j=0;j<Nchi;j++) {
			Fk1_ar[i][j] = 0.;
			Fk2_ar[i][j] = 0.;
			Fk1_Mag_ar[i][j] = 0.;
			Fk2_Mag_ar[i][j] = 0.;
		}
	}

	double chi_ar[Nchi], f1_chi_ar[Nchi], f2_chi_ar[Nchi];
	double f1_chi_RSD_ar[Nchi], f2_chi_RSD_ar[Nchi];
	double f1_chi_Mag_ar[Nchi], f2_chi_Mag_ar[Nchi];

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

	f_chi_for_Psi_cl_Mag(chi_ar, Nchi, f1_chi_Mag_ar, ni);
	if(ni != nj) {f_chi_for_Psi_cl_Mag(chi_ar, Nchi, f2_chi_Mag_ar, nj);}

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
	// char outfilename[] = "c_cl4_rsd.txt";
	// FILE *OUT = fopen(outfilename, "w");


	i_block = 0;
	double cl_temp;

	config my_config, my_config_RSD, my_config_Mag;
	my_config.nu = 1.;
	my_config.c_window_width = 0.25;
	my_config.derivative = 0;
	my_config.N_pad = 200;
	my_config_RSD.nu = 1.01;
	my_config_RSD.c_window_width = 0.25;
	my_config_RSD.derivative = 2;
	my_config_RSD.N_pad = 200;

	my_config_Mag.nu = 1.;
	my_config_Mag.c_window_width = 0.25;
	my_config_Mag.derivative = 0;
	my_config_Mag.N_pad = 200;
	double ell_prefactor;

	double real_coverH0 = cosmology.coverH0 / cosmology.h0;
	double k1_cH0;

	while (fabs(dev) > tolerance){
		//Cl[L] = C_cl_RSD(L,nz,nz);
		for(i=0;i<Nell_block;i++) {ell_ar[i]=i+i_block*Nell_block;}

		cfftlog_ells(chi_ar, f1_chi_ar, Nchi, &my_config, ell_ar, Nell_block, k1_ar, Fk1_ar);
		if(ni != nj) {cfftlog_ells(chi_ar, f2_chi_ar, Nchi, &my_config, ell_ar, Nell_block, k2_ar, Fk2_ar);}

		cfftlog_ells_increment(chi_ar, f1_chi_RSD_ar, Nchi, &my_config_RSD, ell_ar, Nell_block, k1_ar, Fk1_ar);
		if(ni != nj) {cfftlog_ells_increment(chi_ar, f2_chi_RSD_ar, Nchi, &my_config_RSD, ell_ar, Nell_block, k2_ar, Fk2_ar);}

		// Add in lensing magnification contribution
		cfftlog_ells(chi_ar, f1_chi_Mag_ar, Nchi, &my_config_Mag, ell_ar, Nell_block, k1_ar, Fk1_Mag_ar);
		if(ni != nj) {cfftlog_ells(chi_ar, f2_chi_Mag_ar, Nchi, &my_config_Mag, ell_ar, Nell_block, k2_ar, Fk2_Mag_ar);}
		for(i=0;i<Nell_block;i++) {
			ell_prefactor = ell_ar[i]*(ell_ar[i]+1);
			for(j=0;j<Nchi;j++) {
				Fk1_ar[i][j]+= (ell_prefactor / (k1_ar[i][j]*k1_ar[i][j])* Fk1_Mag_ar[i][j]);
				if(ni != nj) {Fk2_ar[i][j]+= (ell_prefactor / (k2_ar[i][j]*k2_ar[i][j])* Fk2_Mag_ar[i][j]);}
			}
		}

		for(i=0;i<Nell_block;i++) {
			cl_temp = 0.;
			for(j=0;j<Nchi;j++) {
				// printf("k,Fk: %d,%d, %lf,%lf\n", i,j, k1_ar[i][j], Fk1_ar[i][j]);
				k1_cH0 = k1_ar[i][j] * real_coverH0;
				if(ni == nj) {
					cl_temp += (Fk1_ar[i][j]) * (Fk1_ar[i][j]) *k1_cH0*k1_cH0*k1_cH0 *p_lin(k1_cH0,1.0)*G_taper(k1_cH0);
					// printf("plin,%lg, %lg\n", k1_ar[i][j],p_lin(k1_cH0,1.0));
				}
				else {
					cl_temp += (Fk1_ar[i][j])*(Fk2_ar[i][j]) *k1_cH0*k1_cH0*k1_cH0 *p_lin(k1_cH0,1.0)*G_taper(k1_cH0);
				}
			}
			Cl[ell_ar[i]] = cl_temp * dlnk * 2./M_PI + C_cl_tomo_nointerp(1.*ell_ar[i],ni,nj) - C_cl_lin_nointerp(1.*ell_ar[i],ni,nj);
			// printf("cl_t/emp: %d, %lg\n", i, cl_temp);
			// fprintf(OUT, "%d %lg %lg %lg", ell_ar[i], Cl[ell_ar[i]], C_cl_tomo_nointerp(1.*ell_ar[i],ni,nj), C_cl_lin_nointerp(1.*ell_ar[i],ni,nj));
			// fprintf(OUT, "\n");
		}

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
	printf("finished bin %d\n", ni);
	free(k1_ar);free(k2_ar);
	free(Fk1_ar);free(Fk2_ar);
	free(Fk1_Mag_ar);free(Fk2_Mag_ar);
	// fclose(OUT);
	// exit(0);
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
			C_cl_mixed(L, LMAX, nz,nz, Cl, dev, tolerance);
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

///////////////////////// 
///////// galaxy-galaxy lensing

double int_for_C_gl_lin(double a, void *params)
{
	double res,ell, fK, k;
	double *ar = (double *) params;
	ell       = ar[2]+0.5;
	fK     = f_K(chi(a));
	k      = ell/fK;
	
	res=W_gal(a,ar[0])*W_kappa(a,fK, ar[1])*dchi_da(a)/fK/fK;
	res= res*p_lin(k,a)*G_taper(k);
	return res;
}

double C_gl_lin_nointerp(double l, int ni, int nj)  //galaxy clustering power spectrum of galaxy bins ni, nj
{
	double array[3] = {1.0*ni,1.0*nj,l};
	return int_gsl_integrate_medium_precision(int_for_C_gl_lin,(void*)array,fmax(amin_lens(ni),amin_lens(nj)),fmin(amax_lens(ni),amax_lens(nj)),NULL,1000);
}


void f_chi_for_Psi_sh(double* chi_ar, int Nchi, double* f_chi_ar, int ns) {
	double g0 =1./growfac(1.);
	double a, z, fK;
	int i;
	double real_coverH0 = cosmology.coverH0 / cosmology.h0;
	double window_L, wkappa;
	for(i=0;i<Nchi;i++) {
		a = a_chi(chi_ar[i] / real_coverH0) ; // first convert unit of chi from Mpc to c/H0
		z = 1./a - 1.;
		fK = f_K(chi_ar[i]/real_coverH0);
		// printf("Here! a, fK, ns: %lg,%lg,%d\n", a, fK, ns);
		wkappa = W_kappa(a, fK, (double)ns);
		window_L = wkappa/ fK / (real_coverH0*real_coverH0);
		// printf("win_L, wkappa, f_K, real_coverH0, %lg %lg %lg %lg\n", window_L, wkappa, fK,real_coverH0);
		f_chi_ar[i] = window_L*growfac(a)*g0; // unit [Mpc^-2]
		// printf("fchi, %lg\n", f_chi_ar[i]);
	}
}


void f_chi_for_Psi_sh_IA(double* chi_ar, int Nchi, double* f_chi_IA_ar, int ns) {
	double g0 =1./growfac(1.);
	double a, z, fK;
	int i;
	double real_coverH0 = cosmology.coverH0 / cosmology.h0;
	double window_L, wkappa;
	for(i=0;i<Nchi;i++) {
		a = a_chi(chi_ar[i] / real_coverH0) ; // first convert unit of chi from Mpc to c/H0
		z = 1./a - 1.;
		fK = f_K(chi_ar[i]/real_coverH0);
		// printf("Here! a, fK, ni: %lg,%lg,%d\n", a, fK, ni);
		wkappa = W_kappa(a, fK, (double)ns);
		window_L = wkappa/ fK / (real_coverH0*real_coverH0);
		// printf("bmag, wkappa, f_K, real_coverH0, %lg %lg %lg %lg\n", gbias.b_mag[ni], wkappa, fK,real_coverH0);
		// pf = (pf_photoz(z,ni)<0.)? 0:pf_photoz(z,ni); // get rid of unphysical negatives
		// f_chi_Mag_ar[i] = chi_ar[i]/a * window_M*growfac(a)*g0;
		f_chi_IA_ar[i] = window_L*growfac(a)*g0; // unit [Mpc^-2]
	}
}
// Mixture of non-Limber and Limber of C_cl (G-G lensing)
void C_gl_mixed(int L, int LMAX, int nl, int ns, double *Cl, double dev, double tolerance) {
	// ni = 4;
	// nj = 4;
	int i,j,i_block;
	long l;
	// run 100 ells at a time, and see if switching to Limber is needed.
	// Save runtime for Limber, and save re-creation time of fftw_plan.
	int Nell_block = 100, Nchi = 1000;
	int ell_ar[Nell_block];
	double **k1_ar, **k2_ar, **Fk1_ar, **Fk2_ar;
	double **Fk1_Mag_ar;

	k1_ar = malloc(Nell_block * sizeof(double *));
	k2_ar = malloc(Nell_block * sizeof(double *));
	Fk1_ar = malloc(Nell_block * sizeof(double *));
	Fk2_ar = malloc(Nell_block * sizeof(double *));

	Fk1_Mag_ar = malloc(Nell_block * sizeof(double *));
	for(i=0;i<Nell_block;i++) {
		k1_ar[i] = malloc(Nchi * sizeof(double));
		k2_ar[i] = malloc(Nchi * sizeof(double));
		Fk1_ar[i] = malloc(Nchi * sizeof(double));
		Fk2_ar[i] = malloc(Nchi * sizeof(double));
		Fk1_Mag_ar[i] = malloc(Nchi * sizeof(double));
		for(j=0;j<Nchi;j++) {
			Fk1_ar[i][j] = 0.;
			Fk2_ar[i][j] = 0.;
			Fk1_Mag_ar[i][j] = 0.;
		}
	}

	double chi_ar[Nchi];
	double f1_chi_ar[Nchi], f1_chi_RSD_ar[Nchi], f1_chi_Mag_ar[Nchi];
	double f2_chi_ar[Nchi], f2_chi_IA_ar[Nchi];

	double chi_min = 60., chi_max = 6000.;
	double dlnchi = log(chi_max/chi_min) / (Nchi - 1.);
	double dlnk = dlnchi;

	for(i=0; i<Nchi; i++) {
		chi_ar[i] = chi_min * exp(dlnchi*i);
	}

	f_chi_for_Psi_cl(chi_ar, Nchi, f1_chi_ar, nl);
	f_chi_for_Psi_cl_RSD(chi_ar, Nchi, f1_chi_RSD_ar, nl);
	f_chi_for_Psi_cl_Mag(chi_ar, Nchi, f1_chi_Mag_ar, nl);

	f_chi_for_Psi_sh(chi_ar, Nchi, f2_chi_ar, ns);
	// f_chi_for_Psi_sh_IA(chi_ar, Nchi, f2_chi_ar_IA, nj);

	// char outfilename[] = "f_chi_gl.txt";
	// FILE *OUT = fopen(outfilename, "w");
	// for(i=0; i<Nchi; i++) {
	// 	// fprintf(OUT, "%lg %lg", chi_ar[i], f1_chi_ar[i]);
	// 	fprintf(OUT, "%lg %lg", chi_ar[i], f2_chi_ar[i]);
	// 	fprintf(OUT, "\n");
	// }
	// for(i=0; i<Nchi; i++) {
	// 	printf("f_chi_ar: %d, %lg\n", i, f2_chi_ar[i]);
	// }
	// exit(0);


	i_block = 0;
	double cl_temp;

	config my_config, my_config_RSD, my_config_Mag;
	my_config.nu = 1.;
	my_config.c_window_width = 0.25;
	my_config.derivative = 0;
	my_config.N_pad = 200;
	my_config_RSD.nu = 1.01;
	my_config_RSD.c_window_width = 0.25;
	my_config_RSD.derivative = 2;
	my_config_RSD.N_pad = 200;

	my_config_Mag.nu = 1.;
	my_config_Mag.c_window_width = 0.25;
	my_config_Mag.derivative = 0;
	my_config_Mag.N_pad = 200;
	double ell_prefactor, ell_prefactor2;

	double real_coverH0 = cosmology.coverH0 / cosmology.h0;
	double k1_cH0;

	while (fabs(dev) > tolerance){
		//Cl[L] = C_cl_RSD(L,nz,nz);
		for(i=0;i<Nell_block;i++) {ell_ar[i]=i+i_block*Nell_block;}

		// galaxy density part
		cfftlog_ells(chi_ar, f1_chi_ar, Nchi, &my_config, ell_ar, Nell_block, k1_ar, Fk1_ar);
		cfftlog_ells_increment(chi_ar, f1_chi_RSD_ar, Nchi, &my_config_RSD, ell_ar, Nell_block, k1_ar, Fk1_ar);
		// Add in lensing magnification contribution
		cfftlog_ells(chi_ar, f1_chi_Mag_ar, Nchi, &my_config_Mag, ell_ar, Nell_block, k1_ar, Fk1_Mag_ar);
		for(i=0;i<Nell_block;i++) {
			ell_prefactor = ell_ar[i]*(ell_ar[i]+1);
			for(j=0;j<Nchi;j++) {
				Fk1_ar[i][j]+= (ell_prefactor / (k1_ar[i][j]*k1_ar[i][j])* Fk1_Mag_ar[i][j]);
				// printf("Fk1: %d,%d, %lg\n", i,j, Fk1_ar[i][j]);
			}
		}

		// shear part
		cfftlog_ells(chi_ar, f2_chi_ar, Nchi, &my_config, ell_ar, Nell_block, k2_ar, Fk2_ar);
		for(i=0;i<Nell_block;i++) {
			ell_prefactor2 =(ell_ar[i]-1)*ell_ar[i]*(ell_ar[i]+1)*(ell_ar[i]+2);
			if(ell_prefactor2<0.) {ell_prefactor2=0.;}
			else {ell_prefactor2 = sqrt(ell_prefactor2);}

			for(j=0;j<Nchi;j++) {
				Fk2_ar[i][j]*= (ell_prefactor2 / (k2_ar[i][j]*k2_ar[i][j]));
				// printf("Fk2: %d,%d, %lg\n", i,j, Fk2_ar[i][j]);
			}
		}
		// exit(0);
		for(i=0;i<Nell_block;i++) {
			cl_temp = 0.;
			for(j=0;j<Nchi;j++) {
				// printf("k,Fk: %d,%d, %lf,%lf\n", i,j, k1_ar[i][j], Fk1_ar[i][j]);
				k1_cH0 = k1_ar[i][j] * real_coverH0;
				cl_temp += (Fk1_ar[i][j])*(Fk2_ar[i][j]) *k1_cH0*k1_cH0*k1_cH0 *p_lin(k1_cH0,1.0)*G_taper(k1_cH0);
			}
			Cl[ell_ar[i]] = cl_temp * dlnk * 2./M_PI + C_gl_tomo_nointerp(1.*ell_ar[i],nl,ns) - C_gl_lin_nointerp(1.*ell_ar[i],nl,ns);
			// printf("cl_t/emp: %d, %lg\n", i, cl_temp);
			// fprintf(OUT, "%d %lg %lg %lg", ell_ar[i], Cl[ell_ar[i]], C_cl_tomo_nointerp(1.*ell_ar[i],ni,nj), C_cl_lin_nointerp(1.*ell_ar[i],ni,nj));
			// fprintf(OUT, "\n");
		}

		i_block++;
		L = i_block*Nell_block -1 ;
		dev = Cl[L]/C_gl_tomo_nointerp(1.0*L,nl,ns)-1.;
	   // printf("ni,L,Cl[L],dev=%d %d %e %e\n",ni,L,Cl[L],dev);
		// printf("i_block: %d\n", i_block);
	}
	L++;
	printf("switching to Limber calculation at l = %d\n",L);
	for (l = 1; l < LMAX; l++){
		Cl[l]=C_gl_tomo((double)l,nl,ns);
	}
	// for (l = L; l < LMAX; l++){
	// 	Cl[l]=C_gl_tomo((double)l,nl,ns);
	// }
	printf("finished bin %d %d\n", nl,ns);
	free(k1_ar);free(k2_ar);
	free(Fk1_ar);free(Fk2_ar);
	free(Fk1_Mag_ar);
	// fclose(OUT);
	// exit(0);
}

double w_gamma_t_nonLimber(int nt, int ni, int nj){
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
		printf("cosmo2D_fullsky.c:w_gamma_t_tomo: like.Ntheta not initialized\nEXIT\n"); exit(1);
	}
	if (Pl ==0){
		Pl =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
		Cl = create_double_vector(0,LMAX-1);
		w_vec = create_double_vector(0,tomo.ggl_Npowerspectra*like.Ntheta-1);
		NTHETA = like.Ntheta;
		double *xmin, *xmax, *Pmin, *Pmax, *dP;
		xmin= create_double_vector(0, like.Ntheta-1);
		xmax= create_double_vector(0, like.Ntheta-1);
		double logdt=(log(like.vtmax)-log(like.vtmin))/like.Ntheta;
		for(i=0; i<like.Ntheta ; i++){
			xmin[i]=cos(exp(log(like.vtmin)+(i+0.0)*logdt));
			xmax[i]=cos(exp(log(like.vtmin)+(i+1.0)*logdt));
		}
		Pmin= create_double_vector(0, LMAX+1);
		Pmax= create_double_vector(0, LMAX+1);

		for (i = 0; i<NTHETA; i ++){
			printf("Tabulating Legendre coefficients %d/%d\n",i+1, NTHETA);
			gsl_sf_legendre_Pl_array(LMAX, xmin[i],Pmin);
			gsl_sf_legendre_Pl_array(LMAX, xmax[i],Pmax);
			for (int l = 2; l < LMAX; l ++){
				//Pl[i][l] = (2.*l+1)/(4.*M_PI*l*(l+1))*gsl_sf_legendre_Plm(l,2,cos(like.theta[i]));	
				Pl[i][l] = (2.*l+1)/(4.*M_PI*l*(l+1)*(xmin[i]-xmax[i]))
				*((l+2./(2*l+1.))*(Pmin[l-1]-Pmax[l-1])
				+(2-l)*(xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
				-2./(2*l+1.)*(Pmin[l+1]-Pmax[l]));
			}
		}
		free_double_vector(xmin,0,like.Ntheta-1);
		free_double_vector(xmax,0,like.Ntheta-1);
		free_double_vector(Pmin,0,LMAX+1);
		free_double_vector(Pmax,0,LMAX+1);
	}
	if (recompute_ggl(C,G,N,ni)){
		//required fractional accuracy in C(l)
		double tolerance= 0.01;
		//dev will be the actual difference between exact and Limber calcuation
		double dev;

		// if (like.IA != 0 && like.IA != 3 && like.IA != 4){printf("cosmo2D_real.c: w_gamma_t_tomo does not support like.IA = %d yet\nEXIT!\n",like.IA); exit(1);}
		// C_tomo_pointer C_gl_pointer = &C_gl_tomo;
		// if (like.IA ==3 || like.IA ==4) C_gl_pointer = &C_ggl_IA_tab;

		for (nz = 0; nz <tomo.ggl_Npowerspectra; nz ++){
			int L = 1;
			// initialize to large value in order to start while loop
			dev=10.*tolerance;
			C_gl_mixed(L, LMAX, ZL(nz),ZS(nz), Cl, dev, tolerance);
			for (i = 0; i < NTHETA; i++){
				w_vec[nz*like.Ntheta+i] =0;
				for (l = 1; l < LMAX; l++){
					w_vec[nz*like.Ntheta+i]+=Pl[i][l]*Cl[l];
				}
			}
		}
		// for (nz = 0; nz <tomo.ggl_Npowerspectra; nz ++){
		// 	for (l = 1; l < LMAX; l++){
		// 		Cl[l]=C_ggl_IA_tab(1.0*l,ZL(nz),ZS(nz));
		// 	}
		// 	for (i = 0; i < NTHETA; i++){
		// 		w_vec[nz*like.Ntheta+i] =0;
		// 		for (l = 1; l < LMAX; l++){
		// 			w_vec[nz*like.Ntheta+i]+=Pl[i][l]*Cl[l];
		// 		}
		// 	}
		// }

		update_cosmopara(&C);
		update_galpara(&G);
		update_nuisance(&N);
	}
	return w_vec[N_ggl(ni,nj)*like.Ntheta+nt];  
}