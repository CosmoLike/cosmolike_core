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
		if( (z<tomo.clustering_zmin[ni]) || (z>tomo.clustering_zmax[ni]) )
		{
			f_chi_ar[i] = 0.;
		}
		else
		{
			pf = (pf_photoz(z,ni)<0.)? 0:pf_photoz(z,ni); // get rid of unphysical negatives
			f_chi_ar[i] = chi_ar[i] * pf*growfac(a)*g0*gbias.b1_function(z,ni)*hoverh0(a)/real_coverH0;
		}
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
		if( (z<tomo.clustering_zmin[ni]) || (z>tomo.clustering_zmax[ni]) )
		{
			f_chi_RSD_ar[i] = 0.;
		}
		else
		{
			pf = (pf_photoz(z,ni)<0.)? 0:pf_photoz(z,ni); // get rid of unphysical negatives
			f_chi_RSD_ar[i] = -chi_ar[i] * pf*growfac(a)*g0*f_growth(z)*hoverh0(a)/real_coverH0;
		}
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
		if( (z<tomo.clustering_zmin[ni]) || (z>tomo.clustering_zmax[ni]) )
		{
			f_chi_Mag_ar[i] = 0.;
		}
		else
		{
			// printf("Here! a, fK, ni: %lg,%lg,%d\n", a, fK, ni);
			wmag = W_mag(a, fK, (double)ni);
			window_M = wmag/ fK / (real_coverH0*real_coverH0);
			// printf("bmag, wkappa, f_K, real_coverH0, %lg %lg %lg %lg\n", gbias.b_mag[ni], wkappa, fK,real_coverH0);
			// pf = (pf_photoz(z,ni)<0.)? 0:pf_photoz(z,ni); // get rid of unphysical negatives
			// f_chi_Mag_ar[i] = chi_ar[i]/a * window_M*growfac(a)*g0;
			f_chi_Mag_ar[i] = window_M*growfac(a)*g0; // unit [Mpc^-2]
		}
	}
}

// Mixture of non-Limber and Limber of C_cl (galaxy clustering)
void C_cl_mixed(int L, int LMAX, int ni, int nj, double *Cl, double dev, double tolerance) {
	// ni = 4;
	// nj = 4;

	extern double S_integrands_cl[5][200][1000];
	extern int S_integrands_cl_flag[5];

	int i,j,i_block;
	long l;
	// run 100 ells at a time, and see if switching to Limber is needed.
	// Save runtime for Limber, and save re-creation time of fftw_plan.
	int Nell_block = 200, Nchi = 1000;
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

	double ell_prefactor, ell_prefactor2;

	double real_coverH0 = cosmology.coverH0 / cosmology.h0;
	double k1_cH0;

	double cl_temp;

	i_block = 0;
	for(i=0;i<Nell_block;i++) {ell_ar[i]=i+i_block*Nell_block;}


	if(S_integrands_cl_flag[ni]) {
		for(i=0;i<Nell_block;i++) {
			for(j=0;j<Nchi;j++) {
				k1_ar[i][j] = (i+1.)/chi_ar[Nchi-1-j];
				k2_ar[i][j] = k1_ar[i][j];
				Fk1_ar[i][j] = S_integrands_cl[ni][i][j];
				Fk2_ar[i][j] = Fk1_ar[i][j];
			}
		}
	}
	else
	{
		f_chi_for_Psi_cl(chi_ar, Nchi, f1_chi_ar, ni);
		if(ni != nj) {f_chi_for_Psi_cl(chi_ar, Nchi, f2_chi_ar, nj);}

		f_chi_for_Psi_cl_RSD(chi_ar, Nchi, f1_chi_RSD_ar, ni);
		if(ni != nj) {f_chi_for_Psi_cl_RSD(chi_ar, Nchi, f2_chi_RSD_ar, nj);}

		f_chi_for_Psi_cl_Mag(chi_ar, Nchi, f1_chi_Mag_ar, ni);
		if(ni != nj) {f_chi_for_Psi_cl_Mag(chi_ar, Nchi, f2_chi_Mag_ar, nj);}

		// char outfilename[] = "f1_chi4.txt";
		// char outfilename[] = "f1_chi4_rsd.txt";
		// char *outfilename = (char*)malloc(40 * sizeof(char));;
		// sprintf(outfilename, "f1_chi_Mag_%d.txt", ni);
		// FILE *OUT = fopen(outfilename, "w");
		
		// for(i=0; i<Nchi; i++) {
		// 	// fprintf(OUT, "%lg %lg", chi_ar[i], f1_chi_ar[i]);
		// 	// fprintf(OUT, "%lg %lg", chi_ar[i], f1_chi_RSD_ar[i]);
		// 	fprintf(OUT, "%lg %lg", chi_ar[i], f1_chi_Mag_ar[i]);
		// 	fprintf(OUT, "\n");
		// }
		// for(i=0; i<Nchi; i++) {
		// 	printf("f_chi_ar: %d, %lg, %lg, %lg, %lg\n", i,chi_ar[i], f1_chi_ar[i],f1_chi_RSD_ar[i], f1_chi_Mag_ar[i]);
		// }
		// exit(0);
		// char *outfilename = (char*)malloc(40 * sizeof(char));;
		// sprintf(outfilename, "c_cl_%d_%d_mag_fft.txt", ni,nj);
		// FILE *OUT = fopen(outfilename, "w");


		config my_config, my_config_RSD, my_config_Mag;
		my_config.nu = 1.;
		my_config.c_window_width = 0.25;
		my_config.derivative = 0;
		my_config.N_pad = 200;
		my_config.N_extrap_low = 0;
		my_config.N_extrap_high = 0;

		my_config_RSD.nu = 1.01;
		my_config_RSD.c_window_width = 0.25;
		my_config_RSD.derivative = 2;
		my_config_RSD.N_pad = 500;
		my_config_RSD.N_extrap_low = 0;
		my_config_RSD.N_extrap_high = 0;

		my_config_Mag.nu = 1.;
		my_config_Mag.c_window_width = 0.25;
		my_config_Mag.derivative = 0;
		my_config_Mag.N_pad = 500;
		my_config_Mag.N_extrap_low = 0;
		my_config_Mag.N_extrap_high = 0;

		cfftlog_ells(chi_ar, f1_chi_ar, Nchi, &my_config, ell_ar, Nell_block, k1_ar, Fk1_ar);
		if(ni != nj) {cfftlog_ells(chi_ar, f2_chi_ar, Nchi, &my_config, ell_ar, Nell_block, k2_ar, Fk2_ar);}

		cfftlog_ells_increment(chi_ar, f1_chi_RSD_ar, Nchi, &my_config_RSD, ell_ar, Nell_block, k1_ar, Fk1_ar);
		if(ni != nj) {cfftlog_ells_increment(chi_ar, f2_chi_RSD_ar, Nchi, &my_config_RSD, ell_ar, Nell_block, k2_ar, Fk2_ar);}

		// Add in lensing magnification contribution
		cfftlog_ells(chi_ar, f1_chi_Mag_ar, Nchi, &my_config_Mag, ell_ar, Nell_block, k1_ar, Fk1_Mag_ar);
		if(ni != nj) {cfftlog_ells(chi_ar, f2_chi_Mag_ar, Nchi, &my_config_Mag, ell_ar, Nell_block, k2_ar, Fk2_Mag_ar);}
		for(i=0;i<Nell_block;i++) {
			ell_prefactor = ell_ar[i]*(ell_ar[i]+1.);
			for(j=0;j<Nchi;j++) {
				Fk1_ar[i][j]+= (ell_prefactor / (k1_ar[i][j]*k1_ar[i][j]) * (gbias.b_mag[ni]) *  Fk1_Mag_ar[i][j]);
				if(ni != nj) {Fk2_ar[i][j]+= (ell_prefactor / (k2_ar[i][j]*k2_ar[i][j])* (gbias.b_mag[nj]) *  Fk2_Mag_ar[i][j]);}

				S_integrands_cl[ni][i][j] = Fk1_ar[i][j];
			}
		}
		S_integrands_cl_flag[ni] = 1;
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
		// Cl[ell_ar[i]] = cl_temp * dlnk * 2./M_PI;
		// printf("cl_t/emp: %d, %lg\n", i, cl_temp);
		// fprintf(OUT, "%d %lg %lg %lg\n", ell_ar[i], Cl[ell_ar[i]], C_cl_tomo_nointerp(1.*ell_ar[i],ni,nj), C_cl_lin_nointerp(1.*ell_ar[i],ni,nj));
	}

	i_block++;
	L = i_block*Nell_block -1 ;
	dev = Cl[L]/C_cl_tomo_nointerp((double)L,ni,nj)-1.;
   // printf("ni,L,Cl[L],dev=%d %d %e %e\n",ni,L,Cl[L],dev);
	// printf("i_block: %d\n", i_block);



	L++;
	// printf("switching to Limber calculation at l = %d\n",L);
	for (l = L; l < LMAX; l++){
		Cl[l]=C_cl_tomo((double)l,ni,nj);
	}
	// printf("finished bin %d\n", ni);
	for(i=0;i<Nell_block;i++) {
		free(k1_ar[i]);free(k2_ar[i]);
		free(Fk1_ar[i]);free(Fk2_ar[i]);
		free(Fk1_Mag_ar[i]);free(Fk2_Mag_ar[i]);
	}
	free(k1_ar);free(k2_ar);
	free(Fk1_ar);free(Fk2_ar);
	free(Fk1_Mag_ar);free(Fk2_Mag_ar);
	// fclose(OUT);
	// exit(0);
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

	double ell_prefactor1 = (ar[2])*(ar[2]+1.);
	double ell_prefactor2 = (ar[2]-1.)*ell_prefactor1*(ar[2]+2.);
	if(ell_prefactor2<=0.) 
		ell_prefactor2=0.;
	else
		ell_prefactor2=sqrt(ell_prefactor2);

	double chi_0,chi_1,a_0,a_1;
	chi_0 = f_K(ell/k);
	chi_1 = f_K((ell+1.)/k);
	if (chi_1 > chi(limits.a_min)){
	return 0;}
	a_0 = a_chi(chi_0);
	a_1 = a_chi(chi_1);

	double wgal = W_gal(a,ar[0]);
	wgal += W_mag(a,fK,ar[0])*(ell_prefactor1/ell/ell -1.) ;
  	res=(wgal + W_RSD(ell, a_0, a_1, ar[0]))*W_kappa(a,fK, ar[1])*dchi_da(a)/fK/fK * ell_prefactor2/ell/ell;
	res= res*p_lin(k,a)*G_taper(k);
	return res;
}

double C_gl_lin_nointerp(double l, int ni, int nj)  //galaxy clustering power spectrum of galaxy bins ni, nj
{
	double array[3] = {1.0*ni,1.0*nj,l};
	return int_gsl_integrate_medium_precision(int_for_C_gl_lin,(void*)array,amin_lens(ni),amax_lens(ni),NULL,1000);
}


double int_for_C_gl_lin_IA(double a, void *params)
{
	double res,ell, fK, k;
	double *ar = (double *) params;
	ell       = ar[2]+0.5;
	fK     = f_K(chi(a));
	k      = ell/fK;
	
	double norm = cosmology.Omega_m*nuisance.c1rhocrit_ia*growfac(0.9999)/growfac(a)*nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia),nuisance.eta_ia);
	res=W_gal(a,ar[0])*(W_kappa(a,fK,ar[1])-W_source(a,ar[1])*norm)*dchi_da(a)/fK/fK;
	res= res*p_lin(k,a)*G_taper(k);
	return res;
}

double C_gl_lin_nointerp_IA(double l, int ni, int nj)  //galaxy clustering power spectrum of galaxy bins ni, nj
{
	double array[3] = {1.0*ni,1.0*nj,l};
	return int_gsl_integrate_medium_precision(int_for_C_gl_lin_IA,(void*)array,amin_lens(ni),amax_lens(ni),NULL,1000);
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
		if( (a<amin_source(ns)) || (a>amax_source(ns)) )
		{
			f_chi_ar[i] = 0.;
		}
		else
		{
			// printf("Here! a, fK, ns: %lg,%lg,%d\n", a, fK, ns);
			wkappa = W_kappa(a, fK, (double)ns);
			window_L = wkappa/ fK / (real_coverH0*real_coverH0);
			// printf("win_L, wkappa, f_K, real_coverH0, %lg %lg %lg %lg\n", window_L, wkappa, fK,real_coverH0);
			f_chi_ar[i] = window_L*growfac(a)*g0; // unit [Mpc^-2]
			// printf("fchi, %lg\n", f_chi_ar[i]);
			// f_chi_ar[i] = f_chi_ar[i];
		}
	}
}

void f_chi_for_Psi_sh_IA(double* chi_ar, int Nchi, double* f_chi_IA_ar, int ns) {
	double g0 =1./growfac(1.);
	double a, z, fK;
	int i;
	double real_coverH0 = cosmology.coverH0 / cosmology.h0;
	double window_L, wsource;
	double norm = cosmology.Omega_m*nuisance.c1rhocrit_ia*growfac(0.9999)/growfac(a)*nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia),nuisance.eta_ia);

	for(i=0;i<Nchi;i++) {
		a = a_chi(chi_ar[i] / real_coverH0) ; // first convert unit of chi from Mpc to c/H0
		z = 1./a - 1.;
		fK = f_K(chi_ar[i]/real_coverH0);
		if( (a<amin_source(ns)) || (a>amax_source(ns)) )
		{
			f_chi_IA_ar[i] = 0.;
		}
		else
		{
			// printf("Here! a, fK, ni: %lg,%lg,%d\n", a, fK, ni);
			wsource = W_source(a, (double)ns);
			window_L = -wsource * norm / fK / (real_coverH0*real_coverH0);
			// printf("bmag, wkappa, f_K, real_coverH0, %lg %lg %lg %lg\n", gbias.b_mag[ni], wkappa, fK,real_coverH0);
			f_chi_IA_ar[i] = window_L*growfac(a)*g0; // unit [Mpc^-2]
		}
	}
}

// Mixture of non-Limber and Limber of C_cl (G-G lensing)
void C_gl_mixed(int L, int LMAX, int nl, int ns, double *Cl, double dev, double tolerance) {
	// nl = 4;
	// ns = 0;

	extern double S_integrands_cl[5][200][1000];
	extern double S_integrands_sh[4][200][1000];
	extern int S_integrands_cl_flag[5];
	extern int S_integrands_sh_flag[4];


	int i,j,i_block;
	long l;
	// run 100 ells at a time, and see if switching to Limber is needed.
	// Save runtime for Limber, and save re-creation time of fftw_plan.
	int Nell_block = 200, Nchi = 1000;
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

	// double f2_chi_temp[Nchi];

	double chi_min = 60., chi_max = 6000.;
	// double chi_min = 6., chi_max = 6000.;

	double dlnchi = log(chi_max/chi_min) / (Nchi - 1.);
	double dlnk = dlnchi;

	for(i=0; i<Nchi; i++) {
		chi_ar[i] = chi_min * exp(dlnchi*i);
	}

	double ell_prefactor, ell_prefactor2;

	double real_coverH0 = cosmology.coverH0 / cosmology.h0;
	double k1_cH0;

	double cl_temp;

	i_block = 0;
	for(i=0;i<Nell_block;i++) {ell_ar[i]=i+i_block*Nell_block;}

	if(S_integrands_cl_flag[nl]) {
		for(i=0;i<Nell_block;i++) {
			for(j=0;j<Nchi;j++) {
				k1_ar[i][j] = (i+1.)/chi_ar[Nchi-1-j];
				Fk1_ar[i][j] = S_integrands_cl[nl][i][j];
			}
		}
	}
	else
	{
		f_chi_for_Psi_cl(chi_ar, Nchi, f1_chi_ar, nl);
		f_chi_for_Psi_cl_RSD(chi_ar, Nchi, f1_chi_RSD_ar, nl);
		f_chi_for_Psi_cl_Mag(chi_ar, Nchi, f1_chi_Mag_ar, nl);
		// for(j=0;j<Nchi;j++) {
		// 	f1_chi_ar[j] += f1_chi_Mag_ar[j];
		// }


		// char outfilename[] = "f_chi_gl1.txt";
		// char outfilename[] = "f1_chi_gl1.txt";
		// char *outfilename = (char*)malloc(40 * sizeof(char));;
		// sprintf(outfilename, "f2_chi_L_%d.txt", ns);
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
		// char *outfilename = (char*)malloc(40 * sizeof(char));;
		// sprintf(outfilename, "c_gl_%d_%d_mag_noIA.txt", nl,ns);
		// FILE *OUT = fopen(outfilename, "w");


		// char *outfilename = (char*)malloc(40 * sizeof(char));;
		// sprintf(outfilename, "pk.txt");
		// FILE *OUT = fopen(outfilename, "w");


		config my_config, my_config_RSD, my_config_Mag;
		my_config.nu = 1.;
		my_config.c_window_width = 0.25;
		my_config.derivative = 0;
		my_config.N_pad = 200;
		my_config.N_extrap_low = 0;
		my_config.N_extrap_high = 0;

		my_config_RSD.nu = 1.01;
		my_config_RSD.c_window_width = 0.25;
		my_config_RSD.derivative = 2;
		my_config_RSD.N_pad = 200;
		my_config_RSD.N_extrap_low = 0;
		my_config_RSD.N_extrap_high = 0;

		my_config_Mag.nu = 1.;
		my_config_Mag.c_window_width = 0.25;
		my_config_Mag.derivative = 0;
		my_config_Mag.N_pad = 500;
		my_config_Mag.N_extrap_low = 0;
		my_config_Mag.N_extrap_high = 0;

		//Cl[L] = C_cl_RSD(L,nz,nz);
		// galaxy density part
		cfftlog_ells(chi_ar, f1_chi_ar, Nchi, &my_config, ell_ar, Nell_block, k1_ar, Fk1_ar);
		cfftlog_ells_increment(chi_ar, f1_chi_RSD_ar, Nchi, &my_config_RSD, ell_ar, Nell_block, k1_ar, Fk1_ar);

		// Add in lensing magnification contribution
		cfftlog_ells(chi_ar, f1_chi_Mag_ar, Nchi, &my_config_Mag, ell_ar, Nell_block, k1_ar, Fk1_Mag_ar);
		for(i=0;i<Nell_block;i++) {
			ell_prefactor = ell_ar[i]*(ell_ar[i]+1);
			for(j=0;j<Nchi;j++) {
				Fk1_ar[i][j]+= (ell_prefactor / (k1_ar[i][j]*k1_ar[i][j])* (gbias.b_mag[nl]) * Fk1_Mag_ar[i][j]);
				// printf("Fk1: %d,%d, %lg\n", i,j, Fk1_ar[i][j]);
				S_integrands_cl[nl][i][j] = Fk1_ar[i][j];
			}
		}
		S_integrands_cl_flag[nl] = 1;
	}

	if(S_integrands_sh_flag[ns]) {
		for(i=0;i<Nell_block;i++) {
			for(j=0;j<Nchi;j++) {
				k2_ar[i][j] = (i+1.)/chi_ar[Nchi-1-j];
				Fk2_ar[i][j] = S_integrands_sh[ns][i][j];
			}
		}
	}
	else
	{
		config my_config_L;
		my_config_L.nu = 1.;
		my_config_L.c_window_width = 0.25;
		my_config_L.derivative = 0;
		my_config_L.N_pad = 500.;
		my_config_L.N_extrap_low = 0;
		my_config_L.N_extrap_high = 0;

		f_chi_for_Psi_sh(chi_ar, Nchi, f2_chi_ar, ns);

		// turn on IA
		if (like.IA ==3 || like.IA ==4) {
			f_chi_for_Psi_sh_IA(chi_ar, Nchi, f2_chi_IA_ar, ns);
			for(i=0;i<Nchi;i++) {
				f2_chi_ar[i] += f2_chi_IA_ar[i];
			}
		}

		cfftlog_ells(chi_ar, f2_chi_ar, Nchi, &my_config_L, ell_ar, Nell_block, k2_ar, Fk2_ar);	

		for(i=0;i<Nell_block;i++) {
			ell_prefactor2 =(ell_ar[i]-1.)*ell_ar[i]*(ell_ar[i]+1.)*(ell_ar[i]+2.);
			if(ell_prefactor2<=0.) {ell_prefactor2=0.;}
			else {ell_prefactor2 = sqrt(ell_prefactor2);}

			for(j=0;j<Nchi;j++) {
				// printf("preFk2: %d,%d, %lg,%d, %lg,%lg\n", i,j, Fk2_ar[i][j], ell_ar[i],ell_prefactor2,k2_ar[i][j]);
				// Fk2_ar[i][j] = -Fk2_ar[i][j]+f2_chi_ar[0]*sqrt(M_PI)/4.* exp(lngamma_lanczos_real(ell_ar[i]/2.)-lngamma_lanczos_real((ell_ar[i]+3.)/2.));
				Fk2_ar[i][j]*= (ell_prefactor2 / (k1_ar[i][j]*k1_ar[i][j]));
				// printf("Fk2: %d,%d, %lg\n", i,j, Fk2_ar[i][j]);
				S_integrands_sh[ns][i][j] = Fk2_ar[i][j];
			}
		}
		S_integrands_sh_flag[ns] = 1;
	}
	// exit(0);
	double (*C_gl_tomo_nointerp_pointer)(double,int,int);
	double (*C_gl_lin_nointerp_pointer)(double,int,int);
	double (*C_gl_tomo_pointer)(double,int,int);
	if (like.IA ==3 || like.IA ==4) {
		C_gl_tomo_nointerp_pointer = &C_ggl_IA_tab;
		C_gl_lin_nointerp_pointer = &C_gl_lin_nointerp_IA;
		C_gl_tomo_pointer = &C_ggl_IA_tab;
	}
	else {
		C_gl_tomo_nointerp_pointer = &C_gl_tomo_nointerp;
		C_gl_lin_nointerp_pointer = &C_gl_lin_nointerp;
		C_gl_tomo_pointer = &C_gl_tomo;
	}

	for(i=0;i<Nell_block;i++) {
		cl_temp = 0.;
		for(j=0;j<Nchi;j++) {
			// printf("k,Fk: %d,%d, %lf,%lf\n", i,j, k1_ar[i][j], Fk1_ar[i][j]);
			k1_cH0 = k1_ar[i][j] * real_coverH0;
			cl_temp += (Fk1_ar[i][j])*(Fk2_ar[i][j]) *k1_cH0*k1_cH0*k1_cH0 *p_lin(k1_cH0,1.0)*G_taper(k1_cH0);
		}
		Cl[ell_ar[i]] = cl_temp * dlnk * 2./M_PI + C_gl_tomo_nointerp_pointer(1.*ell_ar[i],nl,ns) - C_gl_lin_nointerp_pointer(1.*ell_ar[i],nl,ns);
		// Cl[ell_ar[i]] = cl_temp * dlnk * 2./M_PI + C_ggl_IA_tab(1.*ell_ar[i],nl,ns) - C_gl_lin_nointerp_IA(1.*ell_ar[i],nl,ns);
		// printf("cl_temp: %d, %lg\n", i, cl_temp);
		// fprintf(OUT, "%d %lg %lg %lg\n", ell_ar[i], Cl[ell_ar[i]], C_gl_tomo_nointerp(1.*ell_ar[i],nl,ns), C_gl_lin_nointerp(1.*ell_ar[i],nl,ns));
		dev = Cl[ell_ar[i]]/C_gl_tomo_nointerp_pointer(1.0*ell_ar[i],nl,ns)-1.;
		// dev = Cl[L]/C_ggl_IA_tab(1.0*L,nl,ns)-1.;

	   // printf("nl,ns,L,Cl[L],dev=%d %d %d %e %e\n",nl,ns,ell_ar[i],Cl[ell_ar[i]],dev);
	}

	i_block++;
	L = i_block*Nell_block -1 ;
	dev = Cl[L]/C_gl_tomo_nointerp(1.0*L,nl,ns)-1.;
	// dev = Cl[L]/C_ggl_IA_tab(1.0*L,nl,ns)-1.;

 //   printf("ni,L,Cl[L],dev=%d %d %e %e\n",ni,L,Cl[L],dev);
	// printf("i_block: %d\n", i_block);
	// exit(0);
	L++;
	// printf("switching to Limber calculation at l = %d\n",L);
	// for (l = 1; l < 20; l++){
	// 	Cl[l]=C_gl_tomo_nointerp((double)l,nl,ns);
	// }
	// for (l = 20; l < LMAX; l++){
	// 	Cl[l]=C_gl_tomo((double)l,nl,ns);
	// }

	for (l = L; l < LMAX; l++){
		Cl[l]=C_gl_tomo_pointer((double)l,nl,ns);
		// Cl[l]=C_ggl_IA_tab((double)l,nl,ns);
	}
	// printf("finished bin %d %d\n", nl,ns);
	for(i=0;i<Nell_block;i++) {
		free(k1_ar[i]);free(k2_ar[i]);
		free(Fk1_ar[i]);free(Fk2_ar[i]);
		free(Fk1_Mag_ar[i]);
	}
	free(k1_ar);free(k2_ar);
	free(Fk1_ar);free(Fk2_ar);
	free(Fk1_Mag_ar);
	// fclose(OUT);
	// exit(0);
}

/////// Real Space NonLimber Data Vectors

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
			// printf("Tabulating Legendre coefficients %d/%d\n",i+1, NTHETA);
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
			// printf("Tabulating Legendre coefficients %d/%d\n",i+1, NTHETA);
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
