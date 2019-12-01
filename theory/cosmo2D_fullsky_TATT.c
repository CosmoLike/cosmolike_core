double C_EE_tab(double l, int ni, int nj);
double C_BB_tab(double l, int ni, int nj);
double C_ggl_TATT_tab(double l, int ni, int nj);
double w_gamma_t_TATT(int nt,int ni, int nj); //G-G lensing, lens bin ni, source bin nj, including IA contamination if like.IA = 3
double xi_pm_TATT(int pm, int nt, int ni, int nj); //shear tomography correlation functions, including IA contamination if like.IA = 3


/* NLA/TA amplitude C1, nz argument only need if per-bin amplitude*/
double C1_TA(double a, double nz){
	if (like.IA ==3 || like.IA ==5){
		return -nuisance.A_z[(int)nz]*cosmology.Omega_m*nuisance.c1rhocrit_ia*growfac(0.9999)/growfac(a);
	}
	return -cosmology.Omega_m*nuisance.c1rhocrit_ia*growfac(0.9999)/growfac(a)*nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia),nuisance.eta_ia);
}
/* TA source bias parameter, nz argument only need if per-bin amplitude*/
double b_TA(double a, double nz){
	if (like.IA ==5){
		return nuisance.b_ta_z[(int)nz];
	}	
	return nuisance.b_ta_z[0];
}
/* TT amplitude C2, nz argument only need if per-bin amplitude*/
double C2_TT(double a, double nz){
	if (like.IA == 5){
		return nuisance.A2_z[(int)nz]*cosmology.Omega_m*nuisance.c1rhocrit_ia*pow(growfac(0.9999)/growfac(a),2.0);
	}
	return cosmology.Omega_m*nuisance.c1rhocrit_ia*pow(growfac(0.9999)/growfac(a),2.0)*nuisance.A2_ia*pow(1./(a*nuisance.oneplusz0_ia),nuisance.eta_ia_tt);
}


double int_for_C_shear_shear_IA_EE(double a, void *params){
  double res, ell, fK, k,ws1,ws2,wk1,wk2, norm,C1,C1_2,C2,C2_2,b_ta,b_ta_2;
  double *ar = (double *) params;
  if (a >= 1.0) error("a>=1 in int_for_C_II");
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  ws1 = W_source(a,ar[0]);
  ws2 = W_source(a,ar[1]);
  wk1 = W_kappa(a,fK,ar[0]);
  wk2 = W_kappa(a,fK,ar[1]);

  C1 = C1_TA(a,ar[0]); b_ta = b_TA(a,ar[0]); C2 = C2_TT(a,ar[0]);
  /*calculate C1_1 to allow for per-bin NLA*/
  C1_2 = C1_TA(a,ar[1]); 
  /* per-bin TATT not supported yet*/
  //b_ta_2 = b_TA(a,ar[1]); C2_2 = C2_TT(a,ar[1]);

  /*GG cosmic shear */
  res = wk1*wk2*Pdelta(k,a);
  if (C1 || C2){
  	/*II contribution */
  	res += ws1*ws2*TATT_II_EE(k,a,-sqrt(C1*C1_2),C2,b_ta);
  	/*GI contribution */
  	res += ws1*wk2*TATT_GI_E(k,a,C1,C2,b_ta)+ws2*wk1*TATT_GI_E(k,a,C1_2,C2,b_ta);
  }
  return res*dchi_da(a)/fK/fK;
}

double int_for_C_shear_shear_IA_BB(double a, void *params){
  double res = 0., ell, fK, k,ws1,ws2,wk1,wk2, norm,C1,C1_2,C2,C2_2,b_ta,b_ta_2;
  double *ar = (double *) params;
  if (a >= 1.0) error("a>=1 in int_for_C_II");
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  ws1 = W_source(a,ar[0]);
  ws2 = W_source(a,ar[1]);
  wk1 = W_kappa(a,fK,ar[0]);
  wk2 = W_kappa(a,fK,ar[1]);

  C1 = C1_TA(a,ar[0]); b_ta = b_TA(a,ar[0]); C2 = C2_TT(a,ar[0]);
  /* per-bin TATT not supported yet */
  //  C1_2 = C1_TA(a,ar[1]); b_ta_2 = b_TA(a,ar[1]); C2_2 = C2_TA(a,ar[1]);
  if (b_ta || C2){res = ws1*ws2*TATT_II_BB(k,a,C1,C2,b_ta);}
  return res*dchi_da(a)/fK/fK;
}


double int_for_C_ggl_IA_TATT(double a, void *params){
  double res, ell, fK, k,w_density,w_mag, ws,wk, b1,b2,bs2,C1,C2,b_ta;
  double *ar = (double *) params;
  if (a >= 1.0) error("a>=1 in int_for_C_II");
  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;

  ws = W_source(a,ar[1]);
  wk = W_kappa(a,fK,ar[1]);
  C1 = C1_TA(a,ar[1]); b_ta = b_TA(a,ar[1]); C2 = C2_TT(a,ar[1]);

  w_density = W_HOD(a,ar[0]);
  w_mag = W_mag(a,fK,ar[0])*gbias.b_mag[(int)ar[0]];

  b1 = gbias.b1_function(1./a-1.,(int)ar[0]);
  b2 = gbias.b2[(int)ar[0]];
  bs2 = gbias.bs2[(int)ar[0]];
  double g4 = pow(growfac(a)/growfac(1.0),4.);
  double Pnl = Pdelta(k,a);

  double P_1loop =b1*Pnl;
  if (w_density*b2 !=0){
  	P_1loop += g4*(0.5*b2*PT_d1d2(k)+0.5*bs2*PT_d1s2(k)+0.5*b3nl_from_b1(b1)*PT_d1d3(k));
  }

  /*1-loop P_gm ggl terms*/
  res = w_density*wk*P_1loop;
  /* lens magnification x G term*/
  res += w_mag*wk*Pnl;
  /* (linear bias lens density + lens magnification) x TATT_GI terms*/
  if (C1 || C2) res += (b1*w_density+w_mag)*ws*TATT_GI_E(k,a,C1,C2,b_ta);
  return res*dchi_da(a)/fK/fK;
}

double C_EE_TATT(double l, int ni,int  nj){
  double array[3] = {(double) ni, (double) nj, l};
  double EE = int_gsl_integrate_low_precision(int_for_C_shear_shear_IA_EE,(void*)array,fmax(amin_source(ni),amin_source(nj)),amax_source(ni),NULL,1000);
  return EE;
}


double C_BB_TATT(double l, int ni, int nj){
  double array[3] = {(double) ni, (double) nj, l};
  return int_gsl_integrate_low_precision(int_for_C_shear_shear_IA_BB,(void*)array,fmax(amin_source(ni),amin_source(nj)),fmin(amax_source_IA(ni),amax_source_IA(nj)),NULL,1000);
}

double C_ggl_TATT(double l, int nl, int ns)
{
  double array[3] = {(double) nl, (double) ns, l};
  double gE = int_gsl_integrate_medium_precision(int_for_C_ggl_IA_TATT,(void*)array,amin_lens(nl),amax_lens(nl),NULL,1000);
  return gE;
}

/*************** look-up tables for angular correlation functions ***************/
/******************** all angles in radian!   ***********************************/
/******************** full-sky, bin-avergage  ***********************************/

double w_gamma_t_TATT(int nt, int ni, int nj){
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
		printf("cosmo2D_fullsky_TATT.c:w_gamma_t_TATT: like.Ntheta not initialized\nEXIT\n"); exit(1);
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

		for (nz = 0; nz <tomo.ggl_Npowerspectra; nz ++){
			for (l = 2; l < LMAX; l++){
				Cl[l]=C_ggl_TATT_tab(1.0*l,ZL(nz),ZS(nz));
			}
			for (i = 0; i < NTHETA; i++){
				w_vec[nz*like.Ntheta+i] =0;
				for (l = 2; l < LMAX; l++){
					w_vec[nz*like.Ntheta+i]+=Pl[i][l]*Cl[l];
				}
			}
		}
		update_cosmopara(&C);
		update_galpara(&G);
		update_nuisance(&N);
	}
	return w_vec[N_ggl(ni,nj)*like.Ntheta+nt];  
}


double xi_pm_TATT(int pm, int nt, int ni, int nj) //shear tomography correlation functions
{
	static int LMAX = 100000;
	static double **Glplus =0;
	static double **Glminus =0;
	static double *Cl_EE =0,*Cl_BB =0;
	static double *xi_vec_plus =0;
	static double *xi_vec_minus =0;
	static cosmopara C;
	static nuisancepara N;

	int i,l,nz;
	if (like.Ntheta == 0){
		printf("cosmo2D_fullsky_TATT.c:xi_pm_TATT: like.theta not initialized\nEXIT\n"); exit(1);
	}
	if (Glplus ==0){
		Glplus =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
		Glminus =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
		Cl_EE = create_double_vector(0,LMAX-1);
		Cl_BB = create_double_vector(0,LMAX-1);
		xi_vec_plus = create_double_vector(0,tomo.shear_Npowerspectra*like.Ntheta-1);
		xi_vec_minus = create_double_vector(0,tomo.shear_Npowerspectra*like.Ntheta-1);
		double *xmin, *xmax, *Pmin, *Pmax, *dPmin, *dPmax;
		xmin= create_double_vector(0, like.Ntheta-1);
		xmax= create_double_vector(0, like.Ntheta-1);
		double logdt=(log(like.vtmax)-log(like.vtmin))/like.Ntheta;
		for(i=0; i<like.Ntheta ; i++){
			xmin[i]=cos(exp(log(like.vtmin)+(i+0.0)*logdt));
			xmax[i]=cos(exp(log(like.vtmin)+(i+1.0)*logdt));
		}
		Pmin= create_double_vector(0, LMAX+1);
		Pmax= create_double_vector(0, LMAX+1);
		dPmin= create_double_vector(0, LMAX+1);
		dPmax= create_double_vector(0, LMAX+1);
		for (i = 0; i<like.Ntheta; i ++){
			double x = cos(like.theta[i]);
			gsl_sf_legendre_Pl_deriv_array(LMAX, xmin[i],Pmin,dPmin);
			gsl_sf_legendre_Pl_deriv_array(LMAX, xmax[i],Pmax,dPmax);
			for (int l = 3; l < LMAX; l ++){
				/*double plm = gsl_sf_legendre_Plm(l,2,x);
				double plm_1 = gsl_sf_legendre_Plm(l-1,2,x);
				Glplus[i][l] = (2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))
				*(plm*((4-l+2.*x*(l-1))/(1-x*x)-l*(l+1)/2)
				+plm_1*(l-1,2,x)*(l+2)*(x-2)/(1-x*x));


				Glminus[i][l] = (2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))
				*(plm*(l,2,x)*((4-l-2.*x*(l-1))/(1-x*x)-l*(l+1)/2)
				+plm_1*(l-1,2,x)*(l+2)*(x+2)/(1-x*x));*/

				Glplus[i][l] =(2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(

				-l*(l-1.)/2*(l+2./(2*l+1)) * (Pmin[l-1]-Pmax[l-1])
				-l*(l-1.)*(2.-l)/2         * (xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
				+l*(l-1)/(2*l+1)           * (Pmin[l+1]-Pmax[l+1])

				+(4-l)   * (dPmin[l]-dPmax[l])
				+(l+2)   * (xmin[i]*dPmin[l-1] - xmax[i]*dPmax[l-1] - Pmin[l-1] + Pmax[l-1])

				+2*(l-1) * (xmin[i]*dPmin[l]   - xmax[i]*dPmax[l]   - Pmin[l] + Pmax[l])
				-2*(l+2) * (dPmin[l-1]-dPmax[l-1])

				)/(xmin[i]-xmax[i]);           

				Glminus[i][l] =(2.*l+1)/(2.*M_PI*l*l*(l+1)*(l+1))*(

				-l*(l-1.)/2*(l+2./(2*l+1)) * (Pmin[l-1]-Pmax[l-1])
				-l*(l-1.)*(2.-l)/2         * (xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
				+l*(l-1)/(2*l+1)           * (Pmin[l+1]-Pmax[l+1])

				+(4-l)   * (dPmin[l]-dPmax[l])
				+(l+2)   * (xmin[i]*dPmin[l-1] - xmax[i]*dPmax[l-1] - Pmin[l-1] + Pmax[l-1])

				-2*(l-1) * (xmin[i]*dPmin[l]   - xmax[i]*dPmax[l]   - Pmin[l] + Pmax[l])
				+2*(l+2) * (dPmin[l-1]-dPmax[l-1])

				)/(xmin[i]-xmax[i]);

			}
		}
		free_double_vector(xmin,0,like.Ntheta-1);
		free_double_vector(xmax,0,like.Ntheta-1);
		free_double_vector(Pmin,0,LMAX+1);
		free_double_vector(Pmax,0,LMAX+1);
		free_double_vector(dPmin,0,LMAX+1);
		free_double_vector(dPmax,0,LMAX+1);

	}
	if (recompute_shear(C,N)){

		update_cosmopara(&C); update_nuisance(&N);
		for (nz = 0; nz <tomo.shear_Npowerspectra; nz ++){
			for (l = 2; l < LMAX; l++){
				Cl_EE[l]=C_EE_tab(1.0*l,Z1(nz),Z2(nz));
				Cl_BB[l]=C_BB_tab(1.0*l,Z1(nz),Z2(nz));
			}
			for (i = 0; i < like.Ntheta; i++){
				xi_vec_plus[nz*like.Ntheta+i] =0;
				xi_vec_minus[nz*like.Ntheta+i] =0;
				for (l = 2; l < LMAX; l++){
					xi_vec_plus[nz*like.Ntheta+i]+=Glplus[i][l]*(Cl_EE[l]+Cl_BB[l]);
					xi_vec_minus[nz*like.Ntheta+i]+=Glminus[i][l]*(Cl_EE[l]-Cl_BB[l]);
				}
			}
		}
	}
	if (pm> 0) return xi_vec_plus[N_shear(ni,nj)*like.Ntheta + nt];
	return xi_vec_minus[N_shear(ni,nj)*like.Ntheta + nt];
}


/************* tabulated angular power spectra *************/

double C_EE_tab(double l, int ni, int nj)  //shear power spectrum of source galaxies in bins ni, nj
{
  static cosmopara C;
  static nuisancepara N;
  
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  if (ni < 0 || ni >= tomo.shear_Nbin ||nj < 0 || nj >= tomo.shear_Nbin){
    printf("C_shear_shear_EE(l,%d,%d) outside tomo.shear_Nbin range\nEXIT\n",ni,nj); exit(1);
  }
  
  if (recompute_shear(C,N)){
    //printf("calculating C_shear_shear_IA_tab  %e %e %e %e %e\n", nuisance.A_z[0], nuisance.A_z[1], nuisance.A_z[2], nuisance.A_z[3], nuisance.A_z[4]);
    if (table==0) {
      table   = create_double_matrix(0, tomo.shear_Npowerspectra-1, 0, Ntable.N_ell-1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell);
    }
    
    double llog;
    int i,k;
    
    for (k=0; k<tomo.shear_Npowerspectra; k++) {
      llog = logsmin;
      for (i=0; i<Ntable.N_ell; i++, llog+=ds) {
        table[k][i]= log(C_EE_TATT(exp(llog),Z1(k),Z2(k)));
      }
    }
    update_cosmopara(&C); update_nuisance(&N); 
  }
  double f1 = exp(interpol_fitslope(table[N_shear(ni,nj)], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1.));
  if (isnan(f1)){f1 = 0.;}
  // printf("%le %d %d %le\n", l, ni, nj, f1);
  return f1;
}

double C_BB_tab(double l, int ni, int nj)  //shear power spectrum of source galaxies in bins ni, nj
{
  static cosmopara C;
  static nuisancepara N;
  
  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  if (ni < 0 || ni >= tomo.shear_Nbin ||nj < 0 || nj >= tomo.shear_Nbin){
    printf("C_shear_shear_BB(l,%d,%d) outside tomo.shear_Nbin range\nEXIT\n",ni,nj); exit(1);
  }
  
  if (recompute_shear(C,N)){
    //printf("calculating C_shear_shear_IA_tab  %e %e %e %e %e\n", nuisance.A_z[0], nuisance.A_z[1], nuisance.A_z[2], nuisance.A_z[3], nuisance.A_z[4]);
    if (table==0) {
      table   = create_double_matrix(0, tomo.shear_Npowerspectra-1, 0, Ntable.N_ell-1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell);
    }
    
    double llog;
    int i,k;
    
    for (k=0; k<tomo.shear_Npowerspectra; k++) {
      llog = logsmin;
      for (i=0; i<Ntable.N_ell; i++, llog+=ds) {
        table[k][i]= log(C_BB_TATT(exp(llog),Z1(k),Z2(k)));
      }
    }
    update_cosmopara(&C); update_nuisance(&N); 
  }
  double f1 = exp(interpol_fitslope(table[N_shear(ni,nj)], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1.));
  if (isnan(f1)){f1 = 0.;}
  // printf("%le %d %d %le\n", l, ni, nj, f1);
  return f1;
}


double C_ggl_TATT_tab(double l, int ni, int nj)  //G-G lensing power spectrum, lens bin ni, source bin nj
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  
  static double **table, *sig;
  static int osc[100];
  static double ds = .0, logsmin = .0, logsmax = .0;
  
  if (ni < 0 || ni >= tomo.clustering_Nbin ||nj < 0 || nj >= tomo.shear_Nbin){
    printf("C_ggl_TATT_tab(l,%d,%d) outside tomo.X_Nbin range\nEXIT\n",ni,nj); exit(1);
  }
  
  if (recompute_ggl(C,G,N,ni)){
	    //printf("calculating C_ggl_IA_tab  %e %e %e %e %e\n", nuisance.A_z[0], nuisance.A_z[1], nuisance.A_z[2], nuisance.A_z[3], nuisance.A_z[4]);
    if (table==0){
      table   = create_double_matrix(0, tomo.ggl_Npowerspectra-1, 0, Ntable.N_ell-1);
      sig = create_double_vector(0,tomo.ggl_Npowerspectra-1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell);
    }
    int i,k;
    double res,llog;
    
    for (k=0; k<tomo.ggl_Npowerspectra; k++) {
      llog = logsmin;

      sig[k] = 1.;
      osc[k] = 0;
      res = C_ggl_TATT(500.,ZL(k),ZS(k));
      if (res < 0){sig[k] = -1.;}
      for (i=0; i<Ntable.N_ell; i++, llog+=ds) {
        table[k][i] = C_ggl_TATT(exp(llog),ZL(k),ZS(k));
        if (res*sig[k] <0.) {
          osc[k] = 1;
        }
      }
      if (osc[k] == 0){
        for(i = 0; i < Ntable.N_ell; i++){
            res = table[k][i];
            table[k][i] = log(sig[k]*res);}
      }
      
    }

    update_cosmopara(&C); update_nuisance(&N); update_galpara(&G);
    
  }
  int k = N_ggl(ni,nj);
  double f1 = 0.;
  if(test_zoverlap(ni,nj) && osc[k] ==0 ){f1 = sig[k]*exp(interpol_fitslope(table[k], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1.));}
  if(test_zoverlap(ni,nj) && osc[k] ==1 ){f1 = interpol_fitslope(table[k], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1.);}
  if (isnan(f1)){f1 = 0;}
  return f1;
}
