/* bias evolution with redshift*/
double b1_per_bin(double z, int ni); //model b1 using one non-evoling parameter per redshift bin
double b1_per_bin_evolv(double z, int ni); //model b1 using one parameter per redshift bin, power-law evolution within each bin
double b1_per_bin_pass_evolv(double z, int ni); //model b1 using one parameter per redshift bin + passive evolution
double b1_growth_scaling(double z, int ni); //model b1 assuming b1(z) = b_{1,0}*G(z)
double bgal_z(double z, int ni); //bias evolution within redshift bin, used by clustering/G-G-lensing routines without HOD modeling

double b2_from_b1(double b1); //fitting function for b_2(b_1)
double bs2_from_b1(double b1); //theory prediction for b_s2(b_1)
double b3nl_from_b1 (double b1); // theory prediction b3nl(b1) = b1-1
/* HOD model routines */
double n_c(double mh, double a, int nz); //central galaxy occupation as a function of halo mass (in Msun/h) in bin nz at scale factor a
double n_s(double mh, double a, int nz); //satellite galaxy occupation function
double ngal(int nz,double a); //galaxy number density of HOD model  in bin nz at scale factor a (in c/H0 = 1 units)
double bgal(int nz,double a); //mean galaxy bias of HOD model in bin nz at scale factor a
double mmean(int nz,double a);//mean host halo mass of HOD model in bin nz at scale factor a
double fsat(int nz,double a);//satellite fraction of HOD model

/* HOD model galaxy power spectrum + galaxy-matter cross spectrum */
double u_g(double k, double m, double a,int nz);//Fourier transform of normalized galaxy density profile
double int_for_G02 (double logm, void *para);
double G02 (double k, double a, int nz);
double int_GM02 (double logm, void *para);
double GM02 (double k, double a, int nz);
double int_G11 (double logm, void *para);
double G11 (double k, double a, int nz);
double P_gg (double k, double a,int nz);//galaxy-galaxy power spectrum based on HOD model in bin nz
double P_gm (double k, double a,int nz);//galaxy-matter power spectrum based on HOD model in bin nz
void set_HOD(int n); //set HOD parameters


/*********************** galaxy bias & HOD routines *********************************/
//fitting formula from Lazeyras et al. 2016 (Eq. 5.2)
//https://arxiv.org/abs/1511.01096
double b2_from_b1(double b1){
  return 0.412-2.143*b1+0.929*b1*b1+0.008*b1*b1*b1;
}
//e.g. https://arxiv.org/pdf/1405.1447v4.pdf
double bs2_from_b1 (double b1){
  return -4./7.*(b1-1.0);
}

double b3nl_from_b1 (double b1){
  return (b1-1.0);
}

double b1_per_bin_evolv(double z, int ni){
  return gbias.b[ni]*pow((1.+z)/(1.+zmean(ni)),1.0);
}

double b1_per_bin(double z, int ni){
  return gbias.b[ni];
}
double b1_per_bin_pass_evolv(double z, int ni){
   double z_evolv_passiv = growfac (1./(z+1.))/growfac (1./(1.+zmean(ni)));
   return (gbias.b[ni]-1.)/z_evolv_passiv +1.;
}
double b1_growth_scaling(double z, int ni){
  return gbias.b[0]/(growfac (1./(z+1.))/growfac (1.));
}
double b1_powerlaw(double z, int ni){
  return gbias.b[0]*pow(1+z,gbias.b[1]);
}
double bgal_z(double z, int ni){ //bias evolution within redshift bin, used by clustering/G-G-lensing routines without HOD modeling
  //change this into desired redshift evolution as function of (z, z_pivot = gbias[ni][1]), with z_evolv(z_pivot) =1
  //e.g. z_evolv = pow((1+z)/(1+gbias[ni][1]),0.5);
  //passive evolution  following Fry 1996: b(z) = (b(z_pivot)-1)*D(z_pivot)/D(z) + 1

  //no tomography ->  bias = 1
  if (ni==-1) return 1.0;
  // new default: use b1_function specified in init - if different from bgal_z, to avoid circular call
  if ((gbias.b1_function) && (gbias.b1_function != &bgal_z)){return gbias.b1_function(z,ni);}
  //**********old options, for backward compatability********
  // very broad redshift distribution -> no bias evolution model
  if (redshift.clustering_photoz == 4) return gbias.b[ni];
  //HOD parameters specified -> use HOD-derived bias for large-scale bias
  else if (gbias.hod[ni][0] > 10 && gbias.hod[ni][0] < 16){ return bgal(ni,1./(1+z));}
  // old default: passive evolution
  else{
   double z_evolv_passiv = growfac (1./(z+1.))/growfac (1./(1.+0.5*(tomo.clustering_zmin[ni]+tomo.clustering_zmax[ni])));
   return (gbias.b[ni]-1.)/z_evolv_passiv +1.;
  }
}

/******************** standard luminosity threshold HOD routines *******************/
//HOD parameterization as in Eq. (7) of Zehavi et al. 2011

double n_c(double mh, double a, int nz){
  if (gbias.hod[nz][0] < 10 || gbias.hod[nz][0] >16){
    printf("HOD.c:n_c(): HOD parameters in redshift bin %d not set\n",nz);
    exit(EXIT_FAILURE);
  }
	return 0.5*(1.0+gsl_sf_erf((log10(mh)-gbias.hod[nz][0])/gbias.hod[nz][1]));
}
double n_s(double mh, double a, int nz){
	double ns = n_c(mh,a,nz)*pow((mh-pow(10.0,gbias.hod[nz][3]))/pow(10.0,gbias.hod[nz][2]),gbias.hod[nz][4]);
	if (ns >0){return ns;}
	else {return 1.e-15;}
}
double f_c(int nz){
  if (gbias.hod[nz][5]){
    return gbias.hod[nz][5];
  }
  return 1.0;
}
/****  derived HOD quantities *****/

double int_bgal (double m, void *params){
	double *array  = (double*) params;
	double a= array[0];
	return massfunc(exp(m),a)*B1(exp(m),a)*(f_c((int) array[1])*n_c(exp(m), a, (int) array[1])+n_s(exp(m),a, (int) array[1]))*exp(m);
}
double int_fsat (double m, void *params){
	double *array  = (double*) params;
	double a= array[0];
	return massfunc(exp(m),a)*(n_s(exp(m), a,(int) array[1]))*exp(m);
}
double int_ngal (double m, void *params){
	double *array  = (double*) params;
	double a= array[0];
	return massfunc(exp(m),a)*(f_c((int) array[1])*n_c(exp(m),a, (int) array[1])+n_s(exp(m), a,(int) array[1]))*exp(m);
}
double int_mmean (double m, void *params){
	double *array  = (double*) params;
	double a= array[0];
	return massfunc(exp(m),a)*(f_c((int) array[1])*n_c(exp(m),a, (int) array[1])+n_s(exp(m),a, (int) array[1]))*exp(m)*exp(m);
}
double ngal(int nz,double a){
  static cosmopara C;
  static galpara G;


  static double **table;

  static double da = 0.0;
  double aa;
  int i,j;
  double array[2];

  if (recompute_cosmo3D(C) || recompute_galaxies(G,nz))
    {
    if (table==0) {
      table   = create_double_matrix(0, tomo.clustering_Nbin-1, 0, Ntable.N_a_halo-1);
      da = (1./(redshift.clustering_zdistrpar_zmin+1.)-1./(redshift.clustering_zdistrpar_zmax+1.))/(Ntable.N_a_halo-1);
    }

    for (j=0;j<=tomo.clustering_Nbin-1;j++) {
      array[1]=(double) j;
      aa = 1./(redshift.clustering_zdistrpar_zmax+1.);
      for (i=0;i<Ntable.N_a_halo;i++,aa+=da) {
        if (aa >= 1./(1+tomo.clustering_zmax[j])-da && aa <=1./(1+tomo.clustering_zmin[j])+da){
          array[0] = aa;
          table[j][i] = int_gsl_integrate_medium_precision(int_ngal, (void*)array, log(10.)*(gbias.hod[nz][0]-2.),log(limits.M_max),NULL, 5000);
        }
      }
    }
    update_cosmopara(&C); update_galpara(&G);
    }
  if (a<1./(redshift.clustering_zdistrpar_zmax+1.) || nz >= tomo.clustering_Nbin) return 0.0;
  return interpol(table[nz], Ntable.N_a_halo, 1./(redshift.clustering_zdistrpar_zmax+1.), 1./(redshift.clustering_zdistrpar_zmin+1.), da, a, 1.0, 1.0);
}
double bgal(int nz,double a){
  static cosmopara C;
  static galpara G;


  static double **table;

  static double da = 0.0;
  double aa;
  int i,j;
  double array[2];

  if (recompute_cosmo3D(C) || recompute_galaxies(G,nz))
    {
    if (table==0) {
      table   = create_double_matrix(0, tomo.clustering_Nbin-1, 0, Ntable.N_a_halo-1);
      da = (1./(redshift.clustering_zdistrpar_zmin+1.)-1./(redshift.clustering_zdistrpar_zmax+1.))/(Ntable.N_a_halo-1);
    }

    for (j=0;j<=tomo.clustering_Nbin-1;j++) {
      array[1]=(double) j;
      aa = 1./(redshift.clustering_zdistrpar_zmax+1.);
      for (i=0;i<Ntable.N_a_halo;i++,aa+=da) {
        if (aa >= 1./(1+tomo.clustering_zmax[j])-2.*da && aa <=1./(1+tomo.clustering_zmin[j])+2.*da){
          array[0] = aa;
          table[j][i] = int_gsl_integrate_medium_precision(int_bgal, (void*)array, log(10.)*(gbias.hod[nz][0]-2.),log(limits.M_max),NULL, 5000)/int_gsl_integrate_medium_precision(int_ngal, (void*)array, log(10.)*(gbias.hod[nz][0]-2.),log(limits.M_max),NULL, 5000);
        }
      }
    }
    update_cosmopara(&C); update_galpara(&G);
    }
  if (a<1./(redshift.clustering_zdistrpar_zmax+1.) || nz >= tomo.clustering_Nbin) return 0.0;
  return interpol(table[nz], Ntable.N_a_halo, 1./(redshift.clustering_zdistrpar_zmax+1.), 1./(redshift.clustering_zdistrpar_zmin+1.), da, a, 1.0, 1.0);
}
double mmean(int nz,double a){
	double array[2]  ={a, (double) nz};
	return int_gsl_integrate_medium_precision(int_mmean, (void*)array, log(10.)*(gbias.hod[nz][0]-2.),log(limits.M_max),NULL, 5000)/ngal(nz,a);
}
double fsat(int nz,double a){
	double array[2]  ={a, (double) nz};
	return int_gsl_integrate_medium_precision(int_fsat, (void*)array, log(10.)*(gbias.hod[nz][0]-2.),log(limits.M_max),NULL, 5000)/ngal(nz,a);;
}



/********************* galaxy power spectrum from halo model + HOD ***************************/
/* ******************** no halo exclusion or radial biasing so far ***************************/

double u_g(double k, double m, double a,int nz)// Fourier transformed normalized galaxy density profile, NFW with rescaled concentraion so far
{
  // analytic FT of NFW profile, from Cooray & Sheth 01
  if (gbias.cg[nz] <=0){printf("galaxy concentration parameter ill-defined in z-bin %d, exiting.\n",nz); exit(EXIT_FAILURE);}
  return u_nfw_c(conc(m,a)*gbias.cg[nz],k,m,a);
}


/**** 1-halo galaxy-galaxy spectrum ******/
double int_for_G02 (double logm, void *para){
	double *array = (double *) para;
	double m = exp(logm);
	double u,ns,k = array[0],a = array[1];
	int nz = (int) array[2];
	u = u_g(k,m,a,nz);
	ns = n_s(m,a,nz);
	return massfunc(m,a)*m*(u*u*ns*ns+2.*u*ns*n_c(m,a,nz)*f_c(nz));
}

double G02 (double k, double a, int nz){//needs to be devided by ngal(nz, a)^2
	double array[3]  ={k,a, (double) nz};
	return int_gsl_integrate_medium_precision(int_for_G02, (void*)array, log(1.e+8),log(limits.M_max),NULL, 5000);
}
/**** 1-halo galaxy-matter spectrum ******/
double int_GM02 (double logm, void *para){
	double *array = (double *) para;
	double m = exp(logm);
	double k = array[0],a = array[1];
	int nz = (int) array[2];
	return massfunc(m,a)*m*m/(cosmology.rho_crit*cosmology.Omega_m)*u_nfw_c(conc(m,a),k,m,a)*(u_g(k,m,a,nz)*n_s(m,a,nz)+n_c(m,a,nz)*f_c(nz));//u_nfw_c(conc(m,a),k,m,a)
}

double GM02 (double k, double a, int nz){//needs to be devided by ngal(nz, a)
	double array[3]  ={k,a, (double) nz};
	return int_gsl_integrate_medium_precision(int_GM02, (void*)array, log(10.)*(gbias.hod[nz][0]-1.),log(limits.M_max),NULL, 5000);
}
/**** 2-halo term routines *****/
double int_G11 (double logm, void *para){
	double *array = (double *) para;
	double u, m = exp(logm),k = array[0],a = array[1];
  int nz = (int) array[2];
	u = u_g(k,m,a,nz);
	return m*massfunc(m,a)*(u*n_s(m,a,nz)+n_c(m,a,nz)*f_c(nz))*B1(m,a);
}
double G11 (double k, double a, int nz){ //needs to be devided by ngal(nz, a)
	double array[3]  ={k,a, (double) nz};
	return int_gsl_integrate_medium_precision(int_G11, (void*)array, log(10.)*(gbias.hod[nz][0]-1.),log(limits.M_max),NULL, 5000);
}
/**** power spectra ****/
double P_gg (double k, double a,int nz){ //galaxy-galaxy power spectrum based on HOD model in bin nz

  static double logkmin = 0., logkmax = 0., dk = 0., da = 0.,amin =0., amax = 0.;
  static int N_a = 15, N_k_nlin = 50, NZ = -1;
  static double M0 = 0.0;
  static cosmopara C;

  static double **table_P_gg=0;

  double klog,val,aa,kk;
  int i,j;

  if (NZ != nz || M0 != gbias.hod[nz][0]|| recompute_cosmo3D(C)){ //extend this by halo model parameters if these are sampled independently of cosmology parameters
    M0 = gbias.hod[nz][0];
    NZ = nz;
    update_cosmopara(&C);
    if (table_P_gg == 0){
      table_P_gg = create_double_matrix(0, N_a-1, 0, N_k_nlin-1);
    }
    amax = amax_lens(nz);
    amin = amin_lens(nz);
    da = (amax - amin)/(N_a-1.);
    logkmin = log(limits.k_min_cH0);
    logkmax = log(limits.k_max_cH0/10.);
//    printf("tabulating P_gg %d, %e %e, %e %e\n", NZ,exp(logkmin),exp(logkmax), amin,amax);
    dk = (logkmax - logkmin)/(N_k_nlin-1.);
    aa= amin;
    for (i=0; i<N_a; i++, aa +=da) {
      klog  = logkmin;
      for (j=0; j<N_k_nlin; j++, klog += dk) {
        kk = exp(klog);
        table_P_gg[i][j] = log(Pdelta(kk,aa)*pow(bgal(nz,aa),2.0)+G02(kk,aa,nz)/pow(ngal(nz,aa),2.0));
      }
    }
  }
  klog = log(k);
  if (klog <= logkmin || klog >logkmax-dk){return 0.;}
  if (a < amin || a > amax){return 0.;}
  val = interpol2d(table_P_gg, N_a, amin, amax, da, a, N_k_nlin, logkmin, logkmax, dk, klog, 0.0, 0.0);
  return exp(val);
}
//	return (G02(k,a,nz) +  Pdelta(k,a)*pow(G11(k,a,nz),2.0))/pow(ngal(nz,a),2.0);

double P_gm (double k, double a,int nz){//galaxy-matter power spectrum based on HOD model in bin nz
  static double logkmin = 0., logkmax = 0., dk = 0., da = 0.,amin =0., amax = 0.;
  static int N_a = 15, N_k_nlin = 50, NZ = -1;
  static double M0 = 0.0;
  static cosmopara C;

  static double **table_P_gg=0;

  double klog,val,aa,kk;
  int i,j;

  if (NZ != nz || M0 != gbias.hod[nz][0]|| recompute_cosmo3D(C)){ //extend this by halo model parameters if these are sampled independently of cosmology parameters
    M0 = gbias.hod[nz][0];
    NZ = nz;
    update_cosmopara(&C);
    if (table_P_gg == 0){
      table_P_gg = create_double_matrix(0, N_a-1, 0, N_k_nlin-1);
    }
    amax = amax_lens(nz);
    amin = amin_lens(nz);
    da = (amax - amin)/(N_a-1.);
    logkmin = log(limits.k_min_cH0);
    logkmax = log(limits.k_max_cH0/10.);
//    printf("tabulating P_gg %d, %e %e, %e %e\n", NZ,exp(logkmin),exp(logkmax), amin,amax);
    dk = (logkmax - logkmin)/(N_k_nlin-1.);
    aa= amin;
    for (i=0; i<N_a; i++, aa +=da) {
      klog  = logkmin;
      for (j=0; j<N_k_nlin; j++, klog += dk) {
        kk = exp(klog);
        table_P_gg[i][j] = log(Pdelta(kk,aa)*bgal(nz,aa)+GM02(kk,aa,nz)/ngal(nz,aa));
      }
    }
  }
  klog = log(k);
  if (klog <= logkmin || klog >logkmax-dk){return 0.;}
  if (a < amin || a > amax){return 0.;}
  val = interpol2d(table_P_gg, N_a, amin, amax, da, a, N_k_nlin, logkmin, logkmax, dk, klog, 0.0, 0.0);
  return exp(val);
} //	return (GM02(k,a,nz) +  Pdelta(k,a)*G11(k,a,nz)*I1j(1,k,0,0,a))/ngal(nz,a);



void set_HOD(int n){ //n >=0: set HOD parameters in redshift bin n; n = -1: unset HOD parameters (code then uses linear bias + non-linear matter power spectrum instead of halo model)
  int i;
  double z = zmean(n);
  /*set HOD parameters, parameterization of Zehavi et al. */
  /*these example values from Coupon et al. 2012 for red galaxies with M_r < -21.8 (Table B.2)
   hod[zi][] ={lg(M_min), sigma_{lg M}, lg M_1, lg M_0, alpha, f_c}
   gbias.cg[] ={f_g} (shift of concentration parameter: c_g(M) = f_g c(M))*/
  switch (n)
  {
    case -1:
    for (i = 0; i <10; i++){
      gbias.hod[i][0] = 0.;
      gbias.hod[i][1] = 0.;
      gbias.hod[i][2] = 0.;
      gbias.hod[i][3] = 0.;
      gbias.hod[i][4] = 0.;
      gbias.hod[i][5] = 0.;
    }
    printf("unset HOD parameters; code uses non-linear matter power spectrum + linear bias\n");
    break;

    case 0:
    gbias.hod[0][0] =13.17;
    gbias.hod[0][1] = 0.39;
    gbias.hod[0][2] =14.53;
    gbias.hod[0][3] =11.09;
    gbias.hod[0][4] = 1.27;
    gbias.hod[0][5] = 1.00;
    gbias.b[n] = bgal(n,1./(z+1));
    printf("HOD derived quantities in bin %d at <z> = %.2f: <n_g> = %e(h/Mpc)^3 <b_g> = %.2f lg(<M_h> h/M_sun) = %.3f f_sat = %.3f\n",n, z,ngal(n,1./(z+1))*pow(cosmology.coverH0,-3.0),gbias.b[n], log10(mmean(n,1./(z+1))), fsat(n,1./(z+1)));
    break;

    case 1:
    gbias.hod[1][0] =13.18;
    gbias.hod[1][1] = 0.30;
    gbias.hod[1][2] =14.47;
    gbias.hod[1][3] =10.93;
    gbias.hod[1][4] = 1.36;
    gbias.hod[1][5] = 1.00;
    gbias.b[n] = bgal(n,1./(z+1));
    printf("HOD derived quantities in bin %d at <z> = %.2f: <n_g> = %e(h/Mpc)^3 <b_g> = %.2f lg(<M_h> h/M_sun) = %.3f f_sat = %.3f\n",n, z,ngal(n,1./(z+1))*pow(cosmology.coverH0,-3.0),gbias.b[n], log10(mmean(n,1./(z+1))), fsat(n,1./(z+1)));
    break;

    case 2:
    gbias.hod[2][0] =12.96;
    gbias.hod[2][1] = 0.38;
    gbias.hod[2][2] =14.10;
    gbias.hod[2][3] =12.47;
    gbias.hod[2][4] = 1.28;
    gbias.hod[2][5] = 1.00;
    gbias.b[n] = bgal(n,1./(z+1));
    printf("HOD derived quantities in bin %d at <z> = %.2f: <n_g> = %e(h/Mpc)^3 <b_g> = %.2f lg(<M_h> h/M_sun) = %.3f f_sat = %.3f\n",n, z,ngal(n,1./(z+1))*pow(cosmology.coverH0,-3.0),gbias.b[n], log10(mmean(n,1./(z+1))), fsat(n,1./(z+1)));
    break;
    case 3:
    gbias.hod[3][0] =12.80;
    gbias.hod[3][1] = 0.35;
    gbias.hod[3][2] =13.94;
    gbias.hod[3][3] =12.15;
    gbias.hod[3][4] = 1.52;
    gbias.hod[3][5] = 1.00;
    gbias.b[n] = bgal(n,1./(z+1));
    printf("HOD derived quantities in bin %d at <z> = %.2f: <n_g> = %e(h/Mpc)^3 <b_g> = %.2f lg(<M_h> h/M_sun) = %.3f f_sat = %.3f\n",n, z,ngal(n,1./(z+1))*pow(cosmology.coverH0,-3.0),gbias.b[n], log10(mmean(n,1./(z+1))), fsat(n,1./(z+1)));

    break;

    case 4:
    //no information for higher redshift populations - copy 1<z<1.2 values
    gbias.hod[4][0] =12.80;
    gbias.hod[4][1] = 0.35;
    gbias.hod[4][2] =13.94;
    gbias.hod[4][3] =12.15;
    gbias.hod[4][4] = 1.52;
    gbias.hod[4][5] = 1.00;
    gbias.b[n] = bgal(n,1./(z+1));
    printf("HOD derived quantities in bin %d at <z> = %.2f: <n_g> = %e(h/Mpc)^3 <b_g> = %.2f lg(<M_h> h/M_sun) = %.3f f_sat = %.3f\n",n, z,ngal(n,1./(z+1))*pow(cosmology.coverH0,-3.0),gbias.b[n], log10(mmean(n,1./(z+1))), fsat(n,1./(z+1)));

    break;

    default:
    printf("no HOD parameters specified to initialize bin %d\n", n);
  }

  /*for sparse lens populations, update angular galaxy density*/
  survey.n_lens = 2.0;

  printf("\nUsing HOD model for lens population with n_lens = %.2f/arcmin^2.\nLens redshift distribution from redshift.clustering_REDSHIFT_FILE =%s\nSource distribution from redshift.shear_REDSHIFT_FILE =%s\n\n",survey.n_lens,redshift.clustering_REDSHIFT_FILE, redshift.shear_REDSHIFT_FILE);
}
