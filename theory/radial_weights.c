double W_kappa(double a, double fK, double nz);//complete lens efficiency weight
double W2_kappa(double a, double fK, double nz); //complete (lens efficiency)^2 weight for source clustering
double W_source(double a, double nz); //source redshift distribution (radial weight for IA,source clustering)
double W_gal(double a, double nz); //complete weight for galaxy statistics
double W_HOD(double a, double nz); //galaxy weigth without bias factor (for projecting P_gg instead of P_nl)

/* bias evolution with redshift*/
double b1_per_bin(double z, int ni); //model b1 using one non-evoling parameter per redshift bin
double b1_per_bin_evolv(double z, int ni); //model b1 using one parameter per redshift bin, power-law evolution within each bin
double b1_per_bin_pass_evolv(double z, int ni); //model b1 using one parameter per redshift bin + passive evolution
double b1_growth_scaling(double z, int ni); //model b1 assuming b1(z) = b_{1,0}*G(z)
double bgal_z(double z, int ni); //bias evolution within redshift bin, used by clustering/G-G-lensing routines without HOD modeling

double b2_from_b1(double b1); //fitting function for b_2(b_1)
double bs2_from_b1(double b1); //theory prediction for b_s2(b_1)
double b3nl_from_b1 (double b1); // theory prediction b3nl(b1) = b1-1

double C_cl_tomo_nointerp(double l, int ni, int nj);
double C_cl_tomo(double l, int ni, int nj);  //galaxy clustering power spectrum of galaxies in bins ni, nj

double MG_Sigma(double a)
{
//  double aa=a*a;
//  double omegam=cosmology.Omega_m/(aa*a);
  double omegav=omv_vareos(a);
  double hub=hoverh0(a);
  hub = hub*hub;

  return cosmology.MGSigma*omegav/hub/cosmology.Omega_v;
}

double dchi_da(double a){
  return 1./(a*a*hoverh0(a));
}
double W_kappa(double a, double fK, double nz){
  double wkappa = 1.5*cosmology.Omega_m*fK/a*g_tomo(a,(int)nz);
  if(cosmology.MGSigma != 0.){
    wkappa *= (1.+MG_Sigma(a));
  }
  return wkappa;
}
double W2_kappa(double a, double fK, double nz){
  double wkappa = pow(1.5*cosmology.Omega_m*fK/a,2.0)*g2_tomo(a,(int)nz);
  if(cosmology.MGSigma != 0.){
    wkappa *= pow((1.+MG_Sigma(a)),2.);
  }
  return wkappa;
}

double W_mag(double a, double fK, double nz){
  double wmag = 1.5*cosmology.Omega_m*fK/a*g_lens(a,(int)nz);
  if(cosmology.MGSigma != 0.){
    wmag *= (1.+MG_Sigma(a));
  }
  return wmag;
}

double W_gal(double a, double nz){
  double wgal = gbias.b1_function(1./a-1.,(int)nz)*pf_photoz(1./a-1.,(int)nz)*hoverh0(a);
  double wmag = gbias.b_mag[(int)nz]*1.5*cosmology.Omega_m*f_K(chi(a))/a*g_lens(a,(int)nz);
  if(cosmology.MGSigma != 0.){
    wmag *= (1.+MG_Sigma(a));
  }
  return wgal + wmag;
  // return wgal;
}

double W_source(double a, double nz){
  return zdistr_photoz(1./a-1.,(int)nz)*hoverh0(a);
}

double f_rsd (double aa){
  double gamma = 0.55;
  return pow(cosmology.Omega_m /(cosmology.Omega_m +omv_vareos(aa) *aa*aa*aa),gamma);
}
double W_RSD(double l, double a0, double a1, double nz){
  double w = (1+8.*l)/((2*l+1)*(2*l+1))*pf_photoz(1./a0-1.,(int)nz)*hoverh0(a0)*f_rsd(a0);
  w -= 4./(2*l+3)*sqrt((2*l+1.)/(2*l+3.))*pf_photoz(1./a1-1.,(int)nz)*hoverh0(a1)*f_rsd(a1);
  return w;
}

double W_HOD(double a, double nz){
  return pf_photoz(1./a-1.,(int)nz)*hoverh0(a);
}

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
  double z_evolv_passiv = growfac (1./(z+1.))/growfac (1./(1.+0.5*(tomo.clustering_zmin[ni]+tomo.clustering_zmax[ni])));
  return (gbias.b[ni]-1.)/z_evolv_passiv +1.;
}

/*************** other routines from cosmo2D_fourier.c*****************/
/***************** used in non-Limber evaluations *********************/
double int_for_C_cl_tomo(double a, void *params)
{
  double res,ell, fK, k;
  double *ar = (double *) params;
  if (a >= 1.0) error("a>=1 in int_for_C_cl_tomo");

  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;

  res=W_gal(a,ar[0])*W_gal(a,ar[1])*dchi_da(a)/fK/fK;
  res= res*Pdelta(k,a);
  return res;
}

double int_for_C_cl_tomo_b2(double a, void *params)
{
  double res,ell, fK, k;
  double *ar = (double *) params;
  double b1 = gbias.b1_function(1./a-1.,(int)ar[0]);
  double b2 = gbias.b2[(int)ar[0]];
  double bs2 = gbias.bs2[(int)ar[0]];
  double g4 = pow(growfac(a)/growfac(1.0),4.);
  if (a >= 1.0) error("a>=1 in int_for_C_cl_tomo");

  ell       = ar[2]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;

  double s4 = 0.;//PT_sigma4(k);
  res=W_HOD(a,ar[0])*W_HOD(a,ar[1])*dchi_da(a)/fK/fK;
  if(res){
    res= res*(b1*b1*Pdelta(k,a)+g4*(b1*b2*PT_d1d2(k)+0.25*b2*b2*(PT_d2d2(k)-2.*s4)+b1*bs2*PT_d1s2(k)+0.5*b2*bs2*(PT_d2s2(k)-4./3.*s4)+.25*bs2*bs2*(PT_s2s2(k)-8./9.*s4)+b1*b3nl_from_b1(b1)*PT_d1d3(k)));
  }
  res += (W_gal( a, ar[0])*W_mag(a, fK,ar[1])+ W_gal( a, ar[1])*W_mag(a, fK,ar[0]))*dchi_da(a)/fK/fK*Pdelta(k,a);
  return res;
}


double C_cl_tomo_nointerp(double l, int ni, int nj)  //galaxy clustering power spectrum of galaxy bins ni, nj
{ static int init =-1;
  double array[3] = {1.0*ni,1.0*nj,l};
  if (gbias.b2[ni] || gbias.b2[nj]){
    if (ni != nj){
        if (init ==-1){
          printf("\nCalled C_cl(l,z1=%d,z2=%d) with non-linear bias parameters set.\n",ni,nj);
          printf("Cross-clustering beyond linear bias for cross-tomography bins not yet supported.\n");
          printf("Use linear bias only for z1!=z2 clustering\n\n");
          init = 1;}
          return int_gsl_integrate_medium_precision(int_for_C_cl_tomo,(void*)array,fmax(amin_lens(ni),amin_lens(nj)),fmin(amax_lens(ni),amax_lens(nj)),NULL,1000);
        }
    return int_gsl_integrate_medium_precision(int_for_C_cl_tomo_b2,(void*)array,amin_lens(ni),amax_lens(ni),NULL,1000);
  }
  return int_gsl_integrate_medium_precision(int_for_C_cl_tomo,(void*)array,amin_lens(nj),0.99999,NULL,1000); // zi<=zj
}

/*********** angular power spectra - with look-up tables ******************/

double C_cl_tomo(double l, int ni, int nj)  //galaxy clustering power spectrum of galaxies in bins ni, nj
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;

  static double **table;
  static double ds = .0, logsmin = .0, logsmax = .0;
  if (ni < 0 || ni >= tomo.clustering_Nbin ||nj < 0 || nj >= tomo.clustering_Nbin){
    printf("C_cl_tomo(l,%d,%d) outside tomo.clustering_Nbin range\nEXIT\n",ni,nj); exit(1);
  }

  if (recompute_clustering(C,G,N,ni,nj))
  {
    if (table==0) {
      table   = create_double_matrix(0, tomo.clustering_Nbin*tomo.clustering_Nbin-1, 0, Ntable.N_ell-1);
      logsmin = log(limits.P_2_s_min);
      logsmax = log(limits.P_2_s_max);
      ds = (logsmax - logsmin)/(Ntable.N_ell);
    }
    for (int i = 0; i < tomo.clustering_Nbin*tomo.clustering_Nbin; i++){table[i][0] = 123456789.0;}
    update_cosmopara(&C); update_nuisance(&N); update_galpara(&G);
  }
  int j = ni*tomo.clustering_Nbin+nj;
  if(table[j][0] > 123456780.0){ //still need to recompute this tomography bin combination
    //printf("recompute C_cl_tomo, %d, %d\n", ni, nj) ;
    double llog = logsmin;
    double result;
    for (int i=0; i<Ntable.N_ell; i++, llog+=ds) {

//      table[j][i]= log(C_cl_RSD_nointerp(exp(llog),ni,nj));
      result = C_cl_tomo_nointerp(exp(llog),ni,nj);
      if(result<=0) table[j][i] = -100;
      else table[j][i] = log(result);
      table[nj*tomo.clustering_Nbin+ni][i]=table[j][i];
    }
  }
  double f1 = exp(interpol_fitslope(table[j], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1.));
  if (isnan(f1)){f1 = 0.;}
  return f1;
}
