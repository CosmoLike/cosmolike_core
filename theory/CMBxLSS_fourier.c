//================================================================================================
// Naming convention:
// g = galaxy positions ("g" as in "galaxy")
// k = kappa CMB ("k" as in "kappa")
// s = kappa from source galaxies ("s" as in "shear")
// y = Compton-y field from CMB
// And alphabetical order

//================================================================================================

// integrands for power spectra
double W_k(double a, double fK); //lensing efficiency weightfunction for CMB lensing
double beam_SPT(double l);
double beam_cmb(double l, double fwhm_arcmin);
double int_for_C_gk(double a, void *params);
double int_for_C_ks(double a, void *params);
double int_for_C_kk(double a, void *params);

// power spectra - no look-up table
double C_gk_nointerp(double l, int ni);  //CMB kappa x galaxy position power spectrum, lens bin ni
double C_ks_nointerp(double l, int ns); // CMB kappa x galaxy kappa power spectrum, for source z-bin ns
double C_kk_nointerp(double l);
double C_gk(double l, int ni);//CMB kappa x galaxy positions in clustering bin ni
double C_ks(double l, int ni);
double C_kk(double l);

//======= y related power spectra
double W_y(double a, double fK); // efficiency weight function for Compton-y
double int_for_C_gy(double a, void *params);
double int_for_C_ky(double a, void *params);
double int_for_C_sy(double a, void *params);
double int_for_C_yy(double a, void *params);

// power spectra - no look-up table
double C_gy_nointerp(double l, int ni);  //CMB y x galaxy position power spectrum, lens bin ni
double C_sy_nointerp(double l, int ns); // CMB y x galaxy shear power spectrum, for source z-bin ns
double C_ky_nointerp(double l); // CMB y x CMB kappa power spectrum
double C_yy_nointerp(double l); // CMB y x y power spectrum
double C_gy(double l, int ni);
double C_sy(double l, int ns);
double C_ky(double l);
double C_yy(double l);

//================================================================================================
double beam_SPT(double l){
  double fwhm_arcmin =5.4;
  double sigma = fwhm_arcmin/sqrt(8.*log(2.0))*constants.arcmin;
  return exp(-0.5*l*l*sigma*sigma);
}
// integrands for power spectra

double beam_cmb(double l, double fwhm_arcmin){
  double sigma = fwhm_arcmin/sqrt(8.*log(2.0))*constants.arcmin;
  return exp(-0.5*l*l*sigma*sigma);
}

//lensing efficiency weightfunction for CMB lensing
double W_k(double a, double fK){
  double wk = 1.5*cosmology.Omega_m*fK/a*g_cmb(a);
  if(cosmology.MGSigma != 0.){
    wk *= (1.+MG_Sigma(a));
  }
  return wk;
}

// galaxy position x kappa CMB
double int_for_C_gk(double a, void *params)
{
  double *ar = (double *) params;
  double res,ell, fK, k;
  if (a >= 1.0) error("a >1 in int_for_C_gk");

  double ell_prefactor1 = (ar[1])*(ar[1]+1.);
  
  ell       = ar[1]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  res= W_gal(a,ar[0])*W_k(a,fK)*dchi_da(a)/fK/fK * ell_prefactor1/ell/ell ; //W_gal radial weight function for clustering, defined in cosmo2D_fourier.c
  res= res*Pdelta(k,a);
  return res;
}
double int_for_C_gk_b2(double a, void *params)
{
  double res,ell, fK, k;
  double *ar = (double *) params;
  double b1 = gbias.b1_function(1./a-1.,(int)ar[0]);
  double b2 = gbias.b2[(int)ar[0]];
  double bs2 = gbias.bs2[(int)ar[0]];
  double g4 = pow(growfac(a)/growfac(1.0),4.);
  
  ell       = ar[1]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  
  res= W_HOD(a,ar[0])*W_k(a,fK)*dchi_da(a)/fK/fK;
  res= res*(b1*Pdelta(k,a)+g4*(0.5*b2*PT_d1d2(k)+0.5*bs2*PT_d1s2(k)));
  return res;
}


// shear x kappa CMB
double int_for_C_ks(double a, void *params)
{
  double *ar = (double *) params;
  double res,ell, fK, k;
  if (a >= 1.0) error("a >1 in int_for_C_ks");

  double ell_prefactor1 = (ar[1])*(ar[1]+1.);
  double ell_prefactor2 = (ar[1]-1.)*ell_prefactor1*(ar[1]+2.);
  if(ell_prefactor2<=0.) 
    ell_prefactor2=0.;
  else
    ell_prefactor2=sqrt(ell_prefactor2);
  
  ell       = ar[1]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  res= W_kappa(a,fK,ar[0])*W_k(a,fK)*dchi_da(a)/fK/fK * ell_prefactor1*ell_prefactor2/ pow(ell,4);
  res= res*Pdelta(k,a);
  return res;
}

// kappa CMB x kappaCMB
double int_for_C_kk(double a, void *params){
   double *ar = (double *) params;
   double res,ell, fK, k;
   if (a > 1.0) error("a >1 in int_for_C_kk");

   double ell_prefactor1 = (ar[0])*(ar[0]+1.);
   
   ell       = ar[0]+0.5;
   fK     = f_K(chi(a));
   k      = ell/fK;
   res = pow(W_k(a,fK), 2)*dchi_da(a)/fK/fK * ell_prefactor1*ell_prefactor1/pow(ell,4);
   res = res*Pdelta(k,a);
   return res;
}


// IA for kappa CMB x shear: gI + true (fast)
double int_for_C_ks_IA(double a, void *params)
{
   double res, ell, fK, k,ws1,ws2,wk1,wk2, norm;
   double *ar = (double *) params;
   ell       = ar[1]+0.5;
   fK     = f_K(chi(a));
   k      = ell/fK;
   ws1 = W_source(a,ar[0]);
   wk1 = W_kappa(a,fK,ar[0]);
   wk2 = W_k(a,fK);
   norm = A_IA_Joachimi(a)*cosmology.Omega_m*nuisance.c1rhocrit_ia*growfac(1.)/growfac(a);
   res= -ws1*wk2*norm + wk1*wk2;
   return res*Pdelta(k,a)*dchi_da(a)/fK/fK;
}

double int_for_C_ks_IA_Az(double a, void *params)
{
  double res, ell, fK, k,ws1,ws2,wk1,wk2, norm;
  double *ar = (double *) params;
  ell       = ar[1]+0.5;
  fK     = f_K(chi(a));
  k      = ell/fK;
  ws1 = W_source(a,ar[0])*nuisance.A_z[(int)ar[0]];
  wk1 = W_kappa(a,fK,ar[0]);
  wk2 = W_k(a,fK);
  norm = cosmology.Omega_m*nuisance.c1rhocrit_ia*growfac(1.)/growfac(a);
  res= -ws1*wk2*norm + wk1*wk2;
  
  return res*Pdelta(k,a)*dchi_da(a)/fK/fK;
}

double int_for_C_ks_IA_mpp(double a, void *params)
{ // for like.IA==4
   double res, ell, fK, k,ws1,ws2,wk1,wk2, norm;
   double *ar = (double *) params;
   ell       = ar[1]+0.5;
   fK     = f_K(chi(a));
   k      = ell/fK;
   ws1 = W_source(a,ar[0]);
   wk1 = W_kappa(a,fK,ar[0]);
   wk2 = W_k(a,fK);
   norm = cosmology.Omega_m*nuisance.c1rhocrit_ia*growfac(1.)/growfac(a)*nuisance.A_ia*pow(1./(a*nuisance.oneplusz0_ia),nuisance.eta_ia);
   res= -ws1*wk2*norm + wk1*wk2;
   return res*Pdelta(k,a)*dchi_da(a)/fK/fK;
}

double C_ks_IA(double s, int ni)
{
   double array[2] = {(double) ni,s};
   if (like.IA==1) return int_gsl_integrate_medium_precision(int_for_C_ks_IA,(void*)array,amin_source(ni),amax_source(ni),NULL,1000);
   if (like.IA==3) return int_gsl_integrate_medium_precision(int_for_C_ks_IA_Az,(void*)array,amin_source(ni),amax_source(ni),NULL,1000);
   if (like.IA==4) return int_gsl_integrate_medium_precision(int_for_C_ks_IA_mpp,(void*)array,amin_source(ni),0.99999,NULL,1000);
   printf("CMBxLSS.c: C_ks_IA does not support like.IA = %d\nEXIT\n", like.IA);
  exit(1);
}

//================================================================================================
// power spectra - no look-up tables

// galaxy position x kappa CMB, lens z-bin nl
double C_gk_nointerp(double l, int nl)
{
   double array[2] = {(double)nl,l};
   // MANUWARNING
//   printf("amin_lens, amax_lens, %e %e\n", amin_lens(nl), amax_lens(nl));
   if (gbias.b2[nl] || gbias.b2[nl]) return int_gsl_integrate_medium_precision(int_for_C_gk_b2 ,(void*)array,amin_lens(nl),amax_lens(nl),NULL,1000);

   // return int_gsl_integrate_medium_precision(int_for_C_gk,(void*)array,amin_lens(nl),amax_lens(nl),NULL,1000);
   return int_gsl_integrate_medium_precision(int_for_C_gk,(void*)array,amin_lens(nl),0.99999,NULL,1000);
}

// shear x kappa CMB, for source z-bin ns
double C_ks_nointerp(double l, int ns) {
  if (like.IA) return C_ks_IA(l,ns);
  double array[2] = {(double) ns, l};
  // return int_gsl_integrate_medium_precision(int_for_C_ks,(void*)array,amin_source(ns),amax_source(ns),NULL,1000);
  return int_gsl_integrate_medium_precision(int_for_C_ks,(void*)array,amin_source(ns),0.99999,NULL,1000);
}

// kappa CMB x kappa CMB
double C_kk_nointerp(double l){
   double array[1] = {l};
   return int_gsl_integrate_medium_precision(int_for_C_kk, (void*)array, limits.a_min*(1.+1.e-5), 1.-1.e-5, NULL, 1000);
}



//================================================================================================
// power spectra - look-up tables


// galaxy position x kappa CMB, lens bin ni
double C_gk(double l, int ni)
{
   static cosmopara C;
   static nuisancepara N;
   static galpara G;
   
   static double **table;
   static double ds = .0, logsmin = .0, logsmax = .0;
   if (ni < 0 || ni >= tomo.clustering_Nbin){
      printf("Bin %d outside tomo.clustering_Nbin range\nEXIT\n",ni); exit(1);
   }
   
   if (recompute_clustering(C,G,N,ni,ni))
   {
      if (table==0) {
         table   = create_double_matrix(0, tomo.clustering_Nbin-1, 0, Ntable.N_ell-1);
         logsmin = log(limits.P_2_s_min);
         logsmax = log(limits.P_2_s_max);
         ds = (logsmax - logsmin)/(Ntable.N_ell);
      }
      
      double llog;
      int i,l1,k;
      
      for (k=0; k<tomo.clustering_Nbin; k++) {
            llog = logsmin;
            for (i=0; i<Ntable.N_ell; i++, llog+=ds) {
               table[k][i]= log(C_gk_nointerp(exp(llog),k));
            }
         }
      update_cosmopara(&C); update_nuisance(&N); update_galpara(&G);
   }

   double f1 = exp(interpol(table[ni], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1., 1.));
   if (isnan(f1)){f1 = 0.;}
   return f1;
}


// shear x kappa CMB, source bin ni
double C_ks(double l, int ni)
{
   static cosmopara C;
   static nuisancepara N;
   
   static double **table, *sig;
   static int osc[100];
   static double ds = .0, logsmin = .0, logsmax = .0;
   if (ni < 0 || ni >= tomo.clustering_Nbin){
      printf("Bin %d outside tomo.clustering_Nbin range\nEXIT\n",ni); exit(1);
   }

   if (recompute_shear(C,N))
   {
      if (table==0) {
         table   = create_double_matrix(0, tomo.shear_Nbin-1, 0, Ntable.N_ell-1);
         sig = create_double_vector(0,tomo.shear_Nbin-1);
         logsmin = log(limits.P_2_s_min);
         logsmax = log(limits.P_2_s_max);
         ds = (logsmax - logsmin)/(Ntable.N_ell);
      }
      
      double res,llog;
      int i,l1,k;
      
      for (k=0; k<tomo.shear_Nbin; k++) {
        llog = logsmin;

        sig[k] = 1.;
        osc[k] = 0;
        res = C_ks_nointerp(500.,k);
        if (res < 0){sig[k] = -1.;}
        for (i=0; i<Ntable.N_ell; i++, llog+=ds) {
          table[k][i]= C_ks_nointerp(exp(llog),k);
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
      update_cosmopara(&C); update_nuisance(&N);
   }

   double f1 = 0.;
   if(osc[ni] ==0 ){f1 = sig[ni]*exp(interpol_fitslope(table[ni], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1.));}
   if(osc[ni] ==1 ){f1 = interpol_fitslope(table[ni], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1.);}

   // if(osc[k] ==0) {f1 = exp(interpol(table[ni], Ntable.N_ell, logsmin, logsmax, ds, log(l), 1., 1.));}
   if (isnan(f1)){f1 = 0.;}
   return f1;
}


// kappa CMB x kappa CMB
double C_kk(double l)
{
   static cosmopara C;
   
   static double *table;
   static double ds = .0, logsmin = .0, logsmax = .0;
   
   if (recompute_cosmo3D(C))
   {
      if (table==0) {
         table   = create_double_vector(0, Ntable.N_ell-1);
         logsmin = log(limits.P_2_s_min);
         logsmax = log(limits.P_2_s_max);
         ds = (logsmax - logsmin)/(Ntable.N_ell);
      }
      
      double llog;
      int i;
      
      llog = logsmin;
      for (i=0; i<Ntable.N_ell; i++, llog+=ds) {
         table[i]= log(C_kk_nointerp(exp(llog)));
      }

      update_cosmopara(&C); 
   }
   
   double f1 = exp(interpol(table, Ntable.N_ell, logsmin, logsmax, ds, log(l), 1., 1.));
   if (isnan(f1)){f1 = 0.;}
   return f1;
}



//================================================================================================
// y related functions

// efficiency weight function for Compton-y
double W_y(double a, double fK){
  return 0;
}
double int_for_C_gy(double a, void *params){
  return 0;
}
double int_for_C_ky(double a, void *params){
  return 0;
}
double int_for_C_sy(double a, void *params){
  return 0;
}
double int_for_C_yy(double a, void *params){
  return 0;
}

// power spectra - no look-up table
//CMB y x galaxy position power spectrum, lens bin ni
double C_gy_nointerp(double l, int ni){
  return 0;
}

// CMB y x galaxy shear, for source z-bin ns
double C_sy_nointerp(double l, int ns);
{
  return 0;
}

// CMB y x CMB kappa
double C_ky_nointerp(double l);
{
  return 0;
}

double C_yy_nointerp(double l);
{
  return 0;
}

double C_gy(double l, int ni){
  return 0;
}

double C_sy(double l, int ns){
  return 0;
}

double C_ky(double l){
  return 0;
}

double C_yy(double l){
  return 0;
}
