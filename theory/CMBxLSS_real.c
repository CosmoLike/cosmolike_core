// correlation functions
double w_gk_SPT(double theta,int ni); //angular CMB lensing x positions correlation function in tomography bin ni
double w_ks_SPT(double theta, int ni);//angular CMB lensing x galaxy shear correlation function in tomography bin ni


//================================================================================================
// angular correlation functions - no look-up tables
double C_gk_wrapper(double l,int ni, int nj){
  return C_gk(l,ni)*beam_SPT(l);
}

double C_ks_wrapper(double l,int ni, int nj){
  return C_ks(l,ni)*beam_SPT(l);
}

double w_gk_SPT(double theta, int ni) // galaxies (ni) x CMB kappa
{
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;
  if (recompute_clustering(C,G,N,ni,ni)){
    double **tab;
    int i, j,k;
    tab   = create_double_matrix(0, 1, 0, Ntable.N_thetaH-1);
    if (table==0) table   = create_double_matrix(0, tomo.clustering_Nbin-1, 0, Ntable.N_thetaH-1);    
    for (i = 0; i < tomo.clustering_Nbin; i++){
        twopoint_via_hankel(tab, &logthetamin, &logthetamax,&C_gk_wrapper, i,i,0);
        dlogtheta = (logthetamax-logthetamin)/((double)Ntable.N_thetaH);
        for (k = 0; k < Ntable.N_thetaH; k++){
          table[i][k] = tab[0][k];
      }
    }
    free_double_matrix(tab,0, 1, 0, Ntable.N_thetaH-1);   
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }
  return interpol(table[ni], Ntable.N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 1.0, 1.0);
}
double w_ks_SPT(double theta, int ni){ // CMB kappa x shear (ni)
  static cosmopara C;
  static nuisancepara N;
  
  static double **table;
  static double dlogtheta, logthetamin, logthetamax;
  if (recompute_shear(C,N)){
    double **tab;
    int i, j,k;
    tab   = create_double_matrix(0, 1, 0, Ntable.N_thetaH-1);
    if (table==0) table   = create_double_matrix(0, tomo.shear_Nbin-1, 0, Ntable.N_thetaH-1);    
    for (i = 0; i < tomo.shear_Nbin; i++){
        twopoint_via_hankel(tab, &logthetamin, &logthetamax,&C_ks_wrapper, i,i,2);
        dlogtheta = (logthetamax-logthetamin)/((double)Ntable.N_thetaH);
        for (k = 0; k < Ntable.N_thetaH; k++){
          table[i][k] = tab[0][k];
        }
    }
    free_double_matrix(tab,0, 1, 0, Ntable.N_thetaH-1);   
    update_cosmopara(&C);
    update_nuisance(&N);
  }
  return interpol(table[ni], Ntable.N_thetaH, logthetamin, logthetamax,dlogtheta, log(theta), 1.0, 1.0);
}