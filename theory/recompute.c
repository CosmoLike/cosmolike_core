void update_cosmopara (cosmopara *C);
void update_cosmopara_lowz (cosmopara_lowz *C_low);
void update_gal (galpara *G);
void update_nuisance (nuisancepara *N);
//conditions for recomputing look-up tables
int recompute_expansion(cosmopara C);
int recompute_Delta(cosmopara C);
int recompute_cosmo3D(cosmopara C);
int recompute_cosmo3D_lowz(cosmopara_lowz C_low);
int recompute_cosmo3D_CLASS(cosmopara C);
int recompute_zphot_shear(nuisancepara N);
int recompute_zphot_clustering(nuisancepara N);
int recompute_zphot_outlierfrac(nuisancepara N);
int recompute_zphot_magnification(nuisancepara N);
int recompute_shear(cosmopara C, nuisancepara N); //for shear 2-pt statics
int recompute_ii(cosmopara C, nuisancepara N); //for shear 2-pt statics
int recompute_ggl(cosmopara C, galpara G, nuisancepara N,int i);//for gg-lensing statistics
int recompute_clustering(cosmopara C, galpara G, nuisancepara N, int i, int j);//clustering
int recompute_clusters(cosmopara C, nuisancepara N); //recompute criteria
int recompute_PkRatio(barypara B);
void update_PkRatio(barypara *B);
int recompute_DESclusters(cosmopara C, nuisancepara N); //recompute criteria

int recompute_gk(cosmopara C, galpara G, nuisancepara N,int i);//for gk statistics
int recompute_ks(cosmopara C, nuisancepara N);//ks
int recompute_kk(cosmopara C, nuisancepara N);

int recompute_gy(cosmopara C, galpara G, nuisancepara N, int i);
int recompute_sy(cosmopara C, nuisancepara N);
int recompute_ky(cosmopara C, nuisancepara N);
int recompute_yy(cosmopara C, nuisancepara N);

int recompute_halo(nuisancepara N);

void update_cosmopara (cosmopara *C){
  C->Omega_m = cosmology.Omega_m;
  C->Omega_v = cosmology.Omega_v;
  C->Omega_nu = cosmology.Omega_nu;
  C->M_nu = cosmology.M_nu;
  C->sigma_8 = cosmology.sigma_8;
  C->A_s = cosmology.A_s;
  C->n_spec = cosmology.n_spec;
  C->alpha_s = cosmology.alpha_s;
  C->w0 = cosmology.w0;
  C->wa = cosmology.wa;
  C->omb = cosmology.omb;
  C->h0 = cosmology.h0;
  C->f_NL = cosmology.f_NL;
  C->MGSigma = cosmology.MGSigma;
  C->MGmu = cosmology.MGmu;
  C->M_nu = cosmology.M_nu;
  C->theta_s = cosmology.theta_s;
}
void update_cosmopara_lowz (cosmopara_lowz *C_low){
  C_low->sigma_8 = cosmology_lowz.sigma_8;
  C_low->z_low = cosmology_lowz.z_low;
}

void update_galpara (galpara *G){
  int i,j;
  for (i = 0; i < tomo.clustering_Nbin; i++){
    if ((gbias.hod[i][0] > 10 && gbias.hod[i][0]<16) || (gbias.b[i]> 0.2 && gbias.b[i] < 20)){
      G->b[i] = gbias.b[i];
      G->b2[i] = gbias.b2[i];
      G->bs2[i] = gbias.bs2[i];
      for(j = 0; j < 6; j++){
        G->hod[i][j] = gbias.hod[i][j];
      }
      G->cg[i]= gbias.cg[i];
    }
    else{ printf("lens bin %d: neither HOD nor linear bias set, exit\n",i); exit(EXIT_FAILURE);}
  }
}

void update_nuisance (nuisancepara *N){
  int i;
  N->A_ia = nuisance.A_ia;
  N->beta_ia = nuisance.beta_ia;
  N->eta_ia = nuisance.eta_ia;
  N->eta_ia_highz = nuisance.eta_ia_highz;
  N->A2_ia = nuisance.A2_ia;
  N->eta_ia_tt = nuisance.eta_ia_tt;
  N->LF_alpha = nuisance.LF_alpha;
  N->LF_P = nuisance.LF_P;
  N->LF_Q = nuisance.LF_Q;
  N->LF_red_alpha = nuisance.LF_red_alpha;
  N->LF_red_P = nuisance.LF_red_P;
  N->LF_red_Q = nuisance.LF_red_Q;
  for(i = 0; i < tomo.clustering_Nbin; i++){
    N-> fred[i] = nuisance.fred[i];
    N->sigma_zphot_clustering[i] = nuisance.sigma_zphot_clustering[i];
    N->bias_zphot_clustering[i] = nuisance.bias_zphot_clustering[i];
  }
  for(i = 0; i < tomo.shear_Nbin; i++){
    N->sigma_zphot_shear[i] = nuisance.sigma_zphot_shear[i];
    N->bias_zphot_shear[i] = nuisance.bias_zphot_shear[i];
    N->A_z[i] = nuisance.A_z[i];
    N->A2_z[i] = nuisance.A2_z[i];
    N->b_ta_z[i] = nuisance.b_ta_z[i];
  }
  for(i = 0; i < tomo.magnification_Nbin; i++){
    N->sigma_zphot_magnification[i] = nuisance.sigma_zphot_magnification[i];
    N->bias_zphot_magnification[i] = nuisance.bias_zphot_magnification[i];
  }
  N->cluster_Mobs_lgM0 = nuisance.cluster_Mobs_lgM0;
  N->cluster_Mobs_sigma = nuisance.cluster_Mobs_sigma;
  N->cluster_Mobs_alpha = nuisance.cluster_Mobs_alpha;
  N->cluster_Mobs_beta = nuisance.cluster_Mobs_beta;
  N->cluster_Mobs_lgN0 = nuisance.cluster_Mobs_lgN0;
  N->cluster_Mobs_sigma0 = nuisance.cluster_Mobs_sigma0;
  N->cluster_Mobs_sigma_qm = nuisance.cluster_Mobs_sigma_qm;
  N->cluster_Mobs_sigma_qz = nuisance.cluster_Mobs_sigma_qz;

  for(i = 0; i < tomo.cluster_Nbin; i++){
    N->cluster_completeness[i] = nuisance.cluster_completeness[i];
  }
  N->cluster_centering_f0 = nuisance.cluster_centering_f0;
  N->cluster_centering_alpha = nuisance.cluster_centering_alpha;
  N->cluster_centering_sigma = nuisance.cluster_centering_sigma;
  for (int _i = 0; _i < nuisance.N_cluster_MOR; ++_i){
    N->cluster_MOR[_i] = nuisance.cluster_MOR[_i];
  }
  for (int _i = 0; _i < nuisance.N_cluster_selection; ++_i){
    N->cluster_selection[_i] = nuisance.cluster_selection[_i];
  }

  N->frac_lowz = nuisance.frac_lowz;
  N->frac_highz= nuisance.frac_highz;

  if (strcmp(pdeltaparams.runmode,"halomodel") ==0){
    N->gas_beta = nuisance.gas_beta;
    N->gas_lgM0 = nuisance.gas_lgM0; 
    N->gas_eps1 = nuisance.gas_eps1;
    N->gas_eps2 = nuisance.gas_eps2;
    N->gas_beta_v2 = nuisance.gas_beta_v2;
    N->gas_lgM0_v2 = nuisance.gas_lgM0_v2;
    N->gas_eps1_v2 = nuisance.gas_eps1_v2;
    N->gas_eps2_v2 = nuisance.gas_eps2_v2;
    N->gas_alpha = nuisance.gas_alpha; 
    N->gas_A_star = nuisance.gas_A_star; 
    N->gas_lgM_star = nuisance.gas_lgM_star; 
    N->gas_sigma_star = nuisance.gas_sigma_star;
    N->gas_lgT_w = nuisance.gas_lgT_w;
    N->gas_f_H = nuisance.gas_f_H;
    N->gas_Gamma_KS = nuisance.gas_Gamma_KS;
  }
}
int recompute_expansion(cosmopara C){ //rules for recomputing growth factor & comoving distance
  if (C.Omega_m != cosmology.Omega_m || C.Omega_v != cosmology.Omega_v || C.w0 != cosmology.w0 || C.wa != cosmology.wa || C.MGmu != cosmology.MGmu || C.M_nu != cosmology.M_nu){return 1;}
  if (cosmology.theta_s > 0 && C.theta_s != cosmology.theta_s){return 1;}
  else{return 0;}
}

int recompute_Delta(cosmopara C){ //rules for recomputing early time power spectrum Delta_L
  if (C.Omega_m != cosmology.Omega_m || C.Omega_v != cosmology.Omega_v || C.Omega_nu != cosmology.Omega_nu || C.M_nu != cosmology.M_nu || C.h0 != cosmology.h0 || C.omb != cosmology.omb || C.n_spec != cosmology.n_spec|| C.alpha_s != cosmology.alpha_s){return 1;}
  if (cosmology.A_s){
    if(C.A_s != cosmology.A_s){return 1;}
  }
  else{
    if (C.sigma_8 != cosmology.sigma_8){return 1;}
  }
  return 0;
}

int recompute_cosmo3D(cosmopara C){
  if (C.Omega_m != cosmology.Omega_m || C.Omega_v != cosmology.Omega_v || C.Omega_nu != cosmology.Omega_nu || C.M_nu != cosmology.M_nu || C.h0 != cosmology.h0 || C.omb != cosmology.omb || C.n_spec != cosmology.n_spec|| C.alpha_s != cosmology.alpha_s ||  C.w0 != cosmology.w0 || C.wa != cosmology.wa || C.MGSigma != cosmology.MGSigma || C.MGmu != cosmology.MGmu || C.M_nu != cosmology.M_nu){return 1;}
  if (cosmology.A_s){
     if(C.A_s != cosmology.A_s){return 1;}
  }
  else{
     if (C.sigma_8 != cosmology.sigma_8){return 1;}
  }
  if (cosmology.theta_s > 0 && C.theta_s != cosmology.theta_s){return 1;}
  return 0;
}
int recompute_cosmo3D_lowz(cosmopara_lowz C_low){
  if (C_low.sigma_8 != cosmology_lowz.sigma_8 || C_low.z_low != cosmology_lowz.z_low){return 1;}
  return 0;
}
int recompute_cosmo3D_CLASS(cosmopara C){
  if (C.Omega_m != cosmology.Omega_m || C.Omega_v != cosmology.Omega_v || C.Omega_nu != cosmology.Omega_nu || C.M_nu != cosmology.M_nu || C.h0 != cosmology.h0 || C.omb != cosmology.omb || C.n_spec != cosmology.n_spec|| C.alpha_s != cosmology.alpha_s ||  C.w0 != cosmology.w0 || C.wa != cosmology.wa || C.MGSigma != cosmology.MGSigma || C.MGmu != cosmology.MGmu || C.M_nu != cosmology.M_nu){return 1;}
  if (cosmology.A_s > 0){
     if(C.A_s != cosmology.A_s){return 1;}
  }
  else{
     if (C.sigma_8 != cosmology.sigma_8){return 1;}
  }
  if (cosmology.theta_s > 0 && C.theta_s != cosmology.theta_s){return 1;}
  return 0;
}

int recompute_zphot_shear(nuisancepara N){
  if(recompute_zphot_outlierfrac(N)) {return 1;}
  static int photoz = -1;
  if (photoz != redshift.shear_photoz){photoz = redshift.shear_photoz; return 1;}
  if (redshift.shear_photoz != 3 && redshift.shear_photoz != 4){return 0;}
  int i, res = 0;
  for(i = 0; i < tomo.shear_Nbin; i++){
    if (N.sigma_zphot_shear[i]!= nuisance.sigma_zphot_shear[i] || N.bias_zphot_shear[i]!= nuisance.bias_zphot_shear[i]){ res = 1;}
  }
  return res;
}
int recompute_zphot_clustering(nuisancepara N){
  if(recompute_zphot_outlierfrac(N)) {return 1;}
  static int photoz = -1;
  if (photoz != redshift.clustering_photoz){photoz = redshift.clustering_photoz; return 1;}
  if (redshift.clustering_photoz != 3 && redshift.clustering_photoz != 4){return 0;}
  int i, res = 0;
  for(i = 0; i < tomo.clustering_Nbin; i++){
    if (N.sigma_zphot_clustering[i]!= nuisance.sigma_zphot_clustering[i] || N.bias_zphot_clustering[i]!= nuisance.bias_zphot_clustering[i]){ res = 1;}
  }
  return res;
}
int recompute_zphot_outlierfrac(nuisancepara N){
  int i, res = 0;
    if (N.frac_lowz!= nuisance.frac_lowz || N.frac_highz!= nuisance.frac_highz){ return 1;}
  return 0;
}
int recompute_zphot_magnification(nuisancepara N){
  if (redshift.magnification_photoz != 3){return 0;}
  int i, res = 0;
  for(i = 0; i < tomo.magnification_Nbin; i++){
    if (N.sigma_zphot_magnification[i]!= nuisance.sigma_zphot_magnification[i] || N.bias_zphot_magnification[i]!= nuisance.bias_zphot_magnification[i]){ res = 1;}
  }
  return res;
}
int recompute_IA(nuisancepara N){
  if (N.A_ia != nuisance.A_ia || N.eta_ia != nuisance.eta_ia || N.A2_ia != nuisance.A2_ia || N.eta_ia_tt != nuisance.eta_ia_tt) return 1;
  for(int i = 0; i < tomo.shear_Nbin; i++){
    if (N.A_z[i]!= nuisance.A_z[i]){ return 1;}
    if (N.A2_z[i]!= nuisance.A2_z[i]){ return 1;}
    if (N.b_ta_z[i]!= nuisance.b_ta_z[i]){ return 1;}
  }
  return 0;
}

int recompute_shear(cosmopara C, nuisancepara N){
  if (recompute_cosmo3D(C) || recompute_zphot_shear(N) || recompute_zphot_outlierfrac(N)||recompute_IA(N) || recompute_halo(N)){return 1;}
  else{return 0;}
}
int recompute_clusters(cosmopara C, nuisancepara N){
  if (recompute_cosmo3D(C) || N.cluster_Mobs_lgM0 != nuisance.cluster_Mobs_lgM0 || N.cluster_Mobs_sigma != nuisance.cluster_Mobs_sigma || N.cluster_Mobs_alpha != nuisance.cluster_Mobs_alpha || N.cluster_Mobs_beta != nuisance.cluster_Mobs_beta ||   N.cluster_Mobs_lgN0 != nuisance.cluster_Mobs_lgN0 ||  N.cluster_Mobs_sigma0 != nuisance.cluster_Mobs_sigma0 ||  N.cluster_Mobs_sigma_qm != nuisance.cluster_Mobs_sigma_qm ||  N.cluster_Mobs_sigma_qz != nuisance.cluster_Mobs_sigma_qz || N.cluster_completeness[0] != nuisance.cluster_completeness[0] || N.cluster_centering_f0 != nuisance.cluster_centering_f0 || N.cluster_centering_alpha != nuisance.cluster_centering_alpha || N.cluster_centering_sigma != nuisance.cluster_centering_sigma){return 1;}
  else{return 0;}
}

int recompute_DESclusters(cosmopara C, nuisancepara N){
   if (recompute_cosmo3D(C)) return 1;
   for (int _i = 0; _i < nuisance.N_cluster_MOR; ++_i){
      if (N.cluster_MOR[_i] !=nuisance.cluster_MOR[_i] ) return 1; 
   }
   for (int _i = 0; _i < nuisance.N_cluster_selection; ++_i){
      if (N.cluster_selection[_i] !=nuisance.cluster_selection[_i] ) return 1; 
   }
   return 0;
}

int recompute_ii(cosmopara C, nuisancepara N){
  if (recompute_cosmo3D(C) || recompute_halo(N) || recompute_zphot_clustering(N)|| recompute_zphot_shear(N) || recompute_zphot_outlierfrac(N)|| N.A_ia != nuisance.A_ia || N.beta_ia != nuisance.beta_ia ||N.eta_ia != nuisance.eta_ia || N.eta_ia_highz != nuisance.eta_ia_highz || N.LF_alpha != nuisance.LF_alpha || N.LF_red_alpha != nuisance.LF_red_alpha || N.LF_P != nuisance.LF_P || N.LF_Q != nuisance.LF_Q || N.LF_red_P != nuisance.LF_red_P || N.LF_red_Q != nuisance.LF_red_Q){return 1;}
  else{return 0;} 
}

int recompute_galaxies(galpara G, int i){
 int j;
  if (i == -1){return 0;}
  if(G.hod[i][0] > 10 && G.hod[i][0] < 16){
    for(j = 0; j <6; j++){
       if (G.hod[i][j] != gbias.hod[i][j]){return 1;}
    }
  }
  if(G.b[i] != gbias.b[i] || G.b2[i] != gbias.b2[i] || G.bs2[i] != gbias.bs2[i] ||G.cg[i]!= gbias.cg[i]){return 1;}
  return 0;
}

int recompute_ggl(cosmopara C, galpara G, nuisancepara N, int i){
  if (recompute_cosmo3D(C) || recompute_zphot_clustering(N) || recompute_zphot_outlierfrac(N) || recompute_zphot_shear(N) || recompute_galaxies(G,i) ||recompute_IA(N) || recompute_halo(N) ){return 1;}
  else{return 0;}
}

int recompute_clustering(cosmopara C, galpara G, nuisancepara N, int i, int j){
  if (recompute_cosmo3D(C) || recompute_zphot_clustering(N) || recompute_zphot_outlierfrac(N) || recompute_galaxies(G,i)|| recompute_galaxies(G,j) || recompute_halo(N)){return 1;}
  else{return 0;}
 
}

int recompute_gk(cosmopara C, galpara G, nuisancepara N, int i){
  if (recompute_cosmo3D(C) || recompute_zphot_clustering(N) || recompute_zphot_outlierfrac(N) || recompute_galaxies(G,i) || recompute_halo(N) ){return 1;}
  else{return 0;}
}

int recompute_ks(cosmopara C, nuisancepara N){
  if (recompute_cosmo3D(C) || recompute_zphot_shear(N) || recompute_zphot_outlierfrac(N) ||recompute_IA(N) || recompute_halo(N) ){return 1;}
  else{return 0;}
}

int recompute_kk(cosmopara C, nuisancepara N){
  if (recompute_cosmo3D(C) || recompute_halo(N) ){return 1;}
  else{return 0;}
}

int recompute_PkRatio(barypara B){
	if (strcmp(B.scenario,bary.scenario)!=0){return 1;}
	return 0;
}

void update_PkRatio(barypara *B){
	sprintf((*B).scenario,"%s", bary.scenario);
}

int recompute_gy(cosmopara C, galpara G, nuisancepara N, int i){
  if (recompute_cosmo3D(C) || recompute_zphot_clustering(N) || recompute_zphot_outlierfrac(N) || recompute_galaxies(G,i) || recompute_halo(N) ){return 1;}
  else{return 0;}
}

int recompute_sy(cosmopara C, nuisancepara N){
  if (recompute_cosmo3D(C) || recompute_zphot_shear(N) || recompute_zphot_outlierfrac(N) ||recompute_IA(N) || recompute_halo(N) ){return 1;}
  else{return 0;}
}

int recompute_ky(cosmopara C, nuisancepara N){
  if (recompute_cosmo3D(C) || recompute_halo(N) ){return 1;}
  else{return 0;}
}

int recompute_yy(cosmopara C, nuisancepara N){
  if (recompute_cosmo3D(C) || recompute_halo(N) ){return 1;}
  else{return 0;}
}

int recompute_halo(nuisancepara N){
  if (strcmp(pdeltaparams.runmode,"halomodel") ==0){
    if(N.gas_beta != nuisance.gas_beta || N.gas_lgM0 != nuisance.gas_lgM0 || N.gas_eps1 != nuisance.gas_eps1 || N.gas_eps2 != nuisance.gas_eps2 || \
       N.gas_beta_v2 != nuisance.gas_beta_v2 || N.gas_lgM0_v2 != nuisance.gas_lgM0_v2 || N.gas_eps1_v2 != nuisance.gas_eps1_v2 || N.gas_eps2_v2 != nuisance.gas_eps2_v2 || \
       N.gas_alpha != nuisance.gas_alpha || N.gas_A_star != nuisance.gas_A_star || N.gas_lgM_star != nuisance.gas_lgM_star || N.gas_sigma_star != nuisance.gas_sigma_star || \
       N.gas_lgT_w != nuisance.gas_lgT_w || N.gas_f_H != nuisance.gas_f_H || N.gas_Gamma_KS != nuisance.gas_Gamma_KS ) {return 1;}
    else{return 0;}
  }else{
    return 0;
  }
}
