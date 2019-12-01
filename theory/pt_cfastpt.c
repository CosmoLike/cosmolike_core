
double PT_d1d2(double k_coverH0);
double PT_d1d3(double k_coverH0);
double PT_d2d2(double k_coverH0);
double PT_d1s2(double k_coverH0);
double PT_d2s2(double k_coverH0);
double PT_s2s2(double k_coverH0);
double PT_sigma4(double k_coverH0);

double K_CH0(double k_mpch){
  return k_mpch*cosmology.coverH0;
}
void FPT_input(double k[FPT.N], double P[FPT.N]){
  int i;
  double dlgk,lgk;
  if ((int)log10(FPT.k_max/FPT.k_min)*FPT.N_per_dec != FPT.N){
    printf("pt_cfastpt.c:FPT_input: inconsistent k-range and number of bins for FPT\n");
    printf("FPT.k_min=%e, FPT.k_max=%e, FPT.N_per_dec=%d; FPT.N=%d\nEXIT\n",FPT.k_min,FPT.k_max,FPT.N_per_dec,FPT.N);
    exit(1);
  }
  if (FPT.k_min < limits.k_min_cH0 || FPT.k_max > limits.k_max_cH0){
    printf("pt_cfastpt.c:FPT_input: k_min/k_max out of range\n");
    exit(1);
  }
  dlgk = log(10.)/(double) FPT.N_per_dec;
  for (lgk = log(FPT.k_min), i = 0; i< FPT.N; i++,lgk+=dlgk){
    k[i] = exp(lgk);
    P[i] = p_lin(exp(lgk),1.0);
    //printf("%e %e\n",k[i],P[i]);
 } 
}
void get_FPT_bias(void){
  static cosmopara C; 
  if (recompute_cosmo3D(C)){
    update_cosmopara(&C); 
    if (FPT.tab_AB ==0){ //if table doesn't exit yet, create it
      FPT.tab_AB = create_double_matrix(0, FPT.N_AB-1, 0, FPT.N-1);
      FPT.k_min = K_CH0(FPT.k_min);
      FPT.k_max = K_CH0(FPT.k_max);
    }

    double k[FPT.N], Pin[FPT.N], Pout[FPT.N];
    FPT_input(k, Pin); 
    int i;
    Pd1d2(k, Pin, FPT.N, Pout);
    for (i =0; i < FPT.N; i++){
      FPT.tab_AB[0][i] = Pout[i]; //Pd1d2
    }
    Pd2d2(k, Pin, FPT.N, Pout);
    for (i =0; i < FPT.N; i++){
      FPT.tab_AB[1][i] = Pout[i]; //Pd2d2
    }
    Pd1s2(k, Pin, FPT.N, Pout);
    for (i =0; i < FPT.N; i++){
      FPT.tab_AB[2][i] = Pout[i]; //Pd1s2
    }
    Pd2s2(k, Pin, FPT.N, Pout);
    for (i =0; i < FPT.N; i++){
      FPT.tab_AB[3][i] = Pout[i]; //Pd2s2
    }
    Ps2s2(k, Pin, FPT.N, Pout);
    for (i =0; i < FPT.N; i++){
      FPT.tab_AB[4][i] = Pout[i]; //Pd2s2
    }
  }
}


double PT_d1d2(double k_coverH0){ //interpolate FPT.tab_AB[0] - Pd1d2
  static cosmopara C;
  static double logkmin = 0.,logkmax = 0.,dlgk = 0.; 
  // only call FASTPT if cosmology changed since last call
  if (recompute_cosmo3D(C)){
    get_FPT_bias();
    update_cosmopara(&C); 
    logkmin =log(FPT.k_min);
    logkmax = log(FPT.k_max);
    dlgk = log(10.)/(double) FPT.N_per_dec;
  }
  double lgk = log(k_coverH0);
  if (lgk < logkmin || lgk >= logkmax){return 0.;}
  return interpol(FPT.tab_AB[0], FPT.N, logkmin, logkmax, dlgk,lgk, 0.,0.);
}

double PT_d1d3(double k_coverH0){ //interpolate FPT.tab_AB[0] - Pd1d3
  static cosmopara C;
  static double logkmin = 0.,logkmax = 0.,dlgk = 0.; 
  // only call FASTPT if cosmology changed since last call
  if (recompute_cosmo3D(C)){
    get_FPT_bias();
    update_cosmopara(&C); 
    logkmin =log(FPT.k_min);
    logkmax = log(FPT.k_max);
    dlgk = log(10.)/(double) FPT.N_per_dec;
  }
  double lgk = log(k_coverH0);
  if (lgk < logkmin || lgk >= logkmax){return 0.;}
  return 0.0;//interpol(FPT.tab_AB[6], FPT.N, logkmin, logkmax, dlgk,lgk, 0.,0.);
}

double PT_d2d2(double k_coverH0){ //interpolate FPT.tab_AB[1]
  static cosmopara C; 
  static double logkmin = 0.,logkmax = 0.,dlgk = 0.;
  if (recompute_cosmo3D(C)){
    get_FPT_bias();
    update_cosmopara(&C); 
    logkmin =log(FPT.k_min);
    logkmax = log(FPT.k_max);
    dlgk = log(10.)/(double) FPT.N_per_dec;
  }
  double lgk = log(k_coverH0);
  return interpol(FPT.tab_AB[1], FPT.N, logkmin, logkmax, dlgk,lgk, 0.,0.);
}

double PT_d1s2(double k_coverH0){ //interpolate FPT.tab_AB[2]
  static cosmopara C; 
  static double logkmin = 0.,logkmax = 0.,dlgk = 0.;
  // only call FASTPT if cosmology changed since last call
  if (recompute_cosmo3D(C)){
    get_FPT_bias();
    update_cosmopara(&C); 
    double k_unit_conversion = 1./cosmology.coverH0;
    logkmin =log(FPT.k_min);
    logkmax = log(FPT.k_max);
    dlgk = log(10.)/(double) FPT.N_per_dec;
  }
  double lgk = log(k_coverH0);
  if (lgk < logkmin || lgk >= logkmax){return 0.;}
  return interpol(FPT.tab_AB[2], FPT.N, logkmin, logkmax, dlgk,lgk, 0.,0.);
}

double PT_d2s2(double k_coverH0){ //interpolate FPT.tab_AB[3]
  static cosmopara C; 
  static double logkmin = 0.,logkmax = 0.,dlgk = 0.;
  // only call FASTPT if cosmology changed since last call
  if (recompute_cosmo3D(C)){
    get_FPT_bias();
    update_cosmopara(&C); 
    logkmin =log(FPT.k_min);
    logkmax = log(FPT.k_max);
    dlgk = log(10.)/(double) FPT.N_per_dec;
  }
  double lgk = log(k_coverH0);
  if (lgk < logkmin || lgk >= logkmax){return 0.;}
  return interpol(FPT.tab_AB[3], FPT.N, logkmin, logkmax, dlgk,lgk, 0.,0.);
}

double PT_s2s2(double k_coverH0){ //interpolate FPT.tab_AB[4]
  static cosmopara C; 
  static double logkmin = 0.,logkmax = 0.,dlgk = 0.;
  // only call FASTPT if cosmology changed since last call
  if (recompute_cosmo3D(C)){
    get_FPT_bias();
    update_cosmopara(&C); 
    logkmin =log(FPT.k_min);
    logkmax = log(FPT.k_max);
    dlgk = log(10.)/(double) FPT.N_per_dec;
  }
  double lgk = log(k_coverH0);
  if (lgk < logkmin || lgk >= logkmax){return 0.;}
  return interpol(FPT.tab_AB[4], FPT.N, logkmin, logkmax, dlgk,lgk, 0.,0.);
}

double PT_sigma4(double k_coverH0){
  static cosmopara C; 
  // only call FASTPT if cosmology changed since last call
  if (recompute_cosmo3D(C)){
    get_FPT_bias();
    update_cosmopara(&C); 
  }
  return 0.0;
}

