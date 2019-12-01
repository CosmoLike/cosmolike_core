double PT_d1d2(double k_coverH0);
double PT_d1d3(double k_coverH0);
double PT_d2d2(double k_coverH0);
double PT_d1s2(double k_coverH0);
double PT_d2s2(double k_coverH0);
double PT_s2s2(double k_coverH0);
double PT_sigma4(double k_coverH0);
double TATT_II_EE (double k_coverH0, double a, double C1, double C2, double b_ta);
double TATT_II_BB (double k_coverH0, double a, double C1, double C2, double b_ta);
double TATT_GI_E(double k_coverH0, double a, double C1, double C2, double b_ta);
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
      if (FPT.k_min < limits.k_min_cH0){
        FPT.k_min = K_CH0(FPT.k_min);
        FPT.k_max = K_CH0(FPT.k_max);
      }
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


void get_FPT_IA(void){
  static cosmopara C; 
  if (recompute_cosmo3D(C)){
    update_cosmopara(&C); 
    if (FPT.tab_IA ==0){ //if table doesn't exit yet, create it
      FPT.tab_IA = create_double_matrix(0, FPT.N_IA-1, 0, FPT.N-1);
      if (FPT.k_min < limits.k_min_cH0){
        FPT.k_min = K_CH0(FPT.k_min);
        FPT.k_max = K_CH0(FPT.k_max);
      }
    }

    double k[FPT.N], Pin[FPT.N];
    double IA_tt_EE[FPT.N],IA_tt_BB[FPT.N];
    double IA_ta_dE1[FPT.N], IA_ta_dE2[FPT.N], IA_ta_0E0E[FPT.N], IA_ta_0B0B[FPT.N];
    double IA_mix_A[FPT.N], IA_mix_B[FPT.N], IA_mix_DEE[FPT.N], IA_mix_DBB[FPT.N];

    FPT_input(k, Pin); 
    int i;
    IA_tt(k, Pin, FPT.N, IA_tt_EE, IA_tt_BB);
    for (i =0; i < FPT.N; i++){
      FPT.tab_IA[0][i] = IA_tt_EE[i]; //tt_EE
      FPT.tab_IA[1][i] = IA_tt_BB[i]; //tt_BB
      printf("tt k=%e, Plin =%e   %e %e\n",k[i],Pin[i], FPT.tab_IA[0][i],FPT.tab_IA[1][i]);
    }

    IA_ta(k, Pin, FPT.N, IA_ta_dE1, IA_ta_dE2, IA_ta_0E0E, IA_ta_0B0B);
    for (i =0; i < FPT.N; i++){
      FPT.tab_IA[2][i] = IA_ta_dE1[i]; //ta_dE1
      FPT.tab_IA[3][i] = IA_ta_dE2[i]; //ta_dE2
      FPT.tab_IA[4][i] = IA_ta_0E0E[i]; //ta_EE
      FPT.tab_IA[5][i] = IA_ta_0B0B[i]; //ta_BB
      printf("ta k=%e, Plin =%e  %e %e %e %e\n",k[i],Pin[i], FPT.tab_IA[2][i],FPT.tab_IA[3][i],FPT.tab_IA[4][i],FPT.tab_IA[5][i]);
    }

    IA_mix(k,Pin, FPT.N, IA_mix_A, IA_mix_B, IA_mix_DEE, IA_mix_DBB);
    for (i =0; i < FPT.N; i++){
      FPT.tab_IA[6][i] = IA_mix_A[i]; //mix_A
      FPT.tab_IA[7][i] = IA_mix_B[i]; //mix_B
      FPT.tab_IA[8][i] = IA_mix_DEE[i]; //mix_D_EE
      FPT.tab_IA[9][i] = IA_mix_DBB[i]; //mix_D_BB
      printf("mix k=%e, Plin =%e  %e %e %e %e\n",k[i],Pin[i], FPT.tab_IA[6][i],FPT.tab_IA[7][i],FPT.tab_IA[8][i],FPT.tab_IA[9][i]);
    }
  }
}

double TATT_II_EE (double k_coverH0, double a, double C1, double C2, double b_ta){
  double nla_EE = 0., ta_EE = 0., tt_EE = 0., mix_EE = 0.;
  nla_EE = C1*C1*Pdelta(k_coverH0,a);

  static cosmopara C; 
  static double logkmin = 0.,logkmax = 0.,dlgk = 0.;
  // if TATT is specified, i.e. C2 !=0 (TT), or b_ta != 0 (TA)
  if (C2 != 0 || b_ta != 0){
    // only call FASTPT if cosmology changed since last call
    if (recompute_cosmo3D(C)){
      get_FPT_IA();
      update_cosmopara(&C); 
      logkmin =log(FPT.k_min);
      logkmax = log(FPT.k_max);
      dlgk = log(10.)/(double) FPT.N_per_dec;
    }
    double lgk = log(k_coverH0);
    if (lgk <= logkmin || lgk >= logkmax){return nla_EE;}
    double g4 = pow(growfac(a)/growfac(0.99),4.0);


    double P_ta_EE, P_ta_dE1, P_ta_dE2;
    P_ta_EE = g4*interpol(FPT.tab_IA[4], FPT.N, logkmin, logkmax, dlgk,lgk, 0.,0.);
    P_ta_dE1 = g4*interpol(FPT.tab_IA[2], FPT.N, logkmin, logkmax, dlgk,lgk, 0.,0.);
    P_ta_dE2 = g4*interpol(FPT.tab_IA[3], FPT.N, logkmin, logkmax, dlgk,lgk, 0.,0.);

    ta_EE = C1*C1*(b_ta*b_ta*P_ta_EE + 2.*b_ta*(P_ta_dE1+P_ta_dE2));

    if (C2){
      double P_mix_EE, P_mix_A, P_mix_B;
      P_mix_EE = g4*interpol(FPT.tab_IA[8], FPT.N, logkmin, logkmax, dlgk,lgk, 0.,0.);
      P_mix_A = g4*interpol(FPT.tab_IA[6], FPT.N, logkmin, logkmax, dlgk,lgk, 0.,0.);
      P_mix_B = g4*interpol(FPT.tab_IA[7], FPT.N, logkmin, logkmax, dlgk,lgk, 0.,0.);

      mix_EE = 2.*C1*C2*(P_mix_A + 2.*P_mix_B + b_ta*P_mix_EE);
 
      tt_EE = C2*C2*g4*interpol(FPT.tab_IA[0], FPT.N, logkmin, logkmax, dlgk,lgk, 0.,0.);
    }
  }

  return nla_EE + ta_EE + mix_EE + tt_EE;
}

double TATT_II_BB (double k_coverH0, double a, double C1, double C2, double b_ta){
  double ta_BB = 0., tt_BB = 0., mix_BB = 0.;

  static cosmopara C; 
  static double logkmin = 0.,logkmax = 0.,dlgk = 0.;
  // if TATT is specified, i.e. C2 !=0 (TT), or b_ta != 0 (TA)
  if (C2 != 0 || b_ta != 0){
    // only call FASTPT if cosmology changed since last call
    if (recompute_cosmo3D(C)){
      get_FPT_IA();
      update_cosmopara(&C); 
      logkmin =log(FPT.k_min);
      logkmax = log(FPT.k_max);
      dlgk = log(10.)/(double) FPT.N_per_dec;
    }
    double lgk = log(k_coverH0);
    if (lgk < logkmin || lgk >= logkmax){return 0.0;}
    double g4 = pow(growfac(a)/growfac(0.99),4.0);


    ta_BB = C1*C1*b_ta*b_ta*g4*interpol(FPT.tab_IA[5], FPT.N, logkmin, logkmax, dlgk,lgk, 0.,0.);

    if (C2){
      mix_BB = 2.*C1*C2*b_ta*g4*interpol(FPT.tab_IA[9], FPT.N, logkmin, logkmax, dlgk,lgk, 0.,0.);

      tt_BB = C2*C2*g4*interpol(FPT.tab_IA[1], FPT.N, logkmin, logkmax, dlgk,lgk, 0.,0.);
    }
  }
  return ta_BB + mix_BB + tt_BB;
}

double TATT_GI_E (double k_coverH0, double a, double C1, double C2, double b_ta){
  double nla_GI = 0., ta_GI = 0., tt_GI = 0., mix_GI = 0.;
 // printf("in GI\n");
  nla_GI = C1*Pdelta(k_coverH0,a);
  static cosmopara C; 
  static double logkmin = 0.,logkmax = 0.,dlgk = 0.;
  // if TATT is specified, i.e. C2 !=0 (TT), or b_ta != 0 (TA)
  if (C2 != 0 || b_ta != 0){
    // only call FASTPT if cosmology changed since last call
    if (recompute_cosmo3D(C)){
      get_FPT_IA();
      update_cosmopara(&C); 
      logkmin =log(FPT.k_min);
      logkmax = log(FPT.k_max);
      dlgk = log(10.)/(double) FPT.N_per_dec;
    }
    double lgk = log(k_coverH0);
    if (lgk < logkmin || lgk >= logkmax){return nla_GI;}
    double g4 = pow(growfac(a)/growfac(0.99),4.0);

    double P_ta_dE1, P_ta_dE2;
    P_ta_dE1 = g4*interpol(FPT.tab_IA[2], FPT.N, logkmin, logkmax, dlgk,lgk, 0.,0.);
    P_ta_dE2 = g4*interpol(FPT.tab_IA[3], FPT.N, logkmin, logkmax, dlgk,lgk, 0.,0.);

    ta_GI = C1*b_ta*(P_ta_dE1 + P_ta_dE2);

    if (C2){
      double P_mix_A, P_mix_B;
      P_mix_A = g4*interpol(FPT.tab_IA[6], FPT.N, logkmin, logkmax, dlgk,lgk, 0.,0.);
      P_mix_B = g4*interpol(FPT.tab_IA[7], FPT.N, logkmin, logkmax, dlgk,lgk, 0.,0.);

      tt_GI = C2*(P_mix_A + 2.*P_mix_B);
    }
  }
 // printf("completed GI\n");
  return nla_GI + ta_GI  + tt_GI;
}