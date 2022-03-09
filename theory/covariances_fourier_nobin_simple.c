/* 
// non-binned fourier space-only 10x2pt cov
// in a unified way
// -- Xiao Fang */

int func_for_Cl_Nl(double *Cij, double *Nij, double l, int ni, int nj, int is_lsi, int is_lsj);
double cov_G_AB_CD(char ABCD[2][4], double l, double delta_l, int z_ar[4], int is_ls[4]);
double bgal_a(double a, double nz);
double inner_project_tri_cov_AB_CD(double a,void *params);
double cov_NG_AB_CD(char ABCD[2][4], double l1,double l2, int z_ar[4], int is_ls[4]);
double tab_cov_NG_AB_CD(char ABCD[2][4], double l1, double l2, int z_ar[4], int is_ls[4]);

// utility routines

// Reads in the noise N_mv for mv quadratic estimator for d,
// without reduction due to average over ell-bin.
// return interpolated N_kk(\ell)
double kappa_reconstruction_noise(double l){
   
   static double *ell;
   static double *noise;
   static int nEll;

   static double ellmin = .0, ellmax = .0;

   if (noise==0){
      printf("run %s\n", cmb.name);
      // count lines
      nEll = line_count(cmb.pathLensRecNoise);
      printf("Reading CMB lensing noise: %s\n", cmb.pathLensRecNoise);

      // allocate ell and Nlkk
      ell = create_double_vector(0, nEll-1);
      noise = create_double_vector(0, nEll-1);
      // read each line
      FILE *file = fopen(cmb.pathLensRecNoise, "r");
      int iEll;
      for (iEll=0; iEll<nEll; iEll++) {
         fscanf(file, "%le %le", &ell[iEll], &noise[iEll]);
         noise[iEll] = log(noise[iEll]);
      }
      fclose(file);
      ellmax = ell[nEll-1];
      ellmin = ell[0];
   }
   // if l is in the range
   double f1;
   if (l<=0.) {
      return 0;
   } else if ((l>=ellmin) &&(l<=ellmax)){
      int iEll = 0;
      while (ell[iEll] < l) {
         iEll ++;
      }
      f1 = exp((noise[iEll]-noise[iEll-1])/log(ell[iEll]/ell[iEll-1])*log(l/ell[iEll-1]) + noise[iEll-1]);
      if (isnan(f1)){f1 = 0.;}
      // evaluate at that ell
      // C_ell^kk = l*(l+1)/4 * C_ell^dd
   } else{
      f1 = 0.;
   }
   // } else if(l<ellmin){
   //    f1 = exp((noise[1]-noise[0])/log(ell[1]/ell[0])*log(l/ell[0]) + noise[0]);
   // } else{
   //    f1 = exp((noise[nEll-1]-noise[nEll-2])/log(ell[nEll-1]/ell[nEll-2])*log(l/ell[nEll-1]) + noise[nEll-1]);
   // }
   return f1;
}

double y_reconstruction_noise(double l){
   
   static double *ell;
   static double *noise;
   static int nEll;

   static double ellmin = .0, ellmax = .0;
   double dummy;

   if (noise==0){
      printf("run %s\n", cmb.name);
      // count lines
      nEll = line_count(cmb.path_yNoise);
      printf("Reading tSZ y noise: %s\n", cmb.path_yNoise);

      // allocate ell and Nlkk
      ell = create_double_vector(0, nEll-1);
      noise = create_double_vector(0, nEll-1);
      // read each line
      FILE *file = fopen(cmb.path_yNoise, "r");
      int iEll;
      for (iEll=0; iEll<nEll; iEll++) {
         fscanf(file, "%le %le %le %le %le", &ell[iEll], &noise[iEll], &dummy, &dummy, &dummy);
         noise[iEll] = log(noise[iEll]);
         // printf("iELL, ell, noise, %d, %le, %le\n", iEll, ell[iEll], noise[iEll]);
      }
      fclose(file);
      ellmax = ell[nEll-1];
      ellmin = ell[0];
   }
   // if l is in the range
   double f1;
    if ((l>=ellmin) &&(l<=ellmax)){
      int iEll = 0;
      while (ell[iEll] < l) {
         iEll ++;
      }
      f1 = exp((noise[iEll]-noise[iEll-1])/log(ell[iEll]/ell[iEll-1])*log(l/ell[iEll-1]) + noise[iEll-1]);
      if (isnan(f1)){f1 = 0.;}
    } else{
      f1 = 0.;
    }
   return f1;
}


int func_for_Cl_Nl(double *Cij, double *Nij, double l, int ni, int nj, int is_lsi, int is_lsj){
  if(is_lsi==1 && is_lsj==1){// case: ll
    *Cij = C_cl_tomo_nointerp(l,ni,nj);
    if(ni==nj) {*Nij = 1./(nlens(ni)*survey.n_gal_conversion_factor);}
    return 0;
  }
  if(is_lsi*is_lsj==-1){ // case: ls, need to determine lens and src bins
    if(is_lsi==1){*Cij = C_gl_tomo_nointerp(l,ni,nj);}
    else{*Cij = C_gl_tomo_nointerp(l,nj,ni);}
    return 0;
  }
  if(is_lsi==-1 && is_lsj==-1){ // case: ss
    *Cij = C_shear_tomo_nointerp(l,ni,nj);
    if(ni==nj) {*Nij= pow(survey.sigma_e,2.0)/(2.0*nsource(ni)*survey.n_gal_conversion_factor);}
    return 0;
  }
  if(ni==-1 || nj==-1){ // at least one kcmb field
    if(is_lsi==1){*Cij = C_gk_nointerp(l,ni);return 0;}
    else if(is_lsj==1){*Cij = C_gk_nointerp(l,nj);return 0;}
    else if(is_lsi==-1){*Cij = C_ks_nointerp(l,ni);return 0;}
    else if(is_lsj==-1){*Cij = C_ks_nointerp(l,nj);return 0;}
    else if(ni==-2 || nj==-2){*Cij = C_ky_nointerp(l);return 0;}
    else if(ni==-1 && nj==-1){*Cij = C_kk_nointerp(l); *Nij = kappa_reconstruction_noise(l); return 0;}
  }
  if(ni==-2 || nj==-2){ // at least one y field, but ky not included
    if(is_lsi==1){*Cij = C_gy_nointerp(l,ni);return 0;}
    else if(is_lsj==1){*Cij = C_gy_nointerp(l,nj);return 0;}
    else if(is_lsi==-1){*Cij = C_sy_nointerp(l,ni);return 0;}
    else if(is_lsj==-1){*Cij = C_sy_nointerp(l,nj);return 0;}
    else if(ni==-2 && nj==-2){*Cij = C_yy_nointerp(l); *Nij = y_reconstruction_noise(l);return 0;}
  }
  printf("func_for_Cl_Nl: field combination not supported!\n");
  exit(1);
}

double cov_G_AB_CD(char ABCD[2][4], double l, double delta_l, int z_ar[4], int is_ls[4]){
  int i,j;
  double C13, C14, C23, C24, N13=0.0, N14=0.0, N23=0.0, N24=0.0;
  int n1,n2,n3,n4;
  int is_ls1, is_ls2, is_ls3, is_ls4;
  double fsky = survey.area*survey.area_conversion_factor/(4.*M_PI);

  n1 = z_ar[0]; n2 = z_ar[1]; n3 = z_ar[2]; n4 = z_ar[3];
  is_ls1 = is_ls[0]; is_ls2 = is_ls[1]; is_ls3 = is_ls[2]; is_ls4 = is_ls[3];

  func_for_Cl_Nl(&C13, &N13, l,n1,n3,is_ls1,is_ls3);
  func_for_Cl_Nl(&C14, &N14, l,n1,n4,is_ls1,is_ls4);
  func_for_Cl_Nl(&C23, &N23, l,n2,n3,is_ls2,is_ls3);
  func_for_Cl_Nl(&C24, &N24, l,n2,n4,is_ls2,is_ls4);
  
  char cmbprobes[3][4] = {"kk", "ky", "yy"};
  for(i=0;i<2;i++){
    for(j=0;j<3;j++){
      if(strcmp(ABCD[i],cmbprobes[j])==0){
        fsky = (fsky<cmb.fsky) ? cmb.fsky : fsky; // if one of the probes is from CMB, then use the cmb area (if it's larger)
        i=2; j=3; break; // break out double for-loop
      }
    }
  }

  return ((C13+N13)*(C24+N24)+(C14+N14)*(C23+N23))/(fsky*(2.*l+1.)*delta_l);
}

double bgal_a(double a, double nz){
  return gbias.b1_function(1./a-1.,(int)nz);
}
double inner_project_tri_cov_AB_CD(double a,void *params)
{
  double k[2],fK,weights,tri,res = 0.;
  double *ar = (double *) params;
  int i,j;
  fK = f_K(chi(a));
  k[0] = (ar[0]+0.5)/fK;
  k[1] = (ar[1]+0.5)/fK;
  int FLAG_y = 0, ni[4];

  weights = dchi_da(a);
  for(i=2;i<6;i++){
    ni[i-2] = ar[i];
    if(ar[i]==-1.){
      weights *= W_k(a,fK);
    }else if(ar[i]==-2.){
      weights *= W_y(a);
      FLAG_y = 1;
    }else if(ar[4+i]==1.){ // i.e. is_ls==1 -> gal density field
      weights *= W_gal(a,ar[i]);
    }else{ // i.e. is_ls==-1 -> gal shear field
      weights *= W_kappa(a,fK,ar[i]);
    }
  }

  double ssc[2], P_AB=0.;
  double sig_b, fsky1, fsky2, fsky_larger;
  fsky1 = ar[10]; fsky2 = ar[11];

  double (*func_for_delP_SSC)(double, double);

  for(i=0;i<2;i++){
    func_for_delP_SSC = &delP_SSC; // set delP_SSC as for matter field by default
    if(strcmp(pdeltaparams.runmode,"halomodel") ==0){func_for_delP_SSC = &delP_SSC_halo;} // use halo model version

    if(ar[2+2*i]==-2 && ar[2+2*i+1]==-2){ // if both A,B are y-fields
      P_AB = P_yy(k[i],a);
      func_for_delP_SSC = &delP_SSC_yy;
    }else if(ar[2+2*i]!=-2 && ar[2+2*i+1]!=-2){// if neither A,B is y-field
      P_AB = Pdelta(k[i],a);
    }else{ // if one of A,B is y-field, but the other is not.
      P_AB = P_my(k[i],a);
      func_for_delP_SSC = &delP_SSC_my;
    }

    ssc[i]=func_for_delP_SSC(k[i],a); // no galaxy density field

    for(j=0;j<2;j++){// if delta_g field, rescale wrt mean, Eq(A12,A13) in CosmoLike paper
      if(ar[6+2*i+j]==1){ssc[i] -= bgal_a(a,ar[2+2*i+j])*P_AB;}
    }
  }

  if(fsky1==fsky2){
    sig_b = survey_variance(a,fsky1);
  }else{
    sig_b = cross_survey_variance(a,fsky1,fsky2);
  }
  fsky_larger = (fsky1>fsky2 ? fsky1 : fsky2);

  if(weights >0.){
    if (covparams.cng) {
      if(FLAG_y==1){tri = tri_1h_y_cov(k[0],k[1],a,ni);}
      else{tri = tri_matter_cov(k[0],k[1],a);}
      res = tri*pow(fK,-6.)/(fsky_larger*4*M_PI);
    }
    // printf("res cNG:%lg , ", res);
    res += ssc[0]*ssc[1]*sig_b*pow(fK,-4.); //SSC
  }
  // printf("res +SSC, weights:%lg, %lg\n", res, weights);
  res *= weights;
  return res;
}

double cov_NG_AB_CD(char ABCD[2][4], double l1,double l2, int z_ar[4], int is_ls[4]){
  double amin[2], amax[2], a1,a2,array[12];
  int zmin, zmax, zlen, zs;
  int i;
  double fsky[2];

  double fsky_gal = survey.area/41253.0;
  for(i=0;i<2;i++){
    if(strcmp(ABCD[i], "ss")==0){
      zmin = (z_ar[2*i] < z_ar[2*i+1] ? z_ar[2*i] : z_ar[2*i+1]);
      amin[i] = amin_source(zmin); amax[i] = amax_source(zmin); fsky[i]=fsky_gal;
    }else if(strcmp(ABCD[i], "ls")==0 || strcmp(ABCD[i], "lk")==0 || strcmp(ABCD[i], "ly")==0){
      if(is_ls[2*i]==1){zlen = z_ar[2*i];}
      else{zlen = z_ar[2*i+1];}
      amin[i] = amin_lens(zlen); amax[i] = amax_lens(zlen); fsky[i]=fsky_gal;
    }else if(strcmp(ABCD[i], "ll")==0){
      zmin = (z_ar[2*i] < z_ar[2*i+1] ? z_ar[2*i] : z_ar[2*i+1]);
      zmax = (z_ar[2*i] > z_ar[2*i+1] ? z_ar[2*i] : z_ar[2*i+1]);
      amin[i] = amin_lens(zmax); amax[i] = amax_lens(zmin); fsky[i]=fsky_gal;
    }else if(strcmp(ABCD[i], "ks")==0 || strcmp(ABCD[i], "sy")==0){
      if(is_ls[2*i]==-1){zs = z_ar[2*i];}
      else{zs = z_ar[2*i+1];}
      amin[i] = amin_source(zs); amax[i] = amax_source(zs); fsky[i]=fsky_gal;
    }else if(strcmp(ABCD[i], "kk")==0 || strcmp(ABCD[i], "ky")==0 || strcmp(ABCD[i], "yy")==0){
      amin[i] = limits.a_min*(1.+1.e-5); amax[i] = 1.-1.e-5; fsky[i]=cmb.fsky;
    }

    if(strcmp(ABCD[i], "kk")!=0 && amin[i] < limits.a_min_hm) {
      amin[i] = limits.a_min_hm; // limit LOS integration up to a_min_hm as NG cov uses halomodel except for kcmb-kcmb
    }
    if(fsky[i]<=0.) {printf("fsky[%d]=%le is not set correctly!\n", i, fsky[i]); exit(1);}
  }
  a1 = (amin[0]>amin[1] ? amin[0] : amin[1]);
  a2 = (amax[0]<amax[1] ? amax[0] : amax[1]);

  array[0] = l1;
  array[1] = l2;
  for(i=0;i<4;i++){
    array[2+i] = (double) z_ar[i];
    array[6+i] = (double) is_ls[i];
  }
  for(i=0;i<2;i++){array[10+i] = fsky[i];}
  return int_gsl_integrate_low_precision(inner_project_tri_cov_AB_CD,(void*)array,a1,a2,NULL,1000);
}


// double tab_cov_NG_AB_CD(char ABCD[2][4], double l1, double l2, int z_ar[4], int is_ls[4]){
//   static int Z[4] = {-42,-42,-42,-42};
//   static int is_ls_local[4] = {-42,-42,-42,-42};
//   static int Ntab = 40;
//   static double ***table=0;
//   static double ds = .0, logsmin = .0, logsmax = .0;
//   logsmin = log(1.); logsmax = log(5.e+4);
//   ds = (logsmax - logsmin)/(Ntab - 1.);

//   const char *probes[10] = {"ss", "ls", "ll", "lk", "ks", "kk",\
//                             "ly", "sy", "ky", "yy"};
//   int N_covtype = 55;
//   int i_covtype = 0;

//   int i,j;
//   double res, llog1,llog2,ll1,ll2;

//   int COMPUTE = 0;
//   for(i=0;i<4;i++){
//     COMPUTE = (Z[i]!=z_ar[i] || COMPUTE);
//     COMPUTE = (is_ls_local[i]!=is_ls[i] || COMPUTE);
//     if(COMPUTE) {break;}
//   }

//   int isAB;
//   // Determine i_covtype: the first index in the table, corresponding to the cov type
//   for(i=0;i<10;i++){
//     isAB = (strcmp(ABCD[0], probes[i])==0);
//     for(j=0;j<=i;j++){
//       if(isAB && strcmp(ABCD[1], probes[j])==0){
//         i=10; j=10; // break out double for-loop
//         break;
//       }
//       i_covtype++;
//     }
//   }

//   if (COMPUTE){
//     if (table==0) {
//       table = (double***)malloc(N_covtype * sizeof(double**));
//       for(i=0;i<N_covtype;i++){
//         table[i] = NULL;
//       }
//     }
//     if(table[i_covtype]==NULL){
//       table[i_covtype] = create_double_matrix(0, Ntab-1, 0, Ntab-1);
//     }

//     llog1 = logsmin;
//     for (i=0; i<Ntab; i++, llog1+=ds) {
//       ll1 = exp(llog1);
//       llog2 = logsmin;
//       for (j=0; j<Ntab; j++, llog2+=ds) {
//         ll2 = exp(llog2);
//         table[i_covtype][i][j]=cov_NG_AB_CD(ABCD, ll1, ll2, z_ar, is_ls);
//       }
//     }
    
//     for(i=0;i<4;i++){
//       Z[i]=z_ar[i]; is_ls_local[i]=is_ls[i];
//     }
//   }
//   res = 0.;
//   if(l1<=1e-10 || l2<=1e-10) {return 0;}
//   llog1=log(l1);
//   llog2=log(l2);
//   if (llog1 > logsmin && llog2 > logsmin && llog1 < logsmax && llog2 < logsmax){
//     res = interpol2d(table[i_covtype], Ntab, logsmin, logsmax, ds, llog1, Ntab, logsmin, logsmax, ds, llog2,0.0,0.0);}
//   return res;
// }