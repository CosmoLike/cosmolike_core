#include "IA.h"

void set_LF_GAMA(void);
void set_LF_DEEP2(void);
int check_LF (void);
double M_abs(double mag, double a); //absolute magnitude corresponding to aparent magnitude at a; h = 1 units, incl k-corrections
double n_all_LF(double mag, double a);
double f_red_LF(double mag, double a);
double A_LF(double mag, double a);
double A_IA_Joachimi(double a);



/******** select LF ************/
void set_LF_GAMA(void){
  int i;
  for (i = 0; i <5; i++){
    LF_coefficients[0][i] =LF_coefficients_GAMA[0][i];
    LF_coefficients[1][i] =LF_coefficients_GAMA[1][i];
  }
}
void set_LF_DEEP2(void){
  int i;
  for (i = 0; i <5; i++){
    LF_coefficients[0][i] =LF_coefficients_DEEP2[0][i];
    LF_coefficients[1][i] =LF_coefficients_DEEP2[1][i];
  }
}
/************** normalization rouintes *********************/
double M_abs(double mag, double a){ //in h = 1 units, incl. Poggianti 1997 k+e-corrections
  static double *table;
  static double dz;
  double ke;
  if (table ==0){
    int i;
    //read in + tabulate k+e corrections for early types, restframe r band
    //interpolated from http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=J/A%2BAS/122/399
    table = create_double_vector(0,30);
    dz = 0.1;
    /*FILE *ein;
    double d1,d2,d3;
    ein = fopen(survey.Kcorrect_File,"r");
    EXIT_MISSING_FILE(ein, "LoverL0",survey.Kcorrect_File)*/
    for (i = 0; i< 31; i++){
   /*   fscanf(ein,"%le %le %le\n",&d1,&d2,&d3);
       table[i] = d2+d3; */
      table[i] = KE[i];
    }
   /* fclose(ein);*/
  }
  double z=1./a-1.0;
  if(z >= 3.0) {z=2.99;}  // no acceptable k-korrection exists for k>3, also no meaningful IA model
  ke = interpol(table, 31, 0., 3.0, 0.1, z, 1.0, 1.0);

  return mag - 5.0*log10(f_K(chi(a))/a*cosmology.coverH0) -25.0 -ke;
}
double n_all_LF(double mag, double a){
   double Q, alpha, Mstar, Phistar,P,Mlim, LF_all[3];
  if (LF_coefficients[0][0]< 0){printf("missing Lum fct."); exit(1);}

  //r-band LF parameters, Tab. 5 in http://arxiv.org/pdf/1111.0166v2.pdf
  //all galaxies
  alpha =LF_coefficients[1][1]+ nuisance.LF_alpha;
  Mstar = LF_coefficients[1][2];
  Q = LF_coefficients[1][3]+nuisance.LF_Q;
  P = LF_coefficients[1][4]+nuisance.LF_P;
  Phistar = LF_coefficients[1][0];

  LF_all[0] = Phistar*pow(10.0,0.4*P*(1./a-1));
  LF_all[1] = Mstar -Q*(1./a-1. - 0.1);
  LF_all[2]= alpha;
  Mlim = M_abs(mag,a); //also in h = 1 units


  return LF_all[0]*gsl_sf_gamma_inc(LF_all[2]+1,pow(10.0,-0.4*(Mlim-LF_all[1])));
}

double f_red_LF(double mag, double a){
  double Q, alpha, Mstar, Phistar,P,Q_red, alpha_red, Mstar_red, Phistar_red,P_red,Mlim, LF_all[3], LF_red[3];
  if (LF_coefficients[0][0]< 0){printf("missing Lum fct."); exit(1);}

  //r-band LF parameters, Tab. 5 in http://arxiv.org/pdf/1111.0166v2.pdf
  //red galaxies
  alpha_red = LF_coefficients[0][1] + nuisance.LF_red_alpha;
  Mstar_red = LF_coefficients[0][2]; //in h = 1 units
  Q_red = LF_coefficients[0][3]+nuisance.LF_red_Q;
  P_red = LF_coefficients[0][4]+nuisance.LF_red_P;
  Phistar_red = LF_coefficients[0][0];

  LF_red[0] = Phistar_red*pow(10.0,0.4*P_red*(1./a-1));
  LF_red[1] = Mstar_red-Q_red*(1./a-1. - 0.1);
  LF_red[2]= alpha_red;

  //all galaxies
  alpha =LF_coefficients[1][1]+ nuisance.LF_alpha;
  Mstar = LF_coefficients[1][2];
  Q = LF_coefficients[1][3]+nuisance.LF_Q;
  P = LF_coefficients[1][4]+nuisance.LF_P;
  Phistar = LF_coefficients[1][0];

  LF_all[0] = Phistar*pow(10.0,0.4*P*(1./a-1));
  LF_all[1] = Mstar -Q*(1./a-1. - 0.1);
  LF_all[2]= alpha;
  Mlim = M_abs(mag,a); //also in h = 1 units

  return LF_red[0]/LF_all[0]*gsl_sf_gamma_inc(LF_red[2]+1,pow(10.0,-0.4*(Mlim-LF_red[1])))/gsl_sf_gamma_inc(LF_all[2]+1,pow(10.0,-0.4*(Mlim-LF_all[1])));
}

double A_LF(double mag, double a){// averaged (L/L_0)^beta over red galaxy LF
    double Q_red, alpha_red, Mstar_red,x,Mlim, LF_red[3], Lstar,L0;
  if (LF_coefficients[0][0]< 0){printf("missing Lum fct."); exit(1);}

    //r-band LF parameters, Tab. 5 in http://arxiv.org/pdf/1111.0166v2.pdf
    //red galaxies
  alpha_red = LF_coefficients[0][1]+ nuisance.LF_red_alpha;
  Mstar_red = LF_coefficients[0][2]; //in h = 1 units
  Q_red = LF_coefficients[0][3]+nuisance.LF_red_Q;

    LF_red[1] = Mstar_red-Q_red*(1./a-1. -.1);
    LF_red[2]= alpha_red;

    Mlim = M_abs(mag,a);
    Lstar = pow(10.,-0.4*LF_red[1]);
    x = pow(10.,-0.4*(Mlim-LF_red[1])); //Llim/Lstar
    L0 =pow(10.,-0.4*(-22.)); //all in h = 1 units
    return pow(Lstar/L0, nuisance.beta_ia)*gsl_sf_gamma_inc(LF_red[2]+nuisance.beta_ia+1,x)/gsl_sf_gamma_inc(LF_red[2]+1,x);
}


double A_LF_all(double mag, double a){// averaged (L/L_0)^beta over red galaxy LF
    double Q_red, alpha_red, Mstar_red,x,Mlim, LF_red[3], Lstar,L0;
  if (LF_coefficients[0][0]< 0){printf("missing Lum fct."); exit(1);}

    //r-band LF parameters, Tab. 5 in http://arxiv.org/pdf/1111.0166v2.pdf
    //all galaxies
  alpha_red = LF_coefficients[1][1]+ nuisance.LF_alpha;
  Mstar_red = LF_coefficients[1][2]; //in h = 1 units
  Q_red = LF_coefficients[1][3]+nuisance.LF_Q;

    LF_red[1] = Mstar_red-Q_red*(1./a-1. -.1);
    LF_red[2]= alpha_red;

    Mlim = M_abs(mag,a);
    Lstar = pow(10.,-0.4*LF_red[1]);
    x = pow(10.,-0.4*(Mlim-LF_red[1])); //Llim/Lstar
    L0 =pow(10.,-0.4*(-22.)); //all in h = 1 units
    return pow(Lstar/L0, nuisance.beta_ia)*gsl_sf_gamma_inc(LF_red[2]+nuisance.beta_ia+1,x)/gsl_sf_gamma_inc(LF_red[2]+1,x);
}


int check_LF(void){ //return 1 if combination of all + red galaxy LF parameters is unphysical, i.e. if f_red > 1 for some z < redshift.shear_zdistrpar_zmax
  double a=1./(1+redshift.shear_zdistrpar_zmax)+0.005;
  while (a < 1.){
    if( M_abs(survey.m_lim,a) < LF_coefficients[1][2]-(LF_coefficients[1][3]+nuisance.LF_Q)*(1./a-1. - 0.1) || M_abs(survey.m_lim,a) < LF_coefficients[0][2]-(LF_coefficients[0][3]+nuisance.LF_red_Q)*(1./a-1. - 0.1)){return 1;}
    if (f_red_LF(survey.m_lim,a) > 1.0){return 1;}
    a+=0.01;
  }
  return 0;
}

/*=========================================================*/
/*=============  Intrinsic Alignment models  ==============*/
/*=========================================================*/


double A_IA_Joachimi(double a){
  double z, A_red;
  z = 1./a-1;
  A_red = 5.92*A_LF(survey.m_lim,a)*f_red_LF(survey.m_lim,a); //A_0*<(L/L_0)^beta>*f_red
  return -A_red*pow((1.+z)/1.3,-0.47);//*cosmology.Omega_m*nuisance.c1rhocrit_ia*growfac(1.)/growfac(a);
}
