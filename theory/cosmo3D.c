#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "modgrav.h"

#ifndef USE_HICLASS
    #include "../class/include/class.h"
#else
    #include "../hi_class/include/class.h"
#endif
#define EXIT_IF_INVALID_FOR_HORNDESKI

#include <gsl/gsl_odeiv.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_errno.h>

//====
//in 1/Mpc 
#define k_max_coyote 8.5691 
#define k_min_coyote 0.0008694
//#define k_min_coyote 0.0008694
//===
#define a_max_coyote 1.0
#define a_min_coyote 0.2

#define omhh_max_coyote 1.550000e-01
#define omhh_min_coyote 1.200000e-01
#define ombhh_max_coyote 2.350000e-02
#define ombhh_min_coyote 2.150000e-02
#define ns_max_coyote 1.050000e+00
#define ns_min_coyote 8.500000e-01
#define h0_max_coyote 8.500000e-01 
#define h0_min_coyote 5.500000e-01
#define w_max_coyote -7.000000e-01
#define w_min_coyote -1.300000e+00
#define s8_max_coyote 9.000000e-01
#define s8_min_coyote 6.000000e-01

//void omega_a(double aa,double *om_m,double *om_v);
double omv_vareos(double a);
static inline double hoverh0(double a);
double growfac(double a);
int func_for_growfac(double a,const double y[],double f[],void *params);
double Tsqr_EH_wiggle(double khoverMPC);
//double int_for_sigma_r_sqr(double k, void * args);
double sigma_r_sqr();
//double Delta_L_wiggle(double k);
//double Delta_lin_wiggle(double k,double a);

double evaluate_class_tables(double k_coverh0, double a, int NL, int mode, int *status);
double p_lin(double k,double a);
double int_sig_R_knl(double logk, void *args);
double int_neff(double lnk, void *args);
double int_cur(double lnk, void *args);
void nonlin_scale(double amp, double *R_NL, double *neff, double *Curv);
double Halofit(double k, double amp, double omm, double omv,double w_z, double R_NL, double neff,double Curv, double P_delta_Lin);
void Delta_halofit(double **table_P_NL,double logkmin, double logkmax, double dk, double da);
double Delta_NL_Halofit(double k_NL, double a); //k in h/Mpc
double Delta_NL_Coyote(double k_NL,double a); //k in h/Mpc
double Delta_NL_Coyote_only(double k_NL,double a); //k in h/Mpc
double Pdelta(double k_NL,double a); //k in coverH0 units

//double int_for_chi(double a,void * args);
double f_K(double chi);
double chi(double a);
extern void emu(double *xstar, double *ystar, int *outtype);

// Errors returned by CLASS (CLASS_SUCCESS is defined to always be 0)
enum { 
  CLASS_SUCCESS = 0,
  CLASS_ERROR_PARSER,
  CLASS_ERROR_INPUT,
  CLASS_ERROR_BACKGROUND, 
  CLASS_ERROR_THERMODYNAMICS,
  CLASS_ERROR_PERTURB,
  CLASS_ERROR_PRIMORDIAL,
  CLASS_ERROR_NONLINEAR,
  CLASS_ERROR_TRANSFER,
  CLASS_ERROR_SPECTRA,
  CLASS_ERROR_HORNDESKI_STABILITY
};

// Types of data that can be returned from CLASS (used by evaluate_class_tables)
enum { 
  CLASS_RETURN_POWSPEC, 
  CLASS_RETURN_HORNDESKI_MU_EFF, 
  CLASS_RETURN_HORNDESKI_MU_LIGHT, 
  CLASS_RETURN_HORNDESKI_OMEGASMG,
  CLASS_RETURN_HORNDESKI_ALPHA_B,
  CLASS_RETURN_HORNDESKI_ALPHA_K,
  CLASS_RETURN_HORNDESKI_ALPHA_M,
  CLASS_RETURN_HORNDESKI_ALPHA_T,
  CLASS_RETURN_HORNDESKI_CS2,
};

//c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//variable Omega_v
//c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double omv_vareos(double a)
{
  // FIXME: This should properly account for Horndeski
  return(cosmology.Omega_v*exp(-3.*((cosmology.w0+cosmology.wa+1.)*log(a)+cosmology.wa*(1.-a))));
}

//c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//c evolution of omega matter and omega lamda with expansion factor

void omega_a(double aa,double *om_m,double *om_v)
{
  double a2,omega_curv;
  a2=aa*aa;
  omega_curv=1.0-cosmology.Omega_m- cosmology.Omega_v;
  *om_m=cosmology.Omega_m /(cosmology.Omega_m +aa*(omv_vareos(aa) *a2 +omega_curv));
  *om_v=omv_vareos(aa)*a2*aa/(cosmology.Omega_m+aa*(a2*omv_vareos(aa) +omega_curv));
}
//c%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//growth factor including Dark energy parameters w0, wa
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

//function for growfac (DGL)
int func_for_growfac(double a,const double y[],double f[],void *params)
{
  //double *p=(double *)params;
  if (a == 0) {
    printf("a=0 in function 'func_for_growfac'!\n");
    exit(1);
  }
  double aa=a*a;
  double omegam=cosmology.Omega_m/(aa*a);
  double omegav=omv_vareos(a);
  double hub = hoverh0(a);
  double one_plus_mg_mu = 1.;
  hub = hub*hub;
  f[0] = y[1];
  
  // Modified gravity correction to growth
  if(cosmology.MGmu != 0){
    // Phenomenological mu-sigma parametrisation
    one_plus_mg_mu += cosmology.MGmu * omegav/hub/cosmology.Omega_v;
  }else if(cosmology.use_horndeski == 1){
    // Horndeski parametrisation (assumes quasi-static approximation)
    one_plus_mg_mu += horndeski_mu(a);
  }
  
  f[1] = y[0] * 3.*cosmology.Omega_m/(2.*hub*aa*aa*a) * one_plus_mg_mu
       - y[1] / a*(2.-(omegam + (3.*(cosmology.w0+cosmology.wa*(1.-a))+1)*omegav)/(2.*hub));
  return GSL_SUCCESS;
}

static inline double hoverh0(double a){
  return sqrt(cosmology.Omega_m /(a*a*a) + (1.-cosmology.Omega_m -cosmology.Omega_v )/(a*a) + omv_vareos(a) );
}

double growfac(double a)
{
  const double MINA=1.e-8;
  static cosmopara C;
  static double *ai;
  static double *table;
  double res;
  
  
  
  gsl_interp *intf=gsl_interp_alloc(gsl_interp_linear,Ntable.N_a);
  gsl_interp_accel *acc=gsl_interp_accel_alloc();
  
  if (recompute_expansion(C))
  {
    
    if(table!=0) free_double_vector(table,0, Ntable.N_a-1);
    if(ai!=0) free_double_vector(ai,0, Ntable.N_a-1);    
    ai=create_double_vector(0, Ntable.N_a-1);
    table=create_double_vector(0, Ntable.N_a-1);
    
    int i;
    const gsl_odeiv_step_type *T=gsl_odeiv_step_rkf45;
    gsl_odeiv_step *s=gsl_odeiv_step_alloc(T,2);
    gsl_odeiv_control *c=gsl_odeiv_control_y_new(1.e-6,0.0);
    gsl_odeiv_evolve *e=gsl_odeiv_evolve_alloc(2);
    
    double t=MINA;            //start a
    double t1=1.1;                //final a
    double h=1.e-6;              //initial step size
    double y[2]={MINA,MINA};   //initial conditions
    double norm;
    double par[0]={};
    gsl_odeiv_system sys={func_for_growfac,NULL,2,&par};
    
    for (i=1;i<=Ntable.N_a;i++) {
      ai[i-1]=i*t1/(1.*Ntable.N_a);
      while(t<ai[i-1]) 
        gsl_odeiv_evolve_apply(e,c,s,&sys,&t,ai[i-1],&h,y);
      if (i==1) norm=y[0]/ai[i-1];
      table[i-1]=y[0]/norm;
    }
    
    gsl_odeiv_evolve_free(e);
    gsl_odeiv_control_free(c);
    gsl_odeiv_step_free(s);
    update_cosmopara(&C);
  }
  gsl_interp_init(intf,ai,table,Ntable.N_a);
  res=gsl_interp_eval(intf,ai,table,a,acc);
  gsl_interp_accel_free(acc);
  gsl_interp_free(intf);
  return(res);
}



// ---------------------------- Transfer Function from EH98 ---------------------- 
//Input: k -- Wavenumber at which to calculate transfer function, in Mpc^-1. Output: Returns the value of the full transfer function fitting formula. This is the form given in Section 3 of Eisenstein & Hu (1997). Notes: Units are Mpc, not h^-1 Mpc. 

double Tsqr_EH_wiggle(double khoverMPC)
{
  static double omhh=-123.;
  static double obhh=-123.;
  static double OMEGA_V = -123.;
  static double f_baryon;
  
  static double k_equality;
  static double sound_horizon;
  static double beta_c;
  static double alpha_c;
  static double beta_node;
  static double alpha_b;
  static double beta_b;
  static double k_silk;
  
  //invalid_for_horndeski(__func__);
  
  //if (omhh != cosmology.Omega_m*cosmology.h0*cosmology.h0 || obhh != cosmology.omb*cosmology.h0*cosmology.h0|| OMEGA_V != cosmology.Omega_v){  
    double theta_cmb,z_equality,z_drag,R_drag,R_equality;
    double z_drag_b1, z_drag_b2;
    double alpha_c_a1, alpha_c_a2, beta_c_b1, beta_c_b2, alpha_b_G, y;
    
    omhh = cosmology.Omega_m*cosmology.h0*cosmology.h0;  
    obhh = cosmology.omb*cosmology.h0*cosmology.h0;
    OMEGA_V = cosmology.Omega_v;
    f_baryon = obhh/omhh; 
    //printf("%le\n",f_baryon);
    
    theta_cmb=2.728/2.7;// Tcmb in units of 2.7 K 
    z_equality= 2.50e4*omhh/POW4(theta_cmb);//Redshift of matter-radiation equality, really 1+z 
    
    k_equality=0.0746*omhh/SQR(theta_cmb);//Scale of equality, in Mpc^-1 
    z_drag_b1 = 0.313*pow(omhh,-0.419)*(1+0.607*pow(omhh,0.674));
    z_drag_b2 = 0.238*pow(omhh,0.223);
    z_drag=1291*pow(omhh,0.251)/(1+0.659*pow(omhh,0.828))*(1+z_drag_b1*pow(obhh,z_drag_b2));//Redshift of drag epoch 
    R_drag=31.5*obhh/POW4(theta_cmb)*(1000/(1+z_drag));//Photon-baryon ratio at drag epoch 
    
    R_equality= 31.5*obhh/POW4(theta_cmb)*(1000/z_equality);//Photon-baryon ratio at equality epoch 
    sound_horizon=2./3./k_equality*sqrt(6./R_equality)*log((sqrt(1+R_drag)+sqrt(R_drag+R_equality))/(1+sqrt(R_equality)));//Sound horizon at drag epoch, in Mpc 
    k_silk= 1.6*pow(obhh,0.52)*pow(omhh,0.73)*(1+pow(10.4*omhh,-0.95));//Silk damping scale, in Mpc^-1 
    alpha_c_a1 = pow(46.9*omhh,0.670)*(1+pow(32.1*omhh,-0.532));
    alpha_c_a2 = pow(12.0*omhh,0.424)*(1+pow(45.0*omhh,-0.582));
    alpha_c=pow(alpha_c_a1,-f_baryon)*pow(alpha_c_a2,-CUBE(f_baryon)); //CDM suppression 
    
    beta_c_b1 = 0.944/(1+pow(458*omhh,-0.708));
    beta_c_b2 = pow(0.395*omhh, -0.0266);
    beta_c=1.0/(1+beta_c_b1*(pow(1-f_baryon, beta_c_b2)-1));//CDM log shift 
    
    y = z_equality/(1+z_drag);
    alpha_b_G = y*(-6.*sqrt(1+y)+(2.+3.*y)*log((sqrt(1+y)+1)/(sqrt(1+y)-1)));
    alpha_b=2.07*k_equality*sound_horizon*pow(1+R_drag,-0.75)*alpha_b_G;//Baryon suppression 
    beta_node = 8.41*pow(omhh, 0.435);//Sound horizon shift 
    beta_b= 0.5+f_baryon+(3.-2.*f_baryon)*sqrt(pow(17.2*omhh,2.0)+1);//Baryon envelope shift 
  //}  
  // Start of TFfit_onek from the original tf_fit.c at http://background.uchicago.edu/~whu/transfer/transferpage.html
  double T_c_ln_beta,T_c_ln_nobeta , T_c_C_alpha, T_c_C_noalpha;
  double q,qsqr, xx, xx_tilde;
  double T_c_f, T_c, s_tilde, T_b_T0, T_b;
  
  double k=khoverMPC*cosmology.h0; //internally this routine uses Mpc^-1 not h/Mpc
  q = k/13.41/k_equality; 
  qsqr =SQR(q);
  xx = k*sound_horizon;
  
  T_c_ln_beta = log(2.718282+1.8*beta_c*q);
  T_c_ln_nobeta = log(2.718282+1.8*q);
  
  T_c_ln_beta = log(2.718282+1.8*beta_c*q);
  T_c_C_alpha = 14.2/alpha_c + 386.0/(1+69.9*pow(q,1.08));
  T_c_C_noalpha = 14.2 + 386.0/(1+69.9*pow(q,1.08));
  
  T_c_f = 1.0/(1.0+POW4(xx/5.4));
  T_c = T_c_f*T_c_ln_beta/(T_c_ln_beta+T_c_C_noalpha*qsqr) +(1-T_c_f)*T_c_ln_beta/(T_c_ln_beta+T_c_C_alpha*qsqr);
  
  s_tilde = sound_horizon*pow(1+CUBE(beta_node/xx),-1./3.);
  xx_tilde = k*s_tilde;
  
  T_b_T0 = T_c_ln_nobeta/(T_c_ln_nobeta+T_c_C_noalpha*qsqr);
  
  T_b = sin(xx_tilde)/(xx_tilde)*(T_b_T0/(1+SQR(xx/5.2))+alpha_b/(1+CUBE(beta_b/xx))*exp(-pow(k/k_silk,1.4)));
  
  return SQR(f_baryon*T_b + (1-f_baryon)*T_c); 
}


//Calculate Normalization see Cosmology Notes 8.105  
double int_for_sigma_r_sqr(double k, void * args)
{
  double kR, res, x;
  kR = k*8.; // r=8 Mpc/h
  x = (sin(kR) - kR*cos(kR))/(kR*kR*kR);
  res = pow(k,2.+cosmology.n_spec)*Tsqr_EH_wiggle(k)*x*x;
  return res;
}

double sigma_r_sqr()   
{
  static double res = -123.;
  static cosmopara C;
  double integral,array[1];
  
  if (recompute_Delta(C)) //strictly speaking, this is recomputed unnecessarily if only sigma_8 changes
  {
    integral = int_gsl_integrate_medium_precision(int_for_sigma_r_sqr,(void*)array,1e-4,1e6,NULL,512);
    res = 9.0*integral;   //see Peackock97, eq. 29
  update_cosmopara(&C);

  }
  if (!(res>0.0)){
    fprintf(stderr,"failed with sigma_r_sqr = %le\n", res);
  }
  assert(res>0.0);
  return res;
}



double Delta_L_wiggle(double k)
{
  static cosmopara C;
  
  static double *table_P;
  static double dk = .0, logkmin = .0, logkmax = .0;
  
  double klog,f1,norm;
  int i;     
  
  if (k < limits.k_min_mpc || k > limits.k_max_mpc){
    norm=cosmology.sigma_8*cosmology.sigma_8/sigma_r_sqr(); 
    
    return norm*pow(k,cosmology.n_spec+3.0)*Tsqr_EH_wiggle(k);
    //printf("outside Delta_L_tab\n");   
  }
  else{  
    if (recompute_Delta(C))
    {
      if (cosmology.M_nu > 0){
        printf("Implementation of EH transfer function does not support massive neutrinos\n EXIT\n");
      }
      update_cosmopara(&C);
      norm=cosmology.sigma_8*cosmology.sigma_8/sigma_r_sqr();
      
      if(table_P!=0) free_double_vector(table_P,0, Ntable.N_k_lin-1);
      table_P=create_double_vector(0, Ntable.N_k_lin-1);
      
      logkmin = log(limits.k_min_mpc);
      logkmax = log(limits.k_max_mpc);
      dk = (logkmax - logkmin)/(Ntable.N_k_lin-1.);
      klog = logkmin;
      
      for (i=0; i<Ntable.N_k_lin; i++, klog += dk) {  
	table_P[i]=log(norm*pow(exp(klog),cosmology.n_spec+3.0)*Tsqr_EH_wiggle(exp(klog)));
      }
      //printf("finished Delta_L_wiggle\n");   
    }
  }
  klog=log(k);
  f1=interpol(table_P, Ntable.N_k_lin, logkmin, logkmax, dk,klog, 1.0,1.0 );
  return exp(f1);  
}


// double Delta_lin_wiggle(double k,double a)
// {
//   static cosmopara C;
//   static double **table_P_Lz = 0;
//   static double logkmin = 0., logkmax = 0., dk = 0., da = 0.;
  
//   double om_m,om_v,amp,ampsqr,grow0,aa,klog,val;
  
//   //      printf("plin test a=%le k=%le\n",a,k_NL);
//   int i,j;
//   if (a >= 0.99999){a =0.99999;}
//   if (recompute_cosmo3D(C))
//   {
//     update_cosmopara(&C);
//     if (table_P_Lz!=0) free_double_matrix(table_P_Lz,0, Ntable.N_a-1, 0, Ntable.N_k_lin-1);
//     table_P_Lz = create_double_matrix(0, Ntable.N_a-1, 0, Ntable.N_k_lin-1);
//     grow0=growfac(1.);
//     da = (1. - limits.a_min)/(Ntable.N_a-1.);
//     aa = limits.a_min;
//     for (i=0; i<Ntable.N_a; i++, aa +=da) {
//       if(aa>1.0) aa=1.0;
//       omega_a(aa,&om_m,&om_v);
//       amp=growfac(aa)/grow0;
//       ampsqr=amp*amp;
      
//       logkmin = log(limits.k_min_mpc);
//       logkmax = log(limits.k_max_mpc);
//       dk = (logkmax - logkmin)/(Ntable.N_k_lin-1.);
//       klog = logkmin;
//       for (j=0; j<Ntable.N_k_lin; j++, klog += dk) {
//         table_P_Lz[i][j] = log(ampsqr*Delta_L_wiggle(exp(klog)));
//       }
//     }
//   }
//   klog = log(k);
//   val = interpol2d(table_P_Lz, Ntable.N_a, limits.a_min, 1., da, a, Ntable.N_k_lin, logkmin, logkmax, dk, klog, cosmology.n_spec, 0.0);
//   return exp(val);     
// }
void free_class_structs(               
               struct background *ba,
               struct thermo *th,
               struct perturbs *pt,
               struct transfers *tr,
               struct primordial *pm,
               struct spectra *sp,
               struct nonlinear *nl,
               struct lensing *le){
  if (lensing_free(le) == _FAILURE_) {
    printf("\n\nError in lensing_free \n=>%s\n",le->error_message);
  }

  if (spectra_free(sp) == _FAILURE_) {
    printf("\n\nError in spectra_free \n=>%s\n",sp->error_message);
  }

  if (transfer_free(tr) == _FAILURE_) {
    printf("\n\nError in transfer_free \n=>%s\n",tr->error_message);
  }

  if (nonlinear_free(nl) == _FAILURE_) {
    printf("\n\nError in nonlinear_free \n=>%s\n",nl->error_message);
  }

  if (primordial_free(pm) == _FAILURE_) {
    printf("\n\nError in primordial_free \n=>%s\n",pm->error_message);
  }

  if (perturb_free(pt) == _FAILURE_) {
    printf("\n\nError in perturb_free \n=>%s\n",pt->error_message);
  }

  if (thermodynamics_free(th) == _FAILURE_) {
    printf("\n\nError in thermodynamics_free \n=>%s\n",th->error_message);
  }

  if (background_free(ba) == _FAILURE_) {
    printf("\n\nError in background_free \n=>%s\n",ba->error_message);
  }
}


int run_class(
			   struct file_content *fc,
	           struct background *ba,
               struct thermo *th,
               struct perturbs *pt,
               struct transfers *tr,
               struct primordial *pm,
               struct spectra *sp,
               struct nonlinear *nl,
               struct lensing *le){
  struct precision pr;        // for precision parameters 
  struct output op;           /* for output files */
  ErrorMsg errmsg; // for error messages 
  
  if(input_init(fc,&pr,ba,th,pt,tr,pm,sp,nl,le,&op,errmsg) == _FAILURE_) {
    fprintf(stderr,"cosmo3D.c: Error running CLASS input:%s\n",errmsg);
    parser_free(fc);
    return CLASS_ERROR_INPUT;
  }
  if (background_init(&pr,ba) == _FAILURE_) {
    
    if(strstr(ba->error_message, "instability") != NULL){
      // Check for instability errors from Hi_CLASS
      printf("%s\n", ba->error_message);
      printf("Horndeski instability; model is not viable.\n");
      return CLASS_ERROR_HORNDESKI_STABILITY;
    }else{
      fprintf(stderr,"cosmo3D.c: Error running CLASS background:%s\n",ba->error_message);
    }
    return CLASS_ERROR_BACKGROUND;
  }
  if (thermodynamics_init(&pr,ba,th) == _FAILURE_) {
    fprintf(stderr,"cosmo3D.c: Error running CLASS thermodynamics:%s\n",th->error_message);
    background_free(ba);  
    return CLASS_ERROR_THERMODYNAMICS;
  }
  if (perturb_init(&pr,ba,th,pt) == _FAILURE_) {
    fprintf(stderr,"cosmo3D.c: Error running CLASS perturb:%s\n",pt->error_message);
    thermodynamics_free(th);
    background_free(ba);  
    return CLASS_ERROR_PERTURB;
  }
  if (primordial_init(&pr,pt,pm) == _FAILURE_) {
    fprintf(stderr,"cosmo3D.c: Error running CLASS primordial:%s\n",pm->error_message);
    perturb_free(pt);
    thermodynamics_free(th);
    background_free(ba);  
    return CLASS_ERROR_PRIMORDIAL;
  }

  if (nonlinear_init(&pr,ba,th,pt,pm,nl) == _FAILURE_) {
    fprintf(stderr,"cosmo3D.c: Error running CLASS nonlinear:%s\n",nl->error_message);
    primordial_free(pm);
    perturb_free(pt);
    thermodynamics_free(th);
    background_free(ba);  
    return CLASS_ERROR_NONLINEAR;
  }

  if (transfer_init(&pr,ba,th,pt,nl,tr) == _FAILURE_) {
    fprintf(stderr,"cosmo3D.c: Error running CLASS transfer:%s\n",tr->error_message);
    nonlinear_free(nl);
    primordial_free(pm);
    perturb_free(pt);
    thermodynamics_free(th);
    background_free(ba);  
    return CLASS_ERROR_TRANSFER;
  }
  if (spectra_init(&pr,ba,pt,pm,nl,tr,sp) == _FAILURE_) {
    fprintf(stderr,"cosmo3D.c: Error running CLASS spectra:%s\n",sp->error_message);
    transfer_free(tr);
    nonlinear_free(nl);
    primordial_free(pm);
    perturb_free(pt);
    thermodynamics_free(th);
    background_free(ba);  
    return CLASS_ERROR_SPECTRA;
  }
  
  // We can assume that the model is stable, as CLASS could compute it
  cosmology.model_stability_flag = 0;
  return CLASS_SUCCESS;
}


double get_class_s8(struct file_content *fc, int *status){
//structures for class test run
    struct background ba;       // for cosmological background 
    struct thermo th;           // for thermodynamics 
    struct perturbs pt;         // for source functions 
    struct transfers tr;        // for transfer functions 
    struct primordial pm;       // for primordial spectra 
    struct spectra sp;          // for output spectra 
    struct nonlinear nl;        // for non-linear spectra 
    struct lensing le;

  //temporarily overwrite P_k_max_1/Mpc to speed up sigma_8 calculation
  double k_max_old = 0.;
  int position_kmax =2;
  double A_s_guess;
  strcpy(fc->name[1],"non linear");
  strcpy(fc->value[1],"none");
  if (strcmp(fc->name[position_kmax],"P_k_max_1/Mpc")){
    k_max_old = strtof(fc->value[position_kmax],NULL);
    sprintf(fc->value[position_kmax],"%e",10.);  
  }
  *status = run_class(fc,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le);
  if (*status == CLASS_SUCCESS){
      free_class_structs(&ba,&th,&pt,&tr,&pm,&sp,&nl,&le);
  }else{
      fprintf(stderr, "get_class_s8: CLASS failed.\n");
  }
  if (k_max_old >0){
    sprintf(fc->value[position_kmax],"%e",k_max_old);      
  }
  return sp.sigma8;
}

double get_class_As(struct file_content *fc, int position_As,double sigma8, int *status){
//structures for class test run
    struct background ba;       // for cosmological background 
    struct thermo th;           // for thermodynamics 
    struct perturbs pt;         // for source functions 
    struct transfers tr;        // for transfer functions 
    struct primordial pm;       // for primordial spectra 
    struct spectra sp;          // for output spectra 
    struct nonlinear nl;        // for non-linear spectra 
    struct lensing le;

  //temporarily overwrite P_k_max_1/Mpc to speed up sigma_8 calculation
  double k_max_old = 0.;
  int position_kmax =2;
  double A_s_guess;
  strcpy(fc->name[1],"non linear");
  strcpy(fc->value[1],"none");
  if (strcmp(fc->name[position_kmax],"P_k_max_1/Mpc")){
    k_max_old = strtof(fc->value[position_kmax],NULL);
    sprintf(fc->value[position_kmax],"%e",10.);  
  }
  A_s_guess = 2.43e-9*pow(sigma8/0.87659,2.0);
  sprintf(fc->value[position_As],"%e",A_s_guess);

  *status = run_class(fc,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le);
  A_s_guess*=pow(sigma8/sp.sigma8,2.);
  if (*status == CLASS_SUCCESS){
      free_class_structs(&ba,&th,&pt,&tr,&pm,&sp,&nl,&le);
  }else{
      fprintf(stderr, "get_class_As: CLASS failed.\n");
  }

  if (k_max_old >0){
    sprintf(fc->value[position_kmax],"%e",k_max_old);      
  }
  return A_s_guess;
}

// Function to fill CLASS parameters struct with chosen cosmological parameters
int fill_class_parameters(struct file_content * fc, int parser_length, int nonlinear){
  int status = 0;
 // basic CLASS configuration parameters
  strcpy(fc->name[0],"output");
  strcpy(fc->value[0],"mPk, tCl,lCl");

  strcpy(fc->name[2],"P_k_max_h/Mpc");
  //higher k_max makes CLASS very slow!
  sprintf(fc->value[2],"%e",limits.k_max_mpc_class);//limits.k_max_mpc/10.);

  strcpy(fc->name[3],"z_max_pk");
  sprintf(fc->value[3],"%e",1./limits.a_min-1.);

  strcpy(fc->name[4],"modes");
  strcpy(fc->value[4],"s");

  strcpy(fc->name[5],"lensing");
  strcpy(fc->value[5],"no");

  // now, copy over cosmology parameters
  strcpy(fc->name[6],"h");
  sprintf(fc->value[6],"%e",cosmology.h0);

  strcpy(fc->name[7],"Omega_cdm");
  sprintf(fc->value[7],"%e",cosmology.Omega_m-cosmology.Omega_nu-cosmology.omb);

  strcpy(fc->name[8],"Omega_b");
  sprintf(fc->value[8],"%e",cosmology.omb);


  strcpy(fc->name[10],"n_s");
  sprintf(fc->value[10],"%e",cosmology.n_spec);

//cosmological constant?
// set Omega_Lambda = 0.0 if w !=-1
  if ((cosmology.w0 !=-1.0) || (cosmology.wa !=0)){
    // FIXME: Figure out what do do with this block!
    strcpy(fc->name[11],"Omega_Lambda");
    sprintf(fc->value[11],"%e", 0.0);

    strcpy(fc->name[12],"w0_fld");
    sprintf(fc->value[12],"%e", cosmology.w0);

    strcpy(fc->name[13],"wa_fld");
    sprintf(fc->value[13],"%e", cosmology.wa);
  }
// pass neutrino parameters
  if (cosmology.M_nu > 1.e-5 || cosmology.Omega_nu >0.){
    strcpy(fc->name[14],"N_ncdm");
    sprintf(fc->value[14],"%d",1);

    if (cosmology.Omega_nu >0.)
    {
      strcpy(fc->name[15],"Omega_ncdm"); 
      sprintf(fc->value[15],"%e", cosmology.Omega_nu);
    }    
    else{
      strcpy(fc->name[15],"m_ncdm"); //\Sigma(m_nu) in eV
      sprintf(fc->value[15],"%e", cosmology.M_nu);
    }
    strcpy(fc->name[16],"N_ur");
    sprintf(fc->value[16],"%e", 2.0328);
  }
  
  // Pass Horndeski parameters
  int is_gr = (   (cosmology.mg_alpha_xk == 0.) 
               && (cosmology.mg_alpha_xb == 0.) 
               && (cosmology.mg_alpha_xm == 0.)
               && (cosmology.mg_alpha_xt == 0.)
               && (cosmology.mg_alpha_M2 == 1.)
              );
  // FIXME: KP if ((cosmology.use_horndeski == 0) && (is_gr == 0)){
  //   cosmology.use_horndeski = 1;
  // }

  if ((cosmology.use_horndeski == 1) && (is_gr == 0))
  {
    // Use the \alpha_X(a) \propto \Omega_DE(a) redshift parametrisation
    strcpy(fc->name[20], "gravity_model");
    strcpy(fc->value[20], "propto_omega");
    
    // Set the alpha parameters using a list; 
    // KP - it seems like a lot of
    // people are against using sprintf (say it's dangerous)? 
    // Should we change this?
    strcpy(fc->name[21], "parameters_smg");
    sprintf(fc->value[21], "%e, %e, %e, %e, %e", cosmology.mg_alpha_xk,
                                                 cosmology.mg_alpha_xb,
                                                 cosmology.mg_alpha_xm,
                                                 cosmology.mg_alpha_xt,
                                                 cosmology.mg_alpha_M2 );
    
    // Initial conditions of scalar field
    strcpy(fc->name[22], "pert_initial_conditions_smg");
    strcpy(fc->value[22], "single_clock");
    
    // Make sure stability tests are switched on; KP - do you mean off here?
    strcpy(fc->name[23], "skip_stability_tests_smg");
    strcpy(fc->value[23], "no");
    
    // Use a w0-wa expansion history
    strcpy(fc->name[24], "expansion_model");
    strcpy(fc->value[24], "wowa");
    
    // Set the expansion parameters (Omega_smg, w0, wa).
    strcpy(fc->name[25], "expansion_smg");
    sprintf(fc->value[25], "%e, %e, %e", 0.5, cosmology.w0, cosmology.wa);
    // KP: Set Omega_smg to some small number and it is overwritten in the next step.

    // Value of Omega_DE today
    strcpy(fc->name[26], "Omega_smg");
    sprintf(fc->value[26], "%e", cosmology.Omega_v);
    // KP: Setting this in this way allows OmegaK to fuflfill the closure relation.


    
  } // end Horndeski section
  else if((cosmology.use_horndeski == 1) && (is_gr == 1)){
    cosmology.use_horndeski = 0;
  }
  
  // Normalization comes last, so that all other parameters are filled in for 
  // determining A_s if sigma_8 is specified
  if (cosmology.A_s){
//     printf("passing A_s=%e directly\n",cosmology.A_s);
     strcpy(fc->name[parser_length-1],"A_s");
     sprintf(fc->value[parser_length-1],"%e",cosmology.A_s);
  }
  else{
    double A_s = get_class_As(fc,parser_length-1,cosmology.sigma_8, &status);
    strcpy(fc->name[parser_length-1],"A_s");
    sprintf(fc->value[parser_length-1],"%e",A_s);
    if (status == 0){
  	  A_s *=pow(cosmology.sigma_8/get_class_s8(fc,&status),2.0);
	    strcpy(fc->name[parser_length-1],"A_s");
  	  sprintf(fc->value[parser_length-1],"%e",A_s);
    }
    //printf("determined A_s(sigma_8=%e) = %e\n", cosmology.sigma_8,A_s);
  }
  
  // Turn off Halofit if Horndeski is enabled, or if nonlinear mode is turned off
  strcpy(fc->name[1], "non linear");
  // if ( (cosmology.use_horndeski == 0) || (nonlinear == 0) ){
  //   strcpy(fc->value[1], "Halofit");
  // }else{
  //   strcpy(fc->value[1], "");
  // }
  strcpy(fc->value[1], ""); //FIX ME: printing halofit when in horndeski?
  
  
  
  // FIXME: Output all parameters that we are passing to CLASS
  FILE *fileforclass;
  fileforclass = fopen("/Users/kpardo/Documents/hi_class_public/classparams.ini","w");
  for(int i=0; i < 27; i++){
    if (fc->read[i] == 1){
      fprintf(fileforclass,"%s = %s\n", fc->name[i], fc->value[i]);
    }
  }
  fprintf(fileforclass, "root = output/classparams00_");
  
  return status;
}


double evaluate_class_tables(double k_coverh0, double a, int NL, int mode, int *status){
  
  static cosmopara C;
  static double **table_P_L = 0;
  static double **table_P_NL = 0;
  
  static double *table_horndeski_cs2 = 0;
  static double *table_horndeski_alpha_b = 0;
  static double *table_horndeski_alpha_k = 0;
  static double *table_horndeski_alpha_m = 0;
  static double *table_horndeski_alpha_t = 0;
  static double *table_horndeski_mu_eff = 0;
  static double *table_horndeski_mu_light = 0;
  
  static double logkmin = 0., logkmax = 0., dk = 0., da = 0.;
  static int class_status = 0;
  double val, klog;
  
  if (recompute_cosmo3D(C)){
    //printf("+++ Recomputing; parameters changed! +++\n");
    
    update_cosmopara(&C);
    cosmology.model_stability_flag = -1; // We don't know if the model is stable yet
    
    // Allocate memory for power spectrum tables
    if (table_P_L == 0){
      table_P_L = create_double_matrix(0, Ntable.N_a-1, 0, Ntable.N_k_nlin-1); 
      if (NL == 1){
        table_P_NL = create_double_matrix(0, Ntable.N_a-1, 0, Ntable.N_k_nlin-1); 
      }
      da = (1. - limits.a_min)/(Ntable.N_a-1.);
      logkmin = log(limits.k_min_mpc*cosmology.coverH0);
      logkmax = log(limits.k_max_mpc_class*cosmology.coverH0);
      dk = (logkmax-logkmin)/(Ntable.N_k_nlin-1.);
    }
    
    // Allocate CLASS structures
    struct background ba;       // for cosmological background 
    struct thermo th;           // for thermodynamics 
    struct perturbs pt;         // for source functions 
    struct transfers tr;        // for transfer functions 
    struct primordial pm;       // for primordial spectra 
    struct spectra sp;          // for output spectra 
    struct nonlinear nl;        // for non-linear spectra 
    struct lensing le;
    struct output op;

  	ErrorMsg errmsg; // for error messages 
    
    // Initialise parser for CLASS parameters
  	struct file_content fc;
  	int parser_length = 30;
  	if (parser_init(&fc,parser_length,"none",errmsg) == _FAILURE_){
    	fprintf(stderr,"cosmo3D.c: CLASS parser init error:%s\n",errmsg);
    	*status = 1;
    	return CLASS_ERROR_PARSER;
  	}
  	for (int i=0; i < parser_length; i++){
    	strcpy(fc.name[i], " ");
    	strcpy(fc.value[i], " ");
  	}
    
    // Collect input parameters, then run CLASS
  	*status = fill_class_parameters(&fc, parser_length, NL);
  	if(*status>0) return 1; 
  	*status = run_class(&fc,&ba,&th,&pt,&tr,&pm,&sp,&nl,&le);
  	parser_free(&fc);
  	if(*status != CLASS_SUCCESS){
      if (*status == CLASS_ERROR_HORNDESKI_STABILITY) cosmology.model_stability_flag = 0; // model is unstable
      return 1;
    }
  	
  	// Determine normalisation of power spectrum
    double aa, norm, k_class,Pk,ic;
    int i,j,s;
    aa = limits.a_min;
    if (cosmology.A_s){
      norm = 3.*log(cosmology.h0/cosmology.coverH0);
      cosmology.sigma_8 = sp.sigma8;
    }
    else{
      norm = log(pow(cosmology.sigma_8/sp.sigma8,2.)*pow(cosmology.h0/cosmology.coverH0,3.));
    }
    
    // Construct interpolation table for power spectrum
    //printf("power spectrum scaling factor %e\n", pow(cosmology.sigma_8/sp.sigma8,2.));
    if (*status == 0){
      for (i=0; i<Ntable.N_a; i++, aa +=da) { 
        klog = logkmin;
        for (j=0; j<Ntable.N_k_nlin; j++, klog += dk) { 
          k_class = exp(klog)*cosmology.h0/cosmology.coverH0;
          s = spectra_pk_at_k_and_z(&ba, &pm, &sp,k_class,fmax(1./aa-1.,0.), &Pk,&ic);
          table_P_L[i][j] = log(Pk) +norm;
          if (NL == 1){
            s = spectra_pk_nl_at_k_and_z(&ba, &pm, &sp,k_class,fmax(1./aa-1.,0.), &Pk);
            table_P_NL[i][j] = log(Pk) +norm;
          }
        }
      }
    } // end status check
    
    // Construct Horndeski interpolation tables (if requested)
    if (cosmology.use_horndeski == 1){
        double z = 0.;
        double tau = 0.;
        double fac1, fac2;
        int last_index = 0;
        double* pvecback = create_double_vector(0, ba.bg_size-1);
        double alph, bxi, bB, alphaB, alphaK, alphaM, alphaT, cs2;
        
        // Allocate memory for Horndeski function tables
        // (uses same sampling in 'a' as the power spectrum)
        table_horndeski_cs2 = create_double_vector(0, Ntable.N_a-1);
        table_horndeski_alpha_b = create_double_vector(0, Ntable.N_a-1);
        table_horndeski_alpha_k = create_double_vector(0, Ntable.N_a-1);
        table_horndeski_alpha_m = create_double_vector(0, Ntable.N_a-1);
        table_horndeski_alpha_t = create_double_vector(0, Ntable.N_a-1);
        table_horndeski_mu_eff = create_double_vector(0, Ntable.N_a-1);
        table_horndeski_mu_light = create_double_vector(0, Ntable.N_a-1);
        
        // Populate tables with Horndeski function values
        aa = limits.a_min;
        for (i=0; i<Ntable.N_a; i++, aa+=da) {
            
            // Get values of all background parameters at this scale factor
            *status = background_tau_of_z(&ba, fmax(1./aa-1.,0.), &tau);
            if (*status != 0){
                printf("ERROR: Failed to get tau(z).\n"); exit(1); 
            }
            *status = background_at_tau(&ba, tau, ba.long_info, ba.inter_normal, &last_index, pvecback);
            if (*status != 0){
                printf("ERROR: Failed to get background values at z.\n"); exit(1);
            }
            
            // Fetch Horndeski quantities
            cs2 = pvecback[ba.index_bg_cs2_smg];
            alphaB = pvecback[ba.index_bg_braiding_smg];
            alphaK = pvecback[ba.index_bg_kineticity_smg];
            alphaM = 0.; // FIXME: Shouldn't be set to zero by default
            alphaT = pvecback[ba.index_bg_tensor_excess_smg];
            
            alph = alphaK + 6.*alphaB*alphaB;
            fac1 = 2. / alph / cs2; // FIXME: Can be inf
            fac2 = alphaB*(1. + alphaT) + alphaT - alphaM;
            
            // Check for nan/inf
            if (isnan(fac1) || isinf(fac1)) fac1 = 0.;
            if (isnan(fac2) || isinf(fac2)) fac2 = 0.;
            
            //printf("Horndeski params: %3.3e // %3.3e, %3.3e, %3.3e, %3.3e // %3.3e, %3.3e, %3.3e\n", cs2, alphaB, alphaK, alphaM, alphaT, alph, fac1, fac2);
            
            // Populate vectors of Horndeski quantities
            table_horndeski_cs2[i] = cs2;
            table_horndeski_alpha_b[i] = alphaB;
            table_horndeski_alpha_k[i] = alphaK;
            table_horndeski_alpha_m[i] = alphaM;
            table_horndeski_alpha_t[i] = alphaT;
            
            // FIXME: Check for consistency with Hi_CLASS notation!
            // mu_eff (Eq. 3 of arXiv:1705.04714)
            table_horndeski_mu_eff[i] = 1. + alphaT + fac1*fac2*fac2;
            
            // mu_light (Eq. 7 of arXiv:1705.04714)
            table_horndeski_mu_light[i] = 1. + 0.5*alphaT + fac1*fac2*(alphaB + fac2);
            
            
        } // end loop over a
        
        // Free temporary vector
        free_double_vector(pvecback, 0, ba.bg_size-1);
        
    } // end Horndeski check
    
    // Free CLASS structures
    if (*status == 0) free_class_structs(&ba,&th,&pt,&tr,&pm,&sp,&nl,&le);
    
  } // end recompute
  
  
  // Return power spectrum if requested
  if (mode == CLASS_RETURN_POWSPEC){
    klog = log(k_coverh0);
    if (isnan(klog) || class_status) return 0.0;
    if (NL==1) val = interpol2d_fitslope(table_P_NL, Ntable.N_a, limits.a_min, 1., da, fmin(a,.99), Ntable.N_k_nlin, logkmin, logkmax, dk, klog, cosmology.n_spec);
    else val = interpol2d_fitslope(table_P_L, Ntable.N_a, limits.a_min, 1., da, fmin(a,.99), Ntable.N_k_nlin, logkmin, logkmax, dk, klog, cosmology.n_spec);
    if(isnan(val)) return 0.0;
    return exp(val);
  }
  
  // Return Horndeski-related functions if requested
  if (cosmology.use_horndeski == 1){
            
    // Identify which Horndeski function has been requested 
    double* func;
    switch (mode){
      case CLASS_RETURN_HORNDESKI_ALPHA_B: func = table_horndeski_alpha_b; break;
      case CLASS_RETURN_HORNDESKI_ALPHA_K: func = table_horndeski_alpha_k; break;
      case CLASS_RETURN_HORNDESKI_ALPHA_M: func = table_horndeski_alpha_m; break;
      case CLASS_RETURN_HORNDESKI_ALPHA_T: func = table_horndeski_alpha_t; break;
      case CLASS_RETURN_HORNDESKI_CS2: func = table_horndeski_cs2; break;
      case CLASS_RETURN_HORNDESKI_MU_EFF: func = table_horndeski_mu_eff; break;
      case CLASS_RETURN_HORNDESKI_MU_LIGHT: func = table_horndeski_mu_light; break;
      default:
        fprintf(stderr, "cosmo3D.c: Horndeski parameter with ID '%d' not found.\n", mode);
        *status = 1;
        return 0.;
    }
    
    // Interpolate value of this function and return
    /* Interpolates f at the value x, where f is a double[n] array,	*
     * representing a function between a and b, stepwidth dx.	*
     * 'lower' and 'upper' are powers of a logarithmic power law	*
     * extrapolation. If no	extrapolation desired, set these to 0	*/
    val = interpol(func, Ntable.N_a, limits.a_min, 1., da, a, 0, 0);
    //if(isnan(val)) return 0.0;
    return val;
  
  } // end Horndeski check
  
}

// Evaluate Horndeski function at a given scale factor (e.g. alpha_X(a), cs^2(a))
double horndeski_function(double a, char* func_name, int *status){
    int mode = -1;
    int NL = 0; // Non-linear mode (must be zero, i.e. NL switched off)
    double k_dummy = 0.; // Dummy k value (will not be used)
    
    // Identify which Horndeski function is being requested
    if (strcmp(func_name, "alpha_b")) mode = CLASS_RETURN_HORNDESKI_ALPHA_B;
    if (strcmp(func_name, "alpha_k")) mode = CLASS_RETURN_HORNDESKI_ALPHA_K;
    if (strcmp(func_name, "alpha_m")) mode = CLASS_RETURN_HORNDESKI_ALPHA_M;
    if (strcmp(func_name, "alpha_t")) mode = CLASS_RETURN_HORNDESKI_ALPHA_T;
    if (strcmp(func_name, "cs2")) mode = CLASS_RETURN_HORNDESKI_CS2;
    if (strcmp(func_name, "mu_eff")) mode = CLASS_RETURN_HORNDESKI_MU_EFF;
    if (strcmp(func_name, "mu_light")) mode = CLASS_RETURN_HORNDESKI_MU_LIGHT;
    
    // Safety check
    if (mode < 0){ 
        fprintf(stderr, 
            "WARNING: cosmo3D.c:horndeski_function, mode '%d' not recognized.\n", 
            mode);
        *status = 1;
        return 0.;
    }
    
    // Evaluate function from table returned by CLASS
    return evaluate_class_tables(k_dummy, a, NL, mode, status);
}

// Evaluate power spectrum from CLASS
double p_class(double k_coverh0, double a, int NL, int *status){
    double pk = 0.;
    pk = evaluate_class_tables(k_coverh0, a, NL, CLASS_RETURN_POWSPEC, status);
    return pk;
}

// Linear power spectrum routine with k in units H_0/c; used in covariances.c for beat coupling and in halo.c
double p_lin(double k, double a)
{
  
  // Use CLASS instead if Horndeski is being used
  if ( cosmology.use_horndeski != 0 ){
      fprintf(stderr, "WARNING: A function called p_lin(), but this should not be used in 'horndeski' mode.\n");
      int status = 0;
      double pk = p_class(k, a, 0, &status); // Switch off non-linear stuff
      return pk;
  }
  
  static cosmopara C;
  static double **table_P_Lz = 0;
  static double logkmin = 0., logkmax = 0., dk = 0., da = 0.;
  int status;
  if (strcmp(pdeltaparams.runmode,"CLASS")==0 || strcmp(pdeltaparams.runmode,"class")==0) return p_class(k,a,0, &status);
  
  double amp,ampsqr,grow0,aa,klog,val;
  
  int i,j;
  if (a >= 0.99999){a =0.99999;}
  if (recompute_cosmo3D(C)){
    update_cosmopara(&C);
    if (table_P_Lz!=0) free_double_matrix(table_P_Lz,0, Ntable.N_a-1, 0, Ntable.N_k_lin-1);
    table_P_Lz = create_double_matrix(0, Ntable.N_a-1, 0, Ntable.N_k_lin-1);
    grow0=growfac(1.);
    da = (1. - limits.a_min)/(Ntable.N_a-1.);
    aa = limits.a_min;
    for (i=0; i<Ntable.N_a; i++, aa +=da) {
      if(aa>1.0) aa=1.0;
      amp=growfac(aa)/grow0;
      ampsqr=amp*amp;
      
      logkmin = log(limits.k_min_mpc);
      logkmax = log(limits.k_max_mpc);
      dk = (logkmax - logkmin)/(Ntable.N_k_lin-1.);
      klog = logkmin;
      for (j=0; j<Ntable.N_k_lin; j++, klog += dk) {
        table_P_Lz[i][j] = log(ampsqr*Delta_L_wiggle(exp(klog)));
        //printf("%le %le",exp(klog),Delta_L_wiggle(exp(klog));
      }
    }
  }
  if (k/cosmology.coverH0< exp(logkmin) || k/cosmology.coverH0 > exp(logkmax)) return 0.0;
  klog = log(k/cosmology.coverH0);
  val = interpol2d(table_P_Lz, Ntable.N_a, limits.a_min, 1., da, a, Ntable.N_k_lin, logkmin, logkmax, dk, klog, 1.0, 1.0);
  if(isnan(val) || (k==0)) return 0.0;
  return 2.0*constants.pi_sqr*exp(val)/k/k/k;
}


double int_sig_R_knl(double lnk, void *args) // tak12 A4
{
  double krsqr;
  
  double *params= (double *) args;
  double Rscale=params[0];
  //printf("Rscale %le k %le\n",Rscale,exp(lnk));
  krsqr= SQR(exp(lnk)*Rscale);
  return Delta_L_wiggle(exp(lnk))*exp(-krsqr);
}


double int_neff(double lnk, void *args) //tak12 A5
{
  double krsqr;
  double *params= (double *) args;
  double Rscale=params[0];
  krsqr= SQR(exp(lnk)*Rscale);
  return Delta_L_wiggle(exp(lnk))*2.0*krsqr*exp(-krsqr); //see S03 eq. 59
}


double int_cur(double lnk, void *args) //tak12 A5
{
  double krsqr;
  double *params= (double *) args;
  double Rscale=params[0];
  krsqr= SQR(exp(lnk)*Rscale);
  return Delta_L_wiggle(exp(lnk))*4.0*krsqr*(1.0-krsqr)*exp(-krsqr); // S03 eq.60
}


//iterative calculation of the nonlinear scale as defined in tak12 A4
void nonlin_scale(double amp, double *R_NL, double *neff, double *Curv)
{
  double sig_R,kmax,logkmax,sig_R_noamp,neffplus3;
  int iterstep; 
  const int itermax  = 40;
  int converged=0;
  double array[1];
  double logRmin = -3.0;
  double logRmax =  4.0; 
   
  iterstep=0;  
  while(converged==0)
  {        
    array[0]=pow(10.,(logRmin+logRmax)/2.0);
    
    //flexible upper limit of integration depending on drop-off of filter function
    kmax  = sqrt(5.*log(10.))/array[0];
    if (kmax<8000.0) logkmax = log(8000.0);
    sig_R_noamp=sqrt(int_gsl_integrate_medium_precision(int_sig_R_knl,(void*)array,-4.5,logkmax,NULL,512)); //integral goes over ln k exponent correspond to k_min~0.011
    
    sig_R=amp*sig_R_noamp; 
    if (sig_R>1.0)  logRmin=log10(array[0]);
    if (sig_R<1.0)  logRmax=log10(array[0]);
    iterstep=iterstep+1;
    if(fabs(sig_R-1.0) < 0.0001 || iterstep>itermax) converged=1;    
  }  
  *R_NL=array[0]; //R where sig_R==1
  neffplus3=int_gsl_integrate_medium_precision(int_neff,(void*)array,-4.5,logkmax,NULL,512)/sig_R_noamp/sig_R_noamp;
  *neff= neffplus3 - 3.0;
  *Curv= int_gsl_integrate_medium_precision(int_cur,(void*)array,-4.5,logkmax,NULL,512)/sig_R_noamp/sig_R_noamp + SQR(neffplus3);
  
  //printf("%d %le\n",iterstep,amp);  
}                  

 
double Halofit(double k, double amp, double omm, double omv,double w_z, double R_NL, double neff,double Curv, double P_delta_Lin)
{
  double y_scale,n2eff,n3eff,n4eff;
  double a_n,b_n,c_n,gamma_n,alpha_n,beta_n,nu_n,f1,f2,f3;
  
  double Delta_H,Delta_H_Prime,Delta_Q;
  
  // This function is not Horndeski-safe
  invalid_for_horndeski(__func__);
  
  //determine nonlinear scale, neff and curvature, see tak12 A4, A5
  y_scale=k*R_NL;
  
  n2eff=neff*neff;
  n3eff=n2eff*neff;
  n4eff=n2eff*n2eff;
  
  //calculate coefficients 
  a_n = pow(10.,1.5222+2.8553*neff + 2.3706*n2eff+0.9903*n3eff+0.2250*n4eff-0.6038*Curv+0.1749*omv*(1.0+w_z));
  b_n = pow(10., -0.5642+0.5864*neff + 0.5716*n2eff-1.5474*Curv +0.2279*omv*(1.0+w_z));
  c_n = pow(10., 0.3698+ 2.0404*neff + 0.8161*n2eff+0.5869*Curv);
  gamma_n = 0.1971-0.0843*neff + 0.8460*Curv;
  alpha_n = fabs(6.0835 + 1.3373*neff - 0.1959*n2eff - 5.5274*Curv);
  beta_n = 2.0379 - 0.7354*neff + 0.3157*n2eff + 1.2490*n3eff + 0.3980*n4eff - 0.1682*Curv;
  nu_n = pow(10,5.2105+3.6902*neff);
  
  f1 = pow(omm,(-0.0307));
  f2 = pow(omm,(-0.0585));
  f3 = pow(omm,(0.0743));  
  
  //TwoHaloTerm
  Delta_Q=P_delta_Lin*(pow((1.0+P_delta_Lin),beta_n)/(1.0+alpha_n*P_delta_Lin))*exp(-(y_scale/4.0+y_scale*y_scale/8.0));
  //OneHaloterm
  Delta_H_Prime=(a_n*pow(y_scale,3.0*f1))/(1.0+b_n*pow(y_scale,f2)+pow(c_n*f3*y_scale,3.0-gamma_n));
  Delta_H=Delta_H_Prime/(1.0+nu_n*pow(y_scale,-2.0)); // using mu=0.0 Tak A12
  //printf("Delta_Q %le Delta_H %le\n",Delta_Q,Delta_H);
 return Delta_H+Delta_Q;
}


void Delta_halofit(double **table_P_NL,double logkmin, double logkmax, double dk, double da)
{       
  double rk,omm,omv,w_z,amp,grow0,aa,klog;
  double R_NL,Curv,neff,P_delta,P_delta_Lin;
  int i,j;
  
  grow0=growfac(1.);
  aa = limits.a_min;
  //binning in k and a must be the same as in Coyote
  for (i=0; i<Ntable.N_a; i++, aa +=da) { 
    if(aa>1.0) aa=1.0;
    omega_a(aa,&omm,&omv);
    w_z=cosmology.w0+cosmology.wa*(1.-aa);
    amp=growfac(aa)/grow0;
    nonlin_scale(amp, &R_NL, &neff, &Curv);
    //printf("%le %le %le %le\n",aa,R_NL,neff,Curv);
    klog = logkmin;
    for (j=0; j<Ntable.N_k_nlin; j++, klog += dk) { 
      rk=exp(klog);
      P_delta_Lin=amp*amp*Delta_L_wiggle(rk);
      P_delta=Halofit(rk, amp, omm, omv, w_z, R_NL, neff, Curv, P_delta_Lin);
      table_P_NL[i][j]=log(P_delta);
    }
  }
}


double Delta_NL_Halofit(double k_NL, double a)
{     
  static cosmopara C;
  static double logkmin = 0., logkmax = 0., dk = 0., da = 0.;
  
  static double **table_P_NL=0;
  double klog,val; 

  if (recompute_cosmo3D(C)){
    update_cosmopara(&C);
    if (table_P_NL!=0) free_double_matrix(table_P_NL,0, Ntable.N_a-1, 0, Ntable.N_k_nlin-1);
    table_P_NL = create_double_matrix(0, Ntable.N_a-1, 0, Ntable.N_k_nlin-1);     
    
    da = (1. - limits.a_min)/(Ntable.N_a-1.);
    logkmin = log(limits.k_min_mpc);
    logkmax = log(limits.k_max_mpc);
    dk = (logkmax - logkmin)/(Ntable.N_k_nlin-1.);
    
    Delta_halofit(table_P_NL,logkmin, logkmax, dk, da);
  }
  klog = log(k_NL);
  val = interpol2d(table_P_NL, Ntable.N_a, limits.a_min, 1., da, a, Ntable.N_k_nlin, logkmin, logkmax, dk, klog, cosmology.n_spec, 0.0);
  return exp(val);
  // returns the dimensionless power spectrum as a function of scale factor a and k  
}

double nonlinear_scale_computation(double a)
{       
  static cosmopara C;
  static double da = 0.;
  static double *table=0;
  
  double omm,omv,amp,grow0,aa,res;
  double R_NL,Curv,neff;
  int i;
  
  if (recompute_cosmo3D(C)){
    update_cosmopara(&C);
    grow0=growfac(1.);
    da = (1.-limits.a_min)/(Ntable.N_a-1.);
    aa = limits.a_min;
    if (table!=0) free_double_vector(table, 0, Ntable.N_a-1);
    table   = create_double_vector(0, Ntable.N_a-1);
    for (i=0; i<Ntable.N_a; i++, aa+=da) {
      if(aa>1.0) aa=1.0;
      omega_a(aa,&omm,&omv);
      amp=growfac(aa)/grow0;
      nonlin_scale(amp, &R_NL, &neff, &Curv);
      table[i] = 1./R_NL;
    }
  }
  res = interpol(table, Ntable.N_a, limits.a_min, 1., da, a, 0.0, 0.0); 
  return res;
}


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

void determine_coyote_cosmo_calib(double *COSMO_emu, int *calibflag)
{
  COSMO_emu[3] = cosmology.h0;
  if( cosmology.h0<= h0_min_coyote){  
    COSMO_emu[3]=h0_min_coyote+0.000001;
    *calibflag=1;
  }
  if(cosmology.h0>= h0_max_coyote){  
    COSMO_emu[3]=h0_max_coyote-0.000001;
    *calibflag=1;}
  
  
  COSMO_emu[1] = cosmology.Omega_m*COSMO_emu[3]*COSMO_emu[3];
  COSMO_emu[0] = cosmology.omb*COSMO_emu[3]*COSMO_emu[3];
  COSMO_emu[2] = cosmology.n_spec;
  COSMO_emu[4] = cosmology.w0;
  COSMO_emu[5] = cosmology.sigma_8;  

  if (COSMO_emu[0] <= ombhh_min_coyote){ 
    COSMO_emu[0]= ombhh_min_coyote+0.000001;
    *calibflag=1;
  }
  if (COSMO_emu[0] >= ombhh_max_coyote){ 
    COSMO_emu[0]=ombhh_max_coyote-0.000001;
    *calibflag=1;   
  }
  if(COSMO_emu[1]  <= omhh_min_coyote){  
    COSMO_emu[1]=omhh_min_coyote+0.000001;
    *calibflag=1;
  }
  if(COSMO_emu[1] >= omhh_max_coyote){  
    COSMO_emu[1]=omhh_max_coyote-0.000001;
    *calibflag=1;
  }
  if(cosmology.n_spec <= ns_min_coyote){ 
    COSMO_emu[2]=ns_min_coyote+0.000001;
    *calibflag=1; 
  }
  if(cosmology.n_spec >= ns_max_coyote){  
    COSMO_emu[2]=ns_max_coyote-0.000001;
    *calibflag=1;   
  }
  if( cosmology.w0<= w_min_coyote){  
    COSMO_emu[4]=w_min_coyote+0.000001;
    *calibflag=1; 
  }
  if( cosmology.w0>= w_max_coyote){  
    COSMO_emu[4]=w_max_coyote-0.000001;
    *calibflag=1; 
  }
  if(cosmology.sigma_8<= s8_min_coyote){  
    COSMO_emu[5]=s8_min_coyote+0.000001;
    *calibflag=1; 
  }
  if( cosmology.sigma_8>= s8_max_coyote){  
    COSMO_emu[5]=s8_max_coyote-0.000001;
    *calibflag=1;    
  }
  if(cosmology.wa != 0.0){  
    *calibflag=1;    
  } 
}
  

double Delta_NL_Coyote(double k_NL,double a)
{     
  static cosmopara C;
  static double logkmin = 0., logkmax = 0., dk = 0., da = 0.;
  
  static double **table_P_NL=0;
  static double **table_P_NL_halofit=0;
  static double **table_P_NL_halofit_calibrate=0;
  
  double aa,klog,val; 
  double COSMO_orig[7],COSMO_emu[7],ystar[2*582],kstar[582],p_coyote[582],coyote_min,coyote_max;
  int type=1,calibflag=0;
  int i,j,k;
  
  if (recompute_cosmo3D(C)){
    update_cosmopara(&C);

    if (table_P_NL!=0) free_double_matrix(table_P_NL,0, Ntable.N_a-1, 0, Ntable.N_k_nlin-1);     
    if (table_P_NL_halofit!=0) free_double_matrix(table_P_NL_halofit,0, Ntable.N_a-1, 0,Ntable.N_k_nlin-1);     
    table_P_NL = create_double_matrix(0, Ntable.N_a-1, 0,Ntable.N_k_nlin-1);     
    table_P_NL_halofit = create_double_matrix(0, Ntable.N_a-1, 0,Ntable.N_k_nlin-1);     
    da = (1. - limits.a_min)/(Ntable.N_a-1.);
    logkmin = log(limits.k_min_mpc);
    logkmax = log(limits.k_max_mpc);
    dk = (logkmax - logkmin)/(Ntable.N_k_nlin-1.);
    
    //printf("Starting P_delta %le %le %le %le %le %le %le\n",cosmology.Omega_m,cosmology.omb,cosmology.n_spec, cosmology.sigma_8,cosmology.w0,cosmology.wa,cosmology.h0);
    
    //compute Halofit; determine whether outside cosmology of emulator -> use recalibration factor
    Delta_halofit(table_P_NL_halofit,logkmin, logkmax, dk, da);
    determine_coyote_cosmo_calib(COSMO_emu, &calibflag);    
    
    if(calibflag==0){
      //printf("INSIDE Emulator cosmology\n");
      COSMO_emu[3] =COSMO_emu[3]*100;
      aa = limits.a_min;
      //binning in k and a must be the same as in  Delta_halofit
      for (i=0; i<Ntable.N_a; i++, aa +=da) {
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *timspline = gsl_spline_alloc (gsl_interp_cspline, 582);
	COSMO_emu[6] = (1.0/aa)-1.0; //emu takes 7 args 6 cosmopara and 7th redshift    
	if(fabs(COSMO_emu[6])<1.e-10) {
	  COSMO_emu[6]=0.0;
	}
	if(aa >= a_min_coyote){
	  emu(COSMO_emu,ystar, &type);
	  for (k=0; k<582; k++){
	    kstar[k]=ystar[k];
	    p_coyote[k]=ystar[k+582];
	    //printf("%le %le\n",kstar[k],p_coyote[k]);
	  }
	  gsl_spline_init (timspline, kstar, p_coyote, 582);

	  coyote_min=log(p_coyote[0]/Delta_NL_Halofit(kstar[0]/cosmology.h0,aa));
	  coyote_max=log(p_coyote[581]/Delta_NL_Halofit(kstar[581]/cosmology.h0,aa));
	  //printf("%le %le\n",kstar[k-1],p_coyote[k-1]);
	  klog = logkmin; // log k in h/MPC
	  for (j=0; j<Ntable.N_k_nlin; j++, klog += dk) {
	    if ((klog >= log(k_min_coyote/cosmology.h0)) && (klog <= log(k_max_coyote/cosmology.h0))){
	      table_P_NL[i][j]=log(gsl_spline_eval(timspline, exp(klog)*cosmology.h0, acc));
	      //printf("Coyote used\n");
	    }
	    if(klog>log(k_max_coyote/cosmology.h0))  table_P_NL[i][j]=coyote_max+table_P_NL_halofit[i][j];
	    if(klog<log(k_min_coyote/cosmology.h0)) table_P_NL[i][j]=coyote_min+table_P_NL_halofit[i][j];
	      //printf("Halofit used: exceeded Coyote k range k=%le k_min=%le k_max=%le\n",exp(klog),k_min_coyote/cosmology.h0,k_max_coyote/cosmology.h0);
	    }
	  }
	if(aa < a_min_coyote){
	  //printf("non coyote %le\n",aa);
	  for (j=0; j<Ntable.N_k_nlin; j++) {
	    table_P_NL[i][j]=table_P_NL_halofit[i][j]; // coyote goes down to z=4, Halofit and Coyote difference small since pdelta is almost linear -> no need for calibration here
	  }
	}
	gsl_spline_free (timspline);
	gsl_interp_accel_free (acc);
      }
    }
    if(calibflag==1){      
      //printf("OUTSIDE Emulator cosmology\n");
      // to restore the cosmology structure later 
      COSMO_orig[1] = cosmology.Omega_m;
      COSMO_orig[0] = cosmology.omb;
      COSMO_orig[2] = cosmology.n_spec;
      COSMO_orig[3] = cosmology.h0;
      COSMO_orig[4] = cosmology.w0;
      COSMO_orig[5] = cosmology.sigma_8;  
      COSMO_orig[6] = cosmology.wa;
      
      if (table_P_NL_halofit_calibrate!=0) free_double_matrix(table_P_NL_halofit_calibrate,0, Ntable.N_a-1, 0,Ntable.N_k_nlin-1);     
      table_P_NL_halofit_calibrate = create_double_matrix(0, Ntable.N_a-1, 0,Ntable.N_k_nlin-1);     
      
      // set cosmology to compute the Halofit calibration power spectrum 
      cosmology.Omega_m=COSMO_emu[1]/COSMO_emu[3]/COSMO_emu[3];
      cosmology.Omega_v=1.0-cosmology.Omega_m;
      cosmology.omb=COSMO_emu[0]/COSMO_emu[3]/COSMO_emu[3];
      cosmology.n_spec=COSMO_emu[2];
      cosmology.h0 =COSMO_emu[3];
      cosmology.w0=COSMO_emu[4];
      cosmology.sigma_8=COSMO_emu[5];  
      cosmology.wa=0.0;  
      
//      Delta_halofit(table_P_NL_halofit_calibrate,logkmin, logkmax, dk, da);
           
      COSMO_emu[3] = COSMO_emu[3]*100.0;
        aa = limits.a_min;
      //printf("COSMO %le %le %le %le %le %le\n",COSMO_emu[0],COSMO_emu[1],COSMO_emu[2],COSMO_emu[3],COSMO_emu[4],COSMO_emu[5]);
      for (i=0; i<Ntable.N_a; i++, aa +=da) {
	gsl_interp_accel *acc = gsl_interp_accel_alloc ();
	gsl_spline *timspline = gsl_spline_alloc (gsl_interp_cspline, 582);
	COSMO_emu[6] = (1.0/aa)-1.0; //emu takes 7 args 6 cosmopara and 7th redshift    
	if(fabs(COSMO_emu[6])<1.e-10) {
	  COSMO_emu[6]=0.0;
	}
	if(aa >= a_min_coyote){
	  //printf("COSMO %le %le %le %le %le %le\n",COSMO_emu[0],COSMO_emu[1],COSMO_emu[2],COSMO_emu[3],COSMO_emu[4],COSMO_emu[5]);
	  emu(COSMO_emu,ystar, &type);
	  for (k=0; k<582; k++){
	    kstar[k]=ystar[k];
	    p_coyote[k]=ystar[k+582]/Delta_NL_Halofit(kstar[k]/cosmology.h0,aa);
	    //printf("%le %le\n",kstar[k],p_coyote[k]);
	  }
	  gsl_spline_init (timspline, kstar, p_coyote, 582);

	  coyote_min=p_coyote[0];
	  coyote_max=p_coyote[581];
	  //printf("%le %le\n",kstar[k-1],p_coyote[k-1]);
	  klog = logkmin; // log k in h/MPC
	  for (j=0; j<Ntable.N_k_nlin; j++, klog += dk) {
	    if ((klog >= log(k_min_coyote/cosmology.h0)) && (klog <= log(k_max_coyote/cosmology.h0))){
	      table_P_NL[i][j]=log(gsl_spline_eval(timspline, exp(klog)*cosmology.h0, acc))+table_P_NL_halofit[i][j];
	      //printf("Coyote used\n");
	    }
	    if(klog>log(k_max_coyote/cosmology.h0))  table_P_NL[i][j]=log(coyote_max)+table_P_NL_halofit[i][j];
	    if(klog<log(k_min_coyote/cosmology.h0)) table_P_NL[i][j]=log(coyote_min)+table_P_NL_halofit[i][j];
	      //printf("Halofit used: exceeded Coyote k range k=%le k_min=%le k_max=%le\n",exp(klog),k_min_coyote/cosmology.h0,k_max_coyote/cosmology.h0);
	    }
	  }
	if(aa < a_min_coyote){
	  //printf("non coyote %le\n",aa);
	  for (j=0; j<Ntable.N_k_nlin; j++) {
	    table_P_NL[i][j]=table_P_NL_halofit[i][j]; // coyote goes down to z=4, Halofit and Coyote difference small since pdelta is almost linear -> no need for calibration here
	  }
	}
	gsl_spline_free (timspline);
	gsl_interp_accel_free (acc);
	cosmology.Omega_m=COSMO_orig[1];
	cosmology.Omega_v=1.0-cosmology.Omega_m;
	cosmology.omb=COSMO_orig[0];
	cosmology.n_spec=COSMO_orig[2];
	cosmology.h0 =COSMO_orig[3];
	cosmology.w0=COSMO_orig[4];
	cosmology.sigma_8=COSMO_orig[5];  	
	cosmology.wa=COSMO_orig[6];
      }      
    }
  }
//  printf("%le\n",k_NL);
  klog = log(k_NL);
//  if(a < a_min_coyote || klog>log(k_max_coyote/cosmology.h0) || klog<log(k_min_coyote/cosmology.h0))printf("Halofit used: exceeded Coyote a or k range a=%le k=%le\n",a,exp(klog));


  val = interpol2d(table_P_NL, Ntable.N_a, limits.a_min, 1., da, a, Ntable.N_k_nlin, logkmin, logkmax, dk, klog, cosmology.n_spec, 0.0);
//  printf("%le %le\n",k_NL,exp(val));
  return exp(val); 
  // returns the dimensionless power spectrum as a function of scale factor a and k in units of h/Mpc 
}

double Delta_NL_Coyote_only(double k_NL,double a)
{     
  static cosmopara C;

  static double logkmin = 0., logkmax = 0., dk = 0., da = 0.;
  static double **table_P_NL=0;
  
  double aa,klog,val; 
  double COSMO_emu[7],ystar[2*582],kstar[582],p_coyote[582],coyote_min,coyote_max;
  int type=1;
  int i,j,k;
  
  if (recompute_cosmo3D(C)){
    update_cosmopara(&C);

    if (table_P_NL!=0) free_double_matrix(table_P_NL,0, Ntable.N_a-1, 0, Ntable.N_k_nlin-1);     
    table_P_NL = create_double_matrix(0, Ntable.N_a-1, 0,Ntable.N_k_nlin-1);     
    da = (1. - limits.a_min)/(Ntable.N_a-1.);
    logkmin = log(limits.k_min_mpc);
    logkmax = log(limits.k_max_mpc);
    dk = (logkmax - logkmin)/(Ntable.N_k_nlin-1.);
    
    //printf("Starting P_delta %le %le %le %le %le %le %le\n",cosmology.Omega_m,cosmology.omb,cosmology.n_spec,cosmology.sigma_8, cosmology.w0,cosmology.wa,cosmology.h0);
    
    //compute Halofit; determine whether outside cosmology of emulator -> break
   
   COSMO_emu[1] = cosmology.Omega_m*cosmology.h0*cosmology.h0;
   COSMO_emu[0] = cosmology.omb*cosmology.h0*cosmology.h0;
   COSMO_emu[2] = cosmology.n_spec;
   COSMO_emu[3] = cosmology.h0*100.;
   COSMO_emu[4] = cosmology.w0;
   COSMO_emu[5] = cosmology.sigma_8;   
   COSMO_emu[6] = 0.0;

    aa = limits.a_min;
      //binning in k and a must be the same as in  Delta_halofit
     
      for (i=0; i<Ntable.N_a; i++, aa +=da) {
  gsl_interp_accel *acc = gsl_interp_accel_alloc ();
  gsl_spline *timspline = gsl_spline_alloc (gsl_interp_cspline, 582);
  
  COSMO_emu[6] = (1.0/aa)-1.0; //emu takes 7 args 6 cosmopara and 7th redshift    
  if(aa < a_min_coyote) COSMO_emu[6] = 4.0; //max redshift
  if(fabs(COSMO_emu[6])<1.e-10)   COSMO_emu[6]=0.0;
   emu(COSMO_emu,ystar, &type);
    for (k=0; k<582; k++){
      kstar[k]=ystar[k];
      p_coyote[k]=ystar[k+582];
      //printf("%le %le\n",kstar[k],p_coyote[k]);
    }
    gsl_spline_init (timspline, kstar, p_coyote, 582);

    coyote_min=p_coyote[0];
    coyote_max=p_coyote[581];
    //printf("%le %le\n",kstar[k-1],p_coyote[k-1]);
    klog = logkmin; // log k in h/MPC
    for (j=0; j<Ntable.N_k_nlin; j++, klog += dk) {
      if ((klog >= log(k_min_coyote/cosmology.h0)) && (klog <= log(k_max_coyote/cosmology.h0))){
        table_P_NL[i][j]=log(gsl_spline_eval(timspline, exp(klog)*cosmology.h0, acc));
        //printf("Coyote used\n");
      }
      if(klog>log(k_max_coyote/cosmology.h0))  table_P_NL[i][j]=log(coyote_max);
      if(klog<log(k_min_coyote/cosmology.h0)) table_P_NL[i][j]=log(coyote_min);
        //printf("Halofit used: exceeded Coyote k range k=%le k_min=%le k_max=%le\n",exp(klog),k_min_coyote/cosmology.h0,k_max_coyote/cosmology.h0);
      }
      gsl_interp_accel_free(acc);
      gsl_spline_free(timspline);    
     }
   }    

  //  printf("%le\n",k_NL);
  klog = log(k_NL);
//  if(a < a_min_coyote || klog>log(k_max_coyote/cosmology.h0) || klog<log(k_min_coyote/cosmology.h0))printf("Halofit used: exceeded Coyote a or k range a=%le k=%le\n",a,exp(klog));

  val = interpol2d(table_P_NL, Ntable.N_a, limits.a_min, 1., da, a, Ntable.N_k_nlin, logkmin, logkmax, dk, klog, cosmology.n_spec, 0.0);
//  printf("%le %le\n",k_NL,exp(val));
  return exp(val); 
  // returns the dimensionless power spectrum as a function of scale factor a and k in units of h/Mpc 
}


// double fraction_baryons(double k,double a)
// {
//   static int READ_TABLE=0;
//   static double **table=0;
//   int i;
//   double val;
//   double z=1./a-1.;
//   int baryon_zbin=11;
//   int baryon_kbin=100;
//   static double dz=(2.0)/(10.0); //2,0 was the max zdistribution when binning the DES baryonic power spectra
//   static double dk=(10.0-.3)/(99.);
//   if (READ_TABLE==0){
//     table=create_double_matrix(0, baryon_kbin-1, 0, baryon_zbin-1);
//     for(i=0;i<baryon_kbin;i++){
//       table[i][0]=AGN_DESdepth[i][0];
//       table[i][1]=AGN_DESdepth[i][1];
//       table[i][2]=AGN_DESdepth[i][2];
//       table[i][3]=AGN_DESdepth[i][3];
//       table[i][4]=AGN_DESdepth[i][4];
//       table[i][5]=AGN_DESdepth[i][5];
//       table[i][6]=AGN_DESdepth[i][6];
//       table[i][7]=AGN_DESdepth[i][7];
//       table[i][8]=AGN_DESdepth[i][8];
//       table[i][9]=AGN_DESdepth[i][9];
//       table[i][10]=AGN_DESdepth[i][10];
//     }
//     READ_TABLE=1;
//   }
//   val = interpol2d(table, baryon_kbin, 0.3, 10.0, dk, k, baryon_zbin, 0.0, 2.0, dz, z, 0.0, 0.0);
//   return val; 
// }


/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
//Pdelta is called with k in units H0/c since the comoving distance chi is in units c/H0. Upstream Pdelta all routines are in h/mpc 
 double Pdelta(double k_NL, double a)
{ 
  static int P_type = -1;
  if (P_type == -1){
    if (strcmp(pdeltaparams.runmode,"Halofit")==0) P_type = 0;
    if (strcmp(pdeltaparams.runmode,"Coyote")==0) P_type = 1;
    if (strcmp(pdeltaparams.runmode,"Coyote_only")==0) P_type = 2;
    if (strcmp(pdeltaparams.runmode,"linear")==0) P_type = 3;
    if (strcmp(pdeltaparams.runmode,"CLASS")==0) P_type = 4;
    if (strcmp(pdeltaparams.runmode,"class")==0) P_type = 4;
    if (strcmp(pdeltaparams.runmode,"cosmo_sim_test") ==0) P_type = 5;
    if (strcmp(pdeltaparams.runmode,"horndeski")==0) P_type = 6;
  }
  
  // Sanity check (Coyote emulator and Halofit not supported in Hi_CLASS)
  if (cosmology.use_horndeski == 1){
    if (P_type != 6){
      printf("Pdelta(): Run mode '%s' not supported! Use 'horndeski' only.\n", 
             pdeltaparams.runmode);
      exit(1);
    }
  }

  //printf("%s set\n",pdeltaparams.runmode);
  double pdelta = 0.,kintern=k_NL/cosmology.coverH0,error,k_nonlin,res;
  int status;
  switch (P_type){
    case 0: 
            pdelta = 2.0*constants.pi_sqr*Delta_NL_Halofit(kintern,a)/k_NL/k_NL/k_NL;
            break;
    case 1: 
            pdelta = 2.0*constants.pi_sqr*Delta_NL_Coyote(kintern,a)/k_NL/k_NL/k_NL;
            break;
    case 2: 
            pdelta = 2.0*constants.pi_sqr*Delta_NL_Coyote_only(kintern,a)/k_NL/k_NL/k_NL; 
            break;
    case 3: pdelta = p_lin(k_NL,a); break;
    case 4: pdelta = p_class(k_NL,a,1, &status); break;
    case 5: 
            k_nonlin = nonlinear_scale_computation(a);
            if (kintern<0.01) pdelta = 2.0*constants.pi_sqr*Delta_NL_Halofit(kintern,a)/k_NL/k_NL/k_NL;
            else{ 
              error = 0.01*pow((pdeltaparams.DIFF_A*kintern/k_nonlin),pdeltaparams.DIFF_n);
              pdelta = 2.0*constants.pi_sqr*Delta_NL_Halofit(kintern,a)*(1.0+error)/k_NL/k_NL/k_NL;
            }
            break;
    case 6: pdelta = p_class(k_NL, a, 0, &status); break; // Horndeski mode requires nonlinear=0
    default: 
            printf("cosmo3D:Pdelta: %s Pdelta runmode not defined\n",pdeltaparams.runmode);
            printf("using Halofit (standard)\n");
            pdelta = 2.0*constants.pi_sqr*Delta_NL_Halofit(kintern,a)/k_NL/k_NL/k_NL;
            break;
    }
    // if(((1./a-1.)<2.0) && (kintern >0.3) && (kintern <10.)) { 
    //   res=sqrt(pdelta*fraction_baryons(kintern, a)*pdelta);
    // }
    // else res=pdelta;
  return pdelta;  
}    

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*============================================================
 *see BS 2.41 bzw Logbook for detailed calculation of chi from a.*/

double int_for_chi(double a, void * args){
  //double res,asqr;
  //asqr=a*a;
  //res= 1./sqrt(a*cosmology.Omega_m + asqr*(1.-cosmology.Omega_m -cosmology.Omega_v ) + asqr*asqr*omv_vareos(a));
  //return res;
  return 1./(a*a*hoverh0(a)); //changed to call of hoverh0 to be ready for other parametrizations
}

/*for the calculation of chi we have to integrate from a(z2)=a up to a(z1)=1, which means todays expansion factor*/	
double chi(double a)
{
  static cosmopara C;
  static double *table;
  static double da = 0.;
  double aa,res;
  int i;
  double array[1];
  
  if (recompute_expansion(C)){
    update_cosmopara(&C);
    da = (1.-limits.a_min)/(Ntable.N_a-1.);
    aa = limits.a_min;
    if (table!=0) free_double_vector(table, 0, Ntable.N_a-1);
    table   = create_double_vector(0, Ntable.N_a-1);
    for (i=0; i<Ntable.N_a-1; i++, aa+=da) {
      table[i] = int_gsl_integrate_medium_precision(int_for_chi,(void*)array, aa, 1.,NULL,1000);
    }
    table[Ntable.N_a-1] =.0;
  }
  res = interpol(table, Ntable.N_a, limits.a_min, 1., da, a, 0.0, 0.0); // comoving distance in c/H_0
  if (res < 0){printf ("interpolation error in chi(%e)\n",a); res=0.01;}
  return res;
}


/*===============================calculating the angular diameter distance f_K BS01 2.4, 2.30: f_K is a radial function that, depending on the curvature of the Universe, is a trigonometric, linear, or hyperbolic function of chi  */
double f_K(double chi)
{  
  double K, K_h, f;
  K = (cosmology.Omega_m   + cosmology.Omega_v  - 1.);
  if (K > precision.medium) {           /* open */
    K_h = sqrt(K); // K in units H0/c see BS eq. 2.30
    f = 1./K_h*sin(K_h*chi);
    //printf("open\n");
  } else if (K < -precision.medium) {   /* closed */
    K_h = sqrt(-K); 
    f = 1./K_h*sinh(K_h*chi);
    //printf("closed K=%le %le %le\n",K,cosmology.Omega_m,cosmology.Omega_v);
  } else {                     /* flat */
    f = chi;
    //printf("flatK=%le %le %le\n",K,cosmology.Omega_m,cosmology.Omega_v);
  }
  return f;
}




