// correlation functions
double w_gk_planck(double theta,int ni); //angular CMB lensing x positions correlation function in tomography bin ni
double w_ks_planck(double theta, int ni);//angular CMB lensing x galaxy shear correlation function in tomography bin ni

double beam_planck(double l){ // TO BE UPDATED!!!!
  double fwhm_arcmin =7;// 5.4 arcmin
  double sigma = fwhm_arcmin/sqrt(8.*log(2.0))*constants.arcmin;
  return exp(-0.5*l*l*sigma*sigma);
}

//================================================================================================
// angular correlation functions - no look-up tables
double C_gk_wrapper(double l,int ni){
  return C_gk(l,ni)*beam_planck(l);
}

double C_ks_wrapper(double l,int ni){
  return C_ks(l,ni)*beam_planck(l);
}

double w_gk_fullsky(int nt, int ni){
  static int LMAX = 100000;
  static int NTHETA = 0;
  static double ** Pl =0;
  static double *Cl =0;
  static double *w_vec =0;
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  int i,l,nz;
  if (like.Ntheta ==0){
    printf("cosmo2D_real.c:w_tomo_exact: like.Ntheta not initialized\nEXIT\n"); exit(1);
  }
  if (Pl ==0){
    Pl =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
    Cl = create_double_vector(0,LMAX-1);
    w_vec = create_double_vector(0,tomo.clustering_Nbin*like.Ntheta-1);
    NTHETA = like.Ntheta;
    double *xmin, *xmax, *Pmin, *Pmax;
    xmin= create_double_vector(0, like.Ntheta-1);
    xmax= create_double_vector(0, like.Ntheta-1);
    double logdt=(log(like.vtmax)-log(like.vtmin))/like.Ntheta;
    Pmin= create_double_vector(0, LMAX+1);
    Pmax= create_double_vector(0, LMAX+1);
    for(i=0; i<like.Ntheta ; i++){
      xmin[i]=cos(exp(log(like.vtmin)+(i+0.0)*logdt));
      xmax[i]=cos(exp(log(like.vtmin)+(i+1.0)*logdt));
    }

    for (i = 0; i<NTHETA; i ++){
      printf("Tabulating Legendre coefficients %d/%d\n",i+1, NTHETA);
      gsl_sf_legendre_Pl_array(LMAX, xmin[i],Pmin);
      gsl_sf_legendre_Pl_array(LMAX, xmax[i],Pmax);
      //gsl_sf_legendre_Pl_array(LMAX, cos(like.theta[i]),Pmin);      
      for (int l = 1; l < LMAX; l ++){
            //Pl[i][l] = (2*l+1.)/(4.*M_PI)*Pmin[l];
        //Pl[i][l] = (2*l+1.)/(4.*M_PI)*gsl_sf_legendre_Pl(l,cos(like.theta[i]));
        Pl[i][l] = 1./(4.*M_PI)*(Pmin[l+1]-Pmax[l+1]-Pmin[l-1]+Pmax[l-1])/(xmin[i]-xmax[i]);
      }
    }
    free_double_vector(xmin,0,like.Ntheta-1);
    free_double_vector(xmax,0,like.Ntheta-1);
    free_double_vector(Pmin,0,LMAX+1);
    free_double_vector(Pmax,0,LMAX+1);
    }
    if (recompute_gk(C,G,N,ni)){
      for (nz = 0; nz <tomo.clustering_Nbin; nz ++){
        for (l = 1; l < LMAX; l++){ 
//          if (l < 20){Cl[l]=C_cl_RSD_nointerp(l,nz,nz);}
          Cl[l]=C_gk_wrapper(1.0*l,nz);
          // if (l < 20){Cl[l]=C_cl_tomo_nointerp(l,nz,nz);}
          // else Cl[l]=C_cl_tomo(1.0*l,nz,nz);
        }
        for (i = 0; i < NTHETA; i++){
          w_vec[nz*like.Ntheta+i] =0;
          for (l = 1; l < LMAX; l++){
            w_vec[nz*like.Ntheta+i]+=Pl[i][l]*Cl[l];
          }
        }
      }
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }
  return w_vec[ni*like.Ntheta+nt];  
}


double w_ks_fullsky(int nt, int ni){
  static int LMAX = 100000;
  static int NTHETA = 0;
  static double ** Pl =0;
  static double *Cl =0;
  static double *w_vec =0;
  static cosmopara C;
  static nuisancepara N;
  static galpara G;
  int i,l,nz;
  if (like.Ntheta ==0){
    printf("cosmo2D_fullsky.c:w_gamma_t_tomo: like.Ntheta not initialized\nEXIT\n"); exit(1);
  }
  if (Pl ==0){
    Pl =create_double_matrix(0, like.Ntheta-1, 0, LMAX-1);
    Cl = create_double_vector(0,LMAX-1);
    w_vec = create_double_vector(0,tomo.shear_Nbin*like.Ntheta-1);
    NTHETA = like.Ntheta;
    double *xmin, *xmax, *Pmin, *Pmax, *dP;
    xmin= create_double_vector(0, like.Ntheta-1);
    xmax= create_double_vector(0, like.Ntheta-1);
    double logdt=(log(like.vtmax)-log(like.vtmin))/like.Ntheta;
    for(i=0; i<like.Ntheta ; i++){
      xmin[i]=cos(exp(log(like.vtmin)+(i+0.0)*logdt));
      xmax[i]=cos(exp(log(like.vtmin)+(i+1.0)*logdt));
    }
    Pmin= create_double_vector(0, LMAX+1);
    Pmax= create_double_vector(0, LMAX+1);

    for (i = 0; i<NTHETA; i ++){
      printf("Tabulating Legendre coefficients %d/%d\n",i+1, NTHETA);
      gsl_sf_legendre_Pl_array(LMAX, xmin[i],Pmin);
      gsl_sf_legendre_Pl_array(LMAX, xmax[i],Pmax);
      for (int l = 1; l < LMAX; l ++){
        //Pl[i][l] = (2.*l+1)/(4.*M_PI*l*(l+1))*gsl_sf_legendre_Plm(l,2,cos(like.theta[i]));  
        Pl[i][l] = (2.*l+1)/(4.*M_PI*l*(l+1)*(xmin[i]-xmax[i]))
        *((l+2./(2*l+1.))*(Pmin[l-1]-Pmax[l-1])
        +(2-l)*(xmin[i]*Pmin[l]-xmax[i]*Pmax[l])
        -2./(2*l+1.)*(Pmin[l+1]-Pmax[l+1]));
      }
    }
    free_double_vector(xmin,0,like.Ntheta-1);
    free_double_vector(xmax,0,like.Ntheta-1);
    free_double_vector(Pmin,0,LMAX+1);
    free_double_vector(Pmax,0,LMAX+1);
  }
  if (recompute_ks(C,G,N,ni)){
    // if (like.IA != 0 && like.IA != 3 && like.IA != 4){printf("cosmo2D_real.c: w_ks_fullsky does not support like.IA = %d yet\nEXIT!\n",like.IA); exit(1);}
    // C_tomo_pointer C_gl_pointer = &C_gl_tomo;
    // if (like.IA ==3 || like.IA ==4) C_gl_pointer = &C_ggl_IA_tab;

    for (nz = 0; nz <tomo.shear_Nbin; nz ++){
      for (l = 1; l < LMAX; l++){
        // Cl[l]=C_ggl_IA_tab(1.0*l,ZL(nz),ZS(nz));
        Cl[l]=C_ks_wrapper(1.0*l,nz);
      }
      for (i = 0; i < NTHETA; i++){
        w_vec[nz*like.Ntheta+i] =0;
        for (l = 2; l < LMAX; l++){
          w_vec[nz*like.Ntheta+i]+=Pl[i][l]*Cl[l];
        }
      }
    }
    update_cosmopara(&C);
    update_galpara(&G);
    update_nuisance(&N);
  }
  return w_vec[ni*like.Ntheta+nt];  
}
