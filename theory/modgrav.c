
// Flag function that throws an error if Horndeski parameters are enabled in 
// the given cosmology. Should be added to all functions that make assumptions 
// that are invalid for Horndeski models
void invalid_for_horndeski(const char* func_name){
    if (cosmology.use_horndeski == 1){
        fprintf(stderr, "Horndeski model is being used, but function '%s' makes assumptions that are invalid for Horndeski.\n", func_name);
        #ifdef EXIT_IF_INVALID_FOR_HORNDESKI
        exit(1);
        #endif
    }
}

double MG_Sigma(double a)
{
//  double aa=a*a;
//  double omegam=cosmology.Omega_m/(aa*a);
  double omegav=omv_vareos(a);
  double hub=hoverh0(a);
  hub = hub*hub;
 
  return cosmology.MGSigma*omegav/hub/cosmology.Omega_v;
}

double horndeski_sigma(double a){
  // Horndeski: Assumes quasi-static limit, as in Eq. 6 of arXiv:1705.04714
  // mu_light = 1 + Sigma (???)
  if (cosmology.use_horndeski != 1) return 0.;
  
  double sigma;
  int status = 0;
  sigma = horndeski_function(a, "mu_light", &status) - 1.; // FIXME: -1?
  
  //printf("Horndeski: sigma = %4.4e\n", sigma);
  //printf("\t%3.3e, %3.3e, %3.3e, %3.3e\n", cosmology.mg_alpha_xk, cosmology.mg_alpha_xb, cosmology.mg_alpha_xm, cosmology.mg_alpha_xt);
  
  if (status != 0){
    fprintf(stderr, "modgrav.c: horndeski_sigma() failed.\n");
    exit(1);
  }
  return sigma;
}

double horndeski_mu(double a){
  // Horndeski: Assumes quasi-static limit, as in Eq. 2 of arXiv:1705.04714
  // horndeski_mu = 1 + mg_mu (???)
  if (cosmology.use_horndeski != 1) return 0.;
  
  double mu;
  int status = 0;
  mu = horndeski_function(a, "mu_eff", &status); // FIXME: -1?
  
  //printf("Horndeski: mu = %4.4e\n", mu);
  //printf("\t%3.3e, %3.3e, %3.3e, %3.3e\n", cosmology.mg_alpha_xk, cosmology.mg_alpha_xb, cosmology.mg_alpha_xm, cosmology.mg_alpha_xt);
  
  if (status != 0){
    fprintf(stderr, "modgrav.c: horndeski_mu() failed.\n");
    exit(1);
  }
  return mu;
}

