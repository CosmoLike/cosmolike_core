double MG_mu(double omegav, double hub);
double MG_sigma(double a, double omegav, double hub);
double horndeski_mu(double a, double omegav, double hub, double omegam);
double horndeski_sigma(double a, double omegav, double hub, double omegam);
double horndeski_cs_sqrd(double a, double fz, double fzprime, double H, double Hprime, double omegam, double hub);
double horndeski_alpha(double a, double fz);
/**********************************************************/

double MG_mu(double omegav, double hub){
  return cosmology.MGmu*omegav/hub/cosmology.Omega_v;
}

double MG_sigma(double a, double omegav, double hub){
  return cosmology.MGSigma*omegav/hub/cosmology.Omega_v;
}
double horndeski_mu(double a, double omegav, double hub, double omegam){
  //this is all in the quasistatic limit!!
  //given by Eqn 43 in Tessa Baker notes
    //!!!!currently assuming Mpl = Mstar!!!!
  double fz = omegav/hub*1./cosmology.Omega_v;
  double wde = cosmology.w0 + cosmology.wa*(1. - a);
  double H = (100. * cosmology.h0 / constants.lightspeed) * sqrt(hub);
  double Hprime = -0.5*H*H*(omegam/hub+omegav/hub*(1.+3.*wde))-H*H;
  double fzprime = -3.*wde * a*H * omegam/hub * fz;
  double alpha = horndeski_alpha(a, fz);
  double cs_sqrd = horndeski_cs_sqrd(a, fz, fzprime, H, Hprime,omegam, hub);
    /////////////
  double alpha_K = cosmology.MGalpha_K*fz;
  double alpha_B = cosmology.MGalpha_B*fz;
  double alpha_T = cosmology.MGalpha_T*fz;
  double alpha_M = cosmology.MGalpha_M*fz;
  ////////////
  double mu_term1 = alpha*cs_sqrd*(1+alpha_T);
  double mu_term2 = (-alpha_B*(1+alpha_T)+ 2*(alpha_T - alpha_M))*(-alpha_B*(1+alpha_T)+ 2*(alpha_T - alpha_M));
  double mu = (mu_term1+mu_term2)/(alpha*cs_sqrd) - 1.;
  return mu;
}

double horndeski_sigma(double a, double omegav, double hub, double omegam){
  //this is all in the quasistatic limit!!
  //also assumes Mpl = Mstar!!!!
  //alpha function redshift scaling
  //i.e. scale with Omega_DE(a)/Omega_DE(a=1)
  //then, the alphas are equal to the constant values at today.
  double fz = omegav/hub*1./cosmology.Omega_v;
  double wde = cosmology.w0 + cosmology.wa*(1. - a);
  double H = (100. * cosmology.h0 / constants.lightspeed) * sqrt(hub);
  double Hprime = -0.5*H*H*(omegam/hub+omegav/hub*(1.+3.*wde))-H*H;
  double fzprime = -3.*wde * a*H * omegam/hub * fz;
  //define actual alpha functions: fz*alpha_X,0
  double alpha_K = cosmology.MGalpha_K*fz;
  double alpha_B = cosmology.MGalpha_B*fz;
  double alpha_T = cosmology.MGalpha_T*fz;
  double alpha_M = cosmology.MGalpha_M*fz;
  //helpful functions (see Tessa Baker's notes)
  double alpha = horndeski_alpha(a, fz);
  double cs_sqrd = horndeski_cs_sqrd(a, fz, fzprime, H, Hprime, omegam, hub);
  //given by gamma Eqn. 44 in Tessa Baker notes
  double num = alpha*cs_sqrd-alpha_B*(-alpha_B/(2.*(1.+alpha_T))+alpha_T-alpha_M);
  double denom = alpha*cs_sqrd*(1.+alpha_T)+(-alpha_B*(1.+alpha_T)+2.*(alpha_T-alpha_M))*(-alpha_B*(1.+alpha_T)+2.*(alpha_T-alpha_M));
  double mu = horndeski_mu(a, omegav, hub, omegam);
  return 0.5*(1.+mu)*(num/denom+1.);
}
double horndeski_cs_sqrd(double a, double fz, double fzprime, double H, double Hprime, double omegam, double hub){
  double alpha_K = cosmology.MGalpha_K*fz;
  double alpha_B = cosmology.MGalpha_B*fz;
  double alpha_T = cosmology.MGalpha_T*fz;
  double alpha_M = cosmology.MGalpha_M*fz;

  double alpha = horndeski_alpha(a, fz);
  double cs_term1 = (2.-alpha_B)*(Hprime - H*H*((alpha_M -alpha_T)+alpha_B/(2.*(1.+alpha_T))));
  //!!!!!!!!Assuming Mpl = Mstar!!!!!!!!!!!!!
  double cs_term2 = -H*cosmology.MGalpha_B*fzprime + 1.5*H*H*omegam/hub;
  double cs_sqrd = -(cs_term1 + cs_term2)/(H*H*alpha);
  return cs_sqrd;
}
double inline horndeski_alpha(double a, double fz){
	//Helpful function -- see Tessa Baker's notes
	double alpha_K = cosmology.MGalpha_K*fz;
	double alpha_B = cosmology.MGalpha_B*fz;
	return alpha_K + 1.5*alpha_B*alpha_B;
}

