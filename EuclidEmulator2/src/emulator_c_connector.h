#ifndef EMULATOR_C_CONNECTOR_H 
#define EMULATOR_C_CONNECTOR_H 

#ifdef __cplusplus
extern "C" {
#endif
     
     void EuclidEmulator_compute_nlc(double Omega_b, double Omega_m, double Sum_m_nu, double n_s, double h, double w_0, double w_a, double A_s, double* redshift, int n_redshift, double *, double[101][613]);

#ifdef __cplusplus
}
#endif
#endif

