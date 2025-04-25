#include <cstdlib>

#include "emulator_c_connector.h"
#include "emulator.h"
#include "cosmo.h"

#ifdef __cplusplus
extern "C" {
#endif

// Inside this "extern C" block, I can implement functions in C++, which will externally 
//   appear as C functions (which means that the function IDs will be their names, unlike
//   the regular C++ behavior, which allows defining multiple functions with the same name
//   (overloading) and hence uses function signature hashing to enforce unique IDs),


static EuclidEmulator *EuclidEmulator_instance = NULL;

void lazyEuclidEmulator() {
    if (EuclidEmulator_instance == NULL) {
        EuclidEmulator_instance = new EuclidEmulator();
    }
}

void EuclidEmulator_compute_nlc(double Omega_b, double Omega_m, double Sum_m_nu, double n_s, double h, double w_0, double w_a, double A_s, double* redshift, int n_redshift, double* kout, double pkout[101][613]) {
    lazyEuclidEmulator();
    Cosmology csm = Cosmology(Omega_b, Omega_m, Sum_m_nu, n_s, h, w_0, w_a, A_s);
    EuclidEmulator_instance->compute_nlc(csm, vector<double>(redshift, redshift+n_redshift), n_redshift);
    int i=0, j=0;
    for(i=0; i<603; ++i){
        kout[i] = EuclidEmulator_instance->kvec[i];
        for(j=0; j<101; ++j)
            pkout[j][i] =  pow(10.0,EuclidEmulator_instance->Bvec[j][i]);
    }
}

#ifdef __cplusplus
}
#endif
