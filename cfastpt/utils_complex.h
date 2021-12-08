#include <complex.h>
#include <fftw3.h>
void f_z(double z_real, double *z_imag, double complex *fz, long N);

void c_window(double complex *out, double c_window_width, long halfN);

// void resample_fourier_gauss(double *k, double *fk, config *config);

double complex gamma_lanczos(double complex z);
double complex lngamma_lanczos(double complex z);

void g_m_vals(double mu, double q_real, double *q_imag, double complex *gm, long N);


void fftconvolve(fftw_complex *in1, fftw_complex *in2, long N, fftw_complex *out);


void fftconvolve_optimize(fftw_complex *in1, fftw_complex *in2, long N, fftw_complex *out, fftw_complex *a, fftw_complex *b, fftw_complex *a1, fftw_complex *b1, fftw_complex *c, fftw_plan pa, fftw_plan pb, fftw_plan pc);

void fftconvolve_real(double *in1, double *in2, long N1, long N2, double *out);