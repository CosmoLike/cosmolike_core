typedef struct config {
	double nu;
	double c_window_width;
	int derivative;
	long N_pad;
	long N_extrap_low;
	long N_extrap_high;
} config;

void cfftlog(double *x, double *fx, long N, config *config, int ell, double *y, double *Fy);

void cfftlog_ells(double *x, double *fx, long N, config *config, int* ell, long Nell, double **y, double **Fy);

void cfftlog_ells_increment(double *x, double *fx, long N, config *config, int* ell, long Nell, double **y, double **Fy);

void cfftlog_multiple(double *x, double **fx, long N, long Nf, config *config, int ell, double *y, double **Fy);