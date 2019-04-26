typedef struct config {
	double nu;
	double c_window_width;
	int derivative;
} config;

void cfftlog(double *x, double *fx, long N, config *config, int ell, double *y, double *Fy);

void cfftlog_ells(double *x, double *fx, long N, config *config, int* ell, long Nell, double **y, double **Fy);