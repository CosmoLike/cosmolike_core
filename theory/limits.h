#ifndef COSMOLIKE_LIMITS_H
#define COSMOLIKE_LIMITS_H

/* maximum tomographic bin number */
#ifndef MAX_TOMO_BINS
#define MAX_TOMO_BINS 30
#define MAX_PAIRS (MAX_TOMO_BINS * (MAX_TOMO_BINS - 1) / 2)
#endif

/* The following are for exact_fft */
#ifndef NELL_BLOCK
#define NELL_BLOCK 100
#endif
#ifndef NCHI
#define NCHI 1000
#endif

#endif