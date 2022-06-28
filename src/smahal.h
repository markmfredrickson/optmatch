#include <R.h>
#include <R_ext/Lapack.h>

#ifndef _SMAHAL_H_
#define _SMAHAL_H_

// Workaround for pre-3.6 R compatibility:
#ifndef FCONE
#define FCONE
#endif

#define SMAHAL_CUTOFF 1e-10

typedef struct dmat {
  int nr, nc;
  double * data;
} DMAT;

DMAT * smahal(int nr, int nc, double * data, int * z);

#endif
