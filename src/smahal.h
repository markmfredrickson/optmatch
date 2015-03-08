#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>

#define SMAHAL_CUTOFF 1e-10

typedef struct dmat {
  int nr, nc;
  double * data;
} DMAT;

DMAT * smahal(int nr, int nc, double * data, int * z);
