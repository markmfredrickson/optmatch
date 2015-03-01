#include "smahal.h"

SEXP r_smahal(SEXP index, SEXP data, SEXP z) {
  DMAT * ans = smahal(nrows(data), ncols(data), REAL(data), LOGICAL(z));
  if(ans == NULL || ans->nr < 1 || ans->nc <1)
    error("smahal_nosexp returned an invalid answer");

  SEXP out;
  PROTECT(out = allocMatrix(REALSXP, ans->nr, ans->nc));
  memcpy(REAL(out), ans->data, ans->nr * ans->nc * sizeof(double));

  Free(ans->data);
  Free(ans);
  UNPROTECT(1);

  return out;
}
