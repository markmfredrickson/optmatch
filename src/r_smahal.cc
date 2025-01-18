#include "r_smahal.h"
#include "smahal.h"

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP r_smahal(SEXP index, SEXP data, SEXP z) {
  Rcpp::NumericMatrix dataMat(data);
  DMAT * ans = smahal(dataMat.nrow(), dataMat.ncol(), REAL(data), LOGICAL(z));
  if(ans == NULL || ans->nr < 1 || ans->nc <1)
    Rf_error("smahal_nosexp returned an invalid answer");

  SEXP out;
  Rf_protect(out = Rf_allocMatrix(REALSXP, ans->nr, ans->nc));
  memcpy(REAL(out), ans->data, ans->nr * ans->nc * sizeof(double));

  R_Free(ans->data);
  R_Free(ans);
  Rf_unprotect(1);

  return out;
}
