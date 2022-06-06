#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .Call calls */
extern SEXP _optmatch_ismOps(SEXP, SEXP, SEXP);
extern SEXP _optmatch_mahalanobisHelper(SEXP, SEXP, SEXP);
extern SEXP _optmatch_r_smahal(SEXP, SEXP, SEXP);
extern SEXP _optmatch_subsetInfSparseMatrix(SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_optmatch_ismOps",                (DL_FUNC) &_optmatch_ismOps,                3},
  {"_optmatch_mahalanobisHelper",     (DL_FUNC) &_optmatch_mahalanobisHelper,     3},
  {"_optmatch_r_smahal",              (DL_FUNC) &_optmatch_r_smahal,              3},
  {"_optmatch_subsetInfSparseMatrix", (DL_FUNC) &_optmatch_subsetInfSparseMatrix, 3},
  {NULL, NULL, 0}
};

void R_init_optmatch(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
