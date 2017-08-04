#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
  Check these declarations against the C/Fortran source code.
*/

  /* .Call calls */
  extern SEXP optmatch_ismOps(SEXP, SEXP, SEXP);
extern SEXP optmatch_mahalanobisHelper(SEXP, SEXP, SEXP);
extern SEXP optmatch_r_smahal(SEXP, SEXP, SEXP);
extern SEXP optmatch_subsetInfSparseMatrix(SEXP, SEXP, SEXP);

/* .Fortran calls */
extern void F77_NAME(relaxalg)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(relaxalgold)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CallMethodDef CallEntries[] = {
  {"optmatch_ismOps",                (DL_FUNC) &optmatch_ismOps,                3},
  {"optmatch_mahalanobisHelper",     (DL_FUNC) &optmatch_mahalanobisHelper,     3},
  {"optmatch_r_smahal",              (DL_FUNC) &optmatch_r_smahal,              3},
  {"optmatch_subsetInfSparseMatrix", (DL_FUNC) &optmatch_subsetInfSparseMatrix, 3},
  {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
  {"relaxalg",    (DL_FUNC) &F77_NAME(relaxalg),    12},
  {"relaxalgold", (DL_FUNC) &F77_NAME(relaxalgold), 11},
  {NULL, NULL, 0}
};

void R_init_optmatch(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
