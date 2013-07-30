#include <R.h>
#include<Rdefines.h>
#include <Rinternals.h>

SEXP mahalanobisHelper(SEXP data, SEXP index, SEXP invScaleMat);
SEXP subsetInfSparseMatrix(SEXP whichRows, SEXP whichCols, SEXP x);
