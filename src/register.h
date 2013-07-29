#include <R.h>
#include<Rdefines.h>
#include <Rinternals.h>

SEXP mahalanobisHelper(SEXP data, SEXP treat_ids, SEXP control_ids,
                       SEXP invScaleMat);
SEXP subsetInfSparseMatrix(SEXP whichRows, SEXP whichCols, SEXP x);
