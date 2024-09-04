#include "subsetInfSparseMatrix.h"

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
SEXP subsetInfSparseMatrix(SEXP whichRows, SEXP whichCols, SEXP x) {
  Rcpp::S4 ismX(x);
  Rcpp::IntegerVector
    rows( ismX.slot("rows") ),
    cols( ismX.slot("cols") );

    int
        nRows = Rf_length(whichRows),
        nCols = Rf_length(whichCols),
        nData = Rf_length(x),
        nNewRows = 0,
        nNewCols = 0;

    int
        * iWhichRows = INTEGER(whichRows),
        * iWhichCols = INTEGER(whichCols);

    int * newRowNumbers = R_Calloc(nRows, int);
    for(int i = 0; i < nRows; i++) {
        if(iWhichRows[i] == 1) {
            nNewRows++;
            newRowNumbers[i] = nNewRows;
        }
    }

    int * newColNumbers = R_Calloc(nCols, int);
    for(int i = 0; i < nCols; i++) {
        if(iWhichCols[i] == 1) {
            nNewCols++;
            newColNumbers[i] = nNewCols;
        }
    }

    int nPoints = 0;
    for(int i = 0; i < nData; i++) {
        if(iWhichRows[rows[i] - 1] && iWhichCols[cols[i] - 1])
            nPoints++;
    }

    SEXP ans;
    PROTECT(ans = Rf_allocMatrix(REALSXP, nPoints, 3));
    double
        * dData = REAL(x),
        * rans = REAL(ans);

    for(int i = 0, j = 0; i < nData; i++) {
        if(iWhichRows[rows[i] - 1] && iWhichCols[cols[i] - 1]) {
            rans[j] = newRowNumbers[rows[i] - 1];
            rans[j + nPoints] = newColNumbers[cols[i] - 1];
            rans[j + 2 * nPoints] = dData[i];
            j++;
        }
    }

    R_Free(newRowNumbers);
    R_Free(newColNumbers);

    UNPROTECT(1);
    return(ans);
}
