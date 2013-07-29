#include "register.h"

SEXP subsetInfSparseMatrix(SEXP whichRows, SEXP whichCols, SEXP x) {
    SEXP
        rows = GET_SLOT(x, Rf_install("rows")),
        cols = GET_SLOT(x, Rf_install("cols"));
    int
        nRows = length(whichRows),
        nCols = length(whichCols),
        nData = length(x),
        nNewRows = 0,
        nNewCols = 0;

    int
        * iWhichRows = INTEGER(whichRows),
        * iWhichCols = INTEGER(whichCols),
        * iRows = INTEGER(rows),
        * iCols = INTEGER(cols);

    int * newRowNumbers = Calloc(nRows, int);
    for(int i = 0; i < nRows; i++) {
        if(iWhichRows[i] == 1) {
            nNewRows++;
            newRowNumbers[i] = nNewRows;
        }
    }

    int * newColNumbers = Calloc(nCols, int);
    for(int i = 0; i < nCols; i++) {
        if(iWhichCols[i] == 1) {
            nNewCols++;
            newColNumbers[i] = nNewCols;
        }
    }

    int nPoints = 0;
    for(int i = 0; i < nData; i++) {
        if(iWhichRows[iRows[i] - 1] && iWhichCols[iCols[i] - 1])
            nPoints++;
    }

    SEXP ans;
    PROTECT(ans = allocMatrix(REALSXP, nPoints, 3));
    double
        * dData = REAL(x),
        * rans = REAL(ans);

    for(int i = 0, j = 0; i < nData; i++) {
        if(iWhichRows[iRows[i] - 1] && iWhichCols[iCols[i] - 1]) {
            rans[j] = newRowNumbers[iRows[i] - 1];
            rans[j + nPoints] = newColNumbers[iCols[i] - 1];
            rans[j + 2 * nPoints] = dData[i];
            j++;
        }
    }

    Free(newRowNumbers);
    Free(newColNumbers);

    UNPROTECT(1);
    return(ans);
}
