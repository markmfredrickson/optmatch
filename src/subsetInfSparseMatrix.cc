#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix subsetInfSparseMatrix(LogicalVector whichRows,
				    LogicalVector whichCols, S4 ismX) {
  int
    nRows = whichRows.size(),
    nNewRows = 0;

  std::vector<int> newRowNumbers(nRows, 0);
  for (int i = 0; i < nRows; i++) {
    if (whichRows[i]) {
      nNewRows++;
      newRowNumbers[i] = nNewRows;
    }
  }

  int
    nCols = whichCols.size(),
    nNewCols = 0;

  std::vector<int> newColNumbers(nCols, 0);
  for (int i = 0; i < nCols; i++) {
    if (whichCols[i]) {
      nNewCols++;
      newColNumbers[i] = nNewCols;
    }
  }

  NumericVector x( ismX );

  int
    nData = x.size(),
    nPoints = 0;

  IntegerVector
    rows( ismX.slot("rows") ),
    cols( ismX.slot("cols") );

  for (int i = 0; i < nData; i++) {
    if (whichRows[rows[i] - 1] && whichCols[cols[i] - 1])
      nPoints++;
  }

  NumericMatrix ans(nPoints, 3);
  for (int i = 0, j = 0; i < nData; i++) {
    if (whichRows[rows[i] - 1] && whichCols[cols[i] - 1]) {
      ans(j, 0) = newRowNumbers[rows[i] - 1];
      ans(j, 1) = newColNumbers[cols[i] - 1];
      ans(j, 2) = x[i];
      j++;
    }
  }
  return(ans);
}
