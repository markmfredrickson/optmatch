#include "ism.h"
#include <numeric>

// compare row col positions
int cmp(int rowA, int colA, int rowB, int colB) {
  if (rowA > rowB) return  1;
  if (rowA < rowB) return -1;
  if (colA > colB) return  1;
  if (colA < colB) return -1;
  return 0;
}

// modified quicksort
void ismSortIndex (const Rcpp::IntegerVector & rows,
                   const Rcpp::IntegerVector & cols,
                   int * index, int n) {
  int
    i, ii, j, ij,
    p[2], swap;

  if (n < 2) return;

  i = index[n / 2];
  p[0] = rows[i];
  p[1] = cols[i];
    
  for (i = 0, j = n - 1;; ++i, --j) {
    ii = index[i];
    ij = index[j];
    while (rows[ii] < p[0] || (rows[ii] == p[0] && cols[ii] < p[1])) {
      ++i;
      ii = index[i];
    }
    while (p[0] < rows[ij] || (p[0] == rows[ij] && p[1] < cols[ij])) {
      --j;
      ij = index[j];
    }

    if (i >= j) break;

    swap = index[i];
    index[i] = index[j];
    index[j] = swap;
  }
  ismSortIndex(rows, cols, index, i);
  ismSortIndex(rows, cols, index + i, n - i);
}

// modified binary search
int ismLubIndex(int findRow, int findCol,
                const Rcpp::IntegerVector & rows,
                const Rcpp::IntegerVector & cols,
                const int * index, unsigned n) {
  unsigned
    left = 0,
    right = n;

  while (left < right) {
    unsigned m = (left  +  right) / 2;
    if(m >= n) return m;
    int ma = index[m];
    if (cmp(rows[ma], cols[ma], findRow, findCol) < 0) 
      left = m + 1;
    else 
      right = m;
  }
  return right;
}

SEXP ismOps(SEXP o, SEXP a, SEXP b) {
  Rcpp::S4
    ismA(a), ismB(b);
  std::string op = Rcpp::as<std::string>(o);

  int i, ia, ib;
  const Rcpp::IntegerVector
    rowsA( ismA.slot("rows") ),
    colsA( ismA.slot("cols") );
  Rcpp::IntegerVector indexA(rowsA.size());
  for(i = 0; i < indexA.size(); i++)
    indexA[i] = i;
  ismSortIndex(rowsA, colsA, indexA.begin(), indexA.size());

  const Rcpp::IntegerVector
    rowsB( ismB.slot("rows") ),
    colsB( ismB.slot("cols") );
  Rcpp::IntegerVector indexB(rowsB.size());
  for(i = 0; i < indexB.size(); i++)
    indexB[i] = i;
  ismSortIndex(rowsB, colsB, indexB.begin(), indexB.size());

  int
    found = 0, lastFound = 0,
    nPairs = 0;
  Rcpp::IntegerMatrix pairs(indexA.size(), 2);
  
  for (i = 0; i < indexA.size() && lastFound < indexB.size(); ++i) {
    ia = indexA[i];
    found = ismLubIndex(rowsA[ia], colsA[ia],
                   rowsB, colsB,
                   indexB.begin() + lastFound,
                   indexB.size() - lastFound);
    found += lastFound;
    if (found >= indexB.size()) continue;

    ib = indexB[found];
    if (rowsB[ib] != rowsA[ia] ||
        colsB[ib] != colsA[ia]) {
      continue;
    }

    pairs(nPairs, 0) = ia;
    pairs(nPairs, 1) = ib;
    lastFound = found + 1;
    ++nPairs;
  }

  Rcpp::NumericVector
    ansData(nPairs),
    dataA( ismA.slot(".Data") ),
    dataB( ismB.slot(".Data") );

  Rcpp::IntegerVector
    ansRows(nPairs),
    ansCols(nPairs);

  // trying to avoid putting op test
  // inside loop for efficiency (?)
  if (op == "+") {
    for (i = 0; i < nPairs; i++) {
      ia = pairs(i, 0);
      ib = pairs(i, 1);
      ansRows[i] = rowsA[ia];
      ansCols[i] = colsA[ia];
      ansData[i] = dataA[ia] + dataB[ib];
    }
  }  else if (op == "-") {
    for (i = 0; i < nPairs; i++) {
      ia = pairs(i, 0);
      ib = pairs(i, 1);
      ansRows[i] = rowsA[ia];
      ansCols[i] = colsA[ia];
      ansData[i] = dataA[ia] - dataB[ib];
    }
  }  else if (op == "*") {
    for (i = 0; i < nPairs; i++) {
      ia = pairs(i, 0);
      ib = pairs(i, 1);
      ansRows[i] = rowsA[ia];
      ansCols[i] = colsA[ia];
      ansData[i] = dataA[ia] * dataB[ib];
    }
  }  else if (op == "/") {
    for (i = 0; i < nPairs; i++) {
      ia = pairs(i, 0);
      ib = pairs(i, 1);
      ansRows[i] = rowsA[ia];
      ansCols[i] = colsA[ia];
      ansData[i] = dataA[ia] / dataB[ib];
    }
  }

  ismA.slot(".Data") = ansData;
  ismA.slot("rows") = ansRows;
  ismA.slot("cols") = ansCols;

  return ismA;
}
