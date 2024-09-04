#include <Rcpp.h>
#include <unordered_map>

using namespace Rcpp;

/* function: mahalanobisHelper
     consumes: SEXP data, index, invScaleMat
       SEXP data: an matrix with rownames
       SEXP index: a 2 col array of rownames from data;
         each row represents a pair of rownames used to access a pair of
         rows from data. The mahal distance between the two rows will be
         calculated
       SEXP invScaleMat: matrix inverse of the covariates of the data;
         used to calc the mahal distance
     returns: SEXP results, a real vector (double precision) of distances
       between particular pairs of rows from data as indicated by index

       used .Call with SEXP arguments to reduce the memory overhead
       caused by, for example as.double, creating new objects to pass to
       .C

       The difficulty with this approach is that rows of data must be
       accessed using index in C instead of in R (where it's easy). A
       hash map (GNU libc extension) is used to hash the rownames of data
       to row indexes.
 */
// [[Rcpp::export]]
NumericVector mahalanobisHelper(NumericMatrix data, StringMatrix index,
				NumericMatrix invScaleMat) {
  int
    va, vb,
    nv = index.nrow();

  StringVector row_names = rownames(data);
  std::unordered_map<std::string, int> strpos;
  for (int i = 0; i < row_names.size(); i++) {
    std::string row_name_i = as< std::string >(row_names(i));
    strpos[row_name_i] = i;
  }

  // alloc space for the result
  NumericVector result(nv);

  // sqrt((va - vb) * invScaleMat * (va - vb))
  for (int i = 0; i < nv; i++) {
    std::string row_name_a = as< std::string >(index(i, 0));
    va = strpos[row_name_a];

    std::string row_name_b = as< std::string >(index(i, 1));
    vb = strpos[row_name_b];

    double sum_i = 0.0;
    for (int j = 0; j < data.ncol(); j++) {
      double innerSum = 0.0;
      for (int k = 0; k < data.ncol(); k++) {
	innerSum += (data(va, k) - data(vb, k)) * invScaleMat(j, k);
      }
      sum_i += innerSum * (data(va, j) - data(vb, j));
    }
    result(i) = sqrt(sum_i);
  }
  return result;
}
