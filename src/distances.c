#include <limits.h> // for INT_MIN

#include "map.h" // hash map functions for rowname to index
#include "distances.h" // to register the mahalanobisHelper with R

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
SEXP mahalanobisHelper(SEXP data, SEXP index, SEXP invScaleMat) {
  int
    j, k,
    va, vb,
    nv = nrows(index),
    n = ncols(data), nr = nrows(data);
  double
    sum, innerSum;

  // get the rownames, access through dimnames only worked for me
  // tried RownamesSymbol and GetRownames, neither worked
  // we also get colnames and headers for the two dimensions
  // which are not used
  SEXP row_names, cl;
  const char * rn, * cn;
  GetMatrixDimnames(data, &row_names, &cl, &rn, &cn); 
  MAP * strpos = create_map(row_names);

  // alloc space for the result
  SEXP result;
  PROTECT(result = allocVector(REALSXP, nv));

  // get pointers to C types; makes the loop easier to read
  double
    * real_data = REAL(data),
    * mj, * real_m = REAL(invScaleMat),
    * real_result = REAL(result);

  // sqrt((va - vb) * invScaleMat * (va - vb))
  for(int i = 0; i < nv; i++) {
    sum = 0.0;
    // hash col_a rowname and col_b rowname into row indexes
    va = get_pos(CHAR(STRING_ELT(index, i)), strpos);
    vb = get_pos(CHAR(STRING_ELT(index, i + nv)), strpos);

    for(j = 0; j < n; j++) {
      innerSum = 0.0;
      mj = real_m + j * n;
      for(k = 0; k < n; k++) {
	innerSum += (real_data[va + k * nr]
		     - real_data[vb + k * nr])
	            * mj[k];
      }
      sum += innerSum * (real_data[va + j * nr] - real_data[vb + j * nr]);
    }
    real_result[i] = sqrt(sum);
  }
  
  // free the storage for the rowname hash
  delete_map(strpos);

  // return R object
  UNPROTECT(1);
  return result;
}
