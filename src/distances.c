#define _GNU_SOURCE // for the hash map gnu extension to libc

#include <limits.h> // for INT_MIN
#include <search.h> // for gnu extension hash map functions

#include"register.h" // to register the mahalanobisHelper with R

/* function: digits
   computes the number characters required to store an integer base 10
   as a string. It will include a character for the negative sign. It
   will not include a character for the string termination symbol.

   This function should be much faster than the usual call to a
   logarithm function.
 */
int digits(int n) {
  if(n == INT_MIN) return 11;
  if(n < 0) return(1 + digits(-n));

  if(n >= 10000) {
    if(n >= 10000000) {
      if(n >= 100000000) {
	if(n >= 1000000000) return 10;
	return 9;
      }
      return 8;
    }
    if(n >= 100000) {
      if(n >= 1000000) return 7;
      return 6;
    }
    return 5;
  }
  if(n >= 100) {
    if(n >= 1000) return 4;
    return 3;
  }
  if(n >= 10) return 2;
  return 1;
}

/* the MAP type

   This struct was implemented to encapsulate the details of a hash map
   from a rowname label to a row index for an R dataframe, array or matrix.

   The GNU extension hash map was used so that the space for the hash map
   could be alloced from R's heap. So one needs to eventually Free the 
   entries as well as the hash table. This is why the struct is needed.
 */
typedef struct map {
  struct hsearch_data * hash_tab;
  ENTRY * entries;
  size_t n_entries;
} MAP;

/* function: create_map
   This function accepts a character vector SEXP and returns a hash map
   of a string in the character vector to it's position in the vector strs.
   The position is converted to a string for storage. The storage for the
   entries and the map are allocated using R's Calloc. The storage must be
   Freed using the delete_map function. The delete map function will Free
   each entry as well as the hash table.

   If SEXP strs is not a character vector, bad things will happen.
*/
MAP * create_map(SEXP strs) {

  // TODO: check to see that strs is a character vector

  int
    n_strs = length(strs),
    n_map = ceil((4.0 * (double) n_strs) / 3.0);
  // hash must be >= 20% empty for efficient look up

  MAP * strpos = Calloc(1, MAP);
  strpos->hash_tab = Calloc(n_map, struct hsearch_data);
  if( 0 == hcreate_r(n_map, strpos->hash_tab) )
    error("In create_strpos: Failed to create hash table. Out of memory?");

  strpos->entries = Calloc(n_strs, ENTRY);
  strpos->n_entries = n_strs;

  ENTRY * inserted;
  for(int i = 0; i < n_strs; i++) {
    strpos->entries[i].key = (char *) CHAR(STRING_ELT(strs, i));
    strpos->entries[i].data = Calloc(1 + digits(i), char);
    sprintf(strpos->entries[i].data, "%d", i);
    if( 0 == hsearch_r(strpos->entries[i], ENTER, &inserted, strpos->hash_tab) )
	error("In create_strpos: Can't insert key. Table full?");
  }
  return strpos;
}

/* function: delete_map
   Free's each entry of the MAP pointer strpos as well as the hash table.
   The storage must have been allocated with R's Calloc or bad things will
   happen.
 */
void delete_map(MAP * strpos) {

  // TODO: check for NULL pointers (strpos members too)

  hdestroy_r(strpos->hash_tab);
  for(int i = 0; i < strpos->n_entries; i++)
    Free(strpos->entries[i].data);

  Free(strpos->entries);
  Free(strpos->hash_tab);
  Free(strpos);
}

/* function: get_pos
   consumes: a string to_find, a MAP pointer strpos allocated by
     create_map; if strpos has not been eaten by create_map first,
     bad things will happen
   returns: an integer resulting from applying the hash map defined
     by strpos to to_find. The integer will be returned from the GNU
     hsearch_r function as a string and will need to be converted to
     an integer before this function returns
 */
int get_pos(const char * to_find, MAP * strpos) {

  // TODO: check for NULL pointers (strpos members too)

  ENTRY to_find_e, * found;

  // ENTRY's key is not const but R character vectors are always const
  // so we need the cast to avoid compiler warnings.
  to_find_e.key = (char *) to_find;
  if( 0 == hsearch_r(to_find_e, FIND, &found, strpos->hash_tab) )
    error("In get_pos: String not found.");

  // convert hashed string to long and return
  return strtol(found->data, NULL, 0);
}

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
