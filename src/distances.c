#define _GNU_SOURCE

#include <limits.h>
#include <search.h>

#include"register.h"

#include <R_ext/Lapack.h>

int digits(int n) {
  if(n == INT_MIN) return 11;
  if(n < 0) 1 + digits(-n);

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

typedef struct map {
  struct hsearch_data * hash_tab;
  ENTRY * entries;
  size_t n_entries;
} MAP;

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

void delete_map(MAP * strpos) {

  // TODO: check for NULL pointers (strpos members too)

  hdestroy_r(strpos->hash_tab);
  for(int i = 0; i < strpos->n_entries; i++)
    Free(strpos->entries[i].data);

  Free(strpos->entries);
  Free(strpos->hash_tab);
  Free(strpos);
}

int get_pos(const char * to_find, MAP * strpos) {

  // TODO: check for NULL pointers (strpos members too)

  ENTRY to_find_e, * found;
  to_find_e.key = (char *) to_find;
  if( 0 == hsearch_r(to_find_e, FIND, &found, strpos->hash_tab) )
    error("In get_pos: String not found.");

  return strtol(found->data, NULL, 0);
}

SEXP mahalanobisHelper(SEXP data, SEXP index, SEXP invScaleMat) {
  int
    j, k,
    va, vb,
    nv = nrows(index),
    n = ncols(data), nr = nrows(data);
  double
    sum, innerSum;

  SEXP row_names, cl;
  const char * rn, * cn;
  GetMatrixDimnames(data, &row_names, &cl, &rn, &cn); 
  MAP * strpos = create_map(row_names);

  SEXP result;
  PROTECT(result = allocVector(REALSXP, nv));
  double
    * real_data = REAL(data),
    * mj, * real_m = REAL(invScaleMat),
    * real_result = REAL(result);

  for(int i = 0; i < nv; i++) {
    sum = 0.0;
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
  delete_map(strpos);
  UNPROTECT(1);
  return result;
}
