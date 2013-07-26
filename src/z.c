#define _GNU_SOURCE

#include <limits.h>
#include <search.h>

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Lapack.h>

void rPrintV(int n, int stride, int start, const double * v) {
  for(int i = start; i < n; i += stride)
    Rprintf("%g ", v[i]);

  Rprintf("\n");
}

void rPrintMat(int rows, int cols, const double * mat) {
  for(int i = 0; i < rows; i++)
    rPrintV(cols, rows , i, mat);

  Rprintf("\n");
}

int get_name_index(const char * name, SEXP names) {
  int n = length(names);
  for(int i = 0; i < n; i++) {
    if( 0 == strcmp(CHAR(STRING_ELT(names, i)), name) )
      return i;
  }
  return -n;
}

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

void names_to_indexes(SEXP names, SEXP index_names, int * out_index) {
  int
    n_index = length(index_names),
    n_names = length(names),
    n_hash = ceil((4.0 * (double) n_names) / 3.0);

  struct hsearch_data * htab = Calloc(n_hash, struct hsearch_data);

  if( 0 == hcreate_r(n_hash, htab) )
    error("in names_to_indexes: failed to create hash to index data row names\npossibly out of memory\n");

  ENTRY
    to_find, * found,
    * entries = Calloc(n_index, ENTRY);

  for(int i = 0; i < n_names; i++) {
    entries[i].key = (char *) CHAR(STRING_ELT(names, i));
    entries[i].data = Calloc(1 + digits(i), char);
    sprintf(entries[i].data, "%d", i);
    if( 0 == hsearch_r(entries[i], ENTER, &found, htab) )
	error("in names_to_indexes: unknown error on inserting key into hash table; table might be full\n");
  }

  for(int i = 0; i < n_index; i++) {
    to_find.key = (char *) CHAR(STRING_ELT(index_names, i));
    if( 0 == hsearch_r(to_find, FIND, &found, htab) )
	error("in names_to_indexes: element of index not found among rownames\n");
    out_index[i] = (int) strtol(found->data, NULL, 0);
  }

  hdestroy_r(htab);
  for(int i = 0; i < n_index; i++)
    Free(entries[i].data);

  Free(entries);
  Free(htab);
}

SEXP z(SEXP data, SEXP data_row_names, SEXP index, SEXP invScaleMat) {
  int
    nv = nrows(index), n = ncols(data),
    data_rows = nrows(data);

  int
    va, vb,
    * index_row_i = Calloc(2 * nv, int);
  double
    * pairDiff = Calloc(nv * n, double),
    * product = Calloc(nv * n, double);

  names_to_indexes(data_row_names, index, index_row_i);
  for(int i = 0; i < nv; i++) {
    va = index_row_i[i];
    vb = index_row_i[i + nv];
    for(int j = 0; j < n; j++) {
      pairDiff[i + j * nv] = REAL(data)[va + j * data_rows]
	    - REAL(data)[vb + j * data_rows];
    }
  }
  Free(index_row_i);

  char
    side = 'R', uplo = 'U';
  double
    alpha = 1.0, beta = 0.0;

  // pairDiff * invScaleMat
  F77_CALL(dsymm)(&side, &uplo, &nv, &n, &alpha, REAL(invScaleMat), &n,
		  pairDiff, &nv, &beta, product, &nv);

  for(int i = 0; i < nv * n; i++)
    product[i] *= pairDiff[i];

  Free(pairDiff);

  SEXP result;
  PROTECT(result = allocVector(REALSXP, nv));
  for(int i = 0; i < nv; i++) {
    REAL(result)[i] = 0.0;
    for(int j = 0; j < n; j++)
      REAL(result)[i] += product[i + j * nv];
    
    REAL(result)[i] = sqrt(REAL(result)[i]);
  }

  Free(product);
  UNPROTECT(1);

  return result;
}
