#include <Rcpp.h>

#ifndef _ISM_H_
#define _ISM_H_

extern "C" {
  SEXP ismOps(SEXP o, SEXP a, SEXP b);
}

#endif
