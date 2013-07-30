#include "register.h"
#include <R_ext/Rdynload.h>

void R_init_optmatch(DllInfo *info) {
  R_CallMethodDef callMethods[]  = {
    {"mahalanobisHelper", (DL_FUNC) &mahalanobisHelper, 4},
    {"subsetInfSparseMatrix", (DL_FUNC) &subsetInfSparseMatrix, 3},
    {NULL, NULL, 0}
  };
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
