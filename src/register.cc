#include <R.h>

#include "register.h"
#include "distances.h"
#include "subsetInfSparseMatrix.h"
#include "r_smahal.h"
#include "ism.h"

void R_init_optmatch(DllInfo *info) {
  R_CallMethodDef callMethods[]  = {
    {"mahalanobisHelper", (DL_FUNC) &mahalanobisHelper, 3},
    {"subsetInfSparseMatrix", (DL_FUNC) &subsetInfSparseMatrix, 3},
    {"r_smahal", (DL_FUNC) &r_smahal, 3},
    {"ismOps", (DL_FUNC) &ismOps, 3},
    {NULL, NULL, 0}
  };
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
}
