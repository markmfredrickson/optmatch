#include <R_ext/Rdynload.h>

#ifndef _REGISTER_H_
#define _REGISTER_H_

extern "C" {
  void R_init_optmatch(DllInfo *info);
}

#endif
