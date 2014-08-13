#ifndef __TWOEL_H
#define __TWOEL_H

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

  int GFInitialize();
  void GFFinalize();

  void twoel(double schwmax, double *etwo, int nproc);
  void twoel_orig(double schwmax, double *etwo, int nproc);

#ifdef __cplusplus
}
#endif

#undef CLOG

#ifdef CLOG
extern FILE *clog;
#endif

#endif /* __TWOEL_H */
