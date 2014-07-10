#ifndef __ONEEL_H
#define __ONEEL_H

#ifdef __cplusplus
extern "C" {
#endif

  void oneel(double schwmax, double *eone, int nproc);
  void oneel_orig(double schwmax, double *eone, int nproc);

#ifdef __cplusplus
}
#endif

#endif /* __ONEEL_H */
