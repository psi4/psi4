#ifndef __SCF_H
#define __SCF_H

#ifdef __cplusplus
extern "C" {
#endif

  int next_4chunk(int *lo, int *hi, int *ilo, int *jlo, int *klo, int *llo, long int *pitask);
  void clean_chunk(double chunk[][ichunk]);
  void g(double *value, int i, int j, int k, int l);
  double contract_matrices(int g_a, int g_b);

  long int acquire_tasks(int numTasks);
  int translate_task(long int itask, int *lo, int *hi, int *ilo, int *jlo, int *klo, int *llo);

  int next_chunk(int *lo, int *hi);
  double h(int i, int j);

#ifdef __cplusplus
}
#endif

#endif /* __SCF_H */
