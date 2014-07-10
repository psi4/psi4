#ifndef __GFUTEX_H
#define __GFUTEX_H

#include <cstddef>

namespace globalFutures {
  typedef void (*RemoteFuncProto) (int g_a, int ndims, int lo[], int hi[], void *arg);
  typedef size_t GFHandle;

  int GFInitialize();
  void GFFinalize();

  GFHandle GFRegister(RemoteFuncProto func, const size_t argSz, const size_t retSz = 0U);

  void GFEnqueue(GFHandle hndl, int g_a, int ndims, int lo[], int hi[], void *arg);
  void GFEnqueue(GFHandle hndl, int g_a, int ndims, int lo[], int hi[], void *iarg, void **&oarg);

  void GFExecute(GFHandle hndl, int g_a, int ndims, int lo[], int hi[], void *arg);
  void GFExecute(GFHandle hndl, int g_a, int ndims, int lo[], int hi[], void *iarg, void **&oarg);

  int GFMaxConcurrency();

  void GFQuiesce();
  void GFAllQuiesce();

  void GFQuiesce(GFHandle hndl);
  void GFAllQuiesce(GFHandle hndl);

  // Overloaded GA functions for thread-safety

  void GF_Access(int g_a, int lo[], int hi[], void *ptr, int ld[]);
  void GF_Release(int g_a, int lo[], int hi[]);
  void GF_Get(int g_a, int lo[], int hi[], void *buf, int ld[]);
  int GF_Locate_region(int g_a, int lo[], int hi[], int map[], int procs[]);
  void GF_Distribution(int g_a, int iproc, int lo[], int hi[]);
  void GF_Acc(int g_a, int lo[], int hi[], void *buf, int ld[], void *alpha);
  void GF_NbAcc(int g_a, int lo[], int hi[], void *buf, int ld[], void *alpha,
		ga_nbhdl_t *nbhandle);
  void GF_NbGet(int g_a, int lo[], int hi[], void *buf, int ld[], ga_nbhdl_t *nbhandle);
  void GF_NbPut(int g_a, int lo[], int hi[], void *buf, int ld[], ga_nbhdl_t *nbhandle);
  void GF_NbWait(ga_nbhdl_t *nbhandle);

  int GF_Ndim(int g_a);

  int GF_Cluster_nodeid();
  int GF_Cluster_proc_nodeid(int proc);

  void GF_Inquire(int g_a, int *type, int *ndim, int dims[]);

  // Caching functions for remote data

  void GF_CachedGet(int g_a, int lo[], int hi[], void *buf, int ld[]);
  void GF_CachedAcc(int g_a, int lo[], int hi[], void *buf, int ld[], void *alpha);
  void GF_CachedNbGet(int g_a, int lo[], int hi[], void *buf, int ld[], ga_nbhdl_t *nbhandle);
  void GF_CachedNbWait(ga_nbhdl_t *nbhandle);

  void GF_CacheReadOnlyEmpty(int g_a);
  void GF_CacheReadWriteFlush(int g_a);
  void GF_CacheAccFlush(int g_a);

  unsigned int GF_CacheGetMisses(int g_a);
  unsigned int GF_CacheGetHits(int g_a);
}

#endif // __GFUTEX_H
