#ifndef _PARMCI_H_
#define _PARMCI_H_

#include "armci.h"

extern int    PARMCI_Acc(int optype, void *scale, void *src, void* dst, int bytes, int proc);
extern int    PARMCI_AccS(int optype, void *scale, void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc);
extern int    PARMCI_AccV(int op, void *scale, armci_giov_t * darr, int len, int proc);
extern void   PARMCI_AllFence();
extern void   PARMCI_Barrier();
extern int    PARMCI_Create_mutexes(int num);
extern int    PARMCI_Destroy_mutexes();
extern void   PARMCI_Fence(int proc);
extern void   PARMCI_Finalize();
extern int    PARMCI_Free_local(void *ptr);
extern int    PARMCI_Free(void *ptr);
extern int    PARMCI_GetS(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc);
extern double PARMCI_GetValueDouble(void *src, int proc);
extern float  PARMCI_GetValueFloat(void *src, int proc);
extern int    PARMCI_GetValueInt(void *src, int proc);
extern long   PARMCI_GetValueLong(void *src, int proc);
extern int    PARMCI_GetV(armci_giov_t * darr, int len, int proc);
extern int    PARMCI_Get(void *src, void *dst, int bytes, int proc);
extern int    PARMCI_Init();
extern int    PARMCI_Init_args(int *argc, char ***argv);
extern int    PARMCI_Initialized();
extern void   PARMCI_Lock(int mutex, int proc);
extern void*  PARMCI_Malloc_local(armci_size_t bytes);
extern int    PARMCI_Malloc(void **ptr_arr, armci_size_t bytes);
extern void*  PARMCI_Memat(armci_meminfo_t * meminfo, long offset);
extern void   PARMCI_Memget(size_t bytes, armci_meminfo_t * meminfo, int memflg);
extern int    PARMCI_NbAccS(int optype, void *scale, void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc, armci_hdl_t * nb_handle);
extern int    PARMCI_NbAccV(int op, void *scale, armci_giov_t * darr, int len, int proc, armci_hdl_t * nb_handle);
extern int    PARMCI_NbGetS(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc, armci_hdl_t * nb_handle);
extern int    PARMCI_NbGetV(armci_giov_t * darr, int len, int proc, armci_hdl_t * nb_handle);
extern int    PARMCI_NbGet(void *src, void *dst, int bytes, int proc, armci_hdl_t * nb_handle);
extern int    PARMCI_NbPutS(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc, armci_hdl_t * nb_handle);
extern int    PARMCI_NbPutValueDouble(double src, void *dst, int proc, armci_hdl_t* nb_handle);
extern int    PARMCI_NbPutValueFloat(float src, void *dst, int proc, armci_hdl_t* nb_handle);
extern int    PARMCI_NbPutValueInt(int src, void *dst, int proc, armci_hdl_t* nb_handle);
extern int    PARMCI_NbPutValueLong(long src, void *dst, int proc, armci_hdl_t* nb_handle);
extern int    PARMCI_NbPutV(armci_giov_t * darr, int len, int proc, armci_hdl_t * nb_handle);
extern int    PARMCI_NbPut(void *src, void *dst, int bytes, int proc, armci_hdl_t * nb_handle);
extern int    PARMCI_Put_flag(void *src, void *dst, int bytes, int *f, int v, int proc);
extern int    PARMCI_PutS_flag_dir(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int *flag, int val, int proc);
extern int    PARMCI_PutS_flag(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int *flag, int val, int proc);
extern int    PARMCI_PutS(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc);
extern int    PARMCI_PutValueDouble(double src, void *dst, int proc);
extern int    PARMCI_PutValueFloat(float src, void *dst, int proc);
extern int    PARMCI_PutValueInt(int src, void *dst, int proc);
extern int    PARMCI_PutValueLong(long src, void *dst, int proc);
extern int    PARMCI_PutV(armci_giov_t * darr, int len, int proc);
extern int    PARMCI_Put(void *src, void *dst, int bytes, int proc);
extern int    PARMCI_Rmw(int op, void *ploc, void *prem, int extra, int proc);
extern int    PARMCI_Test(armci_hdl_t * nb_handle);
extern void   PARMCI_Unlock(int mutex, int proc);
extern int    PARMCI_WaitAll();
extern int    PARMCI_Wait(armci_hdl_t * nb_handle);
extern int    PARMCI_WaitProc(int proc);
extern void   parmci_msg_barrier();
extern void   parmci_msg_group_barrier(ARMCI_Group *group);
extern int    parmci_notify(int proc);
extern int    parmci_notify_wait(int proc, int *pval);

#endif /* _PARMCI_H_ */
