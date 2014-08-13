
#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <stdio.h>

#include "armci.h"
#include "parmci.h"


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_Acc
#endif
int ARMCI_Acc(int optype, void *scale, void *src, void *dst, int bytes, int proc)
{
    return PARMCI_Acc(optype, scale, src, dst, bytes, proc);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_AccS
#endif
int ARMCI_AccS(int optype, void *scale, void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc)
{
    return PARMCI_AccS(optype, scale, src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, proc);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_AccV
#endif
int ARMCI_AccV(int op, void *scale, armci_giov_t *darr, int len, int proc)
{
    return PARMCI_AccV(op, scale, darr, len, proc);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_AllFence
#endif
void ARMCI_AllFence()
{
    PARMCI_AllFence();
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_Barrier
#endif
void ARMCI_Barrier()
{
    PARMCI_Barrier();
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_Create_mutexes
#endif
int ARMCI_Create_mutexes(int num)
{
    return PARMCI_Create_mutexes(num);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_Destroy_mutexes
#endif
int ARMCI_Destroy_mutexes()
{
    return PARMCI_Destroy_mutexes();
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_Fence
#endif
void ARMCI_Fence(int proc)
{
    PARMCI_Fence(proc);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_Finalize
#endif
void ARMCI_Finalize()
{
    PARMCI_Finalize();
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_Free
#endif
int ARMCI_Free(void *ptr)
{
    return PARMCI_Free(ptr);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_Free_local
#endif
int ARMCI_Free_local(void *ptr)
{
    return PARMCI_Free_local(ptr);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_Get
#endif
int ARMCI_Get(void *src, void *dst, int bytes, int proc)
{
    return PARMCI_Get(src, dst, bytes, proc);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_GetS
#endif
int ARMCI_GetS(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc)
{
    return PARMCI_GetS(src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, proc);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_GetV
#endif
int ARMCI_GetV(armci_giov_t *darr, int len, int proc)
{
    return PARMCI_GetV(darr, len, proc);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_GetValueDouble
#endif
double ARMCI_GetValueDouble(void *src, int proc)
{
    return PARMCI_GetValueDouble(src, proc);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_GetValueFloat
#endif
float ARMCI_GetValueFloat(void *src, int proc)
{
    return PARMCI_GetValueFloat(src, proc);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_GetValueInt
#endif
int ARMCI_GetValueInt(void *src, int proc)
{
    return PARMCI_GetValueInt(src, proc);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_GetValueLong
#endif
long ARMCI_GetValueLong(void *src, int proc)
{
    return PARMCI_GetValueLong(src, proc);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_Init
#endif
int ARMCI_Init()
{
    return PARMCI_Init();
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_Init_args
#endif
int ARMCI_Init_args(int *argc, char ***argv)
{
    return PARMCI_Init_args(argc, argv);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_Initialized
#endif
int ARMCI_Initialized()
{
    return PARMCI_Initialized();
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_Lock
#endif
void ARMCI_Lock(int mutex, int proc)
{
    PARMCI_Lock(mutex, proc);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_Malloc
#endif
int ARMCI_Malloc(void **ptr_arr, armci_size_t bytes)
{
    return PARMCI_Malloc(ptr_arr, bytes);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_Malloc_local
#endif
void* ARMCI_Malloc_local(armci_size_t bytes)
{
    return PARMCI_Malloc_local(bytes);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_Memat
#endif
void* ARMCI_Memat(armci_meminfo_t *meminfo, long offset)
{
    return PARMCI_Memat(meminfo, offset);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_Memget
#endif
void ARMCI_Memget(size_t bytes, armci_meminfo_t *meminfo, int memflg)
{
    PARMCI_Memget(bytes, meminfo, memflg);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_NbAccS
#endif
int ARMCI_NbAccS(int optype, void *scale, void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc, armci_hdl_t *nb_handle)
{
    return PARMCI_NbAccS(optype, scale, src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, proc, nb_handle);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_NbAccV
#endif
int ARMCI_NbAccV(int op, void *scale, armci_giov_t *darr, int len, int proc, armci_hdl_t *nb_handle)
{
    return PARMCI_NbAccV(op, scale, darr, len, proc, nb_handle);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_NbGet
#endif
int ARMCI_NbGet(void *src, void *dst, int bytes, int proc, armci_hdl_t *nb_handle)
{
    return PARMCI_NbGet(src, dst, bytes, proc, nb_handle);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_NbGetS
#endif
int ARMCI_NbGetS(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc, armci_hdl_t *nb_handle)
{
    return PARMCI_NbGetS(src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, proc, nb_handle);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_NbGetV
#endif
int ARMCI_NbGetV(armci_giov_t *darr, int len, int proc, armci_hdl_t *nb_handle)
{
    return PARMCI_NbGetV(darr, len, proc, nb_handle);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_NbPut
#endif
int ARMCI_NbPut(void *src, void *dst, int bytes, int proc, armci_hdl_t *nb_handle)
{
    return PARMCI_NbPut(src, dst, bytes, proc, nb_handle);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_NbPutS
#endif
int ARMCI_NbPutS(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc, armci_hdl_t *nb_handle)
{
    return PARMCI_NbPutS(src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, proc, nb_handle);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_NbPutV
#endif
int ARMCI_NbPutV(armci_giov_t *darr, int len, int proc, armci_hdl_t *nb_handle)
{
    return PARMCI_NbPutV(darr, len, proc, nb_handle);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_NbPutValueDouble
#endif
int ARMCI_NbPutValueDouble(double src, void *dst, int proc, armci_hdl_t *nb_handle)
{
    return PARMCI_NbPutValueDouble(src, dst, proc, nb_handle);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_NbPutValueFloat
#endif
int ARMCI_NbPutValueFloat(float src, void *dst, int proc, armci_hdl_t *nb_handle)
{
    return PARMCI_NbPutValueFloat(src, dst, proc, nb_handle);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_NbPutValueInt
#endif
int ARMCI_NbPutValueInt(int src, void *dst, int proc, armci_hdl_t *nb_handle)
{
    return PARMCI_NbPutValueInt(src, dst, proc, nb_handle);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_NbPutValueLong
#endif
int ARMCI_NbPutValueLong(long src, void *dst, int proc, armci_hdl_t *nb_handle)
{
    return PARMCI_NbPutValueLong(src, dst, proc, nb_handle);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_Put
#endif
int ARMCI_Put(void *src, void *dst, int bytes, int proc)
{
    return PARMCI_Put(src, dst, bytes, proc);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_PutS
#endif
int ARMCI_PutS(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc)
{
    return PARMCI_PutS(src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, proc);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_PutS_flag
#endif
int ARMCI_PutS_flag(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int *flag, int val, int proc)
{
    return PARMCI_PutS_flag(src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, flag, val, proc);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_PutS_flag_dir
#endif
int ARMCI_PutS_flag_dir(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int *flag, int val, int proc)
{
    return PARMCI_PutS_flag_dir(src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, flag, val, proc);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_PutV
#endif
int ARMCI_PutV(armci_giov_t *darr, int len, int proc)
{
    return PARMCI_PutV(darr, len, proc);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_PutValueDouble
#endif
int ARMCI_PutValueDouble(double src, void *dst, int proc)
{
    return PARMCI_PutValueDouble(src, dst, proc);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_PutValueFloat
#endif
int ARMCI_PutValueFloat(float src, void *dst, int proc)
{
    return PARMCI_PutValueFloat(src, dst, proc);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_PutValueInt
#endif
int ARMCI_PutValueInt(int src, void *dst, int proc)
{
    return PARMCI_PutValueInt(src, dst, proc);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_PutValueLong
#endif
int ARMCI_PutValueLong(long src, void *dst, int proc)
{
    return PARMCI_PutValueLong(src, dst, proc);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_Put_flag
#endif
int ARMCI_Put_flag(void *src, void *dst, int bytes, int *f, int v, int proc)
{
    return PARMCI_Put_flag(src, dst, bytes, f, v, proc);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_Rmw
#endif
int ARMCI_Rmw(int op, void *ploc, void *prem, int extra, int proc)
{
    return PARMCI_Rmw(op, ploc, prem, extra, proc);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_Test
#endif
int ARMCI_Test(armci_hdl_t *nb_handle)
{
    return PARMCI_Test(nb_handle);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_Unlock
#endif
void ARMCI_Unlock(int mutex, int proc)
{
    PARMCI_Unlock(mutex, proc);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_Wait
#endif
int ARMCI_Wait(armci_hdl_t *nb_handle)
{
    return PARMCI_Wait(nb_handle);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_WaitAll
#endif
int ARMCI_WaitAll()
{
    return PARMCI_WaitAll();
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak ARMCI_WaitProc
#endif
int ARMCI_WaitProc(int proc)
{
    return PARMCI_WaitProc(proc);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak armci_msg_barrier
#endif
void armci_msg_barrier()
{
    parmci_msg_barrier();
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak armci_msg_group_barrier
#endif
void armci_msg_group_barrier(ARMCI_Group *group)
{
    parmci_msg_group_barrier(group);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak armci_notify
#endif
int armci_notify(int proc)
{
    return parmci_notify(proc);
}


#if HAVE_SYS_WEAK_ALIAS_PRAGMA
#   pragma weak armci_notify_wait
#endif
int armci_notify_wait(int proc, int *pval)
{
    return parmci_notify_wait(proc, pval);
}

