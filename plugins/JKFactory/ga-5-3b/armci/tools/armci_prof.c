
#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <stdio.h>

#include "armci.h"
#include "parmci.h"
#include "armci_profile.h"
#include "armci_profile.c"


int ARMCI_Acc(int optype, void *scale, void *src, void *dst, int bytes, int proc)
{
    int ret;
    armci_profile_start_strided(&bytes, 0, proc, ARMCI_PROF_ACC);
    ret = PARMCI_Acc(optype, scale, src, dst, bytes, proc);
    armci_profile_stop_strided(ARMCI_PROF_ACC);
    return ret;
}


int ARMCI_AccS(int optype, void *scale, void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc)
{
    int ret;
    armci_profile_start_strided(count, stride_levels, proc, ARMCI_PROF_ACCS);
    ret = PARMCI_AccS(optype, scale, src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, proc);
    armci_profile_stop_strided(ARMCI_PROF_ACCS);
    return ret;
}


int ARMCI_AccV(int op, void *scale, armci_giov_t *darr, int len, int proc)
{
    int ret;
    armci_profile_start_vector(darr, len, proc, ARMCI_PROF_ACCV);
    ret = PARMCI_AccV(op, scale, darr, len, proc);
    armci_profile_stop_vector(ARMCI_PROF_ACCV);
    return ret;
}


void ARMCI_AllFence()
{
    
    armci_profile_start(ARMCI_PROF_ALLFENCE);
    PARMCI_AllFence();
    armci_profile_stop(ARMCI_PROF_ALLFENCE);
    
}


void ARMCI_Barrier()
{
    
    armci_profile_start(ARMCI_PROF_BARRIER);
    PARMCI_Barrier();
    armci_profile_stop(ARMCI_PROF_BARRIER);
    
}


int ARMCI_Create_mutexes(int num)
{
    int ret;
    ret = PARMCI_Create_mutexes(num);
    return ret;
}


int ARMCI_Destroy_mutexes()
{
    int ret;
    ret = PARMCI_Destroy_mutexes();
    return ret;
}


void ARMCI_Fence(int proc)
{
    if (!SAMECLUSNODE(proc))
    armci_profile_start(ARMCI_PROF_FENCE);
    PARMCI_Fence(proc);
    if (!SAMECLUSNODE(proc))
    armci_profile_stop(ARMCI_PROF_FENCE);
}


void ARMCI_Finalize()
{
    armci_profile_terminate();
    PARMCI_Finalize();
}


int ARMCI_Free(void *ptr)
{
    int ret;
    ret = PARMCI_Free(ptr);
    return ret;
}


int ARMCI_Free_local(void *ptr)
{
    int ret;
    ret = PARMCI_Free_local(ptr);
    return ret;
}


int ARMCI_Get(void *src, void *dst, int bytes, int proc)
{
    int ret;
    armci_profile_start_strided(&bytes, 0, proc, ARMCI_PROF_GET);
    ret = PARMCI_Get(src, dst, bytes, proc);
    armci_profile_stop_strided(ARMCI_PROF_GET);
    return ret;
}


int ARMCI_GetS(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc)
{
    int ret;
    armci_profile_start_strided(count, stride_levels, proc, ARMCI_PROF_GETS);
    ret = PARMCI_GetS(src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, proc);
    armci_profile_stop_strided(ARMCI_PROF_GETS);
    return ret;
}


int ARMCI_GetV(armci_giov_t *darr, int len, int proc)
{
    int ret;
    armci_profile_start_vector(darr, len, proc, ARMCI_PROF_GETV);
    ret = PARMCI_GetV(darr, len, proc);
    armci_profile_stop_vector(ARMCI_PROF_GETV);
    return ret;
}


double ARMCI_GetValueDouble(void *src, int proc)
{
    double ret;
    ret = PARMCI_GetValueDouble(src, proc);
    return ret;
}


float ARMCI_GetValueFloat(void *src, int proc)
{
    float ret;
    ret = PARMCI_GetValueFloat(src, proc);
    return ret;
}


int ARMCI_GetValueInt(void *src, int proc)
{
    int ret;
    ret = PARMCI_GetValueInt(src, proc);
    return ret;
}


long ARMCI_GetValueLong(void *src, int proc)
{
    long ret;
    ret = PARMCI_GetValueLong(src, proc);
    return ret;
}


int ARMCI_Init()
{
    int ret;
    ret = PARMCI_Init();
    armci_profile_init();
    return ret;
}


int ARMCI_Init_args(int *argc, char ***argv)
{
    int ret;
    ret = PARMCI_Init_args(argc, argv);
    armci_profile_init();
    return ret;
}


void ARMCI_Lock(int mutex, int proc)
{
    
    PARMCI_Lock(mutex, proc);
    
}


int ARMCI_Malloc(void **ptr_arr, armci_size_t bytes)
{
    int ret;
    ret = PARMCI_Malloc(ptr_arr, bytes);
    return ret;
}


void* ARMCI_Malloc_local(armci_size_t bytes)
{
    void* ret;
    ret = PARMCI_Malloc_local(bytes);
    return ret;
}


void* ARMCI_Memat(armci_meminfo_t *meminfo, long offset)
{
    void* ret;
    ret = PARMCI_Memat(meminfo, offset);
    return ret;
}


void ARMCI_Memget(size_t bytes, armci_meminfo_t *meminfo, int memflg)
{
    
    PARMCI_Memget(bytes, meminfo, memflg);
    
}


int ARMCI_NbAccS(int optype, void *scale, void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc, armci_hdl_t *nb_handle)
{
    int ret;
    armci_profile_start_strided(count, stride_levels, proc, ARMCI_PROF_NBACCS);
    ret = PARMCI_NbAccS(optype, scale, src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, proc, nb_handle);
    armci_profile_stop_strided(ARMCI_PROF_NBACCS);
    return ret;
}


int ARMCI_NbAccV(int op, void *scale, armci_giov_t *darr, int len, int proc, armci_hdl_t *nb_handle)
{
    int ret;
    armci_profile_start_vector(darr, len, proc, ARMCI_PROF_NBACCV);
    ret = PARMCI_NbAccV(op, scale, darr, len, proc, nb_handle);
    armci_profile_stop_vector(ARMCI_PROF_NBACCV);
    return ret;
}


int ARMCI_NbGet(void *src, void *dst, int bytes, int proc, armci_hdl_t *nb_handle)
{
    int ret;
    armci_profile_start_strided(&bytes, 0, proc, ARMCI_PROF_NBGET);
    ret = PARMCI_NbGet(src, dst, bytes, proc, nb_handle);
    armci_profile_stop_strided(ARMCI_PROF_NBGET);
    return ret;
}


int ARMCI_NbGetS(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc, armci_hdl_t *nb_handle)
{
    int ret;
    armci_profile_start_strided(count, stride_levels, proc, ARMCI_PROF_NBGETS);
    ret = PARMCI_NbGetS(src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, proc, nb_handle);
    armci_profile_stop_strided(ARMCI_PROF_NBGETS);
    return ret;
}


int ARMCI_NbGetV(armci_giov_t *darr, int len, int proc, armci_hdl_t *nb_handle)
{
    int ret;
    armci_profile_start_vector(darr, len, proc, ARMCI_PROF_NBGETV);
    ret = PARMCI_NbGetV(darr, len, proc, nb_handle);
    armci_profile_stop_vector(ARMCI_PROF_NBGETV);
    return ret;
}


int ARMCI_NbPut(void *src, void *dst, int bytes, int proc, armci_hdl_t *nb_handle)
{
    int ret;
    armci_profile_start_strided(&bytes, 0, proc, ARMCI_PROF_NBPUT);
    ret = PARMCI_NbPut(src, dst, bytes, proc, nb_handle);
    armci_profile_stop_strided(ARMCI_PROF_NBPUT);
    return ret;
}


int ARMCI_NbPutS(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc, armci_hdl_t *nb_handle)
{
    int ret;
    armci_profile_start_strided(count, stride_levels, proc, ARMCI_PROF_NBPUTS);
    ret = PARMCI_NbPutS(src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, proc, nb_handle);
    armci_profile_stop_strided(ARMCI_PROF_NBPUTS);
    return ret;
}


int ARMCI_NbPutV(armci_giov_t *darr, int len, int proc, armci_hdl_t *nb_handle)
{
    int ret;
    armci_profile_start_vector(darr, len, proc, ARMCI_PROF_NBPUTV);
    ret = PARMCI_NbPutV(darr, len, proc, nb_handle);
    armci_profile_stop_vector(ARMCI_PROF_NBPUTV);
    return ret;
}


int ARMCI_NbPutValueDouble(double src, void *dst, int proc, armci_hdl_t *nb_handle)
{
    int ret;
    ret = PARMCI_NbPutValueDouble(src, dst, proc, nb_handle);
    return ret;
}


int ARMCI_NbPutValueFloat(float src, void *dst, int proc, armci_hdl_t *nb_handle)
{
    int ret;
    ret = PARMCI_NbPutValueFloat(src, dst, proc, nb_handle);
    return ret;
}


int ARMCI_NbPutValueInt(int src, void *dst, int proc, armci_hdl_t *nb_handle)
{
    int ret;
    ret = PARMCI_NbPutValueInt(src, dst, proc, nb_handle);
    return ret;
}


int ARMCI_NbPutValueLong(long src, void *dst, int proc, armci_hdl_t *nb_handle)
{
    int ret;
    ret = PARMCI_NbPutValueLong(src, dst, proc, nb_handle);
    return ret;
}


int ARMCI_Put(void *src, void *dst, int bytes, int proc)
{
    int ret;
    armci_profile_start_strided(&bytes, 0, proc, ARMCI_PROF_PUT);
    ret = PARMCI_Put(src, dst, bytes, proc);
    armci_profile_stop_strided(ARMCI_PROF_PUT);
    return ret;
}


int ARMCI_PutS(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc)
{
    int ret;
    armci_profile_start_strided(count, stride_levels, proc, ARMCI_PROF_PUTS);
    ret = PARMCI_PutS(src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, proc);
    armci_profile_stop_strided(ARMCI_PROF_PUTS);
    return ret;
}


int ARMCI_PutS_flag(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int *flag, int val, int proc)
{
    int ret;
    ret = PARMCI_PutS_flag(src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, flag, val, proc);
    return ret;
}


int ARMCI_PutS_flag_dir(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int *flag, int val, int proc)
{
    int ret;
    ret = PARMCI_PutS_flag_dir(src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, flag, val, proc);
    return ret;
}


int ARMCI_PutV(armci_giov_t *darr, int len, int proc)
{
    int ret;
    armci_profile_start_vector(darr, len, proc, ARMCI_PROF_PUTV);
    ret = PARMCI_PutV(darr, len, proc);
    armci_profile_stop_vector(ARMCI_PROF_PUTV);
    return ret;
}


int ARMCI_PutValueDouble(double src, void *dst, int proc)
{
    int ret;
    ret = PARMCI_PutValueDouble(src, dst, proc);
    return ret;
}


int ARMCI_PutValueFloat(float src, void *dst, int proc)
{
    int ret;
    ret = PARMCI_PutValueFloat(src, dst, proc);
    return ret;
}


int ARMCI_PutValueInt(int src, void *dst, int proc)
{
    int ret;
    ret = PARMCI_PutValueInt(src, dst, proc);
    return ret;
}


int ARMCI_PutValueLong(long src, void *dst, int proc)
{
    int ret;
    ret = PARMCI_PutValueLong(src, dst, proc);
    return ret;
}


int ARMCI_Put_flag(void *src, void *dst, int bytes, int *f, int v, int proc)
{
    int ret;
    ret = PARMCI_Put_flag(src, dst, bytes, f, v, proc);
    return ret;
}


int ARMCI_Rmw(int op, void *ploc, void *prem, int extra, int proc)
{
    int ret;
    armci_profile_start(ARMCI_PROF_RMW);
    ret = PARMCI_Rmw(op, ploc, prem, extra, proc);
    armci_profile_stop(ARMCI_PROF_RMW);
    return ret;
}


int ARMCI_Test(armci_hdl_t *nb_handle)
{
    int ret;
    ret = PARMCI_Test(nb_handle);
    return ret;
}


void ARMCI_Unlock(int mutex, int proc)
{
    
    PARMCI_Unlock(mutex, proc);
    
}


int ARMCI_Wait(armci_hdl_t *nb_handle)
{
    int ret;
    armci_profile_start(ARMCI_PROF_WAIT);
    ret = PARMCI_Wait(nb_handle);
    armci_profile_stop(ARMCI_PROF_WAIT);
    return ret;
}


int ARMCI_WaitAll()
{
    int ret;
    ret = PARMCI_WaitAll();
    return ret;
}


int ARMCI_WaitProc(int proc)
{
    int ret;
    ret = PARMCI_WaitProc(proc);
    return ret;
}


void armci_msg_barrier()
{
    
    armci_profile_start(ARMCI_PROF_BARRIER);
    parmci_msg_barrier();
    armci_profile_stop(ARMCI_PROF_BARRIER);
    
}


void armci_msg_group_barrier(ARMCI_Group *group)
{
    
    armci_profile_start(ARMCI_PROF_BARRIER);
    parmci_msg_group_barrier(group);
    armci_profile_stop(ARMCI_PROF_BARRIER);
    
}


int armci_notify_wait(int proc, int *pval)
{
    int ret;
    armci_profile_start(ARMCI_PROF_NOTIFY);
    ret = parmci_notify_wait(proc, pval);
    armci_profile_stop(ARMCI_PROF_NOTIFY);
    return ret;
}

