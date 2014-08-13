
#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <stdio.h>

#include <mpi.h>

#include "armci.h"
#include "parmci.h"

static MPI_Comm comm = MPI_COMM_NULL;
static int me = -1;
static int nproc = -1;
static double global_start = 0.0;
static double global_stop = 0.0;

static long count_PARMCI_Acc = 0;
static long count_PARMCI_AccS = 0;
static long count_PARMCI_AccV = 0;
static long count_PARMCI_AllFence = 0;
static long count_PARMCI_Barrier = 0;
static long count_PARMCI_Create_mutexes = 0;
static long count_PARMCI_Destroy_mutexes = 0;
static long count_PARMCI_Fence = 0;
static long count_PARMCI_Finalize = 0;
static long count_PARMCI_Free = 0;
static long count_PARMCI_Free_local = 0;
static long count_PARMCI_Get = 0;
static long count_PARMCI_GetS = 0;
static long count_PARMCI_GetV = 0;
static long count_PARMCI_GetValueDouble = 0;
static long count_PARMCI_GetValueFloat = 0;
static long count_PARMCI_GetValueInt = 0;
static long count_PARMCI_GetValueLong = 0;
static long count_PARMCI_Init = 0;
static long count_PARMCI_Init_args = 0;
static long count_PARMCI_Initialized = 0;
static long count_PARMCI_Lock = 0;
static long count_PARMCI_Malloc = 0;
static long count_PARMCI_Malloc_local = 0;
static long count_PARMCI_Memat = 0;
static long count_PARMCI_Memget = 0;
static long count_PARMCI_NbAccS = 0;
static long count_PARMCI_NbAccV = 0;
static long count_PARMCI_NbGet = 0;
static long count_PARMCI_NbGetS = 0;
static long count_PARMCI_NbGetV = 0;
static long count_PARMCI_NbPut = 0;
static long count_PARMCI_NbPutS = 0;
static long count_PARMCI_NbPutV = 0;
static long count_PARMCI_NbPutValueDouble = 0;
static long count_PARMCI_NbPutValueFloat = 0;
static long count_PARMCI_NbPutValueInt = 0;
static long count_PARMCI_NbPutValueLong = 0;
static long count_PARMCI_Put = 0;
static long count_PARMCI_PutS = 0;
static long count_PARMCI_PutS_flag = 0;
static long count_PARMCI_PutS_flag_dir = 0;
static long count_PARMCI_PutV = 0;
static long count_PARMCI_PutValueDouble = 0;
static long count_PARMCI_PutValueFloat = 0;
static long count_PARMCI_PutValueInt = 0;
static long count_PARMCI_PutValueLong = 0;
static long count_PARMCI_Put_flag = 0;
static long count_PARMCI_Rmw = 0;
static long count_PARMCI_Test = 0;
static long count_PARMCI_Unlock = 0;
static long count_PARMCI_Wait = 0;
static long count_PARMCI_WaitAll = 0;
static long count_PARMCI_WaitProc = 0;
static long count_parmci_msg_barrier = 0;
static long count_parmci_msg_group_barrier = 0;
static long count_parmci_notify = 0;
static long count_parmci_notify_wait = 0;

static double time_PARMCI_Acc = 0;
static double time_PARMCI_AccS = 0;
static double time_PARMCI_AccV = 0;
static double time_PARMCI_AllFence = 0;
static double time_PARMCI_Barrier = 0;
static double time_PARMCI_Create_mutexes = 0;
static double time_PARMCI_Destroy_mutexes = 0;
static double time_PARMCI_Fence = 0;
static double time_PARMCI_Finalize = 0;
static double time_PARMCI_Free = 0;
static double time_PARMCI_Free_local = 0;
static double time_PARMCI_Get = 0;
static double time_PARMCI_GetS = 0;
static double time_PARMCI_GetV = 0;
static double time_PARMCI_GetValueDouble = 0;
static double time_PARMCI_GetValueFloat = 0;
static double time_PARMCI_GetValueInt = 0;
static double time_PARMCI_GetValueLong = 0;
static double time_PARMCI_Init = 0;
static double time_PARMCI_Init_args = 0;
static double time_PARMCI_Initialized = 0;
static double time_PARMCI_Lock = 0;
static double time_PARMCI_Malloc = 0;
static double time_PARMCI_Malloc_local = 0;
static double time_PARMCI_Memat = 0;
static double time_PARMCI_Memget = 0;
static double time_PARMCI_NbAccS = 0;
static double time_PARMCI_NbAccV = 0;
static double time_PARMCI_NbGet = 0;
static double time_PARMCI_NbGetS = 0;
static double time_PARMCI_NbGetV = 0;
static double time_PARMCI_NbPut = 0;
static double time_PARMCI_NbPutS = 0;
static double time_PARMCI_NbPutV = 0;
static double time_PARMCI_NbPutValueDouble = 0;
static double time_PARMCI_NbPutValueFloat = 0;
static double time_PARMCI_NbPutValueInt = 0;
static double time_PARMCI_NbPutValueLong = 0;
static double time_PARMCI_Put = 0;
static double time_PARMCI_PutS = 0;
static double time_PARMCI_PutS_flag = 0;
static double time_PARMCI_PutS_flag_dir = 0;
static double time_PARMCI_PutV = 0;
static double time_PARMCI_PutValueDouble = 0;
static double time_PARMCI_PutValueFloat = 0;
static double time_PARMCI_PutValueInt = 0;
static double time_PARMCI_PutValueLong = 0;
static double time_PARMCI_Put_flag = 0;
static double time_PARMCI_Rmw = 0;
static double time_PARMCI_Test = 0;
static double time_PARMCI_Unlock = 0;
static double time_PARMCI_Wait = 0;
static double time_PARMCI_WaitAll = 0;
static double time_PARMCI_WaitProc = 0;
static double time_parmci_msg_barrier = 0;
static double time_parmci_msg_group_barrier = 0;
static double time_parmci_notify = 0;
static double time_parmci_notify_wait = 0;


int ARMCI_Acc(int optype, void *scale, void *src, void *dst, int bytes, int proc)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_Acc;
    local_start = MPI_Wtime();
    ret = PARMCI_Acc(optype, scale, src, dst, bytes, proc);
    local_stop = MPI_Wtime();
    time_PARMCI_Acc += local_stop - local_start;
    return ret;
}


int ARMCI_AccS(int optype, void *scale, void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_AccS;
    local_start = MPI_Wtime();
    ret = PARMCI_AccS(optype, scale, src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, proc);
    local_stop = MPI_Wtime();
    time_PARMCI_AccS += local_stop - local_start;
    return ret;
}


int ARMCI_AccV(int op, void *scale, armci_giov_t *darr, int len, int proc)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_AccV;
    local_start = MPI_Wtime();
    ret = PARMCI_AccV(op, scale, darr, len, proc);
    local_stop = MPI_Wtime();
    time_PARMCI_AccV += local_stop - local_start;
    return ret;
}


void ARMCI_AllFence()
{
    double local_start;
    double local_stop;
    
    ++count_PARMCI_AllFence;
    local_start = MPI_Wtime();
    PARMCI_AllFence();
    local_stop = MPI_Wtime();
    time_PARMCI_AllFence += local_stop - local_start;
    
}


void ARMCI_Barrier()
{
    double local_start;
    double local_stop;
    
    ++count_PARMCI_Barrier;
    local_start = MPI_Wtime();
    PARMCI_Barrier();
    local_stop = MPI_Wtime();
    time_PARMCI_Barrier += local_stop - local_start;
    
}


int ARMCI_Create_mutexes(int num)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_Create_mutexes;
    local_start = MPI_Wtime();
    ret = PARMCI_Create_mutexes(num);
    local_stop = MPI_Wtime();
    time_PARMCI_Create_mutexes += local_stop - local_start;
    return ret;
}


int ARMCI_Destroy_mutexes()
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_Destroy_mutexes;
    local_start = MPI_Wtime();
    ret = PARMCI_Destroy_mutexes();
    local_stop = MPI_Wtime();
    time_PARMCI_Destroy_mutexes += local_stop - local_start;
    return ret;
}


void ARMCI_Fence(int proc)
{
    double local_start;
    double local_stop;
    
    ++count_PARMCI_Fence;
    local_start = MPI_Wtime();
    PARMCI_Fence(proc);
    local_stop = MPI_Wtime();
    time_PARMCI_Fence += local_stop - local_start;
    
}


int ARMCI_Free(void *ptr)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_Free;
    local_start = MPI_Wtime();
    ret = PARMCI_Free(ptr);
    local_stop = MPI_Wtime();
    time_PARMCI_Free += local_stop - local_start;
    return ret;
}


int ARMCI_Free_local(void *ptr)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_Free_local;
    local_start = MPI_Wtime();
    ret = PARMCI_Free_local(ptr);
    local_stop = MPI_Wtime();
    time_PARMCI_Free_local += local_stop - local_start;
    return ret;
}


int ARMCI_Get(void *src, void *dst, int bytes, int proc)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_Get;
    local_start = MPI_Wtime();
    ret = PARMCI_Get(src, dst, bytes, proc);
    local_stop = MPI_Wtime();
    time_PARMCI_Get += local_stop - local_start;
    return ret;
}


int ARMCI_GetS(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_GetS;
    local_start = MPI_Wtime();
    ret = PARMCI_GetS(src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, proc);
    local_stop = MPI_Wtime();
    time_PARMCI_GetS += local_stop - local_start;
    return ret;
}


int ARMCI_GetV(armci_giov_t *darr, int len, int proc)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_GetV;
    local_start = MPI_Wtime();
    ret = PARMCI_GetV(darr, len, proc);
    local_stop = MPI_Wtime();
    time_PARMCI_GetV += local_stop - local_start;
    return ret;
}


double ARMCI_GetValueDouble(void *src, int proc)
{
    double local_start;
    double local_stop;
    double ret;
    ++count_PARMCI_GetValueDouble;
    local_start = MPI_Wtime();
    ret = PARMCI_GetValueDouble(src, proc);
    local_stop = MPI_Wtime();
    time_PARMCI_GetValueDouble += local_stop - local_start;
    return ret;
}


float ARMCI_GetValueFloat(void *src, int proc)
{
    double local_start;
    double local_stop;
    float ret;
    ++count_PARMCI_GetValueFloat;
    local_start = MPI_Wtime();
    ret = PARMCI_GetValueFloat(src, proc);
    local_stop = MPI_Wtime();
    time_PARMCI_GetValueFloat += local_stop - local_start;
    return ret;
}


int ARMCI_GetValueInt(void *src, int proc)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_GetValueInt;
    local_start = MPI_Wtime();
    ret = PARMCI_GetValueInt(src, proc);
    local_stop = MPI_Wtime();
    time_PARMCI_GetValueInt += local_stop - local_start;
    return ret;
}


long ARMCI_GetValueLong(void *src, int proc)
{
    double local_start;
    double local_stop;
    long ret;
    ++count_PARMCI_GetValueLong;
    local_start = MPI_Wtime();
    ret = PARMCI_GetValueLong(src, proc);
    local_stop = MPI_Wtime();
    time_PARMCI_GetValueLong += local_stop - local_start;
    return ret;
}


int ARMCI_Init()
{
    double local_start;
    double local_stop;
    int ret;
    if (comm == MPI_COMM_NULL) {
        MPI_Comm_dup(MPI_COMM_WORLD, &comm);
        MPI_Comm_rank(comm, &me);
        MPI_Comm_size(comm, &nproc);
    }
    ++count_PARMCI_Init;
    global_start = MPI_Wtime();
    local_start = MPI_Wtime();
    ret = PARMCI_Init();
    local_stop = MPI_Wtime();
    time_PARMCI_Init += local_stop - local_start;
    return ret;
}


int ARMCI_Init_args(int *argc, char ***argv)
{
    double local_start;
    double local_stop;
    int ret;
    if (comm == MPI_COMM_NULL) {
        MPI_Comm_dup(MPI_COMM_WORLD, &comm);
        MPI_Comm_rank(comm, &me);
        MPI_Comm_size(comm, &nproc);
    }
    ++count_PARMCI_Init_args;
    global_start = MPI_Wtime();
    local_start = MPI_Wtime();
    ret = PARMCI_Init_args(argc, argv);
    local_stop = MPI_Wtime();
    time_PARMCI_Init_args += local_stop - local_start;
    return ret;
}


int ARMCI_Initialized()
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_Initialized;
    local_start = MPI_Wtime();
    ret = PARMCI_Initialized();
    local_stop = MPI_Wtime();
    time_PARMCI_Initialized += local_stop - local_start;
    return ret;
}


void ARMCI_Lock(int mutex, int proc)
{
    double local_start;
    double local_stop;
    
    ++count_PARMCI_Lock;
    local_start = MPI_Wtime();
    PARMCI_Lock(mutex, proc);
    local_stop = MPI_Wtime();
    time_PARMCI_Lock += local_stop - local_start;
    
}


int ARMCI_Malloc(void **ptr_arr, armci_size_t bytes)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_Malloc;
    local_start = MPI_Wtime();
    ret = PARMCI_Malloc(ptr_arr, bytes);
    local_stop = MPI_Wtime();
    time_PARMCI_Malloc += local_stop - local_start;
    return ret;
}


void* ARMCI_Malloc_local(armci_size_t bytes)
{
    double local_start;
    double local_stop;
    void* ret;
    ++count_PARMCI_Malloc_local;
    local_start = MPI_Wtime();
    ret = PARMCI_Malloc_local(bytes);
    local_stop = MPI_Wtime();
    time_PARMCI_Malloc_local += local_stop - local_start;
    return ret;
}


void* ARMCI_Memat(armci_meminfo_t *meminfo, long offset)
{
    double local_start;
    double local_stop;
    void* ret;
    ++count_PARMCI_Memat;
    local_start = MPI_Wtime();
    ret = PARMCI_Memat(meminfo, offset);
    local_stop = MPI_Wtime();
    time_PARMCI_Memat += local_stop - local_start;
    return ret;
}


void ARMCI_Memget(size_t bytes, armci_meminfo_t *meminfo, int memflg)
{
    double local_start;
    double local_stop;
    
    ++count_PARMCI_Memget;
    local_start = MPI_Wtime();
    PARMCI_Memget(bytes, meminfo, memflg);
    local_stop = MPI_Wtime();
    time_PARMCI_Memget += local_stop - local_start;
    
}


int ARMCI_NbAccS(int optype, void *scale, void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc, armci_hdl_t *nb_handle)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_NbAccS;
    local_start = MPI_Wtime();
    ret = PARMCI_NbAccS(optype, scale, src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, proc, nb_handle);
    local_stop = MPI_Wtime();
    time_PARMCI_NbAccS += local_stop - local_start;
    return ret;
}


int ARMCI_NbAccV(int op, void *scale, armci_giov_t *darr, int len, int proc, armci_hdl_t *nb_handle)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_NbAccV;
    local_start = MPI_Wtime();
    ret = PARMCI_NbAccV(op, scale, darr, len, proc, nb_handle);
    local_stop = MPI_Wtime();
    time_PARMCI_NbAccV += local_stop - local_start;
    return ret;
}


int ARMCI_NbGet(void *src, void *dst, int bytes, int proc, armci_hdl_t *nb_handle)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_NbGet;
    local_start = MPI_Wtime();
    ret = PARMCI_NbGet(src, dst, bytes, proc, nb_handle);
    local_stop = MPI_Wtime();
    time_PARMCI_NbGet += local_stop - local_start;
    return ret;
}


int ARMCI_NbGetS(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc, armci_hdl_t *nb_handle)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_NbGetS;
    local_start = MPI_Wtime();
    ret = PARMCI_NbGetS(src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, proc, nb_handle);
    local_stop = MPI_Wtime();
    time_PARMCI_NbGetS += local_stop - local_start;
    return ret;
}


int ARMCI_NbGetV(armci_giov_t *darr, int len, int proc, armci_hdl_t *nb_handle)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_NbGetV;
    local_start = MPI_Wtime();
    ret = PARMCI_NbGetV(darr, len, proc, nb_handle);
    local_stop = MPI_Wtime();
    time_PARMCI_NbGetV += local_stop - local_start;
    return ret;
}


int ARMCI_NbPut(void *src, void *dst, int bytes, int proc, armci_hdl_t *nb_handle)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_NbPut;
    local_start = MPI_Wtime();
    ret = PARMCI_NbPut(src, dst, bytes, proc, nb_handle);
    local_stop = MPI_Wtime();
    time_PARMCI_NbPut += local_stop - local_start;
    return ret;
}


int ARMCI_NbPutS(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc, armci_hdl_t *nb_handle)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_NbPutS;
    local_start = MPI_Wtime();
    ret = PARMCI_NbPutS(src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, proc, nb_handle);
    local_stop = MPI_Wtime();
    time_PARMCI_NbPutS += local_stop - local_start;
    return ret;
}


int ARMCI_NbPutV(armci_giov_t *darr, int len, int proc, armci_hdl_t *nb_handle)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_NbPutV;
    local_start = MPI_Wtime();
    ret = PARMCI_NbPutV(darr, len, proc, nb_handle);
    local_stop = MPI_Wtime();
    time_PARMCI_NbPutV += local_stop - local_start;
    return ret;
}


int ARMCI_NbPutValueDouble(double src, void *dst, int proc, armci_hdl_t *nb_handle)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_NbPutValueDouble;
    local_start = MPI_Wtime();
    ret = PARMCI_NbPutValueDouble(src, dst, proc, nb_handle);
    local_stop = MPI_Wtime();
    time_PARMCI_NbPutValueDouble += local_stop - local_start;
    return ret;
}


int ARMCI_NbPutValueFloat(float src, void *dst, int proc, armci_hdl_t *nb_handle)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_NbPutValueFloat;
    local_start = MPI_Wtime();
    ret = PARMCI_NbPutValueFloat(src, dst, proc, nb_handle);
    local_stop = MPI_Wtime();
    time_PARMCI_NbPutValueFloat += local_stop - local_start;
    return ret;
}


int ARMCI_NbPutValueInt(int src, void *dst, int proc, armci_hdl_t *nb_handle)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_NbPutValueInt;
    local_start = MPI_Wtime();
    ret = PARMCI_NbPutValueInt(src, dst, proc, nb_handle);
    local_stop = MPI_Wtime();
    time_PARMCI_NbPutValueInt += local_stop - local_start;
    return ret;
}


int ARMCI_NbPutValueLong(long src, void *dst, int proc, armci_hdl_t *nb_handle)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_NbPutValueLong;
    local_start = MPI_Wtime();
    ret = PARMCI_NbPutValueLong(src, dst, proc, nb_handle);
    local_stop = MPI_Wtime();
    time_PARMCI_NbPutValueLong += local_stop - local_start;
    return ret;
}


int ARMCI_Put(void *src, void *dst, int bytes, int proc)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_Put;
    local_start = MPI_Wtime();
    ret = PARMCI_Put(src, dst, bytes, proc);
    local_stop = MPI_Wtime();
    time_PARMCI_Put += local_stop - local_start;
    return ret;
}


int ARMCI_PutS(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_PutS;
    local_start = MPI_Wtime();
    ret = PARMCI_PutS(src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, proc);
    local_stop = MPI_Wtime();
    time_PARMCI_PutS += local_stop - local_start;
    return ret;
}


int ARMCI_PutS_flag(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int *flag, int val, int proc)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_PutS_flag;
    local_start = MPI_Wtime();
    ret = PARMCI_PutS_flag(src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, flag, val, proc);
    local_stop = MPI_Wtime();
    time_PARMCI_PutS_flag += local_stop - local_start;
    return ret;
}


int ARMCI_PutS_flag_dir(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int *flag, int val, int proc)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_PutS_flag_dir;
    local_start = MPI_Wtime();
    ret = PARMCI_PutS_flag_dir(src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, flag, val, proc);
    local_stop = MPI_Wtime();
    time_PARMCI_PutS_flag_dir += local_stop - local_start;
    return ret;
}


int ARMCI_PutV(armci_giov_t *darr, int len, int proc)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_PutV;
    local_start = MPI_Wtime();
    ret = PARMCI_PutV(darr, len, proc);
    local_stop = MPI_Wtime();
    time_PARMCI_PutV += local_stop - local_start;
    return ret;
}


int ARMCI_PutValueDouble(double src, void *dst, int proc)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_PutValueDouble;
    local_start = MPI_Wtime();
    ret = PARMCI_PutValueDouble(src, dst, proc);
    local_stop = MPI_Wtime();
    time_PARMCI_PutValueDouble += local_stop - local_start;
    return ret;
}


int ARMCI_PutValueFloat(float src, void *dst, int proc)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_PutValueFloat;
    local_start = MPI_Wtime();
    ret = PARMCI_PutValueFloat(src, dst, proc);
    local_stop = MPI_Wtime();
    time_PARMCI_PutValueFloat += local_stop - local_start;
    return ret;
}


int ARMCI_PutValueInt(int src, void *dst, int proc)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_PutValueInt;
    local_start = MPI_Wtime();
    ret = PARMCI_PutValueInt(src, dst, proc);
    local_stop = MPI_Wtime();
    time_PARMCI_PutValueInt += local_stop - local_start;
    return ret;
}


int ARMCI_PutValueLong(long src, void *dst, int proc)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_PutValueLong;
    local_start = MPI_Wtime();
    ret = PARMCI_PutValueLong(src, dst, proc);
    local_stop = MPI_Wtime();
    time_PARMCI_PutValueLong += local_stop - local_start;
    return ret;
}


int ARMCI_Put_flag(void *src, void *dst, int bytes, int *f, int v, int proc)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_Put_flag;
    local_start = MPI_Wtime();
    ret = PARMCI_Put_flag(src, dst, bytes, f, v, proc);
    local_stop = MPI_Wtime();
    time_PARMCI_Put_flag += local_stop - local_start;
    return ret;
}


int ARMCI_Rmw(int op, void *ploc, void *prem, int extra, int proc)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_Rmw;
    local_start = MPI_Wtime();
    ret = PARMCI_Rmw(op, ploc, prem, extra, proc);
    local_stop = MPI_Wtime();
    time_PARMCI_Rmw += local_stop - local_start;
    return ret;
}


int ARMCI_Test(armci_hdl_t *nb_handle)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_Test;
    local_start = MPI_Wtime();
    ret = PARMCI_Test(nb_handle);
    local_stop = MPI_Wtime();
    time_PARMCI_Test += local_stop - local_start;
    return ret;
}


void ARMCI_Unlock(int mutex, int proc)
{
    double local_start;
    double local_stop;
    
    ++count_PARMCI_Unlock;
    local_start = MPI_Wtime();
    PARMCI_Unlock(mutex, proc);
    local_stop = MPI_Wtime();
    time_PARMCI_Unlock += local_stop - local_start;
    
}


int ARMCI_Wait(armci_hdl_t *nb_handle)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_Wait;
    local_start = MPI_Wtime();
    ret = PARMCI_Wait(nb_handle);
    local_stop = MPI_Wtime();
    time_PARMCI_Wait += local_stop - local_start;
    return ret;
}


int ARMCI_WaitAll()
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_WaitAll;
    local_start = MPI_Wtime();
    ret = PARMCI_WaitAll();
    local_stop = MPI_Wtime();
    time_PARMCI_WaitAll += local_stop - local_start;
    return ret;
}


int ARMCI_WaitProc(int proc)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_PARMCI_WaitProc;
    local_start = MPI_Wtime();
    ret = PARMCI_WaitProc(proc);
    local_stop = MPI_Wtime();
    time_PARMCI_WaitProc += local_stop - local_start;
    return ret;
}


void armci_msg_barrier()
{
    double local_start;
    double local_stop;
    
    ++count_parmci_msg_barrier;
    local_start = MPI_Wtime();
    parmci_msg_barrier();
    local_stop = MPI_Wtime();
    time_parmci_msg_barrier += local_stop - local_start;
    
}


void armci_msg_group_barrier(ARMCI_Group *group)
{
    double local_start;
    double local_stop;
    
    ++count_parmci_msg_group_barrier;
    local_start = MPI_Wtime();
    parmci_msg_group_barrier(group);
    local_stop = MPI_Wtime();
    time_parmci_msg_group_barrier += local_stop - local_start;
    
}


int armci_notify(int proc)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_parmci_notify;
    local_start = MPI_Wtime();
    ret = parmci_notify(proc);
    local_stop = MPI_Wtime();
    time_parmci_notify += local_stop - local_start;
    return ret;
}


int armci_notify_wait(int proc, int *pval)
{
    double local_start;
    double local_stop;
    int ret;
    ++count_parmci_notify_wait;
    local_start = MPI_Wtime();
    ret = parmci_notify_wait(proc, pval);
    local_stop = MPI_Wtime();
    time_parmci_notify_wait += local_stop - local_start;
    return ret;
}

void ARMCI_Finalize()
{
    int ret;
    ++count_PARMCI_Finalize;
    PARMCI_Finalize();
    global_stop = MPI_Wtime();
    /* don't dump info if terminate more than once */
    if (1 == count_PARMCI_Finalize) {

        double recvbuf = 0.0;

        MPI_Reduce(&time_PARMCI_Acc, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_Acc,%ld,%lf\n", count_PARMCI_Acc, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_AccS, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_AccS,%ld,%lf\n", count_PARMCI_AccS, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_AccV, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_AccV,%ld,%lf\n", count_PARMCI_AccV, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_AllFence, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_AllFence,%ld,%lf\n", count_PARMCI_AllFence, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_Barrier, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_Barrier,%ld,%lf\n", count_PARMCI_Barrier, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_Create_mutexes, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_Create_mutexes,%ld,%lf\n", count_PARMCI_Create_mutexes, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_Destroy_mutexes, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_Destroy_mutexes,%ld,%lf\n", count_PARMCI_Destroy_mutexes, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_Fence, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_Fence,%ld,%lf\n", count_PARMCI_Fence, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_Finalize, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_Finalize,%ld,%lf\n", count_PARMCI_Finalize, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_Free, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_Free,%ld,%lf\n", count_PARMCI_Free, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_Free_local, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_Free_local,%ld,%lf\n", count_PARMCI_Free_local, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_Get, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_Get,%ld,%lf\n", count_PARMCI_Get, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_GetS, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_GetS,%ld,%lf\n", count_PARMCI_GetS, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_GetV, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_GetV,%ld,%lf\n", count_PARMCI_GetV, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_GetValueDouble, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_GetValueDouble,%ld,%lf\n", count_PARMCI_GetValueDouble, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_GetValueFloat, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_GetValueFloat,%ld,%lf\n", count_PARMCI_GetValueFloat, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_GetValueInt, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_GetValueInt,%ld,%lf\n", count_PARMCI_GetValueInt, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_GetValueLong, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_GetValueLong,%ld,%lf\n", count_PARMCI_GetValueLong, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_Init, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_Init,%ld,%lf\n", count_PARMCI_Init, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_Init_args, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_Init_args,%ld,%lf\n", count_PARMCI_Init_args, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_Initialized, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_Initialized,%ld,%lf\n", count_PARMCI_Initialized, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_Lock, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_Lock,%ld,%lf\n", count_PARMCI_Lock, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_Malloc, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_Malloc,%ld,%lf\n", count_PARMCI_Malloc, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_Malloc_local, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_Malloc_local,%ld,%lf\n", count_PARMCI_Malloc_local, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_Memat, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_Memat,%ld,%lf\n", count_PARMCI_Memat, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_Memget, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_Memget,%ld,%lf\n", count_PARMCI_Memget, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_NbAccS, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_NbAccS,%ld,%lf\n", count_PARMCI_NbAccS, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_NbAccV, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_NbAccV,%ld,%lf\n", count_PARMCI_NbAccV, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_NbGet, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_NbGet,%ld,%lf\n", count_PARMCI_NbGet, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_NbGetS, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_NbGetS,%ld,%lf\n", count_PARMCI_NbGetS, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_NbGetV, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_NbGetV,%ld,%lf\n", count_PARMCI_NbGetV, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_NbPut, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_NbPut,%ld,%lf\n", count_PARMCI_NbPut, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_NbPutS, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_NbPutS,%ld,%lf\n", count_PARMCI_NbPutS, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_NbPutV, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_NbPutV,%ld,%lf\n", count_PARMCI_NbPutV, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_NbPutValueDouble, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_NbPutValueDouble,%ld,%lf\n", count_PARMCI_NbPutValueDouble, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_NbPutValueFloat, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_NbPutValueFloat,%ld,%lf\n", count_PARMCI_NbPutValueFloat, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_NbPutValueInt, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_NbPutValueInt,%ld,%lf\n", count_PARMCI_NbPutValueInt, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_NbPutValueLong, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_NbPutValueLong,%ld,%lf\n", count_PARMCI_NbPutValueLong, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_Put, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_Put,%ld,%lf\n", count_PARMCI_Put, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_PutS, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_PutS,%ld,%lf\n", count_PARMCI_PutS, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_PutS_flag, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_PutS_flag,%ld,%lf\n", count_PARMCI_PutS_flag, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_PutS_flag_dir, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_PutS_flag_dir,%ld,%lf\n", count_PARMCI_PutS_flag_dir, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_PutV, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_PutV,%ld,%lf\n", count_PARMCI_PutV, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_PutValueDouble, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_PutValueDouble,%ld,%lf\n", count_PARMCI_PutValueDouble, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_PutValueFloat, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_PutValueFloat,%ld,%lf\n", count_PARMCI_PutValueFloat, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_PutValueInt, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_PutValueInt,%ld,%lf\n", count_PARMCI_PutValueInt, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_PutValueLong, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_PutValueLong,%ld,%lf\n", count_PARMCI_PutValueLong, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_Put_flag, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_Put_flag,%ld,%lf\n", count_PARMCI_Put_flag, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_Rmw, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_Rmw,%ld,%lf\n", count_PARMCI_Rmw, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_Test, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_Test,%ld,%lf\n", count_PARMCI_Test, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_Unlock, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_Unlock,%ld,%lf\n", count_PARMCI_Unlock, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_Wait, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_Wait,%ld,%lf\n", count_PARMCI_Wait, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_WaitAll, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_WaitAll,%ld,%lf\n", count_PARMCI_WaitAll, recvbuf);
        }

        MPI_Reduce(&time_PARMCI_WaitProc, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("PARMCI_WaitProc,%ld,%lf\n", count_PARMCI_WaitProc, recvbuf);
        }

        MPI_Reduce(&time_parmci_msg_barrier, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("parmci_msg_barrier,%ld,%lf\n", count_parmci_msg_barrier, recvbuf);
        }

        MPI_Reduce(&time_parmci_msg_group_barrier, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("parmci_msg_group_barrier,%ld,%lf\n", count_parmci_msg_group_barrier, recvbuf);
        }

        MPI_Reduce(&time_parmci_notify, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("parmci_notify,%ld,%lf\n", count_parmci_notify, recvbuf);
        }

        MPI_Reduce(&time_parmci_notify_wait, &recvbuf, 1, MPI_DOUBLE, MPI_SUM, 0, comm);
        if (me == 0) {
            printf("parmci_notify_wait,%ld,%lf\n", count_parmci_notify_wait, recvbuf);
        }

        if (me == 0) {
            printf("global_stop-global_start,0,%lf\n",
                    global_stop-global_start);
        }
    }
    MPI_Comm_free(&comm);
}

