#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "armci.h"
#include "parmci.h"
#include "comex.h"


extern int ARMCI_Default_Proc_Group;
MPI_Comm ARMCI_COMM_WORLD;


static void convert_giov(armci_giov_t *a, comex_giov_t *b, int len)
{
    int i;
    for (i=0; i<len; ++i) {
        b[i].src = a[i].src_ptr_array;
        b[i].dst = a[i].dst_ptr_array;
        b[i].count = a[i].ptr_array_len;
        b[i].bytes = a[i].bytes;
    }
}


int PARMCI_Acc(int optype, void *scale, void *src, void *dst, int bytes, int proc)
{
    return comex_acc(optype, scale, src, dst, bytes, proc, COMEX_GROUP_WORLD);
}


int PARMCI_AccS(int optype, void *scale, void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc)
{
    return comex_accs(optype, scale, src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, proc, COMEX_GROUP_WORLD);
}


int PARMCI_AccV(int op, void *scale, armci_giov_t *darr, int len, int proc)
{
    int rc;
    comex_giov_t *adarr = malloc(sizeof(comex_giov_t) * len);
    convert_giov(darr, adarr, len);
    rc = comex_accv(op, scale, adarr, len, proc, COMEX_GROUP_WORLD);
    free(adarr);
    return rc;
}


/* fence is always on the world group */
void PARMCI_AllFence()
{
    assert(COMEX_SUCCESS == comex_fence_all(COMEX_GROUP_WORLD));
}


void PARMCI_Barrier()
{
    assert(COMEX_SUCCESS == comex_barrier(ARMCI_Default_Proc_Group));
}


int PARMCI_Create_mutexes(int num)
{
    return comex_create_mutexes(num);
}


int PARMCI_Destroy_mutexes()
{
    return comex_destroy_mutexes();
}


/* fence is always on the world group */
void PARMCI_Fence(int proc)
{
    comex_fence_proc(proc, COMEX_GROUP_WORLD);
}


void PARMCI_Finalize()
{
    comex_finalize();
}


int PARMCI_Free(void *ptr)
{
    return comex_free(ptr, ARMCI_Default_Proc_Group);
}


int ARMCI_Free_group(void *ptr, ARMCI_Group *group)
{
    return comex_free(ptr, *group);
}


int PARMCI_Free_local(void *ptr)
{
    return comex_free_local(ptr);
}


int PARMCI_Get(void *src, void *dst, int bytes, int proc)
{
    return comex_get(src, dst, bytes, proc, COMEX_GROUP_WORLD);
}


int PARMCI_GetS(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc)
{
    return comex_gets(src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, proc, COMEX_GROUP_WORLD);
}


int PARMCI_GetV(armci_giov_t *darr, int len, int proc)
{
    int rc;
    comex_giov_t *adarr = malloc(sizeof(comex_giov_t) * len);
    convert_giov(darr, adarr, len);
    rc = comex_getv(adarr, len, proc, COMEX_GROUP_WORLD);
    free(adarr);
    return rc;
}


double PARMCI_GetValueDouble(void *src, int proc)
{
    double val;
    assert(COMEX_SUCCESS == 
            comex_get(src, &val, sizeof(double), proc, COMEX_GROUP_WORLD));
    return val;
}


float PARMCI_GetValueFloat(void *src, int proc)
{
    float val;
    assert(COMEX_SUCCESS == 
            comex_get(src, &val, sizeof(float), proc, COMEX_GROUP_WORLD));
    return val;
}


int PARMCI_GetValueInt(void *src, int proc)
{
    int val;
    assert(COMEX_SUCCESS == 
            comex_get(src, &val, sizeof(int), proc, COMEX_GROUP_WORLD));
    return val;
}


long PARMCI_GetValueLong(void *src, int proc)
{
    long val;
    assert(COMEX_SUCCESS == 
            comex_get(src, &val, sizeof(long), proc, COMEX_GROUP_WORLD));
    return val;
}


int PARMCI_Init()
{
    int rc = comex_init();
    assert(COMEX_SUCCESS == comex_group_comm(COMEX_GROUP_WORLD, &ARMCI_COMM_WORLD));
    return rc;
}


int PARMCI_Init_args(int *argc, char ***argv)
{
    int rc = comex_init_args(argc, argv);
    assert(COMEX_SUCCESS == comex_group_comm(COMEX_GROUP_WORLD, &ARMCI_COMM_WORLD));
    return rc;
}


int PARMCI_Initialized()
{
    return comex_initialized();
}


void PARMCI_Lock(int mutex, int proc)
{
    comex_lock(mutex, proc);
}


int PARMCI_Malloc(void **ptr_arr, armci_size_t bytes)
{
    return comex_malloc(ptr_arr, bytes, ARMCI_Default_Proc_Group);
}


int ARMCI_Malloc_group(void **ptr_arr, armci_size_t bytes, ARMCI_Group *group)
{
    return comex_malloc(ptr_arr, bytes, *group);
}


void* PARMCI_Malloc_local(armci_size_t bytes)
{
    return comex_malloc_local(bytes);
}


void* PARMCI_Memat(armci_meminfo_t *meminfo, long offset)
{
    assert(0);
    return NULL;
}


void PARMCI_Memget(size_t bytes, armci_meminfo_t *meminfo, int memflg)
{
    assert(0);
}


void PARMCI_Memdt(armci_meminfo_t *meminfo, long offset)
{
    assert(0);
}


void PARMCI_Memctl(armci_meminfo_t *meminfo)
{
    assert(0);
}



int PARMCI_NbAccS(int optype, void *scale, void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc, armci_hdl_t *nb_handle)
{
    return comex_nbaccs(optype, scale, src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, proc, COMEX_GROUP_WORLD, nb_handle);
}


int PARMCI_NbAccV(int op, void *scale, armci_giov_t *darr, int len, int proc, armci_hdl_t *nb_handle)
{
    int rc;
    comex_giov_t *adarr = malloc(sizeof(comex_giov_t) * len);
    convert_giov(darr, adarr, len);
    rc = comex_nbaccv(op, scale, adarr, len, proc, COMEX_GROUP_WORLD, nb_handle);
    free(adarr);
    return rc;
}


int PARMCI_NbGet(void *src, void *dst, int bytes, int proc, armci_hdl_t *nb_handle)
{
    return comex_nbget(src, dst, bytes, proc, COMEX_GROUP_WORLD, nb_handle);
}


int PARMCI_NbGetS(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc, armci_hdl_t *nb_handle)
{
    return comex_nbgets(src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, proc, COMEX_GROUP_WORLD, nb_handle);
}


int PARMCI_NbGetV(armci_giov_t *darr, int len, int proc, armci_hdl_t *nb_handle)
{
    int rc;
    comex_giov_t *adarr = malloc(sizeof(comex_giov_t) * len);
    convert_giov(darr, adarr, len);
    rc = comex_nbgetv(adarr, len, proc, COMEX_GROUP_WORLD, nb_handle);
    free(adarr);
    return rc;
}


int PARMCI_NbPut(void *src, void *dst, int bytes, int proc, armci_hdl_t *nb_handle)
{
    return comex_nbput(src, dst, bytes, proc, COMEX_GROUP_WORLD, nb_handle);
}


int PARMCI_NbPutS(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc, armci_hdl_t *nb_handle)
{
    return comex_nbputs(src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, proc, COMEX_GROUP_WORLD, nb_handle);
}


int PARMCI_NbPutV(armci_giov_t *darr, int len, int proc, armci_hdl_t *nb_handle)
{
    int rc;
    comex_giov_t *adarr = malloc(sizeof(comex_giov_t) * len);
    convert_giov(darr, adarr, len);
    rc = comex_nbputv(adarr, len, proc, COMEX_GROUP_WORLD, nb_handle);
    free(adarr);
    return rc;
}


int PARMCI_NbPutValueDouble(double src, void *dst, int proc, armci_hdl_t *nb_handle)
{
    return comex_nbput(&src, dst, sizeof(double), proc, COMEX_GROUP_WORLD, nb_handle);
}


int PARMCI_NbPutValueFloat(float src, void *dst, int proc, armci_hdl_t *nb_handle)
{
    return comex_nbput(&src, dst, sizeof(float), proc, COMEX_GROUP_WORLD, nb_handle);
}


int PARMCI_NbPutValueInt(int src, void *dst, int proc, armci_hdl_t *nb_handle)
{
    return comex_nbput(&src, dst, sizeof(int), proc, COMEX_GROUP_WORLD, nb_handle);
}


int PARMCI_NbPutValueLong(long src, void *dst, int proc, armci_hdl_t *nb_handle)
{
    return comex_nbput(&src, dst, sizeof(long), proc, COMEX_GROUP_WORLD, nb_handle);
}


int PARMCI_Put(void *src, void *dst, int bytes, int proc)
{
    return comex_put(src, dst, bytes, proc, COMEX_GROUP_WORLD);
}


int PARMCI_PutS(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int proc)
{
    return comex_puts(src_ptr, src_stride_arr, dst_ptr, dst_stride_arr, count, stride_levels, proc, COMEX_GROUP_WORLD);
}


int PARMCI_PutS_flag(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int *flag, int val, int proc)
{
    assert(0);
    return 0;
}


int PARMCI_PutS_flag_dir(void *src_ptr, int *src_stride_arr, void *dst_ptr, int *dst_stride_arr, int *count, int stride_levels, int *flag, int val, int proc)
{
    assert(0);
    return 0;
}


int PARMCI_PutV(armci_giov_t *darr, int len, int proc)
{
    int rc;
    comex_giov_t *adarr = malloc(sizeof(comex_giov_t) * len);
    convert_giov(darr, adarr, len);
    rc = comex_putv(adarr, len, proc, COMEX_GROUP_WORLD);
    free(adarr);
    return rc;
}


int PARMCI_PutValueDouble(double src, void *dst, int proc)
{
    return comex_put(&src, dst, sizeof(double), proc, COMEX_GROUP_WORLD);
}


int PARMCI_PutValueFloat(float src, void *dst, int proc)
{
    return comex_put(&src, dst, sizeof(float), proc, COMEX_GROUP_WORLD);
}


int PARMCI_PutValueInt(int src, void *dst, int proc)
{
    return comex_put(&src, dst, sizeof(int), proc, COMEX_GROUP_WORLD);
}


int PARMCI_PutValueLong(long src, void *dst, int proc)
{
    return comex_put(&src, dst, sizeof(long), proc, COMEX_GROUP_WORLD);
}


int PARMCI_Put_flag(void *src, void *dst, int bytes, int *f, int v, int proc)
{
    assert(0);
    return 0;
}


int PARMCI_Rmw(int op, void *ploc, void *prem, int extra, int proc)
{
    return comex_rmw(op, ploc, prem, extra, proc, COMEX_GROUP_WORLD);
}


int PARMCI_Test(armci_hdl_t *nb_handle)
{
    int status;
    assert(COMEX_SUCCESS == comex_test(nb_handle, &status));
    return status;
}


void PARMCI_Unlock(int mutex, int proc)
{
    comex_unlock(mutex, proc);
}


int PARMCI_Wait(armci_hdl_t *nb_handle)
{
    return comex_wait(nb_handle);
}


int PARMCI_WaitAll()
{
    return comex_wait_all(COMEX_GROUP_WORLD);
}


int PARMCI_WaitProc(int proc)
{
    return comex_wait_proc(proc, COMEX_GROUP_WORLD);
}


int parmci_notify(int proc)
{
    assert(0);
    return 0;
}


int parmci_notify_wait(int proc, int *pval)
{
    assert(0);
    return 0;
}


int armci_domain_nprocs(armci_domain_t domain, int id)
{
    return 1;
}


int armci_domain_id(armci_domain_t domain, int glob_proc_id)
{
    return glob_proc_id;
}


int armci_domain_glob_proc_id(armci_domain_t domain, int id, int loc_proc_id)
{
    return id;
}


int armci_domain_my_id(armci_domain_t domain)
{
    int rank;

    assert(comex_initialized());
    comex_group_rank(COMEX_GROUP_WORLD, &rank);

    return rank;
}


int armci_domain_count(armci_domain_t domain)
{
    int size;

    assert(comex_initialized());
    comex_group_size(COMEX_GROUP_WORLD, &size);

    return size;
}


int armci_domain_same_id(armci_domain_t domain, int proc)
{
    int rank;

    comex_group_rank(COMEX_GROUP_WORLD, &rank);

    return (proc == rank);
}


void ARMCI_Error(char *msg, int code)
{
    comex_error(msg, code);
}


void ARMCI_Set_shm_limit(unsigned long shmemlimit)
{
    /* ignore */
}


/* Shared memory not implemented */
int ARMCI_Uses_shm()
{
    return 0;
}


int ARMCI_Uses_shm_group()
{
    return 0;
}


/* Is it memory copy? */
void PARMCI_Copy(void *src, void *dst, int n)
{
    assert(0);
    memcpy(dst, src, sizeof(int) * n);
}


/* Group Functions */
int ARMCI_Uses_shm_grp(ARMCI_Group *group)
{
    assert(0);
}


void ARMCI_Cleanup()
{
    comex_finalize();
}


/* JAD technically not an error to have empty impl of aggregate methods */
void ARMCI_SET_AGGREGATE_HANDLE(armci_hdl_t* handle)
{
}


void ARMCI_UNSET_AGGREGATE_HANDLE(armci_hdl_t* handle)
{
}


/* Always return 0, since shared memory not implemented yet */
int ARMCI_Same_node(int proc)
{
    return 0;
}



