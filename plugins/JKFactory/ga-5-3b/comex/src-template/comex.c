#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* C and/or system headers */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <strings.h>
#include <unistd.h>

/* 3rd party headers */
#include <mpi.h>

/* our headers */
#include "comex.h"
#include "comex_impl.h"
#include "groups.h"

#define DEBUG 0


/* exported state */
local_state l_state;

/* static state */
static int  initialized=0;  /* for comex_initialized(), 0=false */
static char skip_lock=0;    /* don't acquire or release lock */

/* static function declarations */
static void acquire_remote_lock(int proc);
static void release_remote_lock(int proc);
static inline void acc(
        int datatype, int count, void *get_buf,
        void *src_ptr, long src_idx, void *scale);

/* needed for complex accumulate */
typedef struct {
    double real;
    double imag;
} DoubleComplex;

typedef struct {
    float real;
    float imag;
} SingleComplex;


int comex_init()
{
    int status;
    
    if (initialized) {
        return 0;
    }
    initialized = 1;

    /* Assert MPI has been initialized */
    int init_flag;
    status = MPI_Initialized(&init_flag);
    assert(MPI_SUCCESS == status);
    assert(init_flag);
    
    /* Duplicate the World Communicator */
    status = MPI_Comm_dup(MPI_COMM_WORLD, &(l_state.world_comm));
    assert(MPI_SUCCESS == status);
    assert(l_state.world_comm); 

    /* My Rank */
    status = MPI_Comm_rank(l_state.world_comm, &(l_state.rank));
    assert(MPI_SUCCESS == status);

    /* World Size */
    status = MPI_Comm_size(l_state.world_comm, &(l_state.size));
    assert(MPI_SUCCESS == status);
    
    /* groups */
    comex_group_init();

    /* Synch - Sanity Check */
    MPI_Barrier(l_state.world_comm);

    return COMEX_SUCCESS;
}


int comex_init_args(int *argc, char ***argv)
{
    int init_flag;
    
    MPI_Initialized(&init_flag);
    
    if(!init_flag) {
        MPI_Init(argc, argv);
    }
    
    return comex_init();
}


int comex_initialized()
{
    return initialized;
}


void comex_error(char *msg, int code)
{
    fprintf(stderr,"[%d] Received an Error in Communication: (%d) %s\n",
            l_state.rank, code, msg);
    
    MPI_Abort(l_state.world_comm, code);
}


int comex_put(
        void *src, void *dst, int bytes,
        int proc, comex_group_t group)
{
    assert(0);
    return COMEX_SUCCESS;
}


int comex_get(
        void *src, void *dst, int bytes,
        int proc, comex_group_t group)
{
    assert(0);
    return COMEX_SUCCESS;
}


int comex_acc(
        int datatype, void *scale,
        void *src_ptr, void *dst_ptr, int bytes,
        int proc, comex_group_t group)
{
    return comex_accs(
            datatype, scale,
            src_ptr, NULL,
            dst_ptr, NULL,
            &bytes, 0,
            proc, group);
}


int comex_puts(
        void *src_ptr, int *src_stride_ar,
        void *dst_ptr, int *dst_stride_ar,
        int *count, int stride_levels,
        int proc, comex_group_t group)
{
    int i, j;
    long src_idx, dst_idx;  /* index offset of current block position to ptr */
    int n1dim;  /* number of 1 dim block */
    int src_bvalue[7], src_bunit[7];
    int dst_bvalue[7], dst_bunit[7];
    int status;

    /* number of n-element of the first dimension */
    n1dim = 1;
    for(i=1; i<=stride_levels; i++) {
        n1dim *= count[i];
    }

    /* calculate the destination indices */
    src_bvalue[0] = 0; src_bvalue[1] = 0; src_bunit[0] = 1; src_bunit[1] = 1;
    dst_bvalue[0] = 0; dst_bvalue[1] = 0; dst_bunit[0] = 1; dst_bunit[1] = 1;

    for(i=2; i<=stride_levels; i++) {
        src_bvalue[i] = 0;
        dst_bvalue[i] = 0;
        src_bunit[i] = src_bunit[i-1] * count[i-1];
        dst_bunit[i] = dst_bunit[i-1] * count[i-1];
    }

    /* index mangling */
    for(i=0; i<n1dim; i++) {
        src_idx = 0;
        dst_idx = 0;
        for(j=1; j<=stride_levels; j++) {
            src_idx += src_bvalue[j] * src_stride_ar[j-1];
            if((i+1) % src_bunit[j] == 0) {
                src_bvalue[j]++;
            }
            if(src_bvalue[j] > (count[j]-1)) {
                src_bvalue[j] = 0;
            }
        }

        for(j=1; j<=stride_levels; j++) {
            dst_idx += dst_bvalue[j] * dst_stride_ar[j-1];
            if((i+1) % dst_bunit[j] == 0) {
                dst_bvalue[j]++;
            }
            if(dst_bvalue[j] > (count[j]-1)) {
                dst_bvalue[j] = 0;
            }
        }
        
        status = comex_put((char *)src_ptr + src_idx, 
                (char *)dst_ptr + dst_idx, count[0], proc, group);
        assert(status == COMEX_SUCCESS);
    }

    return COMEX_SUCCESS;
}


int comex_gets(
        void *src_ptr, int *src_stride_ar,
        void *dst_ptr, int *dst_stride_ar,
        int *count, int stride_levels,
        int proc, comex_group_t group)
{
    int i, j;
    long src_idx, dst_idx;  /* index offset of current block position to ptr */
    int n1dim;  /* number of 1 dim block */
    int src_bvalue[7], src_bunit[7];
    int dst_bvalue[7], dst_bunit[7];
    int status;

    /* number of n-element of the first dimension */
    n1dim = 1;
    for(i=1; i<=stride_levels; i++) {
        n1dim *= count[i];
    }

    /* calculate the destination indices */
    src_bvalue[0] = 0; src_bvalue[1] = 0; src_bunit[0] = 1; src_bunit[1] = 1;
    dst_bvalue[0] = 0; dst_bvalue[1] = 0; dst_bunit[0] = 1; dst_bunit[1] = 1;

    for(i=2; i<=stride_levels; i++) {
        src_bvalue[i] = 0;
        dst_bvalue[i] = 0;
        src_bunit[i] = src_bunit[i-1] * count[i-1];
        dst_bunit[i] = dst_bunit[i-1] * count[i-1];
    }

    for(i=0; i<n1dim; i++) {
        src_idx = 0;
        for(j=1; j<=stride_levels; j++) {
            src_idx += src_bvalue[j] * src_stride_ar[j-1];
            if((i+1) % src_bunit[j] == 0) {
                src_bvalue[j]++;
            }
            if(src_bvalue[j] > (count[j]-1)) {
                src_bvalue[j] = 0;
            }
        }

        dst_idx = 0;
        
        for(j=1; j<=stride_levels; j++) {
            dst_idx += dst_bvalue[j] * dst_stride_ar[j-1];
            if((i+1) % dst_bunit[j] == 0) {
                dst_bvalue[j]++;
            }
            if(dst_bvalue[j] > (count[j]-1)) {
                dst_bvalue[j] = 0;
            }
        }
        
        status = comex_get((char *)src_ptr + src_idx, 
                (char *)dst_ptr + dst_idx, count[0], proc, group);
        assert(status == COMEX_SUCCESS);
    }
    
    return COMEX_SUCCESS;
}


int comex_accs(
        int datatype, void *scale,
        void *src_ptr, int *src_stride_ar,
        void *dst_ptr, int *dst_stride_ar,
        int *count, int stride_levels,
        int proc, comex_group_t group)
{
    int i, j;
    long src_idx, dst_idx;  /* index offset of current block position to ptr */
    int n1dim;  /* number of 1 dim block */
    int src_bvalue[7], src_bunit[7];
    int dst_bvalue[7], dst_bunit[7];
    void *get_buf;

    /* number of n-element of the first dimension */
    n1dim = 1;
    for(i=1; i<=stride_levels; i++)
        n1dim *= count[i];

    /* calculate the destination indices */
    src_bvalue[0] = 0; src_bvalue[1] = 0; src_bunit[0] = 1; src_bunit[1] = 1;
    dst_bvalue[0] = 0; dst_bvalue[1] = 0; dst_bunit[0] = 1; dst_bunit[1] = 1;

    for(i=2; i<=stride_levels; i++)
    {
        src_bvalue[i] = 0;
        dst_bvalue[i] = 0;
        src_bunit[i] = src_bunit[i-1] * count[i-1];
        dst_bunit[i] = dst_bunit[i-1] * count[i-1];
    }

    if (0 == skip_lock) {
        // grab the atomics lock
        acquire_remote_lock(proc);
    }

    get_buf = (char *)malloc(sizeof(char) * count[0]);
    assert(get_buf);

    for(i=0; i<n1dim; i++) {
        src_idx = 0;
        for(j=1; j<=stride_levels; j++) {
            src_idx += src_bvalue[j] * src_stride_ar[j-1];
            if((i+1) % src_bunit[j] == 0) {
                src_bvalue[j]++;
            }
            if(src_bvalue[j] > (count[j]-1)) {
                src_bvalue[j] = 0;
            }
        }

        dst_idx = 0;

        for(j=1; j<=stride_levels; j++) {
            dst_idx += dst_bvalue[j] * dst_stride_ar[j-1];
            if((i+1) % dst_bunit[j] == 0) {
                dst_bvalue[j]++;
            }
            if(dst_bvalue[j] > (count[j]-1)) {
                dst_bvalue[j] = 0;
            }
        }

        // Get the remote data in a temp buffer
        comex_get((char *)dst_ptr + dst_idx, get_buf, count[0], proc, group);

        // Local accumulate
        acc(datatype, count[0], get_buf, src_ptr, src_idx, scale);

        // Write back to remote data
        comex_put(get_buf, (char *)dst_ptr + dst_idx, count[0], proc, group);
    }

    if (0 == skip_lock) {
        // ungrab the lock
        release_remote_lock(proc);
    }

    free(get_buf);

    return COMEX_SUCCESS;
}


int comex_putv(
        comex_giov_t *iov, int iov_len,
        int proc, comex_group_t group)
{
    int status;
    int i;

    for (i=0; i<iov_len; ++i) {
        int j;
        void **src = iov[i].src;
        void **dst = iov[i].dst;
        int bytes = iov[i].bytes;
        int limit = iov[i].count;
        for (j=0; j<limit; ++j) {
            status = comex_put(src[j], dst[j], bytes, proc, group);
            assert(status == COMEX_SUCCESS);
        }
    }

    return COMEX_SUCCESS;
}


int comex_getv(
        comex_giov_t *iov, int iov_len,
        int proc, comex_group_t group)
{
    int status;
    int i;

    for (i=0; i<iov_len; ++i) {
        int j;
        void **src = iov[i].src;
        void **dst = iov[i].dst;
        int bytes = iov[i].bytes;
        int limit = iov[i].count;
        for (j=0; j<limit; ++j) {
            status = comex_get(src[j], dst[j], bytes, proc, group);
            assert(status == COMEX_SUCCESS);
        }
    }

    return COMEX_SUCCESS;
}


int comex_accv(
        int datatype, void *scale,
        comex_giov_t *iov, int iov_len,
        int proc, comex_group_t group)
{
    int i;
    
    skip_lock = 1;
    acquire_remote_lock(proc);

    for (i=0; i<iov_len; ++i) {
        int j;
        void **src = iov[i].src;
        void **dst = iov[i].dst;
        int bytes = iov[i].bytes;
        int limit = iov[i].count;
        for (j=0; j<limit; ++j) {
            comex_acc(datatype, scale, src[j], dst[j], bytes, proc, group);
        }
    }

    skip_lock = 0;
    release_remote_lock(proc);

    return COMEX_SUCCESS;
}


int comex_fence_all(comex_group_t group)
{
    return comex_wait_all(group);
}


int comex_fence_proc(int proc, comex_group_t group)
{
    return comex_wait_all(group);
}


/* comex_barrier is comex_fence_all + MPI_Barrier */
int comex_barrier(comex_group_t group)
{
    MPI_Comm comm;

    comex_fence_all(group);
    assert(COMEX_SUCCESS == comex_group_comm(group, &comm));
    MPI_Barrier(comm);

    return COMEX_SUCCESS;
}


void *comex_malloc_local(size_t size)
{
    assert(0);

    return NULL;
}


int comex_free_local(void *ptr)
{
    assert(0);

    return COMEX_SUCCESS;
}


int comex_finalize()
{
    /* it's okay to call multiple times -- extra calls are no-ops */
    if (!initialized) {
        return;
    }

    initialized = 0;

    /* Make sure that all outstanding operations are done */
    comex_wait_all(COMEX_GROUP_WORLD);
    
    /* groups */
    comex_group_finalize();

    MPI_Barrier(l_state.world_comm);

    // destroy the communicators
    MPI_Comm_free(&l_state.world_comm);

    return COMEX_SUCCESS;
}


int comex_wait_proc(int proc, comex_group_t group)
{
    assert(0);

    return COMEX_SUCCESS;
}


int comex_wait(comex_request_t* hdl)
{
    assert(0);

    return COMEX_SUCCESS;
}


int comex_test(comex_request_t* hdl, int *status)
{
    assert(0);

    *status = 0;
    return COMEX_SUCCESS;
}


int comex_wait_all(comex_group_t group)
{
    assert(0);

    return COMEX_SUCCESS;
}


int comex_nbput(
        void *src, void *dst, int bytes,
        int proc, comex_group_t group,
        comex_request_t *hdl)
{
    return comex_put(src, dst, bytes, proc, group);
}


int comex_nbget(
        void *src, void *dst, int bytes,
        int proc, comex_group_t group,
        comex_request_t *hdl)
{
    return comex_get(src, dst, bytes, proc, group);
}


int comex_nbacc(
        int datatype, void *scale,
        void *src_ptr, void *dst_ptr, int bytes,
        int proc, comex_group_t group,
        comex_request_t *hdl)
{
    return comex_accs(
            datatype, scale,
            src_ptr, NULL,
            dst_ptr, NULL,
            &bytes, 0,
            proc, group);
}


int comex_nbputs(
        void *src, int *src_stride,
        void *dst, int *dst_stride,
        int *count, int stride_levels, 
        int proc, comex_group_t group,
        comex_request_t *hdl)
{
    return comex_puts(src, src_stride, dst, dst_stride,
            count, stride_levels, proc, group);
}


int comex_nbgets(
        void *src, int *src_stride,
        void *dst, int *dst_stride,
        int *count, int stride_levels, 
        int proc, comex_group_t group,
        comex_request_t *hdl) 
{
    return comex_gets(src, src_stride, dst, dst_stride,
            count, stride_levels, proc, group);
}


int comex_nbaccs(
        int datatype, void *scale,
        void *src, int *src_stride,
        void *dst, int *dst_stride,
        int *count, int stride_levels,
        int proc, comex_group_t group,
        comex_request_t *hdl)
{
    return comex_accs(datatype, scale,
            src, src_stride, dst, dst_stride,
            count, stride_levels, proc, group);
}


int comex_nbputv(
        comex_giov_t *iov, int iov_len,
        int proc, comex_group_t group,
        comex_request_t* handle)
{
    return comex_putv(iov, iov_len, proc, group);
}


int comex_nbgetv(
        comex_giov_t *iov, int iov_len,
        int proc, comex_group_t group,
        comex_request_t* handle)
{
    return comex_getv(iov, iov_len, proc, group);
}


int comex_nbaccv(
        int datatype, void *scale,
        comex_giov_t *iov, int iov_len,
        int proc, comex_group_t group,
        comex_request_t* handle)
{
    return comex_accv(datatype, scale, iov, iov_len, proc, group);
}


int comex_rmw(
        int op, void *ploc, void *prem, int extra,
        int proc, comex_group_t group)
{
    int status;
    if (op == COMEX_FETCH_AND_ADD) {
        int tmp;
        acquire_remote_lock(proc);
        comex_get(prem, ploc, sizeof(int), proc, group);
        tmp = *(int*)ploc + extra;
        comex_put(&tmp, prem, sizeof(int), proc, group);
        release_remote_lock(proc);
    }
    else if (op == COMEX_FETCH_AND_ADD_LONG) {
        long tmp;
        acquire_remote_lock(proc);
        comex_get(prem, ploc, sizeof(long), proc, group);
        tmp = *(long*)ploc + extra;
        comex_put(&tmp, prem, sizeof(long), proc, group);
        release_remote_lock(proc);
    }
    else if (op == COMEX_SWAP) {
        int tmp;
        acquire_remote_lock(proc);
        comex_get(prem, &tmp, sizeof(int), proc, group);
        comex_put(ploc, prem, sizeof(int), proc, group);
        release_remote_lock(proc);
        *(int*)ploc = tmp;
    }
    else if (op == COMEX_SWAP_LONG) {
        long tmp;
        acquire_remote_lock(proc);
        comex_get(prem, &tmp, sizeof(long), proc, group);
        comex_put(ploc, prem, sizeof(long), proc, group);
        release_remote_lock(proc);
        *(long*)ploc = tmp;
    }
    else  {
        assert(0);
    }
    
    return COMEX_SUCCESS;
}


/* Mutex Operations */
int comex_create_mutexes(int num)
{
    assert(0);

    return COMEX_SUCCESS;
}


int comex_destroy_mutexes()
{
    assert(0);

    return COMEX_SUCCESS;
}


int comex_lock(int mutex, int proc)
{
    assert(0);

    return COMEX_SUCCESS;
}


int comex_unlock(int mutex, int proc)
{
    assert(0);

    return COMEX_SUCCESS;
}


int comex_malloc(void *ptrs[], size_t size, comex_group_t group)
{
    comex_igroup_t *igroup = NULL;
    MPI_Comm comm = MPI_COMM_NULL;
    int comm_rank = -1;
    int comm_size = -1;
    void *src_buf = NULL;

    /* preconditions */
    assert(ptrs);
   
    igroup = comex_get_igroup_from_group(group);
    comm = igroup->comm;
    assert(comm != MPI_COMM_NULL);
    MPI_Comm_rank(comm, &comm_rank);
    MPI_Comm_size(comm, &comm_size);

    /* allocate and register segment */
    ptrs[comm_rank] = comex_malloc_local(sizeof(char)*size);
  
    /* exchange buffer address */
    /* @TODO: Consider using MPI_IN_PLACE? */
    memcpy(&src_buf, &ptrs[comm_rank], sizeof(void *));
    MPI_Allgather(&src_buf, sizeof(void *), MPI_BYTE, ptrs,
            sizeof(void *), MPI_BYTE, comm);

    MPI_Barrier(comm);

    return COMEX_SUCCESS;
}


int comex_free(void *ptr, comex_group_t group)
{
    comex_igroup_t *igroup = NULL;
    MPI_Comm comm = MPI_COMM_NULL;
    int comm_rank;
    int comm_size;
    long **allgather_ptrs = NULL;

    /* preconditions */
    assert(NULL != ptr);

    igroup = comex_get_igroup_from_group(group);
    comm = igroup->comm;
    assert(comm != MPI_COMM_NULL);
    MPI_Comm_rank(comm, &comm_rank);
    MPI_Comm_size(comm, &comm_size);

    /* allocate receive buffer for exchange of pointers */
    allgather_ptrs = (long **)malloc(sizeof(void *) * comm_size);
    assert(allgather_ptrs);

    /* exchange of pointers */
    MPI_Allgather(&ptr, sizeof(void *), MPI_BYTE,
            allgather_ptrs, sizeof(void *), MPI_BYTE, comm);

    /* TODO do something useful with pointers */

    /* remove my ptr from reg cache and free ptr */
    comex_free_local(ptr);
    free(allgather_ptrs);

    /* Is this needed? */
    MPI_Barrier(comm);

    return COMEX_SUCCESS;
}


static void acquire_remote_lock(int proc)
{
    assert(0);
}


static void release_remote_lock(int proc)
{
    assert(0);
}


static inline void acc(
        int datatype, int count, void *get_buf,
        void *src_ptr, long src_idx, void *scale)
{
#define EQ_ONE_REG(A) ((A) == 1.0)
#define EQ_ONE_CPL(A) ((A).real == 1.0 && (A).imag == 0.0)
#define IADD_REG(A,B) (A) += (B)
#define IADD_CPL(A,B) (A).real += (B).real; (A).imag += (B).imag
#define IADD_SCALE_REG(A,B,C) (A) += (B) * (C)
#define IADD_SCALE_CPL(A,B,C) (A).real += ((B).real*(C).real) - ((B).imag*(C).imag);\
                              (A).imag += ((B).real*(C).imag) + ((B).imag*(C).real);
#define ACC(WHICH, COMEX_TYPE, C_TYPE)                                  \
    if (datatype == COMEX_TYPE) {                                       \
        int m;                                                          \
        int m_lim = count/sizeof(C_TYPE);                               \
        C_TYPE *iterator = (C_TYPE *)get_buf;                           \
        C_TYPE *value = (C_TYPE *)((char *)src_ptr + src_idx);          \
        C_TYPE calc_scale = *(C_TYPE *)scale;                           \
        if (EQ_ONE_##WHICH(calc_scale)) {                               \
            for (m = 0 ; m < m_lim; ++m) {                              \
                IADD_##WHICH(iterator[m], value[m]);                    \
            }                                                           \
        }                                                               \
        else {                                                          \
            for (m = 0 ; m < m_lim; ++m) {                              \
                IADD_SCALE_##WHICH(iterator[m], value[m], calc_scale);  \
            }                                                           \
        }                                                               \
    } else
    ACC(REG, COMEX_ACC_DBL, double)
    ACC(REG, COMEX_ACC_FLT, float)
    ACC(REG, COMEX_ACC_INT, int)
    ACC(REG, COMEX_ACC_LNG, long)
    ACC(CPL, COMEX_ACC_DCP, DoubleComplex)
    ACC(CPL, COMEX_ACC_CPL, SingleComplex)
    {
        assert(0);
    }
#undef ACC
#undef EQ_ONE_REG
#undef EQ_ONE_CPL
#undef IADD_REG
#undef IADD_CPL
#undef IADD_SCALE_REG
#undef IADD_SCALE_CPL
}

