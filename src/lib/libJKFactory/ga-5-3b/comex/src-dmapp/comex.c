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
#include <hugetlbfs.h>

/* 3rd party headers */
#include <mpi.h>
#include <dmapp.h>

/* our headers */
#include "comex.h"
#include "comex_impl.h"
/*#include "clusterinfo.h"*/
#include "groups.h"
#include "reg_cache.h"
#include "acc.h"

/* Cray */
#if HAVE_DMAPP_LOCK_DESC_T && HAVE_DMAPP_LOCK_HANDLE_T
#   define HAVE_DMAPP_LOCK 1
#endif

#define DEBUG 0


#if HAVE_DMAPP_LOCK
// COMEX_MAX_LOCKS mirrors the default DMAPP_MAX_LOCKS limit
// Larger values of COMEX_MAX_LOCKS will require DMAPP_MAX_LOCKS be set at runtime.
// DMAPP_MAX_LOCKS has a maxium value of 1023
#define COMEX_MAX_LOCKS 128  
static dmapp_lock_desc_t   lock_desc[COMEX_MAX_LOCKS];
static dmapp_lock_handle_t lock_handle[COMEX_MAX_LOCKS];
#endif


/* exported state */
local_state l_state;
int comex_me=-1;
int comex_nproc=-1;

/* static state */
static int  initialized=0;                  /* for comex_initialized(), 0=false */
static int  total_outstanding=0;
static int  max_outstanding_nb=MAX_NB_OUTSTANDING;
static int  malloc_is_using_huge_pages=0;   /* from env var, 0=false */
static int  comex_is_using_huge_pages=0;    /* from env var, 0=false */
static long hugetlb_default_page_size=0;    /* from env var, in bytes */
static long sc_page_size=0;                 /* from sysconf, in bytes */
static long hugepagesize=0;                 /* from libhugetlbfs, in bytes */
static long comex_page_size=0;              /* page size consensus, in bytes */
static char skip_lock=0;                    /* don't acquire or release lock */
static char skip_sync=0;                    /* don't sync implicit nb requests */

/* static function declarations */
static void  check_envs(void);
static void  create_dmapp_locks(void);
static void  destroy_dmapp_locks(void);
static void  dmapp_alloc_buf(void);
static void  dmapp_free_buf(void);
static void  dmapp_initialize(void);
static void  dmapp_network_lock(int proc);
static void  dmapp_network_unlock(int proc);
static void  dmapp_terminate(void);
static void  increment_total_outstanding(void);
static void  wait_and_clear_total_outstanding(void);
static void  my_free(void *ptr);
static void* my_malloc(size_t size);
static int   my_memalign(void **memptr, size_t alignment, size_t size);
static int   comex_get_nbi(void *src, void *dst, int bytes, int proc);
static int   comex_get_nb(void *src, void *dst, int bytes, int proc, dmapp_syncid_handle_t *handle);
static int   comex_put_nbi(void *src, void *dst, int bytes, int proc);
static int   comex_put_nb(void *src, void *dst, int bytes, int proc, dmapp_syncid_handle_t *handle);
static void* _comex_malloc_local(size_t size, dmapp_seg_desc_t *seg);


static void* my_malloc(size_t size)
{
    void *memptr=NULL;

#if DEBUG
    if (0 == l_state.rank) {
        printf("my_malloc(%lu)\n", (long unsigned)size);
    }
#endif

#if HAVE_LIBHUGETLBFS
    if (malloc_is_using_huge_pages) {
        memptr = malloc(size);
    }
    else if (comex_is_using_huge_pages) {
        memptr = get_hugepage_region(size, GHR_DEFAULT);
    }
    else {
        memptr = malloc(size);
    }
#else
    memptr = malloc(size);
#endif

    /* postconditions */
    assert(memptr);

    return memptr;
}


static void my_free(void *ptr)
{
#if DEBUG
    if (0 == l_state.rank) {
        printf("my_free(%p)\n", ptr);
    }
#endif

#if HAVE_LIBHUGETLBFS
    if (malloc_is_using_huge_pages) {
        free(ptr);
    }
    else if (comex_is_using_huge_pages) {
        free_hugepage_region(ptr);
    }
    else {
        free(ptr);
    }
#else
    free(ptr);
#endif
}


static int my_memalign(void **memptr, size_t alignment, size_t size)
{
    int status = 0;

#if DEBUG
    if (0 == l_state.rank) {
        printf("my_memalign(%lu)\n", (long unsigned)size);
    }
#endif

    /* preconditions */
    assert(memptr);

#if HAVE_LIBHUGETLBFS
    if (malloc_is_using_huge_pages) {
        status = posix_memalign(memptr, alignment, size);
    }
    else if (comex_is_using_huge_pages) {
        *memptr = get_hugepage_region(size, GHR_DEFAULT);
    }
    else {
        status = posix_memalign(memptr, alignment, size);
    }
#else
    status = posix_memalign(memptr, alignment, size);
#endif

    /* postconditions */
    assert(*memptr);

    return status;
}


static void increment_total_outstanding(void)
{
    ++total_outstanding;

    if (total_outstanding == max_outstanding_nb) {
        wait_and_clear_total_outstanding();
    }
}


static void wait_and_clear_total_outstanding(void)
{
    int status;
    status = dmapp_gsync_wait();
    assert(status == DMAPP_RC_SUCCESS);
    total_outstanding = 0;
}


/* The blocking implementations should use blocking DMAPP calls */
int comex_put(void *src, void *dst, int bytes, int proc, comex_group_t group)
{
    int status;
    status = comex_put_nbi(src, dst, bytes, proc);
    assert(status == DMAPP_RC_SUCCESS);
    comex_wait_proc(proc, group);
    return COMEX_SUCCESS;
}


int comex_get(void *src, void *dst, int bytes, int proc, comex_group_t group)
{
    int status;
    status = comex_get_nbi(src, dst, bytes, proc);
    assert(status == DMAPP_RC_SUCCESS);
    comex_wait_proc(proc, group);
    return COMEX_SUCCESS;
}


/* The blocking implementations should use blocking DMAPP calls */
static int comex_put_nbi(void *src, void *dst, int bytes, int proc)
{
    int status = DMAPP_RC_SUCCESS;
    int nelems = bytes;
    int type = DMAPP_BYTE;
    int failure_observed = 0;
    reg_entry_t *dst_reg = NULL;
    reg_entry_t *src_reg = NULL;

    /* Corner case */
    if (proc == l_state.rank) {
        memcpy(dst, src, bytes);
        return status;
    }

    /* If the number of bytes is even, use Double word datatype,
     * DMAPP_BYTE performance is much worse */
    if (0 == bytes%16) {
        nelems = bytes/16;
        type = DMAPP_DQW;
    }
    else if (0 == bytes%8) {
        nelems = bytes/8;
        type = DMAPP_QW;
    }
    else if (0 == bytes%4) {
        nelems = bytes/4;
        type = DMAPP_DW;
    }

    /* Find the dmapp seg desc */
    dst_reg = reg_cache_find(proc, dst, bytes);
    assert(dst_reg);

    src_reg = reg_cache_find(l_state.rank, src, bytes);

    status = dmapp_put_nbi(dst, &(dst_reg->mr), proc, src, nelems, type);
    increment_total_outstanding();
    if (status != DMAPP_RC_SUCCESS) {
        failure_observed = 1;
    }

    /* Fallback */
    if (failure_observed) {
        comex_wait_all(COMEX_GROUP_WORLD);
        assert(bytes <= l_state.put_buf_len);
        memcpy(l_state.put_buf, src, bytes);
        status = dmapp_put_nbi(dst, &(dst_reg->mr),
                proc, l_state.put_buf, nelems, type);
        increment_total_outstanding();
        comex_wait_all(COMEX_GROUP_WORLD);

        /* Fallback must work correctly */
        assert(status == DMAPP_RC_SUCCESS);
    }

    return status;
}


/* The blocking implementations should use blocking DMAPP calls */
static int comex_put_nb(void *src, void *dst, int bytes, int proc, dmapp_syncid_handle_t *handle)
{
    int status = DMAPP_RC_SUCCESS;
    int nelems = bytes;
    int type = DMAPP_BYTE;
    int failure_observed = 0;
    reg_entry_t *dst_reg = NULL;
    reg_entry_t *src_reg = NULL;

    /* If the number of bytes is even, use Double word datatype,
     * DMAPP_BYTE performance is much worse */
    if (0 == bytes%16) {
        nelems = bytes/16;
        type = DMAPP_DQW;
    }
    else if (0 == bytes%8) {
        nelems = bytes/8;
        type = DMAPP_QW;
    }
    else if (0 == bytes%4) {
        nelems = bytes/4;
        type = DMAPP_DW;
    }

    /* Find the dmapp seg desc */
    dst_reg = reg_cache_find(proc, dst, bytes);
    assert(dst_reg);

    src_reg = reg_cache_find(l_state.rank, src, bytes);

    status = dmapp_put_nb(dst, &(dst_reg->mr), proc, src, nelems, type, handle);
    assert(status == DMAPP_RC_SUCCESS);

    return status;
}


static int comex_get_nb(void *src, void *dst, int bytes, int proc, dmapp_syncid_handle_t *handle)
{
    int status = DMAPP_RC_SUCCESS;
    int nelems = bytes;
    int type = DMAPP_BYTE;
    int failure_observed = 0;
    reg_entry_t *dst_reg = NULL;
    reg_entry_t *src_reg = NULL;

    /* If the number of bytes is even, use Double word datatype,
     * DMAPP_BYTE performance is much worse */
    if (0 == bytes%16) {
        nelems = bytes/16;
        type = DMAPP_DQW;
    }
    else if (0 == bytes%8) {
        nelems = bytes/8;
        type = DMAPP_QW;
    }
    else if (0 == bytes%4) {
        nelems = bytes/4;
        type = DMAPP_DW;
    }

    /* Find the dmapp seg desc */
    dst_reg = reg_cache_find(proc, src, bytes);
    assert(dst_reg);
    
    src_reg = reg_cache_find(l_state.rank, dst, bytes);

    status = dmapp_get_nb(dst, src, &(dst_reg->mr), proc, nelems, type, handle);
    assert(status == DMAPP_RC_SUCCESS);

    return COMEX_SUCCESS;
}


static int comex_get_nbi(void *src, void *dst, int bytes, int proc)
{
    int status = DMAPP_RC_SUCCESS;
    int nelems = bytes;
    int type = DMAPP_BYTE;
    int failure_observed = 0;
    reg_entry_t *dst_reg = NULL;
    reg_entry_t *src_reg = NULL;

    /* Corner case */
    if (proc == l_state.rank) {
        memcpy(dst, src, bytes);
        return status;
    }

    /* If the number of bytes is even, use Double word datatype,
     * DMAPP_BYTE performance is much worse */
    if (0 == bytes%16) {
        nelems = bytes/16;
        type = DMAPP_DQW;
    }
    else if (0 == bytes%8) {
        nelems = bytes/8;
        type = DMAPP_QW;
    }
    else if (0 == bytes%4) {
        nelems = bytes/4;
        type = DMAPP_DW;
    }

    /* Find the dmapp seg desc */
    dst_reg = reg_cache_find(proc, src, bytes);
    assert(dst_reg);
    
    src_reg = reg_cache_find(l_state.rank, dst, bytes);

    status = dmapp_get_nbi(dst, src, &(dst_reg->mr),
            proc, nelems, type);
    increment_total_outstanding();
    if (status != DMAPP_RC_SUCCESS) {
        failure_observed = 1;
    }
    
    /* Fallback */
    if (failure_observed) {    
        comex_wait_all(COMEX_GROUP_WORLD);
        assert(bytes <= l_state.get_buf_len);
        status = dmapp_get(l_state.get_buf, src, &(dst_reg->mr),
                proc, nelems, type);
        memcpy(dst, l_state.get_buf, bytes);
    }

    /* Original or fallback must work correctly */
    assert(status == DMAPP_RC_SUCCESS);

    return status;
}


static void dmapp_network_lock(int proc)
{
    int dmapp_status;

#if HAVE_DMAPP_LOCK
    dmapp_lock_acquire( &lock_desc[0], &(l_state.job.data_seg), proc, 0, &lock_handle[0]);
#else
    reg_entry_t *dst_reg= reg_cache_find(proc, 
            l_state.atomic_lock_buf[proc], sizeof(long));

    assert(dst_reg);

    do {    
        dmapp_status = dmapp_acswap_qw(l_state.local_lock_buf, 
                l_state.atomic_lock_buf[proc],
                &(dst_reg->mr),
                proc, 0, l_state.rank + 1);

        assert(dmapp_status == DMAPP_RC_SUCCESS);
    }
    while(*(l_state.local_lock_buf) != 0);
#endif
}


static void dmapp_network_unlock(int proc)
{
    int dmapp_status;

# if HAVE_DMAPP_LOCK
    dmapp_lock_release( lock_handle[0], 0 );
#else
    reg_entry_t *dst_reg= reg_cache_find(proc, 
            l_state.atomic_lock_buf[proc], sizeof(long));

    assert(dst_reg);

    do {
        dmapp_status = dmapp_acswap_qw(l_state.local_lock_buf, 
                l_state.atomic_lock_buf[proc],
                &(dst_reg->mr),
                proc, l_state.rank + 1, 0);
        assert(dmapp_status == DMAPP_RC_SUCCESS);
    } while (*(l_state.local_lock_buf) != (unsigned long)(l_state.rank + 1));
#endif
}


int comex_puts(void *src_ptr, int src_stride_ar[/*stride_levels*/],
                void *dst_ptr, int dst_stride_ar[/*stride_levels*/],
                int count[/*stride_levels+1*/], int stride_levels, int proc, comex_group_t group)
{
    int i, j;
    long src_idx, dst_idx;  /* index offset of current block position to ptr */
    int n1dim;  /* number of 1 dim block */
    int src_bvalue[7], src_bunit[7];
    int dst_bvalue[7], dst_bunit[7];
    int dmapp_status;

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
        
        dmapp_status = comex_put_nbi((char *)src_ptr + src_idx, 
                (char *)dst_ptr + dst_idx, count[0], proc);
        assert(dmapp_status == DMAPP_RC_SUCCESS);
    }

    if (0 == skip_sync) {
        comex_wait_proc(proc, group);
    }

    return COMEX_SUCCESS;
}


int comex_gets(void *src_ptr, int src_stride_ar[/*stride_levels*/],
                void *dst_ptr, int dst_stride_ar[/*stride_levels*/],
                int count[/*stride_levels+1*/], int stride_levels, int proc, comex_group_t group)
{
    int i, j;
    long src_idx, dst_idx;  /* index offset of current block position to ptr */
    int n1dim;  /* number of 1 dim block */
    int src_bvalue[7], src_bunit[7];
    int dst_bvalue[7], dst_bunit[7];
    int dmapp_status;

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
        
        dmapp_status = comex_get_nbi((char *)src_ptr + src_idx, 
                (char *)dst_ptr + dst_idx, count[0], proc);
        assert(dmapp_status == DMAPP_RC_SUCCESS);
    }
    
    if (0 == skip_sync) {
        comex_wait_proc(proc, group);
    }
    
    return COMEX_SUCCESS;
}


int comex_acc(int datatype, void *scale,
               void *src_ptr, 
               void *dst_ptr, 
               int bytes, int proc, comex_group_t group)
{

    comex_accs(datatype, scale, src_ptr, NULL, dst_ptr, 
            NULL, &bytes, 0, proc, group);
    return COMEX_SUCCESS;
}


int comex_accs(int datatype, void *scale,
                void *src_ptr, int src_stride_ar[/*stride_levels*/],
                void *dst_ptr, int dst_stride_ar[/*stride_levels*/],
                int count[/*stride_levels+1*/], int stride_levels, int proc, comex_group_t group)
{
    int i, j;
    long src_idx, dst_idx;  /* index offset of current block position to ptr */
    int n1dim;  /* number of 1 dim block */
    int src_bvalue[7], src_bunit[7];
    int dst_bvalue[7], dst_bunit[7];
    int sizetogetput;
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

    sizetogetput = count[0];

    if (0 == skip_lock) {
        // grab the atomics lock
        dmapp_network_lock(proc);
    }

#if PIPELINED_ACCUMULATE
    if (sizetogetput > l_state.pipe_acc_buf_len)
#endif
    {
        /* fall back to sequential with newly allocated buffer */
        if (sizetogetput <= l_state.acc_buf_len) {
            get_buf = l_state.acc_buf;
        }
        else {
            get_buf = (char *)my_malloc(sizeof(char) * sizetogetput);
        }
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
            comex_get((char *)dst_ptr + dst_idx, get_buf, sizetogetput, proc, group);

            _acc(datatype, count[0], get_buf, ((char*)src_ptr)+src_idx, scale);

            // Write back
            comex_put(get_buf, (char *)dst_ptr + dst_idx, sizetogetput, proc, group);
        }
        if (sizetogetput > l_state.acc_buf_len) {
            // unregister temp buffer ?  TODO consider keeping temp buf around
            // in case another large request comes along?
            free(get_buf);
        }
    }
#if PIPELINED_ACCUMULATE
    else {
        printf("pipelined acc\n");
        /* pipelined protocol, sort of */
        long pipe_src_idx[PIPELINED_MAX_BUFFERS];
        long pipe_dst_idx[PIPELINED_MAX_BUFFERS];
        dmapp_syncid_handle_t pipe_handle[PIPELINED_MAX_BUFFERS];
        char pipe_state[PIPELINED_MAX_BUFFERS];
        long pipe_index = 0;
        char done = 0;

        for (i=0; i<PIPELINED_MAX_BUFFERS; i++) {
            /* 0 means empty
             * 1 means outstanding get
             * 2 means outstanding put */
            pipe_state[i] = 0;
        }

        /* issue new nbget requests as available */
        for(i=0; i<n1dim; i++) {
            int found = 0;

            /* find the first completed or available slot */
            found = 0;
            while (found == 0) {
                if (pipe_state[pipe_index] == 0) {
                    found = 1;
                }
                else if (pipe_state[pipe_index] == 1) {
                    int test;
                    dmapp_syncid_test(&pipe_handle[pipe_index], &test);
                    assert(test == 0 || test == 1);
                    if (test == 1) {
                        /* get completed, peform acc, then nb put */
                        _acc(datatype, sizetogetput,
                                l_state.pipe_acc_buf[pipe_index],
                                ((char*)src_ptr)+pipe_src_idx[pipe_index],
                                scale);
                        comex_put_nb(l_state.pipe_acc_buf[pipe_index],
                                (char *)dst_ptr + pipe_dst_idx[pipe_index],
                                sizetogetput, proc,
                                &pipe_handle[pipe_index]);
                        pipe_state[pipe_index] = 2;
                    }
                    pipe_index = (pipe_index+1) % PIPELINED_MAX_BUFFERS;
                }
                else if (pipe_state[pipe_index] == 2) {
                    int test;
                    dmapp_syncid_test(&pipe_handle[pipe_index], &test);
                    assert(test == 0 || test == 1);
                    if (test == 1) {
                        /* put completed */
                        found = 1;
                        pipe_state[pipe_index] = 0;
                    }
                    else {
                        pipe_index = (pipe_index+1) % PIPELINED_MAX_BUFFERS;
                    }
                }
                else {
                    assert(0);
                }
            }

            /* calculate the src_idx */
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

            /* calculate the dst_idx */
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
            comex_get_nb((char *)dst_ptr + dst_idx,
                    l_state.pipe_acc_buf[pipe_index],
                    sizetogetput, proc, &pipe_handle[pipe_index]);
            pipe_dst_idx[pipe_index] = dst_idx;
            pipe_src_idx[pipe_index] = src_idx;
            pipe_state[pipe_index] = 1;
            pipe_index = (pipe_index+1) % PIPELINED_MAX_BUFFERS;
        }

        /* no more gets to issue, but we may have outstanding gets/puts */
        do {
            done = 1;
            for (pipe_index=0; pipe_index<PIPELINED_MAX_BUFFERS; ++pipe_index) {
                if (pipe_state[pipe_index] == 0) {
                    // okay
                }
                else if (pipe_state[pipe_index] == 1) {
                    int test;
                    dmapp_syncid_test(&pipe_handle[pipe_index], &test);
                    assert(test == 0 || test == 1);
                    if (test == 1) {
                        /* get completed, peform acc, then nb put */
                        _acc(datatype, sizetogetput,
                                l_state.pipe_acc_buf[pipe_index],
                                ((char*)src_ptr)+pipe_src_idx[pipe_index],
                                scale);
                        comex_put_nb(l_state.pipe_acc_buf[pipe_index],
                                (char *)dst_ptr + pipe_dst_idx[pipe_index],
                                sizetogetput, proc,
                                &pipe_handle[pipe_index]);
                        pipe_state[pipe_index] = 2;
                    }
                    done = 0;
                }
                else if (pipe_state[pipe_index] == 2) {
                    int test;
                    dmapp_syncid_test(&pipe_handle[pipe_index], &test);
                    assert(test == 0 || test == 1);
                    if (test == 1) {
                        /* put completed */
                        pipe_state[pipe_index] = 0;
                    }
                    else {
                        done = 0;
                    }
                }
                else {
                    assert(0);
                }
            }
        } while (done == 0);
    }
#endif

    if (0 == skip_lock) {
        // ungrab the lock
        dmapp_network_unlock(proc);
    }

    if (sizetogetput > l_state.acc_buf_len)
        my_free(get_buf);

    return COMEX_SUCCESS;
}


int comex_fence_all(comex_group_t group)
{
    comex_wait_all(group);
    /* noop for DMAPP */
    return COMEX_SUCCESS;
}


int comex_fence_proc(int proc, comex_group_t group)
{
    comex_wait_all(group);
    /* noop for DMAPP */
    return COMEX_SUCCESS;
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


void comex_error(char *msg, int code)
{
    if (0 == l_state.rank)
        fprintf(stderr,"Received an Error in Communication\n");
    
    MPI_Abort(l_state.world_comm, code);
}


static void* _comex_malloc_local(size_t size, dmapp_seg_desc_t *seg)
{
    void *ptr;
    int rc;
    int status;

    rc = my_memalign(&ptr, comex_page_size, sizeof(char)*size);
    assert(0 == rc);
    assert(ptr);

    status = dmapp_mem_register(ptr, size, seg);
    assert(status == DMAPP_RC_SUCCESS);
#if DEBUG
    printf("[%d] _comex_malloc_local ptr=%p size=%zu\n",
            l_state.rank, ptr, size);
    printf("[%d] _comex_malloc_local seg=%p size=%zu\n",
            l_state.rank, seg->addr, seg->len);
#endif
#if 0
    assert(seg->addr == ptr);
    assert(seg->len == size); /* @TODO this failed! */
#endif
    reg_cache_insert(l_state.rank, ptr, size, *seg);

    return ptr;
}


void *comex_malloc_local(size_t size)
{
    void *ptr = NULL;
    dmapp_seg_desc_t seg;

    ptr = _comex_malloc_local(size, &seg);

    return ptr;
}


int comex_free_local(void *ptr)
{
    reg_return_t status = RR_FAILURE;

    /* preconditions */
    assert(NULL != ptr);

    /* remove from reg cache */
    status = reg_cache_delete(l_state.rank, ptr);
    assert(RR_SUCCESS == status);

    /* free the memory */
    my_free(ptr);

    return COMEX_SUCCESS;
}


static void destroy_dmapp_locks(void)
{
#if DMAPP_LOCK
#else
    if (l_state.local_lock_buf)
            comex_free_local(l_state.local_lock_buf);

    if (l_state.atomic_lock_buf)
            comex_free(l_state.atomic_lock_buf[l_state.rank], COMEX_GROUP_WORLD);
#endif
}


static void create_dmapp_locks(void)
{
#if DMAPP_LOCK
    bzero(lock_desc, sizeof(lock_desc));
#else
    l_state.local_lock_buf = comex_malloc_local(sizeof(long));
    assert(l_state.local_lock_buf);

    l_state.atomic_lock_buf =
        (unsigned long**)my_malloc(l_state.size * sizeof(unsigned long*));
    assert(l_state.atomic_lock_buf);

    comex_malloc((void**)l_state.atomic_lock_buf, sizeof(long), COMEX_GROUP_WORLD);

    *(long *)(l_state.atomic_lock_buf[l_state.rank]) = 0;
    *(long *)(l_state.local_lock_buf) = 0;
#endif

    MPI_Barrier(l_state.world_comm);
}


static void dmapp_alloc_buf(void)
{
    int i;

    // FAILURE_BUFSIZE should be some multiple of our page size?

#if PIPELINED_ACCUMULATE
    l_state.pipe_acc_buf = malloc(PIPELINED_MAX_BUFFERS * sizeof(void*));
    l_state.pipe_acc_buf_len = comex_page_size;
    for (i=0; i<PIPELINED_MAX_BUFFERS; ++i) {
        l_state.pipe_acc_buf[i] = comex_malloc_local(l_state.pipe_acc_buf_len);
        assert(l_state.pipe_acc_buf[i]);
    }
#endif

    //l_state.acc_buf_len = FAILURE_BUFSIZE;
    l_state.acc_buf_len = comex_page_size;
    l_state.acc_buf = comex_malloc_local( l_state.acc_buf_len);
    assert(l_state.acc_buf);

    //l_state.put_buf_len = FAILURE_BUFSIZE;
    l_state.put_buf_len = comex_page_size;
    l_state.put_buf = comex_malloc_local(l_state.put_buf_len);
    assert(l_state.put_buf);

    //l_state.get_buf_len = FAILURE_BUFSIZE;
    l_state.get_buf_len = comex_page_size;
    l_state.get_buf = comex_malloc_local(l_state.get_buf_len);
    assert(l_state.get_buf);
}


static void dmapp_free_buf(void)
{
    comex_free_local(l_state.acc_buf);
    comex_free_local(l_state.put_buf);
    comex_free_local(l_state.get_buf);
}


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
    comex_me = l_state.rank;

    /* World Size */
    status = MPI_Comm_size(l_state.world_comm, &(l_state.size));
    assert(MPI_SUCCESS == status);
    comex_nproc = l_state.size;
    
    /* compile-time sanity check */
#if !(_POSIX_C_SOURCE >= 200112L || _XOPEN_SOURCE >= 600)
#   error posix_memalign *NOT* available
#endif

    /* groups */
    comex_group_init();

    /* Initialize */
    dmapp_initialize();

    /* mutexes */
    l_state.mutexes = NULL;
    l_state.local_mutex = NULL;
    l_state.num_mutexes = NULL;

    /* cluster info */
    /*comex_init_clusinfo();*/

    /* Synch - Sanity Check */
    MPI_Barrier(l_state.world_comm);

    return COMEX_SUCCESS;
}


int comex_init_args(int *argc, char ***argv)
{
    int rc;
    int init_flag;
    
    MPI_Initialized(&init_flag);
    
    if(!init_flag)
        MPI_Init(argc, argv);
    
    rc = comex_init();
    return rc;
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
    
    dmapp_terminate();

    /* groups */
    comex_group_finalize();

    MPI_Barrier(l_state.world_comm);

    // destroy the communicators
    MPI_Comm_free(&l_state.world_comm);

    return COMEX_SUCCESS;
}


int comex_nbput(void *src, void *dst, int bytes, int proc, comex_group_t group, comex_request_t *hdl)
{
    int rc;
    rc = comex_put_nbi(src, dst, bytes, proc);
    return (rc == DMAPP_RC_SUCCESS) ? COMEX_SUCCESS : COMEX_FAILURE;
}


int comex_nbget(void *src, void *dst, int bytes, int proc, comex_group_t group, comex_request_t *hdl)
{
    int rc;
    rc = comex_get_nbi(src, dst, bytes, proc);
    return (rc == DMAPP_RC_SUCCESS) ? COMEX_SUCCESS : COMEX_FAILURE;
}


int comex_wait_proc(int proc, comex_group_t group)
{
    wait_and_clear_total_outstanding();
    return COMEX_SUCCESS;
}


int comex_wait(comex_request_t* hdl)
{
    wait_and_clear_total_outstanding();
    return COMEX_SUCCESS;
}


int comex_test(comex_request_t* hdl, int *status)
{
    wait_and_clear_total_outstanding();
    *status = 0;
    return COMEX_SUCCESS;
}


int comex_wait_all(comex_group_t group)
{
    wait_and_clear_total_outstanding();
    return COMEX_SUCCESS;
}


int comex_nbputs(
        void *src, int *src_stride,
        void *dst, int *dst_stride,
        int *count, int stride_levels, 
        int proc, comex_group_t group,
        comex_request_t *hdl)
{
    int rc;
    assert(0 == skip_sync);
    skip_sync = 1;
    rc = comex_puts(src, src_stride, dst, dst_stride,
            count, stride_levels, proc, group);
    skip_sync = 0;
    return rc;
}


int comex_nbgets(
        void *src, int *src_stride,
        void *dst, int *dst_stride,
        int *count, int stride_levels, 
        int proc, comex_group_t group,
        comex_request_t *hdl) 
{
    int rc;
    assert(0 == skip_sync);
    skip_sync = 1;
    rc = comex_gets(src, src_stride, dst, dst_stride,
            count, stride_levels, proc, group);
    skip_sync = 0;
    return rc;
}


int comex_nbaccs(
        int datatype, void *scale,
        void *src, int *src_stride,
        void *dst, int *dst_stride,
        int *count, int stride_levels,
        int proc, comex_group_t group,
        comex_request_t *hdl)
{
    int rc;
    rc = comex_accs(datatype, scale,
            src, src_stride, dst, dst_stride,
            count, stride_levels, proc, group);
    return rc;
}


/* Vector Calls */


int comex_putv(comex_giov_t *iov, int iov_len, int proc, comex_group_t group)
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
            status = comex_put_nbi(src[j], dst[j], bytes, proc);
            assert(status == DMAPP_RC_SUCCESS);
        }
    }
    if (0 == skip_sync) {
        comex_wait_proc(proc, group);
    }
    return COMEX_SUCCESS;
}


int comex_getv(comex_giov_t *iov, int iov_len, int proc, comex_group_t group)
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
            status = comex_get_nbi(src[j], dst[j], bytes, proc);
            assert(status == DMAPP_RC_SUCCESS);
        }
    }
    if (0 == skip_sync) {
        comex_wait_proc(proc, group);
    }
    return COMEX_SUCCESS;
}


int comex_accv(int datatype, void *scale, comex_giov_t *iov, 
        int iov_len, int proc, comex_group_t group)
{
    int i;
    skip_lock = 1;
    // grab the atomics lock
    dmapp_network_lock(proc);
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
    // ungrab the lock
    dmapp_network_unlock(proc);
    return COMEX_SUCCESS;
}


int comex_nbputv(comex_giov_t *iov, int iov_len, int proc, comex_group_t group, comex_request_t* handle)
{
    int rc;
    assert(0 == skip_sync);
    skip_sync = 1;
    rc = comex_putv(iov, iov_len, proc, group);
    skip_sync = 0;
    return rc;
}


int comex_nbgetv(comex_giov_t *iov, int iov_len, int proc, comex_group_t group, comex_request_t* handle)
{
    int rc;
    assert(0 == skip_sync);
    skip_sync = 1;
    rc = comex_getv(iov, iov_len, proc, group);
    skip_sync = 0;
    return rc;
}


int comex_nbaccv(int datatype, void *scale, comex_giov_t *iov, 
        int iov_len, int proc, comex_group_t group, comex_request_t* handle)
{
    int rc;
    rc = comex_accv(datatype, scale, iov, iov_len, proc, group);
    return rc;
}


int comex_rmw(int op, void *ploc, void *prem, int extra, int proc, comex_group_t group)
{
    int status;
    if (op == COMEX_FETCH_AND_ADD) {
        /* dmapp doesn't have atomic fadd for int */
        int tmp;
        dmapp_network_lock(proc);
        comex_get(prem, ploc, sizeof(int), proc, group);
        tmp = *(int*)ploc + extra;
        comex_put(&tmp, prem, sizeof(int), proc, group);
        dmapp_network_unlock(proc);
    }
    else if (op == COMEX_FETCH_AND_ADD_LONG) {
#if 1
        reg_entry_t *rem_reg = reg_cache_find(proc, prem, sizeof(long));
        assert(rem_reg);
        status = dmapp_afadd_qw(ploc, prem, &(rem_reg->mr), proc, extra);
        assert(status == DMAPP_RC_SUCCESS);
#else
        long tmp;
        dmapp_network_lock(proc);
        comex_get(prem, ploc, sizeof(long), proc, group);
        tmp = *(long*)ploc + extra;
        comex_put(&tmp, prem, sizeof(long), proc, group);
        dmapp_network_unlock(proc);
#endif
    }
    else if (op == COMEX_SWAP) {
        /* dmapp doesn't have atomic swap for int */
        int tmp;
        dmapp_network_lock(proc);
        comex_get(prem, &tmp, sizeof(int), proc, group);
        comex_put(ploc, prem, sizeof(int), proc, group);
        dmapp_network_unlock(proc);
        *(int*)ploc = tmp;
    }
    else if (op == COMEX_SWAP_LONG) {
        /* dmapp has atomic cswap for long, but it's non-blocking */
        long tmp;
        dmapp_network_lock(proc);
        comex_get(prem, &tmp, sizeof(long), proc, group);
        comex_put(ploc, prem, sizeof(long), proc, group);
        dmapp_network_unlock(proc);
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
    int i=0;

    assert(NULL == l_state.mutexes);
    assert(NULL == l_state.local_mutex);
    assert(NULL == l_state.num_mutexes);

    /* every process knows how many mutexes created on every process */
    l_state.num_mutexes = (unsigned int*)my_malloc(l_state.size * sizeof(unsigned int));
    assert(l_state.num_mutexes);
    /* gather the counts */
    MPI_Allgather(&num, 1, MPI_INT,
            l_state.num_mutexes, 1, MPI_UNSIGNED, l_state.world_comm);

    /* create the 1 element buffer to hold a remote mutex */
    l_state.local_mutex = comex_malloc_local(sizeof(unsigned long));
    assert(l_state.local_mutex);
    /* init the local mutex holder to rank+1, indicating no mutex is held */
    *(unsigned long *)(l_state.local_mutex) = l_state.rank+1;
    MPI_Barrier(l_state.world_comm);

    /* create all of the mutexes */
    l_state.mutexes = (unsigned long**)my_malloc(l_state.size * sizeof(unsigned long*));
    assert(l_state.mutexes);
    comex_malloc((void**)l_state.mutexes, num*sizeof(unsigned long), COMEX_GROUP_WORLD);
    /* init all of my mutexes to 0 */
    for (i=0; i<num; ++i) {
        l_state.mutexes[l_state.rank][i] = 0;
    }

    MPI_Barrier(l_state.world_comm);

    return COMEX_SUCCESS;
}


int comex_destroy_mutexes()
{
    MPI_Barrier(l_state.world_comm);

    /* you cannot free mutexes if one is in use */
    assert(*((unsigned long *)l_state.local_mutex)
            == (unsigned long)(l_state.rank+1));
#ifndef NDEBUG
    {
        unsigned int i;
        for (i=0; i<l_state.num_mutexes[l_state.rank]; ++i) {
            unsigned long *mutexes = l_state.mutexes[l_state.rank];
            assert(mutexes[i] == 0);
        }
    }
#endif

    /* destroy mutex counts */
    my_free(l_state.num_mutexes);
    l_state.num_mutexes = NULL;

    /* destroy the 1 element buffer holding a remote mutex */
    comex_free_local(l_state.local_mutex);
    l_state.local_mutex = NULL;

    /* destroy the mutexes */
    comex_free(l_state.mutexes[l_state.rank], COMEX_GROUP_WORLD);
    my_free(l_state.mutexes);
    l_state.mutexes = NULL;

    MPI_Barrier(l_state.world_comm);

    return COMEX_SUCCESS;
}


int comex_lock(int mutex, int proc)
{
    int dmapp_status;
    reg_entry_t *dst_reg = NULL;

    /* preconditions */
    assert(0 <= proc && proc < l_state.size);
    assert(0 <= mutex && ((unsigned int)mutex) < l_state.num_mutexes[proc]);

    /* locate remote lock */
    dst_reg = reg_cache_find(proc, &(l_state.mutexes[proc][mutex]), sizeof(unsigned long));
    assert(dst_reg);

    do {    
        dmapp_status = dmapp_acswap_qw(l_state.local_mutex, 
                &(l_state.mutexes[proc][mutex]),
                &(dst_reg->mr),
                proc, 0, l_state.rank + 1);
        assert(dmapp_status == DMAPP_RC_SUCCESS);
    }
    while(*(l_state.local_mutex) != 0);

    return COMEX_SUCCESS;
}


int comex_unlock(int mutex, int proc)
{
    int dmapp_status;
    reg_entry_t *dst_reg = NULL;

    /* preconditions */
    assert(0 <= proc && proc < l_state.size);
    assert(0 <= mutex && (unsigned int)(mutex) < l_state.num_mutexes[proc]);

    dst_reg = reg_cache_find(proc, &(l_state.mutexes[proc][mutex]), sizeof(unsigned long));
    assert(dst_reg);

    do {
        dmapp_status = dmapp_acswap_qw(l_state.local_mutex, 
                &(l_state.mutexes[proc][mutex]),
                &(dst_reg->mr),
                proc, l_state.rank + 1, 0);
        assert(dmapp_status == DMAPP_RC_SUCCESS);
    }
    while (*(l_state.local_mutex) != (unsigned long)(l_state.rank + 1));

    return COMEX_SUCCESS;
}


int comex_malloc(void *ptrs[], size_t size, comex_group_t group)
{
    comex_igroup_t *igroup = NULL;
    MPI_Comm comm = MPI_COMM_NULL;
    int comm_rank = -1;
    int comm_size = -1;
    int rc = MPI_SUCCESS; 
    void *src_buf = NULL;
    size_t max_size = size;
    dmapp_seg_desc_t heap_seg;
    dmapp_seg_desc_t *allgather_heap_seg = NULL;
    int i = 0;

    /* preconditions */
    assert(ptrs);
   
    igroup = comex_get_igroup_from_group(group);
    comm = igroup->comm;
    assert(comm != MPI_COMM_NULL);
    MPI_Comm_rank(comm, &comm_rank);
    MPI_Comm_size(comm, &comm_size);

    /* achieve consensus on the allocation size */
    rc = MPI_Allreduce(&size, &max_size, 1, MPI_LONG, MPI_MAX, comm);
    assert(rc == MPI_SUCCESS);
    size = max_size; 
    assert(size > 0);

    /* allocate and register segment */
    ptrs[comm_rank] = _comex_malloc_local(sizeof(char)*max_size, &heap_seg);
  
    /* exchange buffer address */
    /* @TODO: Consider using MPI_IN_PLACE? */
    memcpy(&src_buf, &ptrs[comm_rank], sizeof(void *));
    MPI_Allgather(&src_buf, sizeof(void *), MPI_BYTE, ptrs,
            sizeof(void *), MPI_BYTE, comm);

    /* allocate receive buffer for exchange of registration info */
    allgather_heap_seg = (dmapp_seg_desc_t *)my_malloc(
            sizeof(dmapp_seg_desc_t) * comm_size);
    assert(allgather_heap_seg);

    /* exchange registration info */
    MPI_Allgather(&heap_seg, sizeof(dmapp_seg_desc_t), MPI_BYTE,
            allgather_heap_seg, sizeof(dmapp_seg_desc_t), MPI_BYTE, comm); 

    /* insert this info into registration cache */
    for (i = 0; i < comm_size; ++i) {
        int world_rank;
        assert(COMEX_SUCCESS ==
                comex_group_translate_world(group, i, &world_rank));
        if (i == comm_rank)
            continue;
        reg_cache_insert(world_rank, ptrs[i], size, allgather_heap_seg[i]);
    }

    // Free the temporary buffer
    my_free(allgather_heap_seg);

    MPI_Barrier(comm);

    return COMEX_SUCCESS;
}


int comex_free(void *ptr, comex_group_t group)
{
    comex_igroup_t *igroup = NULL;
    MPI_Comm comm = MPI_COMM_NULL;
    int comm_rank;
    int comm_size;
    int i;
    long **allgather_ptrs = NULL;

    /* preconditions */
    assert(NULL != ptr);

    igroup = comex_get_igroup_from_group(group);
    comm = igroup->comm;
    assert(comm != MPI_COMM_NULL);
    MPI_Comm_rank(comm, &comm_rank);
    MPI_Comm_size(comm, &comm_size);

    /* allocate receive buffer for exchange of pointers */
    allgather_ptrs = (long **)my_malloc(sizeof(void *) * comm_size);
    assert(allgather_ptrs);

    /* exchange of pointers */
    MPI_Allgather(&ptr, sizeof(void *), MPI_BYTE,
            allgather_ptrs, sizeof(void *), MPI_BYTE, comm);

    /* remove all ptrs from registration cache */
    for (i = 0; i < comm_size; i++) {
        int world_rank;
        assert(COMEX_SUCCESS ==
                comex_group_translate_world(group, i, &world_rank));
        if (i == comm_rank)
            continue;
        reg_cache_delete(world_rank, allgather_ptrs[i]);
    }

    /* remove my ptr from reg cache and free ptr */
    comex_free_local(ptr);
    my_free(allgather_ptrs);

    /* Is this needed? */
    MPI_Barrier(comm);

    return COMEX_SUCCESS;
}


/* DMAPP Functions */


static void check_envs(void)
{
    char *value;

    /* COMEX_DMAPP_ROUTING
     *
     * TODO description */
    if ((value = getenv("COMEX_DMAPP_ROUTING")) != NULL){
        l_state.dmapp_routing = (atoi(value));
    }
    else {
        l_state.dmapp_routing = DMAPP_ROUTING_ADAPTIVE;
    }

#if HAVE_LIBHUGETLBFS

    /* hugepagesize
     *
     * set the static variable hugepagesize */
    hugepagesize = gethugepagesize();

    /* HUGETLB_MORECORE
     *
     * this variable controls whether malloc() will use hugepage memory which
     * means when we use malloc() and when the user application calls malloc
     * we will be competing for hugepage memory! */
    if ((value = getenv("HUGETLB_MORECORE")) != NULL) {
        if (0 == strncasecmp(value, "y", 1)) {
            malloc_is_using_huge_pages = 1;
        }
        else if (0 == strncasecmp(value, "n", 1)) {
            malloc_is_using_huge_pages = 0;
        }
    }

    /* COMEX_USE_HUGEPAGES
     *
     * COMEX can be built with hugepages and still allow the user to disable
     * their use. We assume that if libhugetlbfs is linked in, the user wants
     * to use it. This env var is then for the user to disable it, for some
     * reason. */
    comex_is_using_huge_pages = 1; /* the default if libhugetlbfs is linked */
    if ((value = getenv("COMEX_USE_HUGEPAGES")) != NULL) {
        if (0 == strncasecmp(value, "y", 1)) {
            comex_is_using_huge_pages = 1;
        }
        else if (0 == strncasecmp(value, "n", 1)) {
            comex_is_using_huge_pages = 0;
        }
    }

    /* HUGETLB_DEFAULT_PAGE_SIZE
     *
     * controls the page size that will be used for hugepages
     * we look for this value in case it is specified and we want to allocate
     * memory aligned to the same page size */
    if ((value = getenv("HUGETLB_DEFAULT_PAGE_SIZE")) != NULL){
        /* must be one of [128K|512K|2M|8M|16M|64M] */
        if (0 == strncasecmp(value, "128K", 4)) {
            hugetlb_default_page_size = 131072;
        }
        else if (0 == strncasecmp(value, "512K", 4)) {
            hugetlb_default_page_size = 524288;
        }
        else if (0 == strncasecmp(value, "2M", 2)) {
            hugetlb_default_page_size = 2097152;
        }
        else if (0 == strncasecmp(value, "8M", 2)) {
            hugetlb_default_page_size = 8388608;
        }
        else if (0 == strncasecmp(value, "16M", 3)) {
            hugetlb_default_page_size = 16777216;
        }
        else if (0 == strncasecmp(value, "64M", 3)) {
            hugetlb_default_page_size = 67108864;
        }
        else {
            assert(0);
        }
    }

    if (malloc_is_using_huge_pages || comex_is_using_huge_pages) {
        comex_page_size = hugepagesize;
    }

#endif /* HAVE_LIBHUGETLBFS */

    /* get page size for memory allocation */
    sc_page_size = sysconf(_SC_PAGESIZE);
    assert(sc_page_size >= 1);
    if (0 == comex_page_size) {
        comex_page_size = sc_page_size;
    }

//#if DEBUG
#if 1
    if (0 == l_state.rank) {
        printf("gethugepagesize()=%ld\n", hugepagesize);
        printf("hugetlb_default_page_size=%ld\n", hugetlb_default_page_size);
        printf("_SC_PAGESIZE=%ld\n", sc_page_size);
        printf("comex_page_size=%ld\n", comex_page_size);
        printf("comex_is_using_huge_pages=%d\n", comex_is_using_huge_pages);
        printf("malloc_is_using_huge_pages=%d\n", malloc_is_using_huge_pages);
    }
#endif
}


static void dmapp_initialize(void)
{
    dmapp_return_t status;
    dmapp_rma_attrs_ext_t requested_attrs;
    dmapp_rma_attrs_ext_t actual_attrs;
  
    memset(&requested_attrs, 0, sizeof(requested_attrs));
    memset(&actual_attrs, 0, sizeof(actual_attrs));

    // Check envs
    check_envs();

    /* The maximum number of outstanding non-blocking requests supported. You
     * can only specify this flag during initialization. The following is the
     * range of valid values to be supplied: [DMAPP_MIN_OUTSTANDING_NB, ..,
     * DMAPP_MAX_OUTSTANDING_NB] Setting the value to one of the extremes may
     * lead to a slowdown. The recommended value is DMAPP_DEF_OUTSTANDING_NB.
     * Users can experiment with the value to find the optimal setting for
     * their application. */
    requested_attrs.max_outstanding_nb = MAX_NB_OUTSTANDING;
    assert(MAX_NB_OUTSTANDING > DMAPP_MIN_OUTSTANDING_NB);
    assert(MAX_NB_OUTSTANDING < DMAPP_MAX_OUTSTANDING_NB);
    if (0 == l_state.rank) {
        if (MAX_NB_OUTSTANDING != DMAPP_DEF_OUTSTANDING_NB) {
            printf("MAX_NB_OUTSTANDING=%u != DMAPP_DEF_OUTSTANDING_NB=%u\n",
                    MAX_NB_OUTSTANDING, DMAPP_DEF_OUTSTANDING_NB);
        }
    }

    /* The threshold, in bytes, for switching between CPU-based
     * mechanisms and CPU offload mechanisms. This value can be
     * specified at any time and can use any value. The default setting is
     * DMAPP_OFFLOAD_THRESHOLD. Very small or very large settings
     * may lead to suboptimal performance. The default value is 4k bytes.
     * Consider how to best set this threshold. While a threshold increase
     * may increase CPU availability, it may also increase transfer latency
     * due to BTE involvement. */
    requested_attrs.offload_threshold = COMEX_DMAPP_OFFLOAD_THRESHOLD;

    /* Specifies the type of routing to be used. Applies to RMA requests with
     * PUT semantics and all AMOs. The default is DMAPP_ROUTING_ADAPTIVE.
     * The value can be specified at any time. Note that
     * DMAPP_ROUTING_IN_ORDER guarantees the requests arrive in order and may
     * result in poor performance.  Valid settings are:
     * - DMAPP_ROUTING_IN_ORDER
     * - DMAPP_ROUTING_DETERMINISTIC
     * - DMAPP_ROUTING_ADAPTIVE */
    requested_attrs.put_relaxed_ordering = l_state.dmapp_routing;

    /* Specifies the type of routing to be used. Applies to RMA requests with
     * GET semantics. The default is DMAPP_ROUTING_ADAPTIVE. The value can be
     * specified at any time. Note that DMAPP_ROUTING_IN_ORDER may result in
     * poor performance. Valid settings are:
     * - DMAPP_ROUTING_IN_ORDER
     * - DMAPP_ROUTING_DETERMINISTIC
     * - DMAPP_ROUTING_ADAPTIVE */
    requested_attrs.get_relaxed_ordering = l_state.dmapp_routing;

    /* The maximum number of threads that can access DMAPP. You can only use
     * this when thread-safety is enabled. The default is 1. You can only
     * specify this during initialization and it must be >= 1. */
    requested_attrs.max_concurrency = 1;

    /* Defines the PI ordering registration flags used by DMAPP when
     * registering all memory regions with GNI. Applies to the data, symmetric
     * heap, and user or dynamically mapped regions. The default is
     * DMAPP_PI_RELAXED_ORDERING.
     *
     * The dmapp_pi_reg_type_t enumeration defines the modes of PI access
     * ordering to be used by DMAPP during memory registration with uGNI;
     * therefore, these modes apply to the data and symmetric heap and any
     * user or dynamically mapped regions. 
     *
     * These modes do not affect GET operations.
     *
     * Strict ordering ensures that posted and non-posted writes arrive at the
     * target in strict order. Default and relaxed ordering impose no ordering
     * constraints, therefore if an application requires the global visibility
     * of data (for example, after a blocking put or gsync/fence), it must
     * perform extra synchronization in the form of a remote GET from the
     * target node in order to ensure that written data is globally visible.
     * - DMAPP_PI_ORDERING_STRICT   Strict PI (P_PASS_PW=0, NP_PASS_PW=0)
     * - DMAPP_PI_ORDERING_DEFAULT  Default GNI PI (P_PASS_PW=0, NP_PASS_PW=1)
     * - DMAPP_PI_ORDERING_RELAXED  Relaxed PI ordering (P_PASS_PW=1, NP_PASS_PW=1) */
    requested_attrs.PI_ordering = DMAPP_PI_ORDERING_RELAXED;

    // initialize    
    status = dmapp_init_ext(&requested_attrs, &actual_attrs);
    assert(status == DMAPP_RC_SUCCESS);
#define sanity(field) assert(actual_attrs.field == requested_attrs.field)
    sanity(max_outstanding_nb);
    sanity(offload_threshold);
    sanity(put_relaxed_ordering);
    sanity(get_relaxed_ordering);
    sanity(max_concurrency);
    sanity(PI_ordering);
#undef sanity

    // TODO is this the correct place to set this?
    max_outstanding_nb = actual_attrs.max_outstanding_nb;

    status = dmapp_get_jobinfo (&(l_state.job));
    assert(status == DMAPP_RC_SUCCESS);

    // Initialize the reg cache
    reg_cache_init(l_state.size);

    // Allocate buffers
    dmapp_alloc_buf();

    // Create locks
    create_dmapp_locks();

    /* Synchronize */
    MPI_Barrier(l_state.world_comm);
}


static void dmapp_terminate(void)
{
    int status;

    destroy_dmapp_locks();
   
    dmapp_free_buf();

    reg_cache_destroy(l_state.size); 
    
    status = dmapp_finalize();
    assert(status == DMAPP_RC_SUCCESS);

    MPI_Barrier(l_state.world_comm);
}


int comex_initialized()
{
    return initialized;
}
