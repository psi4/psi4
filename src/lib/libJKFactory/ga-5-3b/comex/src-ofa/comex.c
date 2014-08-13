#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* C and/or system headers */
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <unistd.h>

/* 3rd party headers */
#include <mpi.h>
#include <infiniband/verbs.h>

/* our headers */
#include "comex.h"
#include "device.h"
#include "comex_impl.h"
#include "groups.h"
#include "reg_cache.h"


/* exported state */
local_state l_state;

/* static state */
static int initialized=0;       /* for comex_initialized(), 0=false */
static long sc_page_size=-1;    /* from sysconf, in bytes */

/* static function declarations */
static void *_comex_malloc_local(size_t size, void **rinfo);


void comex_error(char *msg, int code)
{
    if (0 == l_state.rank) {
        fprintf(stderr,"Received an Error in Communication\n");
    }
    
    MPI_Abort(l_state.world_comm, code);
}


int comex_malloc(void **ptrs, size_t size, comex_group_t group)
{
    comex_igroup_t *igroup = NULL;
    MPI_Comm comm = MPI_COMM_NULL;
    int comm_size = -1;
    int comm_rank = -1;
    int rc = MPI_SUCCESS;
    void *src_buf = NULL;
    size_t max_size = size;
    void *local_rkey_buf = NULL;
    struct ibv_mr *mr = NULL;
    int *allgather_rkey_buf = NULL;
    int i = 0;

    /* preconditions */
    assert(ptrs);

    igroup = comex_get_igroup_from_group(group);
    comm = igroup->comm;
    assert(comm != MPI_COMM_NULL);
    rc = MPI_Comm_rank(comm, &comm_rank);
    assert(rc == MPI_SUCCESS);
    rc = MPI_Comm_size(comm, &comm_size);
    assert(rc == MPI_SUCCESS);
   
    /* achieve consensus on the allocation size */
    rc = MPI_Allreduce(&size, &max_size, 1, MPI_LONG, MPI_MAX, comm);
    assert(rc == MPI_SUCCESS);
    size = max_size; 

    /* 0 Byte registrations are potential problems */
    if (size < MIN_REGISTRATION_SIZE) {
        size = MIN_REGISTRATION_SIZE;
    }

    /* allocate and register the user level buffer */
    ptrs[comm_rank] = _comex_malloc_local(sizeof(char)*max_size, &local_rkey_buf);
    assert(local_rkey_buf);

    /* exchange buffer address */
    /* @TODO: Consider usijng MPI_IN_PLACE? */
    memcpy(&src_buf, &ptrs[comm_rank], sizeof(void *));
    MPI_Allgather(&src_buf, sizeof(void *), MPI_BYTE, ptrs,
            sizeof(void *), MPI_BYTE, comm);
   
    
    /* allocate buffer for collecting remote keys */
    allgather_rkey_buf = (int *)malloc(sizeof(int) * comm_size);
    assert(allgather_rkey_buf);

    /* set the rkey of the local buffer */
    mr = (struct ibv_mr *)local_rkey_buf;
    allgather_rkey_buf[comm_rank] = mr->rkey;
   
    /* exchange rkeys */
    MPI_Allgather(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL,
            allgather_rkey_buf, 1, MPI_INT, comm);

    /* insert in the remote registration cache */
    for (i = 0; i < comm_size; i++) {
        int world_rank;
        comex_group_translate_world(group, i, &world_rank);
        if (i == comm_rank) {
            continue;
        }
        reg_cache_insert(world_rank, ptrs[i], size, -1, allgather_rkey_buf[i], NULL);
    }

    /* free the temporary buffer */
    free(allgather_rkey_buf);

    MPI_Barrier(comm);

    return 0;
}


int comex_free(void *ptr, comex_group_t group)
{
    comex_igroup_t *igroup = NULL;
    MPI_Comm comm = MPI_COMM_NULL;
    int comm_size = -1;
    int comm_rank = -1;
    int i;
    long **allgather_ptrs = NULL;
    int rc = -1;

    /* preconditions */
    assert(NULL != ptr);

    igroup = comex_get_igroup_from_group(group);
    comm = igroup->comm;
    assert(comm != MPI_COMM_NULL);
    rc = MPI_Comm_rank(comm, &comm_rank);
    assert(rc == MPI_SUCCESS);
    rc = MPI_Comm_size(comm, &comm_size);
    assert(rc == MPI_SUCCESS);
   
    /* allocate receive buffer for exhange of pointers */
    allgather_ptrs = (long **)malloc(sizeof(void *) * comm_size); 
    assert(allgather_ptrs);

    /* exchange of pointers */
    rc = MPI_Allgather(&ptr, sizeof(void *), MPI_BYTE, 
            allgather_ptrs, sizeof(void *), MPI_BYTE, comm);
    assert(rc == MPI_SUCCESS);

    /* remove all ptrs from registration cache */
    for (i = 0; i < comm_size; i++) {
        int world_rank;
        comex_group_translate_world(group, i, &world_rank);
        if (i == comm_rank) {
            continue;
        }
        reg_cache_delete(world_rank, allgather_ptrs[i]);
    }
    
    /* remove my ptr from reg cache and free ptr */
    comex_free_local(ptr); 

    /* Synchronize: required by COMEX semantics */
    free(allgather_ptrs);
    MPI_Barrier(comm);        
    
    return 0;
}


static void *_comex_malloc_local(size_t size, void **rinfo)
{
    void *ptr = NULL;
    int rc = 0;

    /* allocate the user level buffer */
    rc = posix_memalign(&ptr, sc_page_size, sizeof(char)*size);
    assert(0 == rc);
    assert(ptr);

    /* register the buffer and check the return info */
    *rinfo = COMEXD_register_memory(ptr, size);

    return ptr;
}

// recursively find individual regions which can be registered
void recursive_reg(void *ptr, size_t size)
{
    assert(sc_page_size >= 0);
    if (size <= (unsigned long)sc_page_size)
        return;

    if (COMEXD_register_memory(ptr, size))
        return;
    else {
        recursive_reg(ptr, size / 2);
        recursive_reg(ptr + size /2, size / 2);
    } 
}

void *comex_malloc_local(size_t size)
{
    void *ptr = NULL;
    void *rinfo = NULL;
    
    ptr = _comex_malloc_local(size, &rinfo);

    if (!rinfo) {
        // registration failed
        recursive_reg(ptr, size);
    }
    return ptr;
}


int comex_free_local(void *ptr)
{
    assert(ptr != NULL);

    COMEXD_deregister_memory(ptr);

    free(ptr);

    return COMEX_SUCCESS;
}


int comex_init()
{
    int init_flag;
    
    /* Test if initialized has been called more than once */ 
    if (initialized) {
        return 0;
    }
    initialized = 1;

    /* Assert MPI has been initialized */
    MPI_Initialized(&init_flag);
    assert(init_flag);
    
    /* Duplicate the World Communicator */
    MPI_Comm_dup(MPI_COMM_WORLD, &(l_state.world_comm));
   
    /* My Rank */
    MPI_Comm_rank(l_state.world_comm, &(l_state.rank));

    /* World Size */
    MPI_Comm_size(l_state.world_comm, &(l_state.size));
  
    // init groups
    comex_group_init();

    // Initialize the number of outstanding messages
    l_state.num_outstanding = 0;

    // Get the page size for malloc
    sc_page_size = sysconf(_SC_PAGESIZE);

    // Initialize the COMEX device 
    COMEXD_initialize(); 

    /* Synch - Sanity Check */
    MPI_Barrier(l_state.world_comm);
   
    return 0;

}

int comex_init_args(int *argc, char ***argv)
{
    int rc;
    int init_flag;
    
    MPI_Initialized(&init_flag);
    
    if(!init_flag) {
        MPI_Init(argc, argv);
    }
    
    rc = comex_init();
   
    return rc;
}

int comex_finalize()
{
    // Make sure that all outstanding operations are done
    comex_wait_all(COMEX_GROUP_WORLD);

    // Call COMEX device specific finalize 
    COMEXD_finalize();

    MPI_Barrier(l_state.world_comm);

    // destroys all groups and their communicators
    comex_group_finalize();

    // destroy the primary communicator
    assert(MPI_SUCCESS == MPI_Comm_free(&l_state.world_comm));

    return COMEX_SUCCESS;
}

void comex_cleanup()
{
    comex_finalize();
}


int comex_rmw(int op, void *ploc, void *prem, int extra, int proc, comex_group_t group)
{
    if (op == COMEX_FETCH_AND_ADD) {
        COMEXD_network_lock(proc);
        comex_get(prem, ploc, sizeof(int), proc, group);
        (*(int *)ploc) += extra;
        comex_put(ploc, prem, sizeof(int), proc, group);
        (*(int *)ploc) -= extra;
        COMEXD_network_unlock(proc);
    }
    else if (op == COMEX_FETCH_AND_ADD_LONG) {
        COMEXD_network_lock(proc);
        comex_get(prem, ploc, sizeof(long), proc, group);
        (*(long *)ploc) += extra;
        comex_put(ploc, prem, sizeof(long), proc, group);
        (*(long *)ploc) -= extra;
        COMEXD_network_unlock(proc);

    }
    else if (op == COMEX_SWAP) {
        int tmp;
        COMEXD_network_lock(proc);
        comex_get(prem, &tmp, sizeof(int), proc, group);
        comex_put(ploc, prem, sizeof(int), proc, group);
        COMEXD_network_unlock(proc);
        *(int*)ploc = tmp;
    }
    else if (op == COMEX_SWAP_LONG) {
        long tmp;
        COMEXD_network_lock(proc);
        comex_get(prem, &tmp, sizeof(long), proc, group);
        comex_put(ploc, prem, sizeof(long), proc, group);
        COMEXD_network_unlock(proc);
        *(long*)ploc = tmp;
    }
    else  {
        assert(0);
    }

    
    return 0;
}


int comex_initialized()
{
    return initialized;
}
