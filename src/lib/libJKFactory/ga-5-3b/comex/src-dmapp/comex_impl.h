#ifndef COMEX_IMPL_H_
#define COMEX_IMPL_H_

#include <dmapp.h>
#include <mpi.h>

#define COMEX_DMAPP_OFFLOAD_THRESHOLD 2048
#define DEFAULT_SYM_HEAP_SIZE 32*1048576
#define DMAPP_ROUTING 1
#define MAX_NB_OUTSTANDING 1024
#define FAILURE_BUFSIZE 1048576
#define PIPELINED_ACCUMULATE 0
#define PIPELINED_MAX_BUFFERS 5
#define PIPELINED_BUFSIZE 1048576

typedef struct {

    MPI_Comm world_comm;
    int rank;
    int size;

    /* DMAPP Specific */
    dmapp_jobinfo_t job;

    /* buffers for locks */
    unsigned long **atomic_lock_buf; /**< internal lock, one per process */
    unsigned long *local_lock_buf; /**< holds value of remote lock locally */
    unsigned long **mutexes; /**< all mutexes */
    unsigned long *local_mutex; /**< store the remote mutex value */
    unsigned int  *num_mutexes; /**< how many mutexes on each process */

    /* fallback buffers when errors are detected */
    void *acc_buf;
    int acc_buf_len;
    void *put_buf;
    int put_buf_len;
    void *get_buf;
    int get_buf_len;

#if PIPELINED_ACCUMULATE
    /* buffers for pipelining */
    void **pipe_acc_buf;
    int pipe_acc_buf_len;
#endif

    /* envs */ 
    int dmapp_routing;
} local_state;

extern local_state l_state;

#endif /* COMEX_IMPL_H_ */
