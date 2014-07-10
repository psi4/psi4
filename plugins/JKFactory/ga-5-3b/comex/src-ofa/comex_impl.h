#ifndef COMEX_IMPL_H_
#define COMEX_IMPL_H_

#include <mpi.h>
#include "device.h"

#define MIN_REGISTRATION_SIZE 64
#define ACC_PIPELINE_THRESHOLD  8192
typedef struct {

    MPI_Comm world_comm;
    int rank;
    int size;

    // buffer for lock
    void **atomic_lock_buf;
    void *local_lock_buf;
   
    void *acc_buf; 
    int acc_buf_len;
    void *put_buf; 
    int put_buf_len;
    void *get_buf; 
    int get_buf_len;

    // Number of outstanding data transfers
    unsigned long num_outstanding;

    int comex_openib_use_dreg;
} local_state;

extern local_state l_state;

#endif
