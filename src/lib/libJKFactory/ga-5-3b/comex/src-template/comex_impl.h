#ifndef COMEX_IMPL_H_
#define COMEX_IMPL_H_

#include <mpi.h>

typedef struct {
    MPI_Comm world_comm;
    int rank;
    int size;
} local_state;

extern local_state l_state;

#endif /* COMEX_IMPL_H_ */
