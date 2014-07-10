#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>
#include <string.h>

#include "comex.h"
#include "comex_impl.h"
#include "device.h"

int comex_fence_all(comex_group_t group)
{
    COMEXD_waitall();
    return COMEX_SUCCESS;
}


int comex_fence_proc(int proc, comex_group_t group)
{
    COMEXD_waitall();
    return COMEX_SUCCESS;
}


int comex_barrier(comex_group_t group)
{
    assert(l_state.world_comm);

    comex_fence_all(group);
    
    MPI_Barrier(l_state.world_comm);

    return COMEX_SUCCESS;
}

