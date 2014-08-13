#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <mpi.h>

void armci_msg_finalize()
{
    MPI_Finalize();
}
