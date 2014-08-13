#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <mpi.h>

void armci_msg_init(int *argc, char ***argv)
{
    int flag=0;
    MPI_Initialized(&flag);
    if (!flag) {
#if 0
        int provided;
        MPI_Init_thread(argc, argv, MPI_THREAD_MULTIPLE, &provided);
#else
        MPI_Init(argc, argv);
#   endif
    }
}
