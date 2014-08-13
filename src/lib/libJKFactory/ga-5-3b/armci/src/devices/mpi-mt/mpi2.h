#ifndef MPI2_H
#define MPI2_H

#include <mpi.h>

#define MPI_SPAWN_DEBUG 0

#undef MPI_SPAWN_ZEROCOPY /* enables zero-copy for large requests */

#ifdef MPI_SPAWN_ZEROCOPY
#  define MPI_USER_DEF_DATATYPE /* Enables MPI userdefined type for non-contig
                                 * data, if MPI_SPAWN_ZEROCOPY is enabled */
#endif

/* uncomment MULTIPLE_BUFS macro definition to disable multiple buffers */
#undef MULTIPLE_BUFS

#define ARMCI_MPI_SPAWN_INIT_TAG   1000
#define ARMCI_MPI_SPAWN_TAG        2000
#define ARMCI_MPI_SPAWN_DATA_TAG   3000
#define ARMCI_MPI_SPAWN_VDATA_TAG  4000
#define ARMCI_MPI_CLIENT2SERVER_TAG 4500
#define ARMCI_MPI_SERVER2CLIENT_TAG 5000

/* In case of multiple buffers, we use tags from 2001 to 2999 (999 tags
 * total) to ensure flow control at the server side */
#define ARMCI_MPI_SPAWN_TAG_BEGIN  2001
#define ARMCI_MPI_SPAWN_TAG_END    2999

#define GET_SEND_BUFFER   _armci_buf_get
#define FREE_SEND_BUFFER  _armci_buf_release

#define COMPLETE_HANDLE   _armci_buf_complete_nb_request
#define TEST_HANDLE       _armci_buf_test_nb_request

#define SEND 0
#define RECV 1

extern void armci_mpi_strided(int op, void *ptr, int stride_levels,
                              int stride_arr[],  int count[], int proc,
                              MPI_Comm comm);

extern void armci_mpi_strided2(int op, void *ptr, int stride_levels,
                               int stride_arr[], int count[], int proc,
                               MPI_Comm comm);


#endif /* MPI2_H */
