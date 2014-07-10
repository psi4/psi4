#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <mpi.h>

#if HAVE_STDIO_H
#   include <stdio.h>
#endif

int main(int argc, char **argv)
{
    int numprocs, myid;
    int ierr;
    char req=0, ack=0;
    int atag=999, rtag=555, to;
    MPI_Status status;
    MPI_Request request;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

    if(myid==0){
        to = numprocs-1;
        printf("Testing nonblocking receive\n\n"); fflush(stdout);
        ierr = MPI_Irecv(&ack, 1, MPI_CHAR, to, atag ,MPI_COMM_WORLD, &request);
        printf(":after nonblocking receive\n"); fflush(stdout);
        ierr = MPI_Send(&req, 1, MPI_CHAR, to, rtag, MPI_COMM_WORLD);
        printf(":sent request\n"); fflush(stdout);
        ierr = MPI_Wait(&request, &status);
        printf(":received response\n"); fflush(stdout);
        printf("\nnonblocking receive is working\n"); fflush(stdout);
    }
    if(myid==numprocs-1){
        to = 0;
        ierr = MPI_Recv(&req, 1, MPI_CHAR, to, rtag, MPI_COMM_WORLD, &status);
        printf("::request received\n"); fflush(stdout);
        ierr = MPI_Send(&ack, 1, MPI_CHAR, to, atag, MPI_COMM_WORLD);
    }

    MPI_Finalize();

    return 0;
}
