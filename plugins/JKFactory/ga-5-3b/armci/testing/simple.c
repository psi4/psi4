#if HAVE_CONFIG_H
#   include "config.h"
#endif

/*$id$*/
#include <stdio.h>
#include <stdlib.h>

#include <mpi.h>

#include "armci.h"
#include "message.h"

int me, nprocs;
int LOOP = 10;

int main(int argc, char **argv)
{
  int k, i;
  double **myptrs[10];
  double t0, t1, tget = 0, tnbget = 0, tput = 0, tnbput = 0, tnbwait = 0, t2 = 0;
  armci_msg_init(&argc, &argv);
  ARMCI_Init_args(&argc, &argv);
  ARMCI_Init_args(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &me);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  for (k = 0; k < 10; k++) {
    myptrs[k] = (double **)malloc(sizeof(double *) * nprocs);
    ARMCI_Malloc((void **)myptrs[k], 400000 * LOOP * sizeof(double));
    for (i = 0; i < LOOP; i++) {
      myptrs[k][me][i] = me + 0.414;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    for (i = 0; i < LOOP; i++) {
      ARMCI_Get(myptrs[k][(me+1)%nprocs] + i, myptrs[k][me] + i, sizeof(double), (me + 1) % nprocs);
      /*if(myptrs[k][me][i]!=0.414+(me+1)%nprocs)ARMCI_Error("errr",myptrs[k][me][i]);*/
    }
    t0 = t1 = tget = tnbget = tput = tnbput = tnbwait = t2 = 0;
    t0 = MPI_Wtime();
    for (i = 0; i < LOOP; i++) {
      ARMCI_Get(myptrs[k][(me+1)%nprocs] + i, myptrs[k][me] + i, sizeof(double), (me + 1) % nprocs);
    }
    t1 = MPI_Wtime();
    printf("\nGet Latency=%f\n", 1e6 *(t1 - t0) / LOOP);
    fflush(stdout);
    t1 = t0 = 0;
    for (i = 0; i < LOOP; i++) {
      armci_hdl_t nbh;
      ARMCI_INIT_HANDLE(&nbh);
      t0 = MPI_Wtime();
      ARMCI_NbGet(myptrs[k][(me+1)%nprocs] + i, myptrs[k][me] + i, sizeof(double), (me + 1) % nprocs, &nbh);
      t1 = MPI_Wtime();
      ARMCI_Wait(&nbh);
      t2 = MPI_Wtime();
      tnbget += (t1 - t0);
      tnbwait += (t2 - t1);
    }
    printf("\nNb Get Latency=%f Nb Wait=%f\n", 1e6 * tnbget / LOOP, 1e6 * tnbwait / LOOP);
    fflush(stdout);
    MPI_Barrier(MPI_COMM_WORLD);
  }
  for (k = 0; k < 10; k++) {
    ARMCI_Free(myptrs[k][me]);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  ARMCI_Finalize();
  ARMCI_Finalize();
  armci_msg_finalize();
  return 0;
}
