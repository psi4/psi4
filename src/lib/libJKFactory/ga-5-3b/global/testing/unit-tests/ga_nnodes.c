/*
 * Test Program for GA
 * This is to test GA_Nnodes (is a local operation)
 * _Nnodes -- returns numbers of processes defined/assigned 
 */

#include<stdio.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"

node_check(int nprocs)
{
  if(GA_Nnodes() != nprocs)
    printf("ERROR: GA_Nnodes didnt return right number of processes \n");
}

int main(int argc, char **argv)
{
  int rank, nprocs;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);
  GA_Initialize();
  
  node_check(nprocs);

  GA_Sync();
  if(rank == 0)
    printf("Test Completed \n");

  GA_Terminate();
  MPI_Finalize();

  return 0;
}
