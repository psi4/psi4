/*
 * Test Program for GA 
 * This is to test GA_Nodeid (is a local operation)
 * verify _Nodeid -- by comparing rank of each processes with corresponding node-ids 
 */

#include<stdio.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"
#include"ga_unit.h"

#define DIM 2

int main(int argc, char **argv)
{
  int rank, nprocs;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);
  GA_Initialize();
  
  if(rank != GA_Nodeid())
    GA_ERROR_MSG();

  GA_Sync();

  if(rank == 0)
    GA_COMPLETE_MSG();

  GA_Terminate();
  MPI_Finalize();

  return 0;
}
