/*
 * Test Program for GA
 * GA_Create_handle --
 */

#include<stdio.h>
#include<stdlib.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"
#include"ga_unit.h"

main(int argc, char **argv)
{
  int rank, nprocs;
  int g_A, g_B;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);
  GA_Initialize();

  g_A=GA_Create_handle();
  g_B=GA_Create_handle();

  if(!g_A)\
    GA_ERROR_MSG();
  if(!g_B)
    GA_ERROR_MSG();
  
  GA_Sync();
  if(rank == 0)
    GA_PRINT_MSG();
 
  GA_Terminate();
  MPI_Finalize();

}
