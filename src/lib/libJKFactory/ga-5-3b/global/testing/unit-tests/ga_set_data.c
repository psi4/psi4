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

#define SIZE 5
#define DIM 2

int main(int argc, char **argv)
{
  int rank, nprocs;
  int g_A;
  int dims[DIM]={SIZE,SIZE}, ndim;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);
  GA_Initialize();

  g_A=GA_Create_handle();

  if(!g_A)
    GA_ERROR_MSG();

  GA_Set_data(g_A, DIM, dims, C_INT);
  //  printf("%d \n", GA_Ndim(g_A));    

  /*
  if(rank == 0)
    {
      ndim=GA_Ndim(g_A);
      if(ndim == 2)
	printf(" An Error \n");
    }
  */

  GA_Sync();
  if(rank == 0)
    GA_PRINT_MSG();
 
  GA_Terminate();
  MPI_Finalize();

  return 0;
}
