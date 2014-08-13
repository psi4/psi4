/*
 * Test Program for GA
 * This is to test GA_Fill(is a collective operation)
 * GA_Fill -- used to fill any constant into global array --- simple to pass values .., 
 */

#include<stdio.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"

#define DIM 2
#define SIZE 5

fillandscale(int rank, int nprocs)
{
  int g_A,  val1=5, val2=5, local_A[SIZE][SIZE], i, j; 
  int dims[DIM]={SIZE,SIZE}, alo[DIM]={1,1}, ahi[DIM]={2,2}, ld=5;

  g_A = NGA_Create(C_INT, DIM, dims, "array_A", NULL);
  GA_Zero(g_A);
  NGA_Fill_patch(g_A, alo, ahi, &val1);
  GA_Print(g_A);

  GA_Scale(g_A, &val2);
  GA_Print(g_A);

  NGA_Get(g_A, alo, ahi, local_A, &ld);

  if(rank == 1)
    {
      for(i=0; i<DIM; i++)
	{
	  for(j=0; j<DIM; j++) if(local_A[i][j]!=val1*val2)
	    printf(" GA ERROR: \n");
	}
    }
  GA_Destroy(g_A);
}

fillonly(int rank, int nprocs)
{
 int g_A,  val1=5, val2=5, local_A[SIZE][SIZE], i, j; 
  int dims[DIM]={SIZE,SIZE}, alo[DIM]={1,1}, ahi[DIM]={2,2}, ld=5;

  g_A = NGA_Create(C_INT, DIM, dims, "array_A", NULL);
  GA_Zero(g_A);
  NGA_Fill_patch(g_A, alo, ahi, &val1);
  GA_Print(g_A);

  NGA_Get(g_A, alo, ahi, local_A, &ld);

  if(rank == 1)
    {
      for(i=0; i<DIM; i++)
	{
	  for(j=0; j<DIM; j++) if(local_A[i][j]!=val1)
	    printf(" GA ERROR: \n");
	}
    }
  GA_Destroy(g_A);
}

main(int argc, char **argv)
{
  int rank, nprocs;
  
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);

  GA_Initialize();
  fillonly(rank, nprocs);
  fillandscale(rank, nprocs);

  GA_Sync();
  if(rank == 0)
    printf("Test Completed \n");

  GA_Terminate();
  MPI_Finalize();
}
