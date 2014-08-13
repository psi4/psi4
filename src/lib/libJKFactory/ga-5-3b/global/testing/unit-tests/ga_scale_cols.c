/*
 * Test Program for GA
 * This is to test GA_Scale_cols (is a collective operation)
 * GA_Scale -- used to scale the global array value with constant- it to manipulate the value -- 
 *
 * GA_scale_cols -- 
 */

#include<stdio.h>
#include<stdlib.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"

#define DIM 2
#define SIZE 5

main(int argc, char **argv)
{
  int rank, nprocs;
  int g_A, g_V, val2=3, local_A[SIZE][SIZE], dims_V=SIZE, local_V[dims_V], local_Ar[SIZE][SIZE]; 
  int dims[DIM]={SIZE,SIZE}, dims2[DIM], lo[DIM]={0,0}, hi[DIM]={4,4}, ld=5, i, j;
  int loV=0, hiV=dims_V-1, arlo[DIM]={0,0}, arhi[DIM]={4,4};

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);

  GA_Initialize();

  for(i=0; i<SIZE; i++)
    {
      for(j=0; j<SIZE; j++)
	local_Ar[i][j]=rand()%5;
    }

  g_A = NGA_Create(C_INT, DIM, dims, "array_A", NULL);
  g_V = NGA_Create(C_INT, 1, &dims_V, "array_A", NULL);

  //  GA_Fill(g_A, &val1);
  NGA_Put(g_A, arlo, arhi, local_Ar, &ld);
  GA_Print(g_A);

  GA_Fill(g_V, &val2);
  GA_Print(g_V);

  GA_Scale_cols(g_A, g_V);
  GA_Print(g_A);
  
  NGA_Get(g_A, lo, hi, local_A, &ld);
  //  NGA_Get(g_A, &loV, &hiV, local_V, &ld);

  if(rank==1)
    {
      for(i=0; i<SIZE; i++)
	{
	  for(j=0; j<SIZE; j++) printf("%d ",local_A[i][j]);	 
	  printf("\n");
	}
      printf("\n");

      for(i=0; i<SIZE; i++)
	{
	  for(j=0; j<SIZE; j++) if(local_A[i][j]!=local_Ar[i][j]*val2)	 
	    printf("GA ERROR: \n");
	}
    }

  GA_Sync();
  if(rank == 0)
    printf("Test Completed \n");

  GA_Terminate();
  MPI_Finalize();
}
