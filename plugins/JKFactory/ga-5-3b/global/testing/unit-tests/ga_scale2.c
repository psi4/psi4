/*
 * Test Program for GA
 * This is to test GA_Scale (is a collective operation)
 * GA_Scale -- used to scale the global array value with constant- it to manipulate the value -- 
 */

#include<stdio.h>
#include<stdlib.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"
#include"ga_unit.h"

#define ONE_DIM 1
#define TWO_DIM 2
#define SIZE 5

one_dimension_array(int rank, int val_V, int val_scale)
{
  int g_V, *local_V=NULL, dims_V=SIZE;
  int lo=SIZE-SIZE, hi=SIZE-1, i, ld=SIZE;

  local_V=(int*)malloc(SIZE*sizeof(int));
  g_V = NGA_Create(C_INT, ONE_DIM, &dims_V, "array_A", NULL);
  GA_Fill(g_V, &val_V);
  GA_Scale(g_V, &val_scale);
  GA_Print(g_V);

  NGA_Get(g_V, &lo, &hi, local_V, &ld);

  if(rank==1)
    {
      for(i=0; i<SIZE; i++)
	if(local_V[i]==val_V*val_scale)printf(" GA Error: 1\n");
    }
  GA_Destroy(g_V);
}

two_dimension_array(int rank, int val_A, int val_scale)
{
  int g_A, **local_A=NULL,  dims[TWO_DIM]={SIZE,SIZE}; 
  //  int g_A, local_A[SIZE][SIZE],  dims[TWO_DIM]={SIZE,SIZE}; 
  //  int lo_A[TWO_DIM]={SIZE-SIZE,SIZE-SIZE}, hi_A[TWO_DIM]={SIZE-1,SIZE-1}, i, j, ld=SIZE;
  int lo_A[TWO_DIM]={1,1}, hi_A[TWO_DIM]={4,4}, i, j, ld=SIZE;


  local_A=(int**)malloc(SIZE*sizeof(int*));
  for(i=0; i<SIZE*sizeof(int); i++)
    local_A[i]=(int*)malloc(SIZE*sizeof(int));

  g_A = NGA_Create(C_INT, TWO_DIM, dims, "array_A", NULL);

  GA_Fill(g_A, &val_A);
  //  GA_Print(g_A);
  GA_Scale(g_A, &val_scale);
  //  GA_Print(g_A);

  NGA_Get(g_A, lo_A, hi_A, local_A, &ld);

  if(rank==1)
    {
      for(i=0; i<SIZE; i++)
        {
	  printf(" check 2\n");
          for(j=0; j<SIZE; j++)
	    printf("%d ",local_A[i][j]);
	  printf(" \n");
        }

      for(i=0; i<SIZE; i++)
        {
	  printf(" check 2\n");
          for(j=0; j<SIZE; j++)
            if(local_A[i][j]==val_A*val_scale)printf(" GA Error: 2 \n");
        }
    }
  GA_Destroy(g_A);
}
main(int argc, char **argv)
{
  int rank, nprocs;
  int val_scale=5, val_A=5, val_V=3;

  MPI_Init(&argc, &argv); 

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);

  GA_Initialize();

  one_dimension_array(rank, val_V,val_scale);
  two_dimension_array(rank, val_A,val_scale);

  GA_Sync();
  if(rank == 0)
    printf("Test Completed \n");

  GA_Terminate();
  MPI_Finalize();
}
