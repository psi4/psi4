/*
 * Test Program for GA
 * This is to test GA_Transpose (is a collective operation)
 * GA_Create -- used to create a global array using handles like 'g_A'
 * _Transpose -- used to transpose the whole array  
 * 
 * not completed --- 
 */

#include<stdio.h>
#include<stdlib.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"
#include"ga_unit.h"

#define DIM 2
#define SIZE 5

validate_transpose(int g_A, int g_B, int* lo, int* hi, int ld)
{
  int i, j;
  int local_A[SIZE][SIZE], local_B[SIZE][SIZE];
  //  int  **local_A=NULL, **local_B=NULL;    

  /*
  local_A=(int**)malloc(SIZE*sizeof(int*));
  for(i=0; i<SIZE; i++)
    local_value[i]=(int*)malloc(SIZE*sizeof(int));

  local_B=(int**)malloc(SIZE*sizeof(int*));
  for(i=0; i<SIZE; i++)
    local_value[i]=(int*)malloc(SIZE*sizeof(int));
  */

  NGA_Get(g_A, lo, hi, local_A, &ld);
  NGA_Get(g_B, lo, hi, local_B, &ld);
  
  for(i=1; i<SIZE; i++)
    {
      for(j=1; j<SIZE; j++)
	{
	      if(local_B[j][i]!=local_A[i][j])
		GA_Error("ERROR : in passing values", DIM);
	}
    }
  
}
main(int argc, char **argv)
{
  int rank, nprocs, i, j;
  int g_A, g_B, **local_value=NULL;

  int dims[DIM]={SIZE,SIZE}, lo[DIM]={SIZE-SIZE,SIZE-SIZE}, hi[DIM]={SIZE-1,SIZE-1}, ld=SIZE;

  local_value=(int**)malloc(SIZE*sizeof(int*));
  for(i=0; i<SIZE; i++)
    local_value[i]=(int*)malloc(SIZE*sizeof(int));

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);

  GA_Initialize();

  g_A = NGA_Create(C_INT, DIM, dims, "array_A", NULL);
  g_B = NGA_Create(C_INT, DIM, dims, "array_B", NULL);

  for(i=0; i<SIZE; i++)
    for(j=0; j<SIZE; j++)
      local_value[i][j]=rand()%10;
	
  if(rank==0) NGA_Put(g_A, lo, hi, local_value, &ld);
  GA_Transpose(g_A, g_B);
  if(rank==0) validate_transpose(g_A, g_B, lo, hi, ld);

  GA_Sync();
  if(rank == 1) GA_PRINT_MSG();
  
  GA_Terminate();
  MPI_Finalize();
}
