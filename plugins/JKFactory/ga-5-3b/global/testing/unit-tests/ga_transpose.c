/*
 * Test Program for GA
 * This is to test GA_Transpose (is a collective operation)
 * GA_Create -- used to create a global array using handles like 'g_A'
 * _Transpose -- used to transpose the whole array  
 * 
 * not completed --- 
 */

#include<stdio.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"

#define DIM 2
#define SIZE 5

main(int argc, char **argv)
{
  int rank, nprocs, i, j;
  int g_A, g_B, value=5, local_GA[SIZE][SIZE], local_A[SIZE][SIZE], local_B[SIZE][SIZE];  
  int dims[DIM]={SIZE,SIZE}, dims2[DIM], lo[DIM]={0,0}, hi[DIM]={4,4}, ld=5;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);

  GA_Initialize();

  g_A = NGA_Create(C_INT, DIM, dims, "array_A", NULL);
  g_B = NGA_Create(C_INT, DIM, dims, "array_B", NULL);

  for(i=0; i<SIZE; i++)
    for(j=0; j<SIZE; j++) local_GA[i][j]=rand()%10;
	
  NGA_Put(g_A, lo, hi, local_GA, &ld);
  GA_Transpose(g_A, g_B);
  GA_Print(g_A);
  GA_Print(g_B);

  if(rank==0)
    {
      NGA_Get(g_A, lo, hi, local_A, &ld);
      NGA_Get(g_B, lo, hi, local_B, &ld);

      for(i=1; i<3; i++)
	{
	  for(j=1; j<3; j++)
	    printf("%d ", local_B[i][j]);
	  printf("\n");
	}

      printf("\n");
      for(i=1; i<3; i++)
	{
	  for(j=1; j<3; j++)
	    printf("%d ", local_A[i][j]);
	  printf("\n");
	}
      
      printf("\n");
      for(i=1; i<3; i++)
	{
	  for(j=1; j<3; j++)
	    {
	      if(local_B[j][i]!=local_A[i][j])
		printf("ERROR : in passing values \n");
	    }
	}
    }
  GA_Sync();
  if(rank == 1)
    printf("Test Completed \n");

  GA_Terminate();
  MPI_Finalize();
}
