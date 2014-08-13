/*
 * Test Program for GA
 * This is to test GA_Zero_patch (is a collective operation)
 * GA_Create -- used to create a global array using handles like 'g_A'
 * GA_Duplicate --used to duplicate and generate one more global array.., handle 'g_A' to 'g_B'
 * GA_Zero_patch -- is used to pass zero-value only certain patch of an array 
 */

#include<stdio.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"

#define DIM 2

main(int argc, char **argv)
{
  int rank, nprocs, i, j;
  int g_A, g_B, local_A[5][5], local_B[5][5];  
  int dims[DIM]={5,5}, alo[DIM]={1,1}, ahi[DIM]={3,2}, blo[DIM]={1,1}, bhi[DIM]={2,3}, ld=5;
  int value=5, val2=4;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);

  GA_Initialize();
  
  g_A = NGA_Create(C_INT, DIM, dims, "array_A", NULL);
  GA_Fill(g_A, &value);

  g_B = GA_Duplicate(g_A, "array_B");
  GA_Print(g_A);
  GA_Fill(g_B, &val2);
  GA_Sync();
  NGA_Zero_patch(g_B, blo, bhi);
  GA_Print(g_B);

  NGA_Get(g_A, alo, ahi, local_A, &ld);
  NGA_Get(g_B, blo, bhi, local_B, &ld);

  if(rank==0)
    {
      for(i=0; i<3; i++)
	{
	  for(j=0; j<2; j++) printf("%d ", local_A[i][j]);
	  printf("\n");
	}
      printf("\n");

      for(i=0; i<2; i++)
	{
	  for(j=0; j<3; j++) printf("%d ", local_B[i][j]);
	  printf("\n");
	}

      for(i=0; i<2; i++)
	{
	  for(j=0; j<3; j++)
	    if( local_A[i][j]==local_B[i][j]) printf("ERROR : \n");
	}

    }
  

  // The process is confirmed and verified by printing the array in OP-scr
  /*
  if(rank==0)
    if(content(g_A) != content(g_B))printf("ERROE : \n");
  */

  GA_Sync();

  if(rank == 1)
    printf("Test Completed \n");
  GA_Terminate();
  MPI_Finalize();

  
}
