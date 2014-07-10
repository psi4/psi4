/*
 * Test Program for GA
 * This is to test GA_Copy_patch (is a collective operation)
 * GA_Create -- used to create a global array using handles like 'g_A'
 * GA_Duplicate --used to duplicate and generate one more global array.., handle 'g_A' to 'g_B'
 * GA_Copy -- is used to transfer the data from one array to other.., ie) from g_A to g_B
 * GA_Copy_patch -- is used  to copy only certain patch of one array to another
 */

#include<stdio.h>
#include"mpi.h"
#include"ga.h"
#include"macdecls.h"
#include"ga_unit.h"

#define DIM 2
#define SIZE 5

main(int argc, char **argv)
{
  int rank, nprocs, i, j;
  int g_A, g_B, local_A[SIZE][SIZE]; 
  int dims[DIM]={SIZE,SIZE}, alo[DIM]={0,0}, ahi[DIM]={3,3}, blo[DIM]={0,0}, bhi[DIM]={3,3}, ld=5;
  int value=5, val2=4;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);

  GA_Initialize();
  
  g_A = NGA_Create(C_INT, DIM, dims, "array_A", NULL);
  GA_Fill(g_A, &value);

  g_B = GA_Duplicate(g_A, "array_B");
  //  GA_Print(g_A);
  GA_Fill(g_B, &val2);
  GA_Sync();
  NGA_Copy_patch('N', g_A, alo, ahi, g_B, blo, bhi);
  // GA_Print(g_B);

  NGA_Get(g_B, blo, bhi, local_A, &ld);

  if(rank==0)
    {
      /*
      for(i=0; i<3; i++)
	{
	  for(j=0; j<3; j++)// printf("%d ", local_A[i][j]);
	     printf("\n");
	}
      //printf("\n");
      */

      for(i=0; i<3; i++)
	{
	  for(j=0; j<3; j++)
	    if( local_A[i][j]!=value) printf("ERROR : \n");
	}

    }


  // The process is confirmed and verified by printing the array in OP-scr
  /*
  if(rank==0)
    if(content(g_A) != content(g_B))printf("ERROE : \n");
  */

  if(rank == 0)
    printf("Test Completed \n");
  GA_Terminate();
  MPI_Finalize();

  
}
