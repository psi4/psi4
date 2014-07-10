/*
 * Test Program for GA
 * This is to test GA_Add_constant_patch (is a collective operation)
 * GA_Create -- used to create a global array using handles like 'g_A'
 * GA_Duplicate --used to duplicate and generate one more global array.., handle 'g_A' to 'g_B'
 * 
 *_add_constant_patch -- helps to 
 */

#include<stdio.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"

#define DIM 2
#define SIZ 5

main(int argc, char **argv)
{
  int rank, nprocs, i, j;
  int g_A, local_A[SIZ][SIZ], dims[DIM]={SIZ,SIZ}, val1=5, alpha=3;
  int alo[DIM]={2,2}, ahi[DIM]={3,3}, ld=SIZ;
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);

  GA_Initialize();
  
  g_A = NGA_Create(C_INT, DIM, dims, "array_A", NULL);

  GA_Fill(g_A, &val1);
  GA_Print(g_A);

  GA_Add_constant_patch(g_A, alo, ahi, &alpha);
  GA_Print(g_A);
 
  NGA_Get(g_A, alo, ahi, local_A, &ld);

  if(rank==1)
    {
      for(i=0; i<DIM; i++)
	{
	  for(j=0; j<DIM; j++)
	    printf(" %d", local_A[i][j]);
	  printf("\n");
	}

      for(i=0; i<DIM; i++)
	for(j=0; j<DIM; j++) if(local_A[i][j]!=val1+alpha)
	  printf(" GA Error: \n");
    }

  GA_Sync();
  if(rank==0)
    printf(" GA TEST: Completed \n");
  GA_Terminate();
  MPI_Finalize();
 
}
