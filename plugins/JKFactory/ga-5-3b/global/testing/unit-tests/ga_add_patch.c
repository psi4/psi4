/*
 * Test Program for GA
 * This is to test GA_Add_patch (is a collective operation)
 * GA_Create -- used to create a global array using handles like 'g_A'
 * GA_Duplicate --used to duplicate and generate one more global array.., handle 'g_A' to 'g_B'
 * 
 *
 */

#include<stdio.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"
#include"ga_unit.h"

#define DIM 2

main(int argc, char **argv)
{
  int rank, nprocs, i, j;
  int g_A, g_B, g_C, local_C[DIM][DIM], dims[DIM]={5,5}, val1=5, val2=4, alpha=3, beta=2, ld=5;
  int alo[DIM]={2,2}, ahi[DIM]={3,3}, blo[DIM]={2,2}, bhi[DIM]={3,3}, clo[DIM]={1,1}, chi[DIM]={2,2};

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);

  GA_Initialize();
  
  g_A = NGA_Create(C_INT, DIM, dims, "array_A", NULL);

  g_B = GA_Duplicate(g_A, "array_B");
  g_C = GA_Duplicate(g_A, "array_C");

  GA_Fill(g_A, &val1);
  GA_Fill(g_B, &val2);
  GA_Zero(g_C);

  NGA_Add_patch(&alpha, g_A, clo, chi, &beta, g_B, blo, bhi, g_C, clo, chi);

  GA_Sync();
  GA_Print(g_A);
  GA_Print(g_B);
  GA_Print(g_C);

  NGA_Get(g_C, clo, chi, local_C, &ld);

  //printf("check 1 \n");

  for(i=0; i<DIM; i++)
    {
      for(j=0; j<DIM; j++)printf("%d ", local_C[i][j]);
      printf("\n");
    }
  
  if(rank == 0)
    {
      printf("check 2\n");
    
      for(i=0; i<DIM; i++)
	{
	  for(j=0; j<DIM; j++)
	    if(local_C[i][j]!=(alpha*val1)+(beta*val2)) printf("GA Error : \n");
	}
    }
  
  if(rank==0)
    GA_PRINT_MSG();

  GA_Sync();

  /*
  GA_Destroy(g_A);
  GA_Destroy(g_B);
  GA_Destroy(g_C);
  */

  //*******************************************************************

  /* what would be the possible reason for GA_destroy to get failed .., 
   * solve this before consolidate the whole
   */

  GA_Terminate();
  MPI_Finalize();
 
}
