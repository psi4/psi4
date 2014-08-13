/*
 * Test Program for GA
 * This is to test GA_Elem_multiply_patch (is a collective operation)
 * GA_Create -- used to create a global array using handles like 'g_A' 
 * GA_Elem_multiply -- used to multiply two array and store it in one
 * Here _elem_multiply patch -- help to does the same function of_elem_multiply ...to an specified patch in array 
 */

#include<stdio.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"

#define DIM 2

main(int argc, char **argv)
{
  int rank, nprocs, i, j;
  int g_A, g_B, g_C, local_C[DIM][DIM]; 
  int dims[DIM]={5,5}, dims2[DIM], ndim, type;
  int val_A=5, val_B=3, ld=DIM; 
  int lo[DIM]={2,2}, hi[DIM]={4,4}, blo[DIM]={0,0}, bhi[DIM]={2,2}, clo[DIM]={1,1}, chi[DIM]={3,3};

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);

  GA_Initialize();
  
  g_A = NGA_Create(C_INT, DIM, dims, "array_A", NULL);
  g_B = NGA_Create(C_INT, DIM, dims, "array_B", NULL);
  g_C = NGA_Create(C_INT, DIM, dims, "array_C", NULL);

  GA_Fill(g_A, &val_A);
  GA_Fill(g_B, &val_B);
  GA_Zero(g_C);
  GA_Elem_multiply_patch(g_A, lo, hi, g_B, blo, bhi, g_C, clo, chi);
  GA_Print(g_C);
  GA_Sync();
  

  NGA_Get(g_C, clo, chi, local_C, &ld);
  if(rank==1)
    {
  
      for(i=0; i<DIM; i++)
	{
	  for(j=0; j<DIM; j++)printf("%d ", local_C[i][j]);
	  printf("\n");
	}
      
      for(i=0; i<DIM; i++)
	{
	  for(j=0; j<DIM; j++)
	    if(local_C[i][j]!=val_A*val_B) printf("GA Error : \n");
	}
      
    }
    
  GA_Sync();
  if(rank == 0)
    printf("Test Completed \n");

  GA_Terminate();
  MPI_Finalize();

}
