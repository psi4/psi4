/*
 * Test Program for GA
 * This is to test GA_Inquire (is a local operation)
 * GA_Create -- used to create a global array using handles like 'g_A' 
 * GA_Inquire --  to verify that g_A handle returns the right spec. of _array created under it.
 * This is one method used to verify the created array. Here to test the specified function _inquire
 */

#include<stdio.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"

#define DIM 2

int main(int argc, char **argv)
{
  int rank, nprocs, i;
  double g_A;
  int dims[DIM]={5,5}, dims2[DIM], ndim, type;
  double value=5;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_DBL, 1000, 1000);

  GA_Initialize();
  
  g_A = NGA_Create(C_DBL, DIM, dims, "array_A", NULL);

  GA_Fill(g_A, &value);
  GA_Print(g_A);

  if(!g_A) printf("ERROR : \n");

  NGA_Inquire(g_A, &type, &ndim, dims2);
  printf(" %d -- %d,,\n", type, ndim);

  for(i=0; i<DIM; i++)
    printf("%d: %d[ %d] ...* \n", rank, i, dims2[i]);
  
  GA_Sync();
  if(rank == 0)
    printf("Test Completed \n");

  GA_Terminate();
  MPI_Finalize();
 
  return 0;
}
