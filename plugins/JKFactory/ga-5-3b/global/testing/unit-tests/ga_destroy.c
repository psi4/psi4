/*
 * Test Program for GA
 * This is to test GA_Destroy (is a collective operation)
 * GA_Create -- used to create a global array using handles like 'g_A' 
 * GA_destroy -- used to trash the created/used arrays
 */

#include<stdio.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"

#define DIM 2

main(int argc, char **argv)
{
  int rank, nprocs, i;
  int g_A,  value=5; 
  int dims[DIM]={5,5}, dims2[DIM], ndim, type;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);

  GA_Initialize();
  
  g_A = NGA_Create(C_INT, DIM, dims, "array_A", NULL);
  GA_Fill(g_A, &value);
  GA_Print(g_A);

  //GA_Sync();
  printf(" -----------%d\n", rank); 
  GA_Destroy(g_A);
  printf(" %d-----------\n", rank); 
  GA_Print(g_A);

  if(rank == 0)
    printf("Test Completed \n");
  GA_Terminate();
  MPI_Finalize();

  
}
