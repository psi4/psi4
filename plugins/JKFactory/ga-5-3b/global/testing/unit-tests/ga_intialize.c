/*
 * Test Program for GA
 * This is to test GA_initialize
 * _intitalize and _terminate ---- are function start and end the space space where other GA-Functions are called     xcfor manipulation
 * Using GA_Ndim -- ehich returns dimension of the g_A 
 */

#include<stdio.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"

#define DIM 2

int main(int argc, char **argv)
{
  int rank, nprocs, n=1;
  int g_A, dims[DIM]={5,5};

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);

  GA_Initialize();
  
  printf("check %d \n", n);  

  g_A = NGA_Create(C_INT, DIM, dims, "array_A", NULL);
  if(GA_Ndim(g_A)!=DIM)
    printf("ERROR: GA_Ndim didnt return nDimension after GA_Initialize\n");
  printf("%d : %d \n", rank, GA_Ndim(g_A));

  GA_Terminate();

  if(rank==0)
    printf(" GA: Test Completed \n");
  MPI_Finalize();

  return 0;
}
