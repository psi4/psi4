/*
 * Test Program for GA
 * This is to test GA_Ndim (is a local operation)
 * verifing GA_Ndim -- which returns dimension of the g_A 
 */

#include<stdio.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"

#define SIZE 5
#define ONE_DIM 1
#define TWO_DIM 2
#define THREE_DIM 3

one_dimension()
{
  int g_A, dims[ONE_DIM]={SIZE};

  g_A = NGA_Create(C_INT, ONE_DIM, dims, "array_A", NULL);

  if(GA_Ndim(g_A) != ONE_DIM)
    printf("ERROR: GA_Ndim didnt return nDimension after GA_Initialize \n");

  GA_Destroy(g_A);
}

two_dimension()
{
  int g_A, dims[TWO_DIM]={SIZE,SIZE};

  g_A = NGA_Create(C_INT, TWO_DIM, dims, "array_A", NULL);

  if(GA_Ndim(g_A) != TWO_DIM)
    printf("ERROR: GA_Ndim didnt return nDimension after GA_Initialize \n");

  GA_Destroy(g_A);
}

three_dimension()
{
  int g_A, dims[THREE_DIM]={SIZE,SIZE,SIZE};

  g_A = NGA_Create(C_INT, THREE_DIM, dims, "array_A", NULL);

  if(GA_Ndim(g_A) != THREE_DIM)
    printf("ERROR: GA_Ndim didnt return nDimension after GA_Initialize \n");

  GA_Destroy(g_A);
}

int main(int argc, char **argv)
{
  int rank, nprocs, n=1, temp=0;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);

  GA_Initialize();
  
  one_dimension();
  two_dimension();
  three_dimension();

  GA_Sync();
  if(rank == 0)
    printf("Test Completed \n");

  GA_Terminate();
  MPI_Finalize();

  return 0;
}
