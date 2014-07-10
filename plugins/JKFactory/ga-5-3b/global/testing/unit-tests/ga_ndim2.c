/*
 * Test Program for GA
 * This is to test GA_Ndim (is a local operation)
 * verifing GA_Ndim -- which returns dimension of the g_A 
 */

#include<stdio.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"

#define SIZE 20
#define NDIM 7

verify_ga_dim(int ndim)
{
  int g_A, dims[ndim], i;

  for(i=0; i<ndim; i++) dims[i]=SIZE;

  g_A = NGA_Create(C_INT, ndim, dims, "array_A", NULL);

  if(GA_Ndim(g_A) != ndim)
    printf("ERROR: GA_Ndim -- %d returned wrong \n", ndim);

  GA_Destroy(g_A);
 }

int main(int argc, char **argv)
{
  int rank, nprocs, i;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);

  GA_Initialize();

  for(i=1; i<=NDIM; i++)
    {
      verify_ga_dim(i);
    }

  GA_Sync();
  if(rank == 0)
    printf("Test Completed \n");

  GA_Terminate();
  MPI_Finalize();

  return 0;
}
