/*
 * Test Program for GA
 * This is to test GA_Scale_rows (is a collective operation)
 * GA_Scale -- used to scale the global array value with constant- it to manipulate the value -- 
 *
 * GA_scale_rows -- 
 */

#include<stdio.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"

#define DIM 2
#define SIZE 5

main(int argc, char **argv)
{
  int rank, nprocs;
  int g_A, g_V,  val1=5, val2=5, local_A[SIZE][SIZE], dims_V=SIZE, local_V[dims_V]; 
  int dims[DIM]={SIZE,SIZE}, dims2[DIM], lo[DIM]={1,1}, hi[DIM]={2,2}, ld=5, i, j;
  int loV=0, hiV=dims_V-1;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);

  GA_Initialize();

  g_A = NGA_Create(C_INT, DIM, dims, "array_A", NULL);
  g_V = NGA_Create(C_INT, 1, &dims_V, "array_A", NULL);
  GA_Fill(g_A, &val1);
  GA_Print(g_A);
  printf("\n");

  GA_Scale(g_A, &val2);
  GA_Print(g_A);

  GA_Get_diag(g_A, g_V);
  GA_Print(g_V);
  
  NGA_Get(g_A, lo, hi, local_A, &ld);
  NGA_Get(g_V, &loV, &hiV, local_V, &ld);

  if(rank==1)
    {
      for(i=0; i<dims_V; i++)
	if(local_V[i]!=val1*val2) printf(" GA Error: \n");
    }

  if(rank == 0)
    printf("Test Completed \n");

  GA_Terminate();
  MPI_Finalize();
}
