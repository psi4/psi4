/*
 * Test Program for GA
 * This is to test GA_Scale_patch (is a collective operation)
 * GA_Scale -- used to scale the global array value with constant- it to manipulate the value -- 
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
  int g_A, g_V,  val1=5, val2=5, local_A[SIZE][SIZE], dims_V=SIZE; 
  int dims[DIM]={SIZE,SIZE}, dims2[DIM], lo[DIM]={1,1}, hi[DIM]={2,2}, ld=5, i, j;

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

  if(rank==1)
    {
      for(i=0; i<2; i++)
        {
          for(j=0; j<2; j++)
            if(local_A[i][j]!=val1*val2)printf(" GA Error: \n");
        }
    }


  if(rank == 0)
    printf("Test Completed \n");

  GA_Terminate();
  MPI_Finalize();
}
