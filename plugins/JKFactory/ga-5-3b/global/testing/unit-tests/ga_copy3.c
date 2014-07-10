/*
 * Test Program for GA
 * This is to test GA_Duplicate (is a collective operation)
 * GA_Create -- used to create a global array using handles like 'g_A'
 * GA_Duplicate --used to duplicate and generate one more global array.., handle 'g_A' to 'g_B'
 * But GA_Copy -- is used to transfer the data from one array to other.., ie) from g_A to g_B
 * Here used GA_Inquire to verify that g_A hanle returns the right values of created_array
 */

#include<stdio.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"
#include"ga_unit.h"

#define DIM 2
#define SIZE 5

main(int argc, char **argv)
{
  int rank, nprocs, i, j;
  int g_A, g_B,  local_A[SIZE][SIZE], local_B[SIZE][SIZE]; 
  int dims[DIM]={SIZE,SIZE}, dims2[DIM], ndim2, type2, dims3[DIM], ndim3, type3, lo[DIM]={0,0}, hi[DIM]={4,4};
  int value=5, val2=4, ld=5;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);

  GA_Initialize();
  
  g_A = NGA_Create(C_INT, DIM, dims, "array_A", NULL);
  GA_Fill(g_A, &value);

  g_B = GA_Duplicate(g_A, "array_B");
  GA_Print(g_A);
  GA_Sync();
  GA_Copy(g_A, g_B);
  GA_Print(g_B);

  // The process is confirmed and verified by printing the array in OP-scr
  /*
  if(rank==0)
    if(content(g_A) != content(g_B))printf("ERROE : \n");
  */

  if(rank==0)
    {
      NGA_Get(g_A, lo, hi, local_A, &ld);
      NGA_Get(g_B, lo, hi, local_B, &ld);

      for(i=1; i<3; i++)
        {
          for(j=1; j<3; j++)
            printf("%d ", local_B[i][j]);
          printf("\n");
        }

      printf("\n");
      for(i=1; i<3; i++)
        {
          for(j=1; j<3; j++)
            printf("%d ", local_A[i][j]);
          printf("\n");
        }

      printf("\n");
      for(i=1; i<3; i++)
        {
          for(j=1; j<3; j++)
            {
              if(local_B[i][j]!=local_A[i][j])
                printf("ERROR : in passing values \n");
            }
        }
    }

  if(rank == 0)
    GA_PRINT_MSG();
  GA_Terminate();
  MPI_Finalize();

}
