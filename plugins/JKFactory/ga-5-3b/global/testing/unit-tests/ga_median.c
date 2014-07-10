#include<stdio.h>
#include<stdlib.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"

#define DIM 2
#define SIZE 5

 main(int argc, char **argv)
{

  int rank, nprocs;
  int g_A, g_B, g_C, g_M, **local_M=NULL;
  int dims[DIM]={SIZE, SIZE}, val_A=5, val_B=4, val_C=3;
  int mlo[DIM]={0,0}, mhi[DIM]={4,4}, ld=SIZE, i, j;
  int local_T[SIZE][SIZE];
  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);
  GA_Initialize();

  g_A = NGA_Create(C_INT, DIM, dims, "array_A", NULL);
  g_B = NGA_Create(C_INT, DIM, dims, "array_B", NULL);
  g_C = NGA_Create(C_INT, DIM, dims, "array_C", NULL);

  g_M = NGA_Create(C_INT, DIM, dims, "array_M", NULL);

  GA_Fill(g_A, &val_A);
  GA_Fill(g_B, &val_B);
  GA_Fill(g_C, &val_C);

  GA_Median(g_A,g_B, g_C, g_M);

  GA_Print(g_M);

  /*
  GA_Print(g_A);
  GA_Print(g_B);
  GA_Print(g_C);
  */

  if(rank==0)
    {
      local_M=(int**)malloc(SIZE*sizeof(int*));
      for(i=0; i<SIZE; i++)
	local_M[i]=(int*)malloc(SIZE*sizeof(int));
      printf("check 1\n");
      
      for(i=0; i<SIZE; i++)
	{
	  for(j=0; j<SIZE; j++) printf("%d ", local_M[i][j]=rand()%5);
	  printf("\n");
	}
      
      for(i=0; i<SIZE; i++)
	{
	  for(j=0; j<SIZE; j++)
	    if((val_A+val_B+val_C)/3==local_M[i][j]) printf("GA ERROR: [%d, %d]\n", i, j);
	}
      
      printf("check 2\n");
      NGA_Get(g_M, mlo, mhi, local_M, &ld);
      NGA_Get(g_M, mlo, mhi, local_T, &ld);

      for(i=0; i<SIZE; i++)
	{
	  for(j=0; j<SIZE; j++) printf("%d ", local_M[i][j]);
	  printf("\n");
	}

      printf("\n");
      for(i=0; i<SIZE; i++)
	{
	  for(j=0; j<SIZE; j++) printf("%d ", local_T[i][j]);
	  printf("\n");
	}

      for(i=0; i<SIZE; i++)
	{
	  for(j=0; j<SIZE; j++)
	    if((val_A+val_B+val_C)/3==local_M[i][j]) printf("GA ERROR: \n");
	}
    }

  GA_Terminate();
  MPI_Finalize();
}
