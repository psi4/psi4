/*
 * Test Program for GA
 * This is to test GA_Add (is a collective operation)
 * GA_Create -- used to create a global array using handles like 'g_A'
 * GA_Duplicate --used to duplicate and generate one more global array.., handle 'g_A' to 'g_B'
 * 
 *
 */

#include<stdio.h>
#include<stdlib.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"

#define DIM 2
#define SIZE 5
main(int argc, char **argv)
{
  int rank, nprocs, i, j;
  int g_A, g_B, g_C, **local_C=NULL, dims[DIM]={SIZE,SIZE}, val1=5, val2=4, alpha=3, beta=2;
  int clo[DIM]={SIZE-SIZE,SIZE-SIZE}, chi[DIM]={SIZE-1,SIZE-1}, ld=SIZE;
  int **local_tm=NULL;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);

  GA_Initialize();

  local_C=(int**)malloc(SIZE*sizeof(int*));
  for(i=0; i<SIZE; i++)
    local_C[i]=(int*)malloc(SIZE*sizeof(int));

  local_tm=(int**)malloc(SIZE*sizeof(int*));
  for(i=0; i<SIZE; i++)
    local_tm[i]=(int*)malloc(SIZE*sizeof(int));
  
  g_A = NGA_Create(C_INT, DIM, dims, "array_A", NULL);

  g_B = GA_Duplicate(g_A, "array_B");
  g_C = GA_Duplicate(g_A, "array_C");

  GA_Fill(g_A, &val1);
  GA_Fill(g_B, &val2);
  
  GA_Add(&alpha, g_A, &beta, g_B, g_C);

  GA_Sync();

  GA_Print(g_A);
  GA_Print(g_B);
  GA_Print(g_C);

  //printf("check 1\n");
  NGA_Get(g_C, clo, chi, local_tm, &ld);
  //printf("check 2\n");
  // GA_Sync();
  if(rank==0)
    {
      for(i=0; i<SIZE; i++)
	{
	  for(j=0; j<SIZE; j++)printf("%d ", local_tm[i][j]);
	  printf("\n");
	}
    }
  /*
  if(rank==0)
    {
      NGA_Get(g_C, clo, chi, local_C, &ld);
      
      printf("check 1 \n");
      
      for(i=0; i<SIZE; i++)
	{
	  for(j=0; j<SIZE; j++)printf("%d ", local_C[i][j]);
	  printf("\n");
	}
      
      printf("check 2\n");
      
      for(i=0; i<SIZE; i++)
	{
	  for(j=0; j<SIZE; j++)
	    if(local_C[i][j]!=(alpha*val1)+(beta*val2)) printf("GA Error : \n");
	}
    }

  */
  //GA_Sync();
  if(rank==0)
    printf("Test Completed \n");

  //  GA_Sync();

  /*
  GA_Destroy(g_A);
  GA_Destroy(g_B);
  GA_Destroy(g_C);
  */

  //*******************************************************************

  /* what would be the possible reason for GA_destroy to get failed .., 
   * solve this before consolidate the whole
   */

  GA_Terminate();
  MPI_Finalize();
 
}
