/*                                                                                                       
 * Test Program for GA                                                                                   
 * This is to test GA_Scatter (is a one-sided operation)                                                 
 * GA_Create -- used to create a global array using handles like 'g_A'                                   
  */

#include<stdio.h>
#include<stdlib.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"

#define N 5
#define D 2

main(int argc, char **argv)
{
  int rank, nprocs;
  int g_A, dims[D]={5,10}, local_A[N], sub_array[N][D], **s_array=NULL;
  int i, j, value=0;
  
  MPI_Init(&argc, &argv);
  GA_Initialize();
  MA_init(C_INT, 1000, 1000);
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  s_array=(int**)malloc(N*sizeof(int*));
  for(i=0; i<N; i++)
    {
      s_array[i]=(int*)malloc(D*sizeof(int));
      for(j=0; j<D; j++) s_array[i][j]=rand()%5;
    }

  for(i=0; i<N; i++)
    for(j=0; j<D; j++)sub_array[i][j]=rand()%5;
  
  for(i=0; i<N; i++) local_A[i]=rand()%5+1;
  /*
   * depends on the value of array ..we can generate the location values in randon 
   * we can also use the if-condition
   */
  if(rank==0)
    {  
      for(i=0; i<N; i++)
	{
	  for(j=0; j<D; j++)printf("%d ",s_array[i][j]);
	  printf("\n");
	} 
      
      printf("\n");
      for(i=0; i<N; i++)
	{
	  for(j=0; j<D; j++)printf("%d ",sub_array[i][j]);
	  printf("\n");
	} 
      
      printf("\n");
      for(i=0; i<N; i++)printf("%d \n", local_A[i]);
    }

  g_A=NGA_Create(C_INT, D, dims, "array_A", NULL);
  GA_Fill(g_A, &value);
  GA_Sync();
                                                   
  //NGA_Scatter(g_A, local_A, sub_array, N);
  NGA_Scatter(g_A, local_A, s_array, N);
  GA_Sync();
  GA_Print(g_A);
  
  if(rank==1)
    printf(" GA Test: Completed\n");
  GA_Destroy(g_A);
  GA_Terminate();
  MPI_Finalize();
  return 0;
}
