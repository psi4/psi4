/*
 * Test Program for GA
 * This is to test GA_Create (is a collective operation)
 * GA_Create -- used to create a global array using handles like 'g_A' 
 * Here used GA_Inquire to verify that g_A hanle returns the right values of created_array
 */

#include<stdio.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"

#define GSIZE 5

one_dimension(int rank, int nprocs)
{
  int g_V; 
  int ndim=1, dims=GSIZE;
  int dims2, ndim2, type, value=5;

  g_V = NGA_Create(C_INT, ndim, &dims, "array_V", NULL);

  GA_Fill(g_V, &value);
  GA_Print(g_V);
  GA_Sync();

  //  NGA_Inquire(g_V, &type, &ndim2, dims2);
  //printf(" %d -- %d,,\n", type, ndim2);

  //for(i=0; i<DIM; i++)
  //printf("%d: %d[ %d] ...* \n", rank, i, dims2[i]);	
  GA_Destroy(g_V);														
}			
				
two_dimension(int rank, int nprocs)			
{							
  int g_A;						
  int ndim=2, dims[2]={GSIZE,GSIZE};			
  int dims2[ndim], ndim2, type, value=5, i;	
							
  g_A = NGA_Create(C_INT, ndim, dims, "array_A", NULL);	
							
  GA_Fill(g_A, &value);					
  GA_Print(g_A);					
  GA_Sync();						
							
  NGA_Inquire(g_A, &type, &ndim2, dims2);		
  printf(" %d -- %d,,\n", type, ndim2);			
							
  for(i=0; i<ndim; i++)					
    printf("%d: %d[ %d] ...* \n", rank, i, dims2[i]);	
  GA_Destroy(g_A);							
}							

three_dimension(int rank, int nprocs)			
{							
  int g_A;						
  int ndim=3, dims[3]={GSIZE,GSIZE,GSIZE};			
  int dims2[ndim], ndim2, type, value=5, i;	
							
  g_A = NGA_Create(C_INT, ndim, dims, "array_A", NULL);	
							
  GA_Fill(g_A, &value);					
  GA_Print(g_A);					
  GA_Sync();						
							
  NGA_Inquire(g_A, &type, &ndim2, dims2);		
  printf(" %d -- %d,,\n", type, ndim2);			
							
  for(i=0; i<ndim; i++)					
    printf("%d: %d[ %d] ...* \n", rank, i, dims2[i]);	
  GA_Destroy(g_A);							
}							

fourth_dimension(int rank, int nprocs)			
{							
  int g_A;						
  int ndim=4, dims[4]={GSIZE,GSIZE,GSIZE,GSIZE};			
  int dims2[ndim], ndim2, type, value=5, i;	
							
  g_A = NGA_Create(C_INT, ndim, dims, "array_A", NULL);	
							
  GA_Fill(g_A, &value);					
  GA_Print(g_A);					
  GA_Sync();						
							
  NGA_Inquire(g_A, &type, &ndim2, dims2);		
  printf(" %d -- %d,,\n", type, ndim2);			
							
  for(i=0; i<ndim; i++)					
    printf("%d: %d[ %d] ...* \n", rank, i, dims2[i]);	
  GA_Destroy(g_A);							
}							

						
main(int argc, char **argv)				
{							
  int rank, nprocs, i;					

  MPI_Init(&argc, &argv);				
							
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);			
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);		
							
  MA_init(C_INT, 1000, 1000);				
							
  GA_Initialize();					
	
  if(rank==1)  printf(" ONE DIMENSION\n");						
  one_dimension(rank, nprocs);
  GA_Sync();
  if(rank==1)  printf(" TWO DIMENSION\n");						
  two_dimension(rank, nprocs);
  GA_Sync();
  if(rank==1)  printf(" THREE DIMENSION\n");						
  three_dimension(rank, nprocs);

  //  fourth_dimension(rank, nprocs);
  GA_Sync();
  if(rank == 0)
    printf("Test Completed \n");
  GA_Terminate();
  MPI_Finalize();
  
}
