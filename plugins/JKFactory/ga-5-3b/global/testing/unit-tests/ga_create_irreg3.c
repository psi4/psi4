/*
 * Test Program for GA
 * This is to test GA_Create_irreg (is a collective operation)
 * GA_Create -- used to create a global array using handles like 'g_A' 
 * Here used GA_Inquire to verify that g_A hanle returns the right values of created_array
 */

#include<stdio.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"

#include"ga_unit.h"

#define SIZE 5
#define MAX_DIM 2

create_irreg_ga(int ndim)
{
  int g_A;						
  int dims[ndim], i;			
  int block[ndim], *map=NULL;

  for(i=0; i<ndim; i++) dims[i]=SIZE;
  
  g_A = NGA_Create(C_INT, ndim, dims, "array_A", NULL);	

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
  
  for(i=2; i<=MAX_DIM; i++)
    create_irreg_ga(i);
  
  GA_Sync();
  if(rank == 0)
    GA_PRINT_MSG();

  GA_Terminate();
  MPI_Finalize();
  
}
