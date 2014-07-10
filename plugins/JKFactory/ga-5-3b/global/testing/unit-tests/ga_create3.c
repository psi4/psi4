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
#include"ga_unit.h"

#define SIZE 10
#define MAX_DIM 7

//#define NUM_TYPES 3

void create_ga(int ndim, int datatypes)
{
  int g_A;						
  int dims[MAX_DIM], i, val=4;			

  for(i=0; i<ndim; i++) dims[i]=SIZE;
   
  g_A = NGA_Create(datatypes, ndim, dims, "array_A", NULL);	
    
  if(!g_A)
    GA_Error("GA Error: no global array exists \n", ndim);
  GA_Destroy(g_A);							
}

main(int argc, char **argv)				
{							
  int rank, nprocs, i, j;					
  //  int datatypes[NUM_TYPES] = {C_INT, C_DBL, C_LONG};

  MPI_Init(&argc, &argv);				
							
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);			
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);		
							
  //MA_init(C_INT, 1000, 1000);				
  GA_Initialize();	
  
  for(i=1; i<=MAX_DIM; i++)
    {
      for(j=0; j<NUM_TYPES; j++)
	create_ga(i, TYPES[j]);
    }
  GA_Sync();
  if(rank == 0)
    GA_PRINT_MSG();
  GA_Terminate();
  MPI_Finalize();
}

/* 
 * TO-DO : assign SIZE to a bigger value
 */ 
