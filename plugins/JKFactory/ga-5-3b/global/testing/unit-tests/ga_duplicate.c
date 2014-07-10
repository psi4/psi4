/*
 * Test Program for GA
 * This is to test GA_Duplicate (is a collective operation)
 * GA_Create -- used to create a global array using handles like 'g_A'
 * GA_Duplicate --helps to duplicate the content from handle 'g_A' 
 * Here used GA_Inquire to verify that g_A hanle returns the right values of created_array
 */

#include<stdio.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"

#define DIM 2

main(int argc, char **argv)
{
  int rank, nprocs, i;
  double g_A, g_B; 
  int dims[DIM]={5,5}, dims2[DIM], ndim2, type2, dims3[DIM], ndim3, type3;
  double value=5, val2=4;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_DBL, 1000, 1000);

  GA_Initialize();
  
  g_A = NGA_Create(C_DBL, DIM, dims, "array_A", NULL);
  GA_Fill(g_A, &value);

  g_B = GA_Duplicate(g_A, "array_B");
  GA_Print(g_A);
  GA_Sync();
  GA_Fill(g_B, &val2);
  GA_Print(g_B);

  if(rank==0)
    {
      NGA_Inquire(g_A, &type2, &ndim2, dims2);
      NGA_Inquire(g_A, &type3, &ndim3, dims3);
      printf(" %d -- %d,,\n", type2, ndim2);
      if(type2!=type3 || ndim2!=ndim3 )
	printf("ERROE : \n");

      if(type2==type3 || ndim2==ndim3 )
	printf("ERROE : Equal \n");

      for(i=0; i<DIM; i++)
	{
	  if(dims2[i]!=dims3[i])
	    printf("ERROE : \n");
	  printf("%d: %d[ %d] ...* \n", rank, i, dims2[i]);
	}
    }

  if(rank == 1)
    printf("Test Completed \n");
  GA_Terminate();
  MPI_Finalize();

  
}
