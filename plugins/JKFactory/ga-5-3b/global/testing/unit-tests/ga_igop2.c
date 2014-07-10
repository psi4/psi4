/*
 * Test Program for GA
 * This is to test GA_Igop (is a collective operation)
 * 
 * 
 */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"
#include"ga_unit.h"

#define NUM 10

addition_operator(int *x, int nprocs)
{
  int i, temp[NUM];
  char *op="min";
  for(i=0; i<NUM; i++) temp[i]=x[i];

  GA_Igop(x, NUM, op);
  for(i=0; i<NUM; i++) printf("%d :: %d\n", x[i], temp[i]);
  printf("\n");
}

checking_operator (int rank, int nprocs)
{
  
  int x[NUM], i, temp[NUM];

  for(i=0; i<NUM; i++) x[i]=rand()%7;
  for(i=0; i<NUM; i++) temp[i]=x[i];
  
  if(rank==0)
    for(i=0; i<NUM; i++) printf("%d \n", x[i]);
    
  GA_Igop(x, NUM, "max");
  printf("\n");

  addition_operator(x, nprocs);
 
  /*
    
  multiply_operator(rank, nprocs, n);
  max_operator(rank, nprocs, n);
  min_operator(rank, nprocs, n);
  absmax_operator(rank, nprocs, n);
  absmin_operator(rank, nprocs, n);
  */

  /*
    
  if(rank==0)
  {
  for(i=0; i<NUM; i++) printf("%d :: %d\n", x[i], temp[i]);
  
  for(i=0; i<n; i++)
  if(x[i]!=temp[i]*sizeof(int)) printf(" GA ERROR: \n");
  
  }
  
  */
}

main(int argc, char **argv)
{

  int rank, nprocs;

  MPI_Init(&argc, &argv);
  GA_Initialize();
  MA_init(C_INT, 1000, 1000);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
 
  checking_operator (rank, nprocs);
  
  GA_Sync();

  if(rank==0)
    GA_PRINT_MSG();
  
  GA_Terminate();

  MPI_Finalize();

}
