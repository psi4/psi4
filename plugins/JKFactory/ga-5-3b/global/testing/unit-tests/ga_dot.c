#include<stdio.h>
#include<stdlib.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"

#define DIM 2
#define SIZE 5 

integer_dot(int rank, int nprocs)
{
  
  int g_A, g_B;
  int dims[DIM]={SIZE,SIZE}, val_A=5, val_B=10, op;
  
  MA_init(C_INT, 1000, 1000);
  
  g_A = NGA_Create(C_INT, DIM, dims, "array_A", NULL);
  g_B = NGA_Create(C_INT, DIM, dims, "array_B", NULL);
  
  GA_Fill(g_A, &val_A);
  GA_Fill(g_B, &val_B);
  
  if(rank==0)
    {
      op=GA_Idot(g_A, g_B);
      
      printf("%d \n", op);
    }
      
}


 /* when i use int_type is working...have find the reason behind the other datatype and how to use it.., */

/*
float_dot(int rank, int nprocs)
{
int dims[DIM]={SIZE,SIZE};
double g_A, g_B; 
double val_A=5, val_B=10, op;

//MA_init(C_DBL, 1000, 1000);

g_A = NGA_Create(C_DBL, DIM, dims, "array_A", NULL);
g_B = NGA_Create(C_DBL, DIM, dims, "array_B", NULL);

GA_Fill(g_A, &val_A);
GA_Fill(g_B, &val_B);

if(rank==0)
op=GA_Ddot(g_A, g_B);

//  printf("%d \n", op);

}
*/

main(int argc, char **argv)
{
  int rank, nprocs;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);
  //  MA_init(C_DBL, 1000, 1000);
  
  GA_Initialize();
  
  integer_dot(rank, nprocs);
  //float_dot(rank, nprocs);
  
  GA_Terminate();
  MPI_Finalize();
}
