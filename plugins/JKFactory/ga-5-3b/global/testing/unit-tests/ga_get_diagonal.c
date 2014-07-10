/*
 * Test Program for GA
 * This is to test GA_Get_diagonal (is a collective operation)
 * GA_Create -- used to create a global array using handles like 'g_A'
 * handle G_V -- represents the array of single dimension  
 *
 */

#include<stdio.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"
#include"ga_unit.h"

#define DIM 2
#define GSIZE 5

main(int argc, char **argv)
{
  int rank, nprocs, i, j;
  int g_A, g_B, g_V, g_V2;
  int dims[DIM]={GSIZE,GSIZE}, alo[DIM]={1,1}, ahi[DIM]={3,2}, blo[DIM]={1,1}, bhi[DIM]={2,3}, ld=5;
  int value=5, val2=4, dims_V=GSIZE;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);

  GA_Initialize();
  
  g_A = NGA_Create(C_INT, DIM, dims, "array_A", NULL);
  g_V =  NGA_Create(C_INT, 1, &dims_V, "array_V", NULL);
  g_V2 =  NGA_Create(C_INT, 1, &dims_V, "array_V2", NULL);
  
  GA_Fill(g_A, &value);
  GA_Zero_diagonal(g_A);
  GA_Print(g_A);

  g_B = GA_Duplicate(g_A, "array_B");
  GA_Fill(g_B, &val2);
  GA_Sync();
  NGA_Zero_patch(g_B, blo, bhi);
  GA_Print(g_B);

  GA_Fill(g_V, &val2);
  GA_Print(g_V);

  GA_Add_diagonal(g_A, g_V);
  GA_Add_diagonal(g_B, g_V);

  GA_Print(g_A);
  GA_Sync();

  GA_Get_diag(g_B, g_V2);
  GA_Print(g_V2);
  GA_Print(g_B);

  if(rank == 1)
    GA_PRINT_MSG();

  GA_Terminate();
  MPI_Finalize();

  
}
