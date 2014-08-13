/*
 * Test Program for GA
 * This is to test GA_Create_irreg (is a collective operation)
 * GA_Create -- used to create a global array of regular size
 * GA_Create_IRREG -- used to create G_array of irregular size -- helps user to define the distribution 
 *
 * GA_Get_Block_info -- used to obtain info (like number of blocks and block dimension) about the global array   
 */

#include<stdio.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"

#define DIM 2

main(int argc, char **argv)
{
  int rank, nprocs, i;
  int g_A, g_B; 
  int dims[DIM]={5,10}, dims2[DIM], ndim, type, value=5, block[DIM]={2,3}, map[5]={0,2,0,4,6}, val=7;
  int n_block[DIM], block_dims[DIM];

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);

  GA_Initialize();
  
  g_A = NGA_Create(C_INT, DIM, dims, "array_A", NULL);
  g_B = NGA_Create_irreg(C_INT, DIM, dims, "array_B", block, map);

  GA_Fill(g_A, &value);
  GA_Print(g_A);

  GA_Fill(g_B, &val);
  GA_Print(g_B);
  GA_Sync();

  GA_Get_block_info(g_B, n_block, block_dims);
  for(i=0; i<DIM; i++)
    printf(" %d:  %d ___ %d --- \n", rank, n_block[i], block_dims[i]);

  if(rank == 0)
    printf("Test Completed \n");
  GA_Terminate();
  MPI_Finalize();
}
