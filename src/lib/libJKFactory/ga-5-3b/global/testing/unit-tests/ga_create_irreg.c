/*
 * Test Program for GA
 * This is to test GA_Create_irreg (is a collective operation)
 * GA_Create -- used to create a global array of regular size
 * GA_Create_IRREG -- used to create G_array of irregular size -- helps user to define the distribution 
 * Here used GA_Inquire to verify that g_A hanle returns the right values of created_array
 */

#include<stdio.h>
#include<stdlib.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"

#define DIM 2

#define GSIZE 10

irregular_array1(int rank)
{

  int g_A, g_B; 
  int dims[DIM]={5,10}, dims2[DIM], ndim, type, value=5, block[DIM]={2,3}, map[5]={0,2,0,4,6}, val=7;
  int n_block[DIM], block_dims[DIM], i;

  g_A = NGA_Create(C_INT, DIM, dims, "array_A", NULL);
  g_B = NGA_Create_irreg(C_INT, DIM, dims, "array_B", block, map);

  GA_Fill(g_A, &value);
  GA_Print(g_A);

  GA_Fill(g_B, &val);
  GA_Print(g_B);
  GA_Sync();

  NGA_Inquire(g_A, &type, &ndim, dims2);
  //printf(" %d -- %d,,\n", type, ndim);

  /*
  GA_Get_block_info(g_B, n_block, block_dims);
  for(i=0; i<DIM; i++)
    printf(" %d:  %d ___ %d --- \n", rank, n_block[i], block_dims[i]);
  */

  GA_Destroy(g_A);
  GA_Destroy(g_B);
}

irregular_array2(int rank)
{

  int g_A, g_B; 
  int dims[DIM]={GSIZE,GSIZE}, dims2[DIM], block[DIM]={3,2}, map[5]={0,2,6,0,4}, val=7;
  int n_block[DIM], block_dims[DIM], i;

  g_A = NGA_Create(C_INT, DIM, dims, "array_A", NULL);
  g_B = NGA_Create_irreg(C_INT, DIM, dims, "array_B", block, map);

  GA_Fill(g_B, &val);
  GA_Print(g_B);
  GA_Sync();

  /*
  GA_Get_block_info(g_B, n_block, block_dims);
  for(i=0; i<DIM; i++)
    printf(" %d:  %d ___ %d --- \n", rank, n_block[i], block_dims[i]);
  */

  GA_Destroy(g_A);
  GA_Destroy(g_B);
}

/* In these function the values of blocks and maps are auto-generated number based on number of processes we use */

auto_number1(int rank, int nprocs)
{

  int g_A, g_B; 
  int dims[DIM]={GSIZE, GSIZE}, dims2[DIM], block[DIM], *map=NULL, val=7;
  int n_block[DIM], block_dims[DIM], b_temp, i;
  int b1, b2, inc=0;

  do{
      b1=DIM+inc;
      b2=nprocs/b1;
      inc++;
    }while(nprocs/b1>=GSIZE);

  block[0]=b1;
  block[1]=b2;

  map=(int*)malloc(nprocs*sizeof(int));

  for(i=0; i<b1; i++)
    map[i]=i*DIM;

  for(i=b1; i<(b2+b1); i++)
    map[i]=i-b1;

  if(rank==0)
    {
      for(i=0; i<(b1+b2); i++)
	printf("map[%d] - %d\n", i, map[i]);
      for(i=0; i<DIM; i++)
	printf("BLOCK[%d] - %d\n", i, block[i]);
    }

  g_B = NGA_Create_irreg(C_INT, DIM, dims, "array_B", block, map);

  GA_Fill(g_B, &val);
  GA_Print(g_B);
  GA_Sync();

  if(rank==1)
    {
      GA_Get_block_info(g_B, n_block, block_dims);
      for(i=0; i<DIM; i++)
	printf(" %d:  %d --- %d ... %d\n", rank, n_block[i], block_dims[i], b_temp);
    }
  GA_Destroy(g_B);
}

auto_number2(int rank, int nprocs)
{

  int g_A, g_B; 
  int dims[DIM]={GSIZE, GSIZE}, dims2[DIM], block[DIM], *map=NULL, val=7;
  int n_block[DIM], block_dims[DIM], b_temp, i;
  int b1, b2, inc=0;

  do{
    
    b2=DIM+inc;
    b1=nprocs/b2;    
      inc++;
    }while(nprocs/b2>=GSIZE);

  block[0]=b1;
  block[1]=b2;

  map=(int*)malloc(nprocs*sizeof(int));

  for(i=0; i<b1; i++)
    map[i]=i;

  for(i=b1; i<(b2+b1); i++)
    map[i]=i-b1;

  if(rank==0)
    {
      for(i=0; i<(b1+b2); i++)
	printf("map[%d] - %d\n", i, map[i]);
      for(i=0; i<DIM; i++)
	printf("BLOCK[%d] - %d\n", i, block[i]);
    }

  g_B = NGA_Create_irreg(C_INT, DIM, dims, "array_B", block, map);

  GA_Fill(g_B, &val);
  GA_Print(g_B);
  GA_Sync();

  if(rank==1)
    {
      GA_Get_block_info(g_B, n_block, block_dims);
      for(i=0; i<DIM; i++)
	printf(" %d:  %d --- %d ... %d\n", rank, n_block[i], block_dims[i], b_temp);
    }
  GA_Destroy(g_B);
}

main(int argc, char **argv)
{
  int rank, nprocs;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);

  GA_Initialize();

  //irregular_array1(rank);  
  //  irregular_array2(rank);  
 
  //  auto_number1(rank, nprocs);
  auto_number2(rank, nprocs);

  if(rank == 0)
    printf("Test Completed \n");

  GA_Terminate();
  MPI_Finalize();
  
}
