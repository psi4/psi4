/*
 * Test Program for GA
 * GA_Set_restricted -- 
 */

#include<stdio.h>
#include<stdlib.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"
#include"ga_unit.h"

#define SIZE 5
#define DIM 2

list_by_two(int nprocs)
{
  int i, l_size, mod, l_size_mod, *list_one=NULL, *list_two=NULL;

  l_size = nprocs / 2;
  mod = nprocs % 2;
  l_size_mod = l_size + mod;

  list_one = (int*)malloc(l_size*sizeof(int));
  list_two = (int*)malloc(l_size*sizeof(int));

  for(i=0; i<l_size; i++)
    list_one[i] = i;

  for(i=0; i<l_size_mod; i++)
    list_two[i] = l_size + i;

  printf("\n");
  printf("list one \n");
  for(i=0; i<l_size; i++)
    printf(" %d : ", list_one[i]);
  printf("\n");

  printf("list two \n");
  for(i=0; i<l_size_mod; i++)
    printf(" %d : ", list_two[i]);
  printf("\n");
}

group_of_four(int nprocs)
{
  int i, l_size, *list_a=NULL, *list_b=NULL, *list_c=NULL, *list_d=NULL;
  l_size = nprocs / 4;
  list_a = (int*)malloc(l_size*sizeof(int));
  list_b = (int*)malloc(l_size*sizeof(int));
  list_c = (int*)malloc(l_size*sizeof(int));
  list_d = (int*)malloc(l_size*sizeof(int));

  for(i=0; i<l_size; i++)
    {
    list_a[i] = i;
    list_b[i] = l_size + i;
    list_c[i] = (l_size * 2) + i;
    list_d[i] = (l_size * 3) + i;
    }

  printf("list_a : list_b: list_c: list_d \n");
  for(i=0; i<l_size; i++)
    printf("-%d\t-%d\t-%d\t-%d\t\n", list_a[i], list_b[i], list_c[i], list_d[i]);
}

main(int argc, char **argv)
{
  int rank, nprocs;
  int g_A, g_B;
  int dims[DIM]={SIZE,SIZE}, *list_even=NULL, *list_odd=NULL, i, j, val=SIZE, l_size, mod, l_size_mod;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);
  GA_Initialize();

  if(rank ==0)
    {
      l_size = nprocs / 2;
      mod = nprocs % 2;
      l_size_mod = l_size + mod;
      list_even = (int*)malloc(l_size*sizeof(int));
      list_odd = (int*)malloc(l_size*sizeof(int));
      
      j=0;
      for(i=0; i<nprocs; i++)
	{
	  if(i%2==0)
	    {
	      list_even[j]=i;
	      j++;
	    }
	}
      j=0;
      for(i=0; i<nprocs; i++)
	{
	  if(i%2==1)
	    {
	      list_odd[j]=i;
	      j++;
	    }
	}

      printf("Even list \n");
      for(i=0; i<l_size_mod; i++)
	printf(" %d : ", list_even[i]);
      
      printf("\nOdd list \n");
      for(i=0; i<l_size; i++)
	printf(" %d : ", list_odd[i]);

      // seperate function for creating two list by splite nprocs in two group of list 
      list_by_two(nprocs);

      // function for group of four :
      if(nprocs % 4 == 0)
	group_of_four(nprocs);
    }


  g_A=GA_Create_handle();

  GA_Set_data(g_A, DIM, dims, C_INT);
  GA_Set_chunk(g_A, NULL);
  GA_Set_array_name(g_A, "array_A");

  GA_Allocate(g_A);
  
  if(!g_A)
    GA_ERROR_MSG();

  GA_Fill(g_A, &val);
  GA_Print(g_A);

  GA_Set_restricted(g_A, list_even, l_size_mod);
  //  GA_Print(g_A);

  GA_Sync();
  if(rank == 0)
    GA_PRINT_MSG();
 
  GA_Terminate();
  MPI_Finalize();

}
