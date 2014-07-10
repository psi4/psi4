/*                                                                                                       
 * Test Program for GA                                                                                   
 * This is to test GA_PGroup_nnodes and GA_Pgroup_nodeid (is a collective operation)                           
 * _pgroup_create -- helps to create different sized group from the given processes and helps to alocate          different job for diffferent group 
 * _pgroup_nnodes -- to count number of nodes in particular processor group
 * _pgroup_nodeid -- to give node-id for nodes present in particular processor group
 * groups are identified using the group handle
 */

#include<stdio.h>
#include<stdlib.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"
#include"ga_unit.h"

even_odd_pgroup(int rank, int nprocs)
{
  int i, j;
  int p_Geven, p_Godd, p_size, mod, p_size_mod, *list_even=NULL, *list_odd=NULL;

  p_size=nprocs/2;
  mod=nprocs%2;
  p_size_mod=p_size+mod;
  list_even = (int*)malloc(p_size*sizeof(int));
  list_odd = (int*)malloc(p_size*sizeof(int));
  
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

  p_Geven=GA_Pgroup_create(list_even, p_size_mod);
  p_Godd=GA_Pgroup_create(list_odd, p_size);

  if(rank%2==0)
    printf("%d: My ID is %d :: %d -- even \n", rank, GA_Pgroup_nodeid(p_Geven), GA_Pgroup_nnodes(p_Geven));
  else
    printf("%d: My ID is %d :: %d --- odd\n", rank, GA_Pgroup_nodeid(p_Godd), GA_Pgroup_nnodes(p_Godd));
}

two_half_pgroup(int rank, int nprocs)
{
  int i, j, k;
  int p_Gone, p_Gtwo, p_size, mod, p_size_mod, *list_one=NULL, *list_two=NULL;

  p_size=nprocs/2;
  mod=nprocs%2;
  p_size_mod=p_size+mod;
  list_one = (int*)malloc(p_size*sizeof(int));
  list_two = (int*)malloc(p_size*sizeof(int));
  
  j=0;
  k=0;
  for(i=0; i<nprocs; i++)
    {
      if(i<p_size)
	{
	  list_one[j]=i;
	  j++;
	}
      else
	{
	  list_two[k]=i;
	  k++;
	}
    }

  p_Gone=GA_Pgroup_create(list_one, p_size);
  p_Gtwo=GA_Pgroup_create(list_two, p_size_mod);

  if(rank<p_size)
    printf("%d: My ID is %d :: %d -- one \n", rank, GA_Pgroup_nodeid(p_Gone), GA_Pgroup_nnodes(p_Gone));
  else
    printf("%d: My ID is %d :: %d --- two\n", rank, GA_Pgroup_nodeid(p_Gtwo), GA_Pgroup_nnodes(p_Gtwo));
}

main(int argc, char **argv)
{
  int rank, nprocs;
  
  MPI_Init(&argc, &argv);
  
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);

  GA_Initialize();
  
  even_odd_pgroup(rank, nprocs);
  printf("\n");
  two_half_pgroup(rank, nprocs);

  GA_Sync();    
  
  if(rank==0)
    GA_PRINT_MSG();
  
  GA_Terminate();
  MPI_Finalize();
}
