/*                                                                                                       
 * Test Program for GA                                                                                   
 * This is to test GA_PGroup_create (is a collective operation)                                                 
 * _pgroup_create -- helps to create different sized group from the given processes and helps to alocate different    job for diffferent group 
 * groups are identified using the group handle
 */

#include<stdio.h>
#include<stdlib.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"

main(int argc, char **argv)
{
  int rank, nprocs, i, j;
  int p_Geven, p_Godd, p_size, *list_even=NULL, *list_odd=NULL;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);

  GA_Initialize();

  p_size=nprocs/2;
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
  if(rank==0)
    {    
      for(i=0; i<p_size; i++)
	printf(" %d--> %d --:: -- %d\n", i, list_even[i], list_odd[i]);
    }
  printf("\n");  
  GA_Sync();
  printf("%d: My ID is %d :: %d \n", rank, GA_Nodeid(), GA_Nnodes());

  p_Geven=GA_Pgroup_create(list_even, p_size);
  p_Godd=GA_Pgroup_create(list_odd, p_size);
 
  GA_Sync();
  if(rank==1)
    printf(" GA: Test Completed \n");

  GA_Terminate();
  MPI_Finalize();
}
