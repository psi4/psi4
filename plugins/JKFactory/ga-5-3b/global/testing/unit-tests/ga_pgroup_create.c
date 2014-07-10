/*                                                                                                       
 * Test Program for GA                                                                                   
 * This is to test GA_PGroup_create (is a collective operation)                                                 
 *
 */

#include<stdio.h>
#include<stdlib.h>

#include"mpi.h"
#include"ga.h"
#include"macdecls.h"

#define PGSIZE 4
main(int argc, char **argv)
{
  int rank, nprocs, i, j;
  int p_Geven, p_Godd, list[PGSIZE], p_size, *list_even=NULL, *list_odd=NULL, even, odd;

  MPI_Init(&argc, &argv);

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MA_init(C_INT, 1000, 1000);

  GA_Initialize();

  for(i=0; i<PGSIZE; i++)
    list[i]=i+4;

  p_size=nprocs/2;

  list_even = (int*)malloc(p_size*sizeof(int));
  list_odd = (int*)malloc(p_size*sizeof(int));
  
  GA_Sync();
 
  
  if(rank==0)
    {
      for(j=0; j<nprocs; j++)
	{

		  if(j%2==0)
		    list_even[i]=j;
		  printf("%d: %d -- %d\n", j, list_even[i], p_size);

		  if(j%2==1)
		    list_odd[i]=j;
		  printf("%d: %d \n", j, list_odd[i]);

	}
      
      for(i=0; i<p_size; i++)
	printf(" %d ::: %d \n", list_even[i], list_odd[i]);
    }
  
  

  //--------------------------------------

  /*
  if(rank==0)
    {
      {
	for(i=0; i<p_size; i++)
	  {
	    for(j=0; j<nprocs; j++)
	      {
		if(j%2==0)
		  {
		    list_even[i]=even;
		    printf("%d: %d -- %d\n", j, list_even[i], p_size);
		  }
		even++;
		
		if(j%2==1)
		  {
		    list_odd[i]=odd;
		    printf("%d: %d \n", j, list_odd[i]);
		  }
		odd++;
	      }
	  }
      }
      
      for(i=0; i<p_size; i++)
	printf(" %d ::: %d \n", list_even[i], list_odd[i]);
    }
  */


  GA_Sync();
  printf("%d: My ID is %d :: %d \n", rank, GA_Nodeid(), GA_Nnodes());

  
  /*
  p_Geven=GA_Pgroup_create(list_even, p_size);
  p_Godd=GA_Pgroup_create(list_odd, p_size);

  
if(rank==0)
    printf("\n");

  GA_Sync();

  if(rank==0 || rank/2==0)
    printf("%d: My ID is %d :: %d \n", rank, GA_Pgroup_nodeid(p_Geven), GA_Pgroup_nnodes(p_Geven));
  else
    printf("%d: My ID is %d :: %d \n", rank, GA_Pgroup_nodeid(p_Godd), GA_Pgroup_nnodes(p_Godd));

    //  p_G=GA_Pgroup_create(list, nprocs);
    */

  if(rank==1)
    printf(" GA: Test Completed \n");

  GA_Terminate();
  MPI_Finalize();
}
