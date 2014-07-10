#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <iostream>
#include <cstdlib>
#include <cstdio>
using namespace std;
#include "testutil.h"
#include "ga++.h"

#define N     10    // First dimension
#define NDIM  4     // Number of dimensions
#define BASE  0
#define PERMUTE_

#define GA_DATA_TYPE MT_F_REAL

void
fillPatch(double *ptr, int dim[], int ld[], int ndim, double val) {
  int i, j, stride=1;
  
  switch (ndim){
  case 0: GA_Error((char *)"fill_patch: error",ndim);
  case 1: for(i=0;i <dim[0];i++)ptr[i]=val;
    break;
  case 2: for(i=0; i< dim[0]; i++){
    double *arr = ptr + i*ld[0];
    for(j=0; j< dim[1]; j++)arr[j]=val;
  }
    break;
  default:
    for(i=0; i<ndim-1; i++)stride *=ld[i];
    for(i=0; i<dim[0]; i++){
      double *arr = ptr + stride*i;
      fillPatch(arr, dim+1, ld+1, ndim-1, val);
    }
  }
}

int 
main(int argc, char *argv[]) {
 
  int me, nproc, proc, loop;
  int dims[NDIM], lo[NDIM], hi[NDIM], block[NDIM], ld[NDIM-1];
  int i,d,*proclist, offset;
  int adims[NDIM], ndim,type;
  typedef struct {
    int lo[NDIM];
    int hi[NDIM];
  } patch_t;
  patch_t *regions;
  int *map;
  double *buf;
  int heap  = 3000000, stack = 3000000;
  
  GA::Initialize(argc, argv, heap, stack, GA_DATA_TYPE, 0);
  me=GA_Nodeid(); nproc=GA_Nnodes();


  /***** create array A with default distribution  *****/
  if(me==0){printf("Creating array A\n"); fflush(stdout);}
  for(i = 0; i<NDIM; i++)
    dims[i] = N*(i+1);
  GA::GlobalArray *g_a = GA::SERVICES.createGA(MT_F_DBL, NDIM, dims, 
					       (char *)"array A", NULL);
  if(me==0)printf("OK\n\n");
  if(me==0)printf("handle = %d\n", g_a->handle());
  /* print info about array we got */
  g_a->inquire(&type, &ndim, adims);
  if(me==0)printf("After Inquire\n");
  g_a->printDistribution();

  GA::SERVICES.sync();
  /* duplicate array A with ga_create irreg rather than ga_duplicate
   * -- want to show distribution control 
   * -- with ga_duplicate it would be g_b=GA_Duplicate(g_a,name) 
   */
  if(me==0)printf("\nReconstructing distribution description for A\n");

  /* get memory for arrays describing distribution */
  proclist = (int*)malloc(nproc*sizeof(int));
  if(!proclist)GA_Error((char *)"malloc failed for proclist",0);
  regions = (patch_t*)malloc(nproc*sizeof(patch_t));
  if(!regions)GA_Error((char *)"malloc failed for regions",0);
  map = (int*)malloc((nproc+ndim)*sizeof(int)); /* ubound= nproc+mdim */
  if(!map)GA_Error((char *)"malloc failed for map",0);

  /* first find out how array g_a is distributed */
  for(i=0;i<ndim;i++)lo[i]=BASE;
  for(i=0;i<ndim;i++)hi[i]=adims[i] -1 + BASE;
  proc = g_a->locateRegion(lo, hi, (int*)regions, proclist);
  if(proc<1) GA_Error((char *)"error in NGA_Locate_region",proc);

  /* determine blocking for each dimension */
  for(i=0;i<ndim;i++)block[i]=0;
  for(i=0;i<ndim;i++)adims[i]=0;
     
  offset =0;
  for(d=0; d<ndim; d++)
    for(i=0;i<proc;i++)
      if( regions[i].hi[d]>adims[d] ){
	map[offset] = regions[i].lo[d];
	offset++;
	block[d]++;
	adims[d]= regions[i].hi[d];
      }
            
  if(me==0){
    printf("Distribution map contains %d elements\n",offset); 
    print_subscript((char *)"number of blocks for each dimension",
		    ndim,block,(char *)"\n");
    print_subscript((char *)"distribution map",offset,map,(char *)"\n\n");
    fflush(stdout);
  }
     
  if(me==0)printf("Creating array B applying distribution of A\n");

#    ifdef USE_DUPLICATE
  GA::GlobalArray *g_b = GA::SERVICES.createGA(g_a, "array B");
#    else
  GA::GlobalArray *g_b = GA::SERVICES.createGA(MT_F_DBL, NDIM, dims, 
					       (char *)"array B", block, map);
#    endif
  if(!g_b) GA_Error((char *)"create failed: B",0); 
  if(me==0)printf("OK\n\n");
  free(proclist); free(regions); free(map);
     
  g_b->printDistribution();

  GA::SERVICES.sync();

  if(me==0){
    printf("\nCompare distributions of A and B\n");
    if(g_a->compareDistr(g_b))
      printf("Failure: distributions NOT identical\n");
    else 
      printf("Success: distributions identical\n");
    fflush(stdout);
  }
       

  if(me==0){
    printf("\nAccessing local elements of A: set them to the owner process id\n");
    fflush(stdout);
  }
  GA::SERVICES.sync();

  g_a->distribution(me, lo, hi);

  if(hi[0]>=0){/* -1 means no elements stored on this processor */
    double *ptr;
    int locdim[NDIM];
    g_a->access(lo, hi, &ptr, ld);
    for(i=0;i<ndim;i++)locdim[i]=hi[i]-lo[i]+1;
    fillPatch(ptr, locdim, ld, ndim,(double)me);
  }

  for(i=0;i<nproc; i++){
    if(me==i && hi[0]>=0){
      char msg[100];
      sprintf(msg,(char *)"%d: leading dimensions",me);
      print_subscript(msg,ndim-1,ld,(char *)"\n");
      fflush(stdout);
    }
    GA::SERVICES.sync();
  }
     
  GA::SERVICES.sync();
  if(me==0)printf("\nRandomly checking the update using ga_get on array sections\n");
  GA::SERVICES.sync();

  /* show ga_get working and verify array updates 
   * every process does N random gets
   * for simplicity get only a single row at a time
   */
  srand(me); /* different seed for every process */
  hi[ndim-1]=adims[ndim-1] -1 + BASE;
  for(i=1;i<ndim-1; i++)ld[i]=1; ld[ndim-2]=adims[ndim-1] -1 + BASE;

  /* get buffer memory */
  buf = (double*)malloc(adims[ndim-1]*sizeof(double));
  if(!buf)GA_Error((char *)"malloc failed for buf",0);

  /* half of the processes check the result */
  if(me<=nproc/2) 
    for(loop = 0; loop< N; loop++){ /* task parallel loop */
      lo[ndim-1]=BASE;
      for (i= 0; i < ndim -1; i ++){
	lo[i] = hi[i] = rand()%adims[i]+BASE; 
      }

      /* print_subscript("getting",ndim,lo,"\n");*/
      g_a->get(lo, hi, buf, ld); 
         
      /* check values */
      for(i=0;i<adims[ndim-1]; i++){
	int p = g_a->locate(lo);
	if((double)p != buf[i]) {
	  char msg[100];
	  sprintf(msg,(char *)"%d: wrong value: %d != %f a",me, p, buf[i]);
	  print_subscript(msg,ndim,lo,(char *)"\n");
	  GA_Error((char *)"Error - bye",i);  
	}
	lo[ndim-1]++;
      }
    }
             
  free(buf);
  GA::SERVICES.sync();
           
  if(me==0)printf("OK\n");
     
  g_a->destroy();
  g_b->destroy();
  if(!me) cout << "Terminating\n";
  GA::Terminate();
}
