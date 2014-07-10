#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_MATH_H
#   include <math.h>
#endif

#include "ga.h"
#include "macdecls.h"
#include "mp3.h"

#define N 4            /* dimension of matrices */


int main( int argc, char **argv ) {
  int g_a, g_b,i;
  int n=N, type=MT_F_DBL;
  int dims[6]={N,N,N,N,N,N};
  int lo[6], hi[6];

  int heap=30000, stack=20000;
  int me, nproc;

  int datatype, elements;
  double *prealloc_mem;

  MP_INIT(argc,argv);

  GA_INIT(argc,argv);                            /* initialize GA */
  me=GA_Nodeid(); 
  nproc=GA_Nnodes();
  if(me==0) {
    if(GA_Uses_fapi())GA_Error("Program runs with C array API only",1);
    printf("Using %ld processes\n",(long)nproc);
    fflush(stdout);
  }

  heap /= nproc;
  stack /= nproc;
  if(! MA_init(MT_F_DBL, stack, heap)) 
    GA_Error("MA_init failed",stack+heap);  /* initialize memory allocator*/ 


  /* This is a regular matrix. */
  if(me==0)printf("Creating matrix A\n");
  g_a = NGA_Create(type, 2, dims, "A", NULL);
  if(!g_a) GA_Error("create failed: A",n); 
  if(me==0)printf("OK\n");
  NGA_Distribution( g_a, me, lo, hi );

#if 1
  /* This is just allocating and freeing memory. */
  datatype = type;
  elements = 1;
  for ( i=0; i<2; i++ ) {
    elements *= (hi[i] - lo[i] + 1);
  }
  prealloc_mem = GA_Getmem(datatype, elements, 0);
  for ( i=0; i<elements; i++ ) prealloc_mem[i] = 3.141592654;
  GA_Freemem(prealloc_mem);
  if(me==0){printf("getmem&freemem OK\n"); fflush(stdout); }
#endif

  /* This is a matrix using preallocated memory. */
  datatype = type;
  elements = 1;
  for ( i=0; i<2; i++ ) {
    elements *= (hi[i] - lo[i] + 1);
  }
  prealloc_mem = GA_Getmem(datatype, elements, 0);
  g_b = GA_Assemble_duplicate(g_a, "Matrix B", (void*)prealloc_mem );
  if(GA_Compare_distr(g_a,g_b)){ 
     if(me==0){printf("GA_Assemble_duplicate failed\n"); fflush(stdout); }
  }else{
     if(me==0){printf("GA_Assemble_duplicate OK\n"); fflush(stdout); }
  }
  GA_Destroy(g_b);

  GA_Destroy(g_a);
  if(me==0)printf("\nSuccess\n");
  GA_Terminate();

  MP_FINALIZE();

 return 0;
}

