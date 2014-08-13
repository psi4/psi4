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

#define N 100            /* dimension of matrices */


int main( int argc, char **argv ) {
  int g_a, g_b, i, j, size, size_me;
  int icnt, idx, jdx, ld;
  int n=N, type=MT_F_INT, one;
  int *values, *ptr;
  int **indices;
  int dims[2]={N,N};
  int lo[2], hi[2];

  int heap=3000000, stack=2000000;
  int me, nproc;

  int datatype, elements;
  double *prealloc_mem;
  MP_INIT(argc,argv);

  GA_INIT(argc,argv);                            /* initialize GA */
  me=GA_Nodeid(); 
  nproc=GA_Nnodes();
  if(me==0) {
    if(GA_Uses_fapi())GA_Error("Program runs with C array API only",1);
    printf("\nUsing %ld processes\n",(long)nproc);
    fflush(stdout);
  }

  heap /= nproc;
  stack /= nproc;
  if(! MA_init(MT_F_DBL, stack, heap)) 
    GA_Error("MA_init failed",stack+heap);  /* initialize memory allocator*/ 

  /* Create a regular matrix. */
  if(me==0)printf("\nCreating matrix A of size %d x %d\n",N,N);
  g_a = NGA_Create(type, 2, dims, "A", NULL);
  if(!g_a) GA_Error("create failed: A",n); 

  /* Fill matrix using scatter routines */
  size = N*N;
  if (size%nproc == 0) {
    size_me = size/nproc;
  } else {
    i = size - size%nproc;
    size_me = i/nproc;
    if (me < size%nproc) size_me++;
  }

  /* Check that sizes are all okay */
  i = size_me;
  GA_Igop(&i,1,"+");
  if (i != size) {
    GA_Error("Sizes don't add up correctly: ",i);
  } else if (me==0) {
    printf("\nSizes add up correctly\n");
  }

  /* Allocate index and value arrays */
  indices = (int**)malloc(size_me*sizeof(int*));
  values = (int*)malloc(size_me*sizeof(int));
  icnt = me;
  for (i=0; i<size_me; i++) {
    values[i] = icnt;
    idx = icnt%N; 
    jdx = (icnt-idx)/N;
    if (idx >= N || idx < 0) {
      printf("p[%d] Bogus index i: %d\n",me,idx);
    }
    if (jdx >= N || jdx < 0) {
      printf("p[%d] Bogus index j: %d\n",me,jdx);
    }
    indices[i] = (int*)malloc(2*sizeof(int));
    (indices[i])[0] = idx;
    (indices[i])[1] = jdx;
    icnt += nproc;
  }

  /* Scatter values into g_a */
  NGA_Scatter(g_a, values, indices, size_me);
  GA_Sync();

  /* Check to see if contents of g_a are correct */
  NGA_Distribution( g_a, me, lo, hi );
  NGA_Access(g_a, lo, hi, &ptr, &ld);
  for (i=lo[0]; i<hi[0]; i++) {
    idx = i-lo[0];
    for (j=lo[1]; j<hi[1]; j++) {
      jdx = j-lo[1];
      if (ptr[idx*ld+jdx] != j*N+i) {
        printf("p[%d] (Scatter) expected: %d actual: %d\n",me,j*N+i,ptr[idx*ld+jdx]);
      }
    }
  }
  if (me==0) printf("\nCompleted test of NGA_Scatter\n");

  for (i=0; i<size_me; i++) {
    values[i] = 0;
  }
  GA_Sync();
  NGA_Gather(g_a, values, indices, size_me);
  icnt = me;
  for (i=0; i<size_me; i++) {
    if (icnt != values[i]) {
      printf("p[%d] (Gather) expected: %d actual: %d\n",me,icnt,values[i]);
    }
    icnt += nproc;
  }
  if (me==0) printf("\nCompleted test of NGA_Gather\n");
  GA_Sync();

  /* Scatter-accumulate values back into GA*/
  one = 1;
  NGA_Scatter_acc(g_a, values, indices, size_me, &one);
  GA_Sync();

  /* Check to see if contents of g_a are correct */
  for (i=lo[0]; i<hi[0]; i++) {
    idx = i-lo[0];
    for (j=lo[1]; j<hi[1]; j++) {
      jdx = j-lo[1];
      if (ptr[idx*ld+jdx] != 2*(j*N+i)) {
        printf("p[%d] (Scatter_acc) expected: %d actual: %d\n",me,2*(j*N+i),ptr[idx*ld+jdx]);
      }
    }
  }
  if (me==0) printf("\nCompleted test of NGA_Scatter_acc\n");
  NGA_Release(g_a, lo, hi);

  /* Test fixed buffer size */
  NGA_Alloc_gatscat_buf(size_me);

  /* Scatter-accumulate values back into GA*/
  GA_Sync();
  NGA_Scatter_acc(g_a, values, indices, size_me, &one);
  GA_Sync();

  /* Check to see if contents of g_a are correct */
  for (i=lo[0]; i<hi[0]; i++) {
    idx = i-lo[0];
    for (j=lo[1]; j<hi[1]; j++) {
      jdx = j-lo[1];
      if (ptr[idx*ld+jdx] != 3*(j*N+i)) {
        printf("p[%d] (Scatter_acc) expected: %d actual: %d\n",me,3*(j*N+i),ptr[idx*ld+jdx]);
      }
    }
  }
  if (me==0) printf("\nCompleted test of NGA_Scatter_acc using fixed buffers\n");
  NGA_Release(g_a, lo, hi);
  NGA_Free_gatscat_buf();

  GA_Destroy(g_a);
  if(me==0)printf("\nSuccess\n");
  GA_Terminate();

  MP_FINALIZE();

 return 0;
}
