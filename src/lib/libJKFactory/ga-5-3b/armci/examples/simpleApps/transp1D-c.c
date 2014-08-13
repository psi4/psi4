#if HAVE_CONFIG_H
#   include "config.h"
#endif

/** @file
 * transpose of 1-d array.
 * E,g: (1 2 3 4 5 6 7 8 9 10) => (10 9 8 7 6 5 4 3 2 1)
 */
#define   TOTALELEMS   1007031

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_MATH_H
#   include <math.h>
#endif
#if HAVE_WINDOWS_H
#   include <windows.h>
#endif
#if HAVE_UNISTD_H
#   include <unistd.h>
#endif
#if HAVE_ASSERT_H
#   include <assert.h>
#endif
#include <stdlib.h>

#include "armci.h"
#include "message.h"

void VERIFY(void **b_ptr, int *dims, int *map) {
    int i, j, length, icnt, ichk, lmin, lmax;
    int *buf, *b;
    void *src_ptr, *dst_ptr;
    int me, nprocs;

    /* Find local processor ID and number of processors */
    me     = armci_msg_me();
    nprocs = armci_msg_nproc();

    /* Process 0 verifies that inversion is correct. Start by allocating
       buffer and guarantee that it is big enough */

    length = (int)(((double)dims[0])/((double)nprocs)) + 1;
    buf = (int*)malloc(length*sizeof(int));
    if (me == 0) {
      icnt = 0;
      ichk = 0;
      for (i=0; i<nprocs; i++) {
        /* Find min and max indices owned by processor i */
        lmin = map[i];
        if (i<nprocs-1) {
          lmax = map[i+1]-1;
        } else {
          lmax = dims[0]-1;
        }
        /* evaluate parameters for get call */
        length = sizeof(int)*(lmax-lmin+1);
        src_ptr = b_ptr[i];
        dst_ptr = (void*)buf;
        ARMCI_Get(src_ptr, dst_ptr, length, i);

        /* check values in buffer */
        length = lmax-lmin+1;
        b = (int*)dst_ptr;
        for (j=0; j<length; j++) {
          /* printf("p[%d] b[%d]: %d\n",me,icnt,b[j]); */
          if (b[j] != dims[0] - icnt) {
            printf("Error found for element %d b: %d != a: %d\n",
                icnt,b[j],dims[0]-icnt);
            ichk = 1;
          }
          icnt++;
        }
      }
      if (ichk == 0) {
        printf("1D transpose successful. No errors found\n");
      } else {
        printf("1D transpose failed\n");
      }
    }
    free(buf);
}

void TRANSPOSE1D() {
    
    int dims[1];
    int nelem, i, ierr, min, max, cmin, cmax, lmin, lmax, pmin, pmax;    
    int src_offset, dst_offset, length;
    int *buf, *map;
    void *src_ptr, *dst_ptr;
    void **a_ptr, **b_ptr;
    int *a, *b;

    /* Find local processor ID and number of processors */
    int me, nprocs;
    me     = armci_msg_me();
    nprocs = armci_msg_nproc();

    /* Allocate pointers to data on all processors */
    a_ptr = (void**)malloc(nprocs*sizeof(int*));
    b_ptr = (void**)malloc(nprocs*sizeof(int*));
    map = (int*)malloc(nprocs*sizeof(int));

    /* Configure array dimensions. Force an unequal data distribution */
    dims[0]  = nprocs*TOTALELEMS + nprocs/2;
    if (me == 0) printf("Size of array: %d\n\n",dims[0]);
    /* Find first (zero-based) index of chunk owned by each processor and
       store it in map array */
    for (i=0; i<nprocs; i++) {
      map[i] = (int)(((double)i)*(((double)dims[0])/((double)nprocs)));
    }

    /* Figure out what size my portion of array is */
    if (me<nprocs-1) {
      nelem = map[me+1]-map[me];
    } else {
      nelem = dims[0]-map[me];
    }

    /* Allocate memory for array A */
    ierr = ARMCI_Malloc(a_ptr, nelem*sizeof(int));
    assert(ierr == 0);
    assert(a_ptr[me]);

    /* Allocate memory for array B */
    ierr = ARMCI_Malloc(b_ptr, nelem*sizeof(int));
    assert(ierr == 0);
    assert(b_ptr[me]);
    
    /* initialize data in array A and zero data in array B */
    a = (int*)a_ptr[me];
    b = (int*)b_ptr[me];
    for (i=0; i<nelem; i++) {
      a[i] = i + map[me] + 1;
      b[i] = 0;
    }

    /* Synchronize all processors to guarantee that everyone has data
       before proceeding to the next step. */
    armci_msg_barrier();

    /* Create local buffer for performing inversion */
    buf = (int*)malloc(nelem*sizeof(int));

    /* Copy inverted data into local buffer */
    a = (int*)a_ptr[me];
    for (i=0; i<nelem; i++) {
      buf[i] = a[nelem-i-1]; 
    }

    /* Find out which blocks of array B inverted block should be copied to.
       Start by finding min and max indices of data in array B*/
    min = dims[0] - (map[me] + nelem);
    max = dims[0] - map[me] - 1;

    /* Locate processors containing the endpoints */
    pmin = 0;
    for (i=0; i<nprocs; i++) {
      if (min >= map[i]) {
        pmin = i;
      } else {
        break;
      }
    }
    pmax = nprocs-1;
    for (i=nprocs-2; i>=0; i--) {
      if (max < map[i+1]) {
        pmax = i;
      } else {
        break;
      }
    }

    /* Loop over processors that will receive data and copy inverted data to
       processors */
    for (i=pmin; i<=pmax; i++) {
      /* Find min and max indices owned by processor i */
      lmin = map[i];
      if (i<nprocs-1) {
        lmax = map[i+1]-1;
      } else {
        lmax = dims[0]-1;
      }

      /* Find min and max indices that should be sent to processor i */
      if (lmin > min) {
        cmin = lmin;
      } else {
        cmin = min;
      }
      if (lmax < max) {
        cmax = lmax;
      } else {
        cmax = max;
      }

      /* Find offsets on source and destination processors */
      src_offset = cmin - min;
      src_ptr = (void*)(buf + src_offset);
      dst_offset = cmin - lmin;
      dst_ptr = ((char*)b_ptr[i]) + sizeof(int)*dst_offset;
      
      /* Find length of data (in bytes) to be sent to processor i */
      length = sizeof(int)*(cmax-cmin+1);

      /* Send data to processor */
      ARMCI_Put(src_ptr, dst_ptr, length, i);
    }
    ARMCI_AllFence();
    armci_msg_barrier();
    
    free(buf);

    VERIFY(b_ptr, dims, map);

    free(map);
    armci_msg_barrier();
    ARMCI_Free(a_ptr[me]);
    ARMCI_Free(b_ptr[me]);
    free(a_ptr);
    free(b_ptr);
}

int main(int argc, char **argv) {
    /* int heap=300000, stack=300000; */
    int me, nprocs;
    
    /* Step1: Initialize Message Passing library */
    armci_msg_init(&argc, &argv);

    /* Step2: Initialize ARMCI */
    ARMCI_Init();
    
    /* Step3: Initialize Memory Allocator (MA) */
    /*bjp
    if(! MA_init(C_DBL, stack, heap) ) ARMCI_Error("MA_init failed",stack+heap);
    */

    me     = armci_msg_me();
    nprocs = armci_msg_nproc();
    if(me==0) {
       printf("\nUsing %d processes\n\n", nprocs); fflush(stdout);
    }
    
       
    TRANSPOSE1D();
    
    if(me==0)printf("\nTerminating ..\n");
    ARMCI_Finalize();
    
    armci_msg_finalize();    
    return(0);
}
