/*
 * As reported by yanait@gmail.com, large (16GB) arrays caused seeg fault on
 * 2 nodes 2 cores per node.
 *
 * The following program (just doing ga_acc to 16GB GA (Array
 * Dimensions:8192x262144) did not work, resulting in segmentation fault. 
 *
 * % mpirun -npernode 2 -np 4 ./a.out 
 *
 * Process=0 (node#0)       owns array section: [0:4095,0:131071]
 * Process=1 (node#0)        owns array section: [0:4095,131072:262143]
 * Process=2 (node#1)        owns array section: [4096:8191,0:131071]
 * Process=3 (node#1)        owns array section: [4096:8191,131072:262143]
 *
 * I found that a couple of subroutines in armci fail to handle large-size
 * arrays. 
 * My debugging showed that the variables of int type overflow to calculate
 * address or pointer (over 2GB). 
 * I made a patch file (see below).
 * This is related to my previous post regarding the kr_shmem_free error in
 * ga_matmul (also using ga_acc).
 */
#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>

#include "mp3.h"
#include "ga.h"
#include "macdecls.h"

/* create N1 x N2 matrix (double precision) */
#define N1 1024*8
#define N2 1024*256

int main(int argc, char **argv)
{
    int me, nproc, g_a, i, j;

    int ndim=2, type=MT_F_DBL, dims[2]={N1,N2};

    double *buf;

    int lo[2], hi[2], ld[1];
    double alpha = 1.0;

    MP_INIT(argc,argv);
    GA_Initialize_ltd(-1);

    me=GA_Nodeid();
    nproc=GA_Nnodes();

    if(me==0) printf("Using %ld processes\n",(long)nproc);
    if(me==0) printf("memory = %ld bytes\n",((long)N1)*((long)N2)*8);

    g_a = NGA_Create(type, ndim, dims, "A", NULL);

    GA_Zero(g_a);   /* zero the matrix */

    GA_Print_distribution(g_a);

    if(me == 0) {
        buf = (double*)(malloc(N1*1024*sizeof(double)));
        for(j = 0; j < N1*1024; ++j) buf[j] = 1.0;
        for(i = 0; i < N2/1024; ++i) {
            lo[0] = 0;
            hi[0] = lo[0] + N1   -1;
            lo[1] = i*1024;
            hi[1] = lo[1] + 1024 -1;
            ld[0] = 1024;
            printf("NGA_Acc.%d:  %d:%d %d:%d\n",i,lo[0],hi[0],lo[1],hi[1]);
            NGA_Init_fence();
            NGA_Acc(g_a, lo, hi, buf, ld, &alpha);
            NGA_Fence();
        }
    }

    GA_Sync();

    GA_Destroy(g_a);

    GA_Terminate();
    MP_FINALIZE();

    return 0;
}
