#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "ga.h"
#include "macdecls.h"
#include "mp3.h"
#include "xgemm.h"

#if defined(FUJITSU) || defined(CRAY_YMP)
#   define THRESH 1.0e-10
#else
#   define THRESH 1.0e-20
#endif
#define ABS(x) ((x) >= 0.0 ? (x) : -(x))
#define MAX(x,y) ((x) >= (y) ? (x) : (y))
#define MISMATCH(x,y) (ABS((x)-(y)) / MAX(1.0,ABS((x)))) > THRESH
#define NMAX 100
#define POW2(x) ((x)*(x))
#define POW4(x) ((x)*(x)*(x)*(x))
#define SPECIFIC_CASE 0
#define VERBOSE 0

static void dpatch_test(
        int size,
        int dist_same,
        int ampos, int akpos,
        int bkpos, int bnpos,
        int cmpos, int cnpos);
static void dpatch_test2();

int main(int argc, char **argv)
{
    int bufsize, gasize;
    time_t t;
    long seed;
    int sizes[] = {10, 50, 100};
    int s, same_dist, ampos, akpos, bkpos, bnpos, cmpos, cnpos;

    MP_INIT(argc,argv);
    GA_Initialize_args(&argc,&argv);

    if (0 == GA_Nodeid()) {
        printf("  GA initialized\n");
        fflush(stdout);
    }

    /* we need srandom seed to be the same on all procs */
    t = time(NULL);
    seed = (long)t;
    GA_Lgop(&seed, 1, "max");
    srandom(seed);
#if VERBOSE
    printf("seed=%ld\n", seed);
    fflush(stdout);
#endif
    GA_Sync();

    /* we want to force distribution of innermost loop in nga_mulmat_patch
     by providing less buffer memory than needed */

    if(GA_Uses_ma()) {
        gasize = (POW4(NMAX) * 3)/GA_Nnodes();
    }
    else {
        gasize = 0;
    }

    bufsize = (NMAX/2 + 1)*(NMAX/3 + 1)*2 + POW2(NMAX/2 + 1);
    bufsize = bufsize*6/7;

    if (!MA_init(MT_DBL, 10, gasize+bufsize+500000)) {
        GA_Error("MA_init failed", -1);
    }
    if (0 == GA_Nodeid()) {
        printf(" \n");
        printf(" CHECKING MATRIX MULTIPLICATION FOR PATCHES \n");
#if VERBOSE
        printf("gasize and bufsize are %d %d\n", gasize, bufsize);
#endif
        printf(" \n");
        fflush(stdout);
    }
#if SPECIFIC_CASE
    dpatch_test(10, 0, 0, 1, 1, 2, 2, 3);
#else
    for (s=0; s<1; ++s) {
        for (same_dist=0; same_dist<2; ++same_dist) {
            for (ampos=0; ampos<4; ++ampos) {
                for (akpos=ampos+1; akpos<4; ++akpos) {
                    for (bkpos=0; bkpos<4; ++bkpos) {
                        for (bnpos=bkpos+1; bnpos<4; ++bnpos) {
                            for (cmpos=0; cmpos<4; ++cmpos) {
                                for (cnpos=cmpos+1; cnpos<4; ++cnpos) {
                                    dpatch_test(sizes[s],
                                            same_dist,
                                            ampos, akpos,
                                            bkpos, bnpos,
                                            cmpos, cnpos);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
#endif
#if 0
    dpatch_test2();
#endif

    if(0 == GA_Nodeid()) {
        printf(" All tests successful \n");
        fflush(stdout);
    }

    GA_Terminate();
    MP_FINALIZE();
    return 0;
}


/**
 * We start with a 4-dimensional array and multiply various 2D patches,
 * comparing results against a locally computed dgemm.
 */
static void dpatch_test(
        int size,
        int dist_same,
        int ampos, int akpos,
        int bkpos, int bnpos,
        int cmpos, int cnpos)
{
    const int ndim=4;
    double *a, *b, *c, *r, *v;
    double alpha, beta;
    int nproc, me;
    int i, j;
    int dims[4], chunk[4], rld[3];
    int alo[4], ahi[4], ald[3];
    int blo[4], bhi[4], bld[3];
    int clo[4], chi[4], cld[3];
    int g_a, g_b, g_c;
    char ta, tb;
    int m, n, k;

    me = GA_Nodeid();
    nproc = GA_Nnodes();

    assert(size <= NMAX);

#if SPECIFIC_CASE
    m = 1;
    n = 1;
    k = 4;
#else
    m = random() % (size/2);
    n = random() % (size/2);
    k = random() % (size/2);
    if (m<=0) m = 1;
    if (n<=0) n = 1;
    if (k<=0) k = 1;
#endif
    if (0 == me) {
        printf("size=%d, dist_same=%d ampos=%d akpos=%d bkpos=%d bnpos=%d cmpos=%d cnpos=%d m=%d n=%d k=%d\n", size, dist_same, ampos, akpos, bkpos, bnpos, cmpos, cnpos, m, n, k);
        fflush(stdout);
    }

    a = malloc(sizeof(double)*size*size);
    b = malloc(sizeof(double)*size*size);
    c = malloc(sizeof(double)*size*size);
    r = malloc(sizeof(double)*size*size);
    memset(a, 0, sizeof(double)*size*size);
    memset(b, 0, sizeof(double)*size*size);
    memset(c, 0, sizeof(double)*size*size);
    memset(r, 0, sizeof(double)*size*size);

    /* establish the shape and default chunking of the global arrays */
    for (i=0; i<ndim; ++i) {
        alo[i] = 0;
        ahi[i] = size-1;
        dims[i] = size;
        chunk[i] = -1;
    }

    /* create g_a */
    g_a = NGA_Create(C_DBL, ndim, dims, "a", chunk);
    if (0 == g_a) {
        printf("NGA_Create failed\n");
        fflush(stdout);
        GA_Error("... exiting", 1);
    }

    /* create g_b and g_c */
    if (dist_same) {
        g_b = GA_Duplicate(g_a, "a_duplicated");
        if(1 == GA_Compare_distr(g_a, g_b)) {
            GA_Error("g_b distribution different",1);
        }
        g_c = GA_Duplicate(g_a, "a_duplicated_again");
        if(1 == GA_Compare_distr(g_a, g_c)) {
            GA_Error("g_c distribution different",1);
        }
    }
    else { 
        chunk[ndim-1] = size;
        g_b = NGA_Create(C_DBL, ndim, dims, "b1", chunk);
        if (0 == g_b) {
            GA_Error("NGA_Create failed:b1",1);
        }
        chunk[ndim-1] = 0;
        chunk[ndim-2] = size;
        g_c = NGA_Create(C_DBL, ndim, dims, "c1", chunk);
        if (0 == g_c) {
            GA_Error("NGA_Create failed:c1",1);
        }
    }

#if 0
    if (0 == me) {
        printf("\n");
        printf("> Checking NGA_Matmul_patch ... \n");
        fflush(stdout);
    }
#endif

    /* fill the g_a and g_b global arrays entirely with data */
    for (i=0; i<ndim; ++i) {
        alo[i] = 0;
        ahi[i] = size-1;
        blo[i] = 0;
        bhi[i] = size-1;
        if (i<ndim-1) {
            ald[i] = size;
            bld[i] = size;
        }
    }
    if (0 == me) {
        /* generate some local data (an enumerated range) */
        double *v = malloc(sizeof(double)*POW4(size));
        memset(v, 0, sizeof(double)*POW4(size));
        for (i=0; i<POW4(size); ++i) {
            v[i] = i;
        }
        NGA_Put(g_a,alo,ahi,v,ald);
        NGA_Put(g_b,blo,bhi,v,bld);
        free(v);
    }
    GA_Zero(g_c);
    GA_Sync();

    /* for g_a, g_b, g_c generate a random starting index for patches */
    for (i=0; i<ndim; ++i) {
        ahi[i] = alo[i] = random() % (size/2);
        bhi[i] = blo[i] = random() % (size/2);
        chi[i] = clo[i] = random() % (size/2);
        if (i<ndim-1) {
            ald[i] = 1;
            bld[i] = 1;
            cld[i] = 1;
        }
    }
    ahi[ampos] += m - 1;
    ahi[akpos] += k - 1;
    if (ampos>0) ald[ampos-1] = m;
    if (akpos>0) ald[akpos-1] = k;
    bhi[bkpos] += k - 1;
    bhi[bnpos] += n - 1;
    if (bkpos>0) bld[bkpos-1] = k;
    if (bnpos>0) bld[bnpos-1] = n;
    chi[cmpos] +=  m - 1;
    chi[cnpos] +=  n - 1;
    if (cmpos>0) cld[cmpos-1] = m;
    if (cnpos>0) cld[cnpos-1] = n;

#if 0
    printf("a[%d:%d,%d:%d,%d:%d,%d:%d] %dx%dx%dx%d (%dx%d)\n",
            alo[0], ahi[0],
            alo[1], ahi[1],
            alo[2], ahi[2],
            alo[3], ahi[3],
            m, 1, k, 1, m, k);
    printf("ald={%d,%d,%d}\n", ald[0], ald[1], ald[2]);
    printf("b[%d:%d,%d:%d,%d:%d,%d:%d] %dx%dx%dx%d (%dx%d)\n",
            blo[0], bhi[0],
            blo[1], bhi[1],
            blo[2], bhi[2],
            blo[3], bhi[3],
            1, k, 1, n, k, n);
    printf("bld={%d,%d,%d}\n", bld[0], bld[1], bld[2]);
    printf("c[%d:%d,%d:%d,%d:%d,%d:%d] %dx%dx%dx%d (%dx%d)\n",
            clo[0], chi[0],
            clo[1], chi[1],
            clo[2], chi[2],
            clo[3], chi[3],
            m, 1, 1, n, m, n);
    printf("cld={%d,%d,%d}\n", cld[0], cld[1], cld[2]);
#endif

    /* reset our buffers, just in case */
    memset(a, 0, sizeof(double)*size*size);
    memset(b, 0, sizeof(double)*size*size);
    memset(c, 0, sizeof(double)*size*size);
    memset(r, 0, sizeof(double)*size*size);

    /* get patches locally and compute locally */
    NGA_Get(g_a, alo, ahi, a, ald);
    NGA_Get(g_b, blo, bhi, b, bld);
    GA_Sync();
    ta = 'n';
    tb = 'n';
    alpha = 1e0;
    beta = 0e0;
    xb_dgemm(&tb, &ta, &n, &m, &k, &alpha, b, &n, a, &k, &beta, c, &n);

    /* perform global computation */
    NGA_Matmul_patch(ta, tb, &alpha, &beta, 
            g_a, alo, ahi,
            g_b, blo, bhi,
            g_c, clo, chi);
    GA_Sync();

    /* get global result into local buf and compare results */
    NGA_Get(g_c, clo, chi, r, cld);
    GA_Sync();
    for (i=0; i<1; ++i) {
        if (MISMATCH(c[i],r[i])) {
            printf("at %d %f != %f\n", i, c[i], r[i]);
            GA_Error("mismatch", 1);
        }
    }

    free(a);
    free(b);
    free(c);
    free(r);
    GA_Destroy(g_a);
    GA_Destroy(g_b);
    GA_Destroy(g_c);
}


static void dpatch_test2()
{
}
