/**
 * testmultrect.c
 * Test rectangular matrix multiplication.
 *
 * Nawab Ali <nawab.ali@pnl.gov>
 */
#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "ga.h"
#include "macdecls.h"
#include "mp3.h"

#define NDIMS 2                 /* 2-dimensional matrices */

#define HEAP  4000000
#define STACK 4000000

int m = 2;
int k = 3;
int n = 4;
int g_a, g_b, g_c;              /* GA handles for matrices A, B, and C */

int mr_create_matrices(void);
int mr_cleanup(void);


int main(int argc, char **argv) {
    int ret, nprocs;
    int me, heap, stack;

    /* Initialize message passing and GA */
    MP_INIT(argc,argv);
    GA_INIT(argc,argv);

    me = GA_Nodeid();
    nprocs = GA_Nnodes();
    heap = HEAP/nprocs;
    stack = STACK/nprocs;

    if (!MA_init(C_DBL, stack, heap)) {
        GA_Error("MA_init failed", stack+heap);
    }

    /* Create the matrices A, B, C */
    ret = mr_create_matrices();

    GA_Print(g_a);
    GA_Print(g_b);
    GA_Print(g_c);

    /* C = A x B */
    GA_Dgemm('N', 'N', m, n, k, 1.0, g_a, g_b, 0.0, g_c);

    /* Clean up */
    ret = mr_cleanup();
    GA_Terminate();
    MP_FINALIZE();

    return 0;
}

int mr_create_matrices(void) {
    double val_a = 10.0;
    double val_b = 20.0;
    int ret, dims[NDIMS];

    /* Create a 2-dimensional matrix A[m][k] */
    dims[0] = m;
    dims[1] = k;

    g_a = GA_Create_handle();
    GA_Set_array_name(g_a, "Matrix A");
    GA_Set_data(g_a, NDIMS, dims, C_DBL);
    ret = GA_Allocate(g_a);
    GA_Fill(g_a, &val_a);

    /* Create a 2-dimensional matrix B[k][n] */
    dims[0] = k;
    dims[1] = n;

    g_b = GA_Create_handle();
    GA_Set_array_name(g_b, "Matrix B");
    GA_Set_data(g_b, NDIMS, dims, C_DBL);
    ret = GA_Allocate(g_b);
    GA_Fill(g_b, &val_b);

    /* Create a 2-dimensional matrix C[m][n] */
    dims[0] = m;
    dims[1] = n;

    g_c = GA_Create_handle();
    GA_Set_array_name(g_c, "Matrix C");
    GA_Set_data(g_c, NDIMS, dims, C_DBL);
    ret = GA_Allocate(g_c);
    GA_Zero(g_c);

    return 0;
}

int mr_cleanup(void) {
    GA_Destroy(g_a);
    GA_Destroy(g_b);
    GA_Destroy(g_c);
    return 0;
}
