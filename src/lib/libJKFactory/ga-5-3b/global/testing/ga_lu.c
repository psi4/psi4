#if HAVE_CONFIG_H
#   include "config.h"
#endif


#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif

#include "ga.h"
#include "macdecls.h"
#include "mp3.h"
#include "galinalg.h"

#define BLOCK_CYCLIC
#define BLOCK_SIZE 500
/*#define USE_SCALAPACK_DISTR*/

#define DEBUG 0

static int nprocs, me;

static void init_array(double *a, int n) 
{
    int i, j;
    double max=0.0; /* for pivoting */
    
    for(i=0; i<n; i++) 
    {
       for(j=0; j<n; j++) 
       {
          a[i*n+j] = (double)rand()/RAND_MAX; /* (i*n+j+1); */
          if(max<a[i*n+j]) max = a[i*n+j];
       }
    }

    /* Forcing the diagonal is always max to avoid pivoting*/
    for(i=0; i<n; i++) 
    {
       a[i*n+i] = max + (double)i;
    }
}

void copy_array(double *a, double *b, int n) 
{
    int i, j;
    
    for(i=0; i<n; i++) 
    {
       for(j=0; j<n; j++) 
       {
          b[i*n+j] = a[i*n+j]; 
       }
    }
}

void print_array(double *a, int n) 
{
    int i, j;
    
    printf("Print Matrix of size %d x %d:\n", n, n);
    for(i=0; i<n; i++) 
    {
       for(j=0; j<n; j++) 
       {
          if(a[i*n+j]<0) printf("%2.4e ",  a[i*n+j]);
          else           printf(" %2.4e ", a[i*n+j]);
       }
       printf("\n");
    }
    printf("\n");
}

/* print U, i.e. upper triangular */
void print_U(double *a, int n) 
{
    int i, j;
    
    printf("Print Upper triangular matrix of size %d x %d:\n", n, n);
    for(i=0; i<n; i++) 
    {
       for(j=0; j<i; j++) 
       {
          printf(" %2.4e ", 0.0);
       }
       
       for(j=i; j<n; j++) 
       {
          if(a[i*n+j]<0) printf("%2.4e ",  a[i*n+j]);
          else           printf(" %2.4e ", a[i*n+j]);
       }
       
       printf("\n");
    }
    printf("\n");
}

/* print L, i.e. unit lower triangular */
void print_L(double *a, int n) 
{
    int i, j;
    
    printf("Print Unit Lower triangular matrix of size %d x %d:\n", n, n);
    for(i=0; i<n; i++) 
    {
       for(j=0; j<i; j++) 
       {
          if(a[i*n+j]<0) printf("%2.4e ",  a[i*n+j]);
          else           printf(" %2.4e ", a[i*n+j]);
       }
       for(j=i; j<n; j++) 
       {
          if(i==j) printf(" %2.4e ", 1.0);
          else     printf(" %2.4e ", 0.0);
       }

       printf("\n");
    }
    printf("\n");
}

void lu_basic(double *a, int n) {

    int i, j, k;
    
    for (i=0; i<n-1; i++) 
    {
       /* scale lower matrix column */
       for (j=i+1; j<n; j++) 
       {
          a[j*n+i] = a[j*n+i]/a[i*n+i]; /* a[j][i] = a[j][i]/a[i][i] */
       }
       
       /* rank-one update of matrix */
       for (j=i+1; j<n; j++) 
       {
          for (k=i+1; k<n; k++) 
          {
             /*a[k][j] = a[k][j] - a[i][j]/a[k][i] */
             a[k*n+j]  = a[k*n+j] - a[i*n+j]/a[k*n+i];
          }   
       }   
    }
}

int lu_lapack(double *a, int n) 
{
    int i, j;
    double *aa=NULL;
    BlasInt *ipiv=NULL, info;
    BlasInt ld=(BlasInt)n;
    BlasInt N=(BlasInt)n;
    
    aa = (double*)malloc(n*n*sizeof(double));
    ipiv = (BlasInt*)malloc(n*sizeof(BlasInt));

    /* row-major to column-major (NOTE: dgetrf_ is a fortran function) */
    for(i=0; i<n; i++) 
    {
       for(j=0; j<n; j++) 
       {
          aa[j*n+i] = a[i*n+j];
       }
    }
    
    LAPACK_DGETRF(&N, &N, aa, &ld, ipiv, &info); /* LAPACK's LU */
    
    /* column-major to row-major */
    for(i=0; i<n; i++) 
    {
       for(j=0; j<n; j++) 
       {
          a[i*n+j] = aa[j*n+i];
       }
    }
    
#if DEBUG
    printf("ipiv[] = [ ");
    for(i=0; i<n; i++) printf("%ld ", ipiv[i]);
    printf("]\n");
#endif
    
    free(aa);
    return (int)info;
}

void lu(double *A, int matrix_size) 
{
    double *a=NULL, *a_verify=NULL;
    int info;
    
    a        = (double*)malloc(matrix_size*matrix_size*sizeof(double));
    a_verify = (double*)malloc(matrix_size*matrix_size*sizeof(double));
    if(a == NULL || a_verify == NULL) 
    {
       GA_Error("lu(): malloc failed", matrix_size*matrix_size*sizeof(double));
    }
    
    copy_array(A, a, matrix_size);
    copy_array(A, a_verify, matrix_size);
    
#if 0
    printf("\nDoing LU Factorization\n\n");
    lu_basic(a, matrix_size);
    printf("LU = \n");
    print_array(a, matrix_size);
    
    printf("\nDoing LAPACK's LU Factorization\n\n");
#endif
    info = lu_lapack(a_verify, matrix_size);

#if DEBUG
    printf("LU = ");
    print_array(a_verify, matrix_size);
#endif
    
    if(info!=0) 
    {
       printf("\nError: dgetrf() of lapack is NOT successful (INFO=%d)\n\n",
              info);
       printf("NOTE:\n INFO=0:  successful exit\n INFO<0:  if INFO = -i, the i-th argument had an illegal value\n INFO>0:  if INFO = i, U(i,i) is exactly zero. The factorization\n has been completed, but the factor U is exactly singular,  and\n division by zero will occur if it is used to solve a system of\n equations.\n");
       exit(0);
    }
    

    free(a_verify);
    free(a);

    
}

void dtrsm_lapack(double *a, double *b, int n,
                  char side, char uplo, char transa, char diag) 
{
    int i, j;
    double *aa=NULL;
    double *bb=NULL;
    BlasInt ld = (BlasInt)n;
    BlasInt N  = (BlasInt)n;

    aa = (double*)malloc(n*n*sizeof(double));
    bb = (double*)malloc(n*n*sizeof(double));
    
    /* row-major to column-major (NOTE: dgetrf_ is a fortran function) */
    for(i=0; i<n; i++) 
    {
       for(j=0; j<n; j++) 
       {
          aa[j*n+i] = a[i*n+j];
          bb[j*n+i] = b[i*n+j];
       }
    }
    
    {
       double alpha= 1.0;
       
       LAPACK_DTRSM(&side, &uplo, &transa, &diag, &N, &N, &alpha, aa, &ld, bb, &ld);
    }   
    
    /* column-major to row-major */
    for(i=0; i<n; i++) 
    {
       for(j=0; j<n; j++) 
       {
          a[i*n+j] = aa[j*n+i];
          b[i*n+j] = bb[j*n+i];
       }
    }
    
    free(aa);
    free(bb);    
}

void transpose(double *a, int n) 
{
    int i, j;
    double *aa=NULL;
    
    aa = (double*)malloc(n*n*sizeof(double));
    copy_array(a, aa, n);
    
    /* row-major to column-major (NOTE: dgetrf_ is a fortran function) */
    for(i=0; i<n; i++) 
    {
       for(j=0; j<n; j++) 
       {
          a[j*n+i] = aa[i*n+j];
       }
    }
}


/* input is matrix size */
void ga_lu(double *A, int matrix_size) 
{
    int g_a, g_b, dims[2], type=C_DBL;
    int lo[2], hi[2], ld;
    int block_size[2];
#ifdef USE_SCALAPACK_DISTR
    int proc_grid[2];
#endif
    double time, gflops;
    
    /* create a 2-d GA (global matrix) */
    dims[0] = matrix_size;
    dims[1] = matrix_size;
    block_size[0] = BLOCK_SIZE;
    block_size[1] = BLOCK_SIZE;
#ifdef USE_SCALAPACK_DISTR
    int proc_grid[2];
    proc_grid[0] = 2;
    proc_grid[1] = nprocs/2;
    if(nprocs%2) GA_Error("For ScaLAPACK stle distribution, nprocs must be "
                         " divisible by 2", 0);
#endif
    
    
#ifndef BLOCK_CYCLIC
    g_a = NGA_Create(type, 2, dims, "A", NULL);
    g_b = GA_Duplicate(g_a, "transposed array B");
#else
    g_a = GA_Create_handle();
    GA_Set_data(g_a, 2, dims, type);
    GA_Set_array_name(g_a,"A");
#  ifdef USE_SCALAPACK_DISTR
    GA_Set_block_cyclic_proc_grid(g_a, block_size, proc_grid);
#  else
    GA_Set_block_cyclic(g_a, block_size);    
#  endif
    GA_Allocate(g_a);
    
    g_b = GA_Create_handle();
    GA_Set_data(g_b, 2, dims, type);
    GA_Set_array_name(g_b,"B");
#  ifdef USE_SCALAPACK_DISTR
    GA_Set_block_cyclic_proc_grid(g_b, block_size, proc_grid);
#  else
    GA_Set_block_cyclic(g_b, block_size);
#  endif
    GA_Allocate(g_b);
    
#endif
    
    /* copy the local matrix into GA */
    if(me==0) 
    {
       lo[0] = 0;
       hi[0] = matrix_size - 1;
       lo[1] = 0;
       hi[1] = matrix_size - 1;
       ld    = matrix_size;
       
       NGA_Put(g_a, lo, hi, A, &ld);
    }
    GA_Sync();

    GA_Transpose(g_a, g_b);
    time = MP_TIMER();
    /* The following function does not exist. Not sure what to replace it
     * with. GA_Lu_solve(char trans, int g_a, int g_b) requiresa an
     * additioanl GA. */
    /* GA_Lu('n', g_b); */
    time = MP_TIMER() - time;

    /* 2/3 N^3 - 1/2 N^2 flops for LU and 2*N^2 for solver */
    gflops = ( (((double)matrix_size) * matrix_size)/(time*1.0e+9) *
               (2.0/3.0 * (double)matrix_size - 0.5) );
    if(me==0) printf("\nGA_Lu: N=%d flops=%2.5e Gflops, time=%2.5e secs\n\n",
                     matrix_size, gflops, time);

#if DEBUG
    GA_Print(g_a);
    GA_Print(g_b);
#endif
    /* if(me==0) lu(A, matrix_size);     */

    GA_Destroy(g_a);
    GA_Destroy(g_b);
}


int main(int argc, char **argv) {
    int i, matrix_size;
    int heap=20000000, stack=200000000;
    double *A=NULL;
    
    if(argc == 1)
    {
      matrix_size = 1234;
    }
    else if (argc == 2)
    {
      matrix_size = atoi(argv[1]);
    }
    else
    {
       printf("Usage Error\n\t Usage: <program> [matrix_size]\n");
       exit(0);
    }

    if(matrix_size <= 0) 
    {
       printf("Error: matrix size (%d) should be > 0\n",
              matrix_size);
       GA_Error("matrix size should be >0", 1);
    }
    
    /* *****************************************************************
     * Initialize MPI/TCGMSG-MPI, GA and MA
     * *****************************************************************/
    MP_INIT(argc,argv);
    
    GA_INIT(argc,argv);        /* initialize GA */

    me     = GA_Nodeid();
    nprocs = GA_Nnodes();

    heap /= nprocs;
    stack /= nprocs;

    if(! MA_init(MT_F_DBL, stack, heap)) /* initialize MA */
    {
       GA_Error("MA_init failed",stack+heap);
    }

    /* create/initialize the matrix */    
    if((A = (double*)malloc(matrix_size*matrix_size*sizeof(double))) == NULL) 
    {
       GA_Error("malloc failed", matrix_size*matrix_size*sizeof(double));
    }

    for(i=0; i<2; i++) /* 5 runs */
    {
       init_array(A, matrix_size);
#if DEBUG
       if(me==0) print_array(A, matrix_size);
#endif
       /* *****************************************************************
        * Perform LU Factorization
        * *****************************************************************/
       ga_lu(A, matrix_size);
    }
    
    free(A);
    
    /* *****************************************************************
     * Terminate MPI/TCGMSG-MPI, GA and MA
     * *****************************************************************/
    if(me==0)printf("Success\n");
    GA_Terminate();

    MP_FINALIZE();
    
    return 0;
}

/**
 * TODO:
 *   - LU for non-square matrix
 *
 */
