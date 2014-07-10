#if HAVE_CONFIG_H
#   include "config.h"
#endif

/**************************************************
 *             LU factorization                   *
 *             Armci Version                      *
 *             Block distribution                 *
 *             Multi-threaded                     *
 **************************************************/
#define DEBUG
#define DEBUG1_
#define DEBUG2_
#define DEBUG3_
#define USE_MUTEX_

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_MATH_H
#   include <math.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_STDARG_H
#   include <stdarg.h>
#endif

#include "armci.h"
#include "utils.h"
#include "message.h"

#define MAXRAND    32767.0
#define DEFAULT_N      8
#define DEFAULT_B      2

#define MAX_THREADS 8

/* global variables */
int n = DEFAULT_N;         /* The size of the matrix */
int block_size = DEFAULT_B;/* Block dimension */
int nblocks;               /* Number of blocks in each dimension */
int num_rows;              /* Number of processors per row of processor grid */
int num_cols;              /* Number of processors per col of processor grid */
double **a;                /* a = lu; l and u both placed back in a */
int nproc, th_per_p = 1, nthreads, me = 0, me_th[MAX_THREADS];
thread_t threads[MAX_THREADS];
int proc_bytes, thread_doubles[MAX_THREADS];
int num;

int doprint = 1;
int d = 0; /* delay */
thread_lock_t mutex;
FILE *rep[MAX_THREADS];
char fname[] = "threadXX.log";
void report(int th_idx, char *fmt, ...)
{
#ifdef DEBUG3
    va_list ap;
    va_start(ap, fmt);
    vfprintf(rep[th_idx], fmt, ap);
    va_end(ap);
#endif
}



/* function declaration */
void *lu(void *);
void lu0(double *,int, int);
void bdiv(double *, double *, int, int, int, int);
void bmodd(double *, double*, int, int, int, int);
void bmod(double *, double *, double *, int, int, int, int, int, int);
void daxpy(double *, double *, int, double);
int block_owner(int, int);
void init_array();
double touch_array(int, int);
void print_array(int);
void print_block_dbg(double *, const char *, int, int, int);

void prefetch(double **A, double *buf, int I, int J, int th_idx, armci_hdl_t **hdlp);
void get_remote(double *buf, int I, int J, armci_hdl_t **hdlp);
int next_block(int th_idx, int bs, int kl, int ci, int cj, int cI, int cJ, int K, int *pI, int *pJ);

/* timing functions */
extern void start_timer(void);
extern double elapsed_time(void);
extern double stop_time(void);

main(int argc, char *argv[])
{
    int i, j, l;
    int ch;
    extern char *optarg;
    int edge;
    int size;
    int lu_arg[MAX_THREADS][3];

    /* ARMCI */
    void **ptr;
    double **ptr_loc;

    THREAD_LOCK_INIT(mutex);

    armci_msg_init(&argc,&argv);
    nproc = armci_msg_nproc();
    me = armci_msg_me();

    while ((ch = getopt(argc, argv, "n:b:p:t:d:h")) != -1) {
        switch(ch) {
            case 'n': n = atoi(optarg); break;
            case 'b': block_size = atoi(optarg); break;
            case 'p': nproc = atoi(optarg); break;
            case 't': th_per_p = atoi(optarg); break;
            case 'd': d = atoi(optarg); break;
            case 'h': {
                printf("Usage: LU, or \n");
        printf("       LU -nMATRIXSIZE -bBLOCKSIZE -pNPROC -tTH_PER_P\n");
                armci_msg_barrier();
                armci_msg_finalize();
                exit(0);
            }
        }
    }

    if(th_per_p>MAX_THREADS) {
        th_per_p=MAX_THREADS;
        if(me==0)printf("Warning: cannot run more than %d threads, adjust MAX_THREADS",MAX_THREADS);
    }

    if (d) {
        fprintf(stderr, "%d: %d\n", me, getpid());
        sleep(d);
    }

    nthreads = th_per_p * nproc;
    if(me == 0) {
        printf("\n Blocked Dense LU Factorization\n");
        printf("     %d by %d Matrix\n", n, n);
        printf("     %d Processors\n", nproc);
        printf("     %d thread(s) per processor, %d threads total\n", th_per_p, nthreads);
        printf("     %d by %d Element Blocks\n", block_size, block_size);
        printf("\n");
    }

    num_rows = (int) sqrt((double) nthreads);
    for (;;) {
        num_cols = nthreads/num_rows;
        if (num_rows*num_cols == nthreads)
            break;
        num_rows--;
    }

    nblocks = n/block_size;
    if (block_size * nblocks != n) {
        nblocks++;
    }

    num = (nblocks * nblocks)/nthreads;
    if((num * nthreads) != (nblocks * nblocks))
        num++;

    edge = n%block_size;
    if (edge == 0) {
        edge = block_size;
    }
#ifdef DEBUG
    if(me == 0)
        for (i=0;i<nblocks;i++) {
            for (j=0;j<nblocks;j++)
                printf("%d ", block_owner(i, j));
            printf("\n");
        }
    armci_msg_barrier();
/*    armci_msg_finalize(); */
/*    exit(0); */
#endif

    for (l = 0; l < th_per_p; l++) {
        me_th[l] = me * th_per_p + l;
        for (i=0;i<nblocks;i++) {
            for (j=0;j<nblocks;j++) {
                if(block_owner(i,j) == me_th[l]) {
                    if ((i == nblocks-1) && (j == nblocks-1)) {
                        size = edge*edge;
                    }
                    else if ((i == nblocks-1) || (j == nblocks-1)) {
                        size = edge*block_size;
                    }
                    else {
                        size = block_size*block_size;
                    }
                    thread_doubles[l] += size;
                }
            }
        }
        proc_bytes += thread_doubles[l] * sizeof(double);
    }

    /* initialize ARMCI */
    ARMCI_Init();
    ptr = (void **)malloc(nproc * sizeof(void *));
    ARMCI_Malloc(ptr, proc_bytes);

    a = (double **)malloc(nblocks*nblocks*sizeof(double *));
    if (a == NULL) {
        fprintf(stderr, "Could not malloc memory for a\n");
        exit(-1);
    }
    ptr_loc = (double **)malloc(nthreads*sizeof(double *));
    for (i = 0; i < nproc; i++) {
        ptr_loc[i * th_per_p] = (double *)ptr[i];
        for (j = 1; j < th_per_p; j++)
            ptr_loc[i * th_per_p + j] = ptr_loc[i * th_per_p + j - 1] + thread_doubles[j - 1];
    }
    for(i=0; i<nblocks;i ++) {
        for(j=0; j<nblocks; j++) {
            a[i+j*nblocks] = ptr_loc[block_owner(i, j)];
            if ((i == nblocks-1) && (j == nblocks-1)) {
                size = edge*edge;
            } else if ((i == nblocks-1) || (j == nblocks-1)) {
                size = edge*block_size;
            } else {
                size = block_size*block_size;
            }
            ptr_loc[block_owner(i, j)] += size;
        }
    }
#if 0
    for(i=0; i<nblocks*nblocks;i ++) printf("%d: a[%d]=%p\n", me, i, a[i]);
    fflush(stdout);
#endif

    /* initialize the array */
    init_array();

    /* barrier to ensure all initialization is done */
    armci_msg_barrier();

    /* to remove cold-start misses, all processors touch their own data */
/*    for (l = 0; l < th_per_p; l++) touch_array(block_size, me_th[l]); */
    armci_msg_barrier();

    if(doprint) {
        if(me == 0) {
            printf("Matrix before LU decomposition\n");
            print_array(me);
        }
        armci_msg_barrier();
    }

#if 1
    for (i = 0; i < nblocks; i++)
        for (j = 0; j < nblocks; j++)
            print_block_dbg(a[i + j * nblocks], "proc %d, a[%d, %d]:\n", me, i, j);
#endif

    TH_INIT(nproc,th_per_p);

    /* Starting the timer */
    if(me == 0) start_timer();

    for (l = 0; l < th_per_p; l++) {
#ifdef DEBUG3
        fname[6] = '0' + me_th[l] / 10;
        fname[7] = '0' + me_th[l] % 10;
        if (!(rep[l] = fopen(fname, "w"))) { perror(fname); abort(); }
        report(l, "started report on process %d thread %d(%d)\n", me, l, me_th[l]);
#endif
        lu_arg[l][0] = n;
        lu_arg[l][1] = block_size;
        lu_arg[l][2] = l;
        THREAD_CREATE(threads + l, lu, lu_arg[l]);
    }

    for (l = 0; l < th_per_p; l++) {
        THREAD_JOIN(threads[l], NULL);
#ifdef DEBUG3
        fclose(rep[l]);
#endif
    }
    armci_msg_barrier();

    /* Timer Stops here */
    if(me == 0)
        printf("\nRunning time = %lf milliseconds.\n\n",  elapsed_time());

    if(doprint) {
        if(me == 0) {
            printf("after LU\n");
            print_array(me);
        }
        armci_msg_barrier();
    }

    /* done */
    ARMCI_Free(ptr[me]);
    ARMCI_Finalize();
    armci_msg_finalize();

    THREAD_LOCK_DESTROY(mutex);
}

void lu0(double *a, int n, int stride)
{
    int j;
    int k;
    int length;
    double alpha;

    for (k=0; k<n; k++) {
        /* modify subsequent columns */
        for (j=k+1; j<n; j++) {
            a[k+j*stride] /= a[k+k*stride];
            alpha = -a[k+j*stride];
            length = n-k-1;
            daxpy(&a[k+1+j*stride], &a[k+1+k*stride], n-k-1, alpha);
        }
    }
}

void bdiv(double *a, double *diag, int stride_a, int stride_diag,
          int dimi, int dimk)
{
    int j;
    int k;
    double alpha;

    for (k=0; k<dimk; k++) {
        for (j=k+1; j<dimk; j++) {
            alpha = -diag[k+j*stride_diag];
            daxpy(&a[j*stride_a], &a[k*stride_a], dimi, alpha);
        }
    }
}

void bmodd(double *a, double *c, int dimi, int dimj,
           int stride_a, int stride_c)
{
    int i;
    int j;
    int k;
    int length;
    double alpha;

    for (k=0; k<dimi; k++) {
        for (j=0; j<dimj; j++) {
            c[k+j*stride_c] /= a[k+k*stride_a];
            alpha = -c[k+j*stride_c];
            length = dimi - k - 1;
            daxpy(&c[k+1+j*stride_c], &a[k+1+k*stride_a], dimi-k-1, alpha);
        }
    }
}

void bmod(double *a, double *b, double *c, int dimi, int dimj, int dimk,
          int stridea, int strideb, int stridec)
{
    int i;
    int j;
    int k;
    double alpha;

    for (k=0; k<dimk; k++) {
        for (j=0; j<dimj; j++) {
            alpha = -b[k+j*strideb];
            daxpy(&c[j*stridec], &a[k*stridea], dimi, alpha);
        }
    }
}

void daxpy(double *a, double *b, int n, double alpha)
{
    int i;
    for (i=0; i<n; i++) {
        a[i] += alpha*b[i];
    }
}

int block_owner(int I, int J)
{
/*    return((J%num_cols) + (I%num_rows)*num_cols); */
    return((I+J*nblocks)/num);
}

void init_array()
{
    int i, j, l;
    int ii, jj;
    int edge;
    int ibs;
    int jbs, skip;

    srand48((long) 1);
    edge = n%block_size;
    for (l = 0; l < th_per_p; l++) {
        for (j=0; j<n; j++) {
            for (i=0; i<n; i++) {
                if(block_owner((i/block_size), (j/block_size)) == me_th[l]) {
                    if ((n - i) <= edge) {
                        ibs = edge;
                        ibs = n-edge;
                        skip = edge;
                    } else {
                        ibs = block_size;
                        skip = block_size;
                    }
                    if ((n - j) <= edge) {
                        jbs = edge;
                        jbs = n-edge;
                    } else {
                        jbs = block_size;
                    }
                    ii = (i/block_size) + (j/block_size)*nblocks;
                    jj = (i%ibs)+(j%jbs)*skip;
                    a[ii][jj] = i + j*6 + 1;
                    if (i == j) {
                        a[ii][jj] *= 10;
                    }
                }
            }
        }
    }
}

double touch_array(int bs, int me)
{
    int i, j, I, J, l;
    double tot = 0.0;
    int ibs;
    int jbs;

    /* touch my portion of A[] */

    for (l = 0; l < th_per_p; l++) {
        for (J=0; J<nblocks; J++) {
            for (I=0; I<nblocks; I++) {
                if (block_owner(I, J) == me_th[l]) {
                    if (J == nblocks-1) {
                        jbs = n%bs;
                        if (jbs == 0) {
                            jbs = bs;
                        }
                    } else {
                        jbs = bs;
                    }
                    if (I == nblocks-1) {
                        ibs = n%bs;
                        if (ibs == 0) {
                            ibs = bs;
                        }
                    } else {
                        ibs = bs;
                    }
                    for (j=0; j<jbs; j++) {
                        for (i=0; i<ibs; i++) {
                            tot += a[I+J*nblocks][i+j*ibs];
                        }
                    }
                }
            }
        }
    }
    return(tot);
}

void print_array(int myid)
{
    int i, j, k;
    double **buf;
    armci_hdl_t *null_hdl = NULL;

    int ii, jj;
    int edge;
    int ibs, jbs, skip;

    buf = (double **)malloc(nblocks*nblocks*sizeof(double *));
    for(i=0; i<nblocks; i++)
        for(j=0; j<nblocks; j++)
            if(block_owner(i, j)/th_per_p == myid)
                buf[i+j*nblocks] = a[i+j*nblocks];
            else {
                buf[i+j*nblocks] = (double *)malloc(block_size*block_size*
                                                   sizeof(double));
                get_remote(buf[i+j*nblocks], i, j, &null_hdl);
            }

    /* copied from lu.C */
    edge = n%block_size;
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            if ((n - i) <= edge) {
                ibs = edge;
                ibs = n-edge;
                skip = edge;
            } else {
                ibs = block_size;
                skip = block_size;
            }
            if ((n - j) <= edge) {
                jbs = edge;
                jbs = n-edge;
            } else {
                jbs = block_size;
            }
            ii = (i/block_size) + (j/block_size)*nblocks;
            jj = (i%ibs)+(j%jbs)*skip;
            printf("%8.1f ", buf[ii][jj]);
        }
        printf("\n");
    }
    fflush(stdout);

    for(i=0; i<nblocks; i++)
        for(j=0; j<nblocks; j++)
            if(block_owner(i, j)/th_per_p != myid) free(buf[i+j*nblocks]);
    free(buf);
}

/* prints only square blocks */
void print_block_dbg(double *block, const char *fmt, int x, int y, int z)
{
#ifdef DEBUG2
    int i, j;

    THREAD_LOCK(mutex);

    printf("# ");
    printf(fmt, x, y, z);
    for (i = 0; i < block_size; i++) {
        printf("# ");
        for (j = 0; j < block_size; j++) printf(" %.1f", block[i + j * block_size]);
        printf("\n");
    }
    printf("#\n");
    fflush(stdout);

    THREAD_UNLOCK(mutex);
#endif
}

void *lu(void *lu_arg)
{
    int n, bs, th_idx;
    int i, il, j, jl, k, kl;
    int I, J, K, pI, pJ, AB;
    double *A, *A1, *A2, *B, *B1, *B2, *C, *D;
    int dimI, dimJ, dimK;
    int strI, strJ, strK;
    unsigned int t1, t2, t3, t4, t11, t22;
    int diagowner;
    double *buf1, *buf2, *buf3, *buf4;
    armci_hdl_t hdl1, hdl2, *hdl1p, *hdl2p;

    n = ((int *)lu_arg)[0];
    bs = ((int *)lu_arg)[1];
    th_idx = ((int *)lu_arg)[2];

#ifdef DEBUG
    printf("DBG: starting thread %d(idx=%d) on node %d\n", me_th[th_idx], th_idx, me);
    fflush(stdout);
#endif

    /* temporary memories */
    buf1 = (double *)malloc(block_size*block_size*sizeof(double));
    buf2 = (double *)malloc(block_size*block_size*sizeof(double));
    buf3 = (double *)malloc(block_size*block_size*sizeof(double));
    buf4 = (double *)malloc(block_size*block_size*sizeof(double));

    for (k=0, K=0; k<n; k+=bs, K++) {
        MT_BARRIER(); /* andriy: I think it should be here */

        kl = k + bs;
        if (kl > n) {
            kl = n;
            strK = kl - k;
        } else {
            strK = bs;
        }

        /* factor diagonal block */
        diagowner = block_owner(K, K);
        if (diagowner == me_th[th_idx]) {
            A = a[K+K*nblocks];
            lu0(A, strK, strK);
        }
        MT_BARRIER();

        /* divide column k by diagonal block */
        if(block_owner(K, K) == me_th[th_idx])
            D = a[K+K*nblocks];
        else {
            D = buf1;
            hdl1p = NULL;
            get_remote(D, K, K, &hdl1p);
        }
        for (i=kl, I=K+1; i<n; i+=bs, I++) {
            if (block_owner(I, K) == me_th[th_idx]) {  /* parcel out blocks */
                il = i + bs;
                if (il > n) {
                    il = n;
                    strI = il - i;
                } else {
                    strI = bs;
                }
                A = a[I+K*nblocks];
                bdiv(A, D, strI, strK, strI, strK);
            }
        }

        /* modify row k by diagonal block */
        for (j=kl, J=K+1; j<n; j+=bs, J++) {
            if (block_owner(K, J) == me_th[th_idx]) { /* parcel out blocks */
                jl = j+bs;
                if (jl > n) {
                    jl = n;
                    strJ = jl - j;
                } else {
                    strJ = bs;
                }
                A = a[K+J*nblocks];
                bmodd(D, A, strK, strJ, strK, strK);
            }
        }

        MT_BARRIER();
        /* prefetch (A1 and B1) */
        AB = 0;  /* 0 if using A1 and B1 (buf1 and buf2), 1 if A2 and B2 (buf3 and buf4) */
        if (next_block(th_idx, bs, kl, kl, kl, K+1, K+1, K, &I, &J)) {
            report(th_idx, "ij: next %d,%d\n", I, J);
            /* next block to be computed (I,J) needs blocks A=(I,K) and J=(K,J) */
            hdl1p = &hdl1; hdl2p = &hdl2;
            prefetch(&A1, buf1, I, K, th_idx, &hdl1p);
            prefetch(&B1, buf2, K, J, th_idx, &hdl2p);
        } else continue;

        /* modify subsequent block columns */
        for (i=kl, I=K+1; i<n; i+=bs, I++) {
           il = i+bs;
            if (il > n) {
                il = n;
                strI = il - i;
            } else {
                strI = bs;
            }

            for (j=kl, J=K+1; j<n; j+=bs, J++) {
                jl = j + bs;
                if (jl > n) {
                    jl = n;
                    strJ= jl - j;
                } else {
                    strJ = bs;
                }
                if (block_owner(I, J) == me_th[th_idx]) {  /* parcel out blocks */
                    report(th_idx, "ij: real %d,%d\n", I, J);
                    /* wait for previously prefetched */
                    A = AB ? A2 : A1;
                    B = AB ? B2 : B1;
                    /* actual wait */
                    if (hdl1p) ARMCI_Wait(hdl1p);
                    if (hdl2p) ARMCI_Wait(hdl2p);
                    /* prefetch next A and B */
                    if (next_block(th_idx, bs, kl, i, j+bs, I, J+1, K, &pI, &pJ)) {
                        report(th_idx, "ij: next %d,%d\n", pI, pJ);
                        hdl1p = &hdl1; hdl2p = &hdl2;
                        if (AB) { /* prefetch into A2 and B2 */
                            prefetch(&A1, buf1, pI, K, th_idx, &hdl1p);
                            prefetch(&B1, buf2, K, pJ, th_idx, &hdl2p);
                        } else {  /* prefetch into A1 and B1 */
                            prefetch(&A2, buf3, pI, K, th_idx, &hdl1p);
                            prefetch(&B2, buf4, K, pJ, th_idx, &hdl2p);
                         }
                        AB = AB ? 0 : 1;
                    }
                    C = a[I+J*nblocks];
                    bmod(A, B, C, strI, strJ, strK, strI, strK, strI);
                }
            }
        }
    }

    free(buf1);
    free(buf2);
    free(buf3);
    free(buf4);
}

void prefetch(double **A, double *buf, int I, int J, int th_idx, armci_hdl_t **hdlp)
{
    if (block_owner(I,J) == me_th[th_idx]) {
        *A = a[I+J*nblocks];
        *hdlp = NULL; /* local should not ARMCI_Wait */
    } else {
        *A = buf;
        get_remote(*A, I, J, hdlp);
    }
}

void get_remote(double *buf, int I, int J, armci_hdl_t **hdlp)
{
    int proc_owner;
    int edge, size;

#ifdef USE_MUTEX
    THREAD_LOCK(mutex);
#endif

    proc_owner = block_owner(I, J) / th_per_p;

    edge = n%block_size;
    if (edge == 0) {
        edge = block_size;
    }

    if ((I == nblocks-1) && (J == nblocks-1)) {
        size = edge*edge;
    }
    else if ((I == nblocks-1) || (J == nblocks-1)) {
        size = edge*block_size;
    }
    else {
        size = block_size*block_size;
    }
    size = size * sizeof(double);

    if (proc_owner == me) {
        memcpy(buf, a[I+J*nblocks], size);
        *hdlp = NULL; /* local should not ARMCI_Wait */
    } else
        if (*hdlp) ARMCI_NbGet(a[I+J*nblocks], buf, size, proc_owner, *hdlp);
        else ARMCI_Get(a[I+J*nblocks], buf, size, proc_owner);

#ifdef USE_MUTEX
    THREAD_UNLOCK(mutex);
#endif
}

/* returns 1 if there is another block on current processor to be computed for some K
 * returns 0 otherwise; location of the block is stored in pI, pJ */
int next_block(int th_idx, int bs, int kl, int ci, int cj, int cI, int cJ, int K, int *pI, int *pJ)
{
    int i, j, I, J;

    j = cj;
    J = cJ;
    for (i=ci, I=cI; i<n; i+=bs, I++) {
        for (; j<n; j+=bs, J++) {
            if (block_owner(I, J) == me_th[th_idx]) {
                *pI = I;
                *pJ = J;
                return 1;
            }
        }
        j=kl;
        J=K+1;
    }
    return 0;
}
