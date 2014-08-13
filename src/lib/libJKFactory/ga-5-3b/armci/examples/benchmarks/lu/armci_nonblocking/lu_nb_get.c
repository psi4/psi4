#if HAVE_CONFIG_H
#   include "config.h"
#endif

/**************************************************
 *             LU factorization                   *
 *             Armci Version                      *
 **************************************************/
 
/*****************
Non-blocking Version 
Pre-GETing
******************/ 

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_MATH_H
#   include <math.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_UNISTD_H
#   include <unistd.h>
#endif
#if HAVE_STRING_H
#   include <string.h>
#endif

#include "armci.h"
#include "message.h"

/* #define DEBUG */
#define MAXRAND                         32767.0
#define DEFAULT_N                       8
#define DEFAULT_B                       2

#define MAXPROC 256 /* Maximum number of processors */
#define MAXBLOCK 2048 /* Maximum number of blocks in a row/column */
#define ANULL (armci_hdl_t *)NULL

/* global variables */
int n = DEFAULT_N;         /* The size of the matrix */
int block_size = DEFAULT_B;/* Block dimension */
int nblocks;               /* Number of blocks in each dimension */
int num_rows;              /* Number of processors per row of processor grid */
int num_cols;              /* Number of processors per col of processor grid */
double **a;                /* a = lu; l and u both placed back in a */
int nproc, me = 0;
int proc_bytes;

int doprint = 1; /* make it 1 to see the LU decomposition output */

/*funnction declaration */
void lu(int, int, int);
void lu0(double *,int, int);
void bdiv(double *, double *, int, int, int, int);
void bmodd(double *, double*, int, int, int, int);
void bmod(double *, double *, double *, int, int, int, int, int, int);
void daxpy(double *, double *, int, double);
int block_owner(int, int);
void init_array();
double touch_array(int, int);
void print_block();
void print_array(int);
void get_remote(double *, int, int, armci_hdl_t *);

/* timing functions */
extern void start_timer(void);
extern double elapsed_time(void);
extern double stop_time(void);

int main(int argc, char *argv[])
{
  int i, j;
  int ch;
  int edge;
  int size;
    
  /* ARMCI */
  void **ptr;
  double **ptr_loc;
  
  armci_msg_init(&argc,&argv);
  nproc = armci_msg_nproc();
  me = armci_msg_me();
    
  while ((ch = getopt(argc, argv, "n:b:p:h")) != -1) {
    switch(ch) {
    case 'n': n = atoi(optarg); break;
    case 'b': block_size = atoi(optarg); break;
    case 'p': nproc = atoi(optarg); break;
    case 'h': {
      printf("Usage: LU, or \n");
      printf("       LU -nMATRIXSIZE -bBLOCKSIZE -pNPROC\n");
      armci_msg_barrier();
      armci_msg_finalize();
      exit(0);
    }            
    }
  }
    
  if(me == 0) {
    printf("\n Using pre_FETCHing \n");
    printf("\n Blocked Dense LU Factorization\n");
    printf("     %d by %d Matrix\n", n, n);
    printf("     %d Processors\n", nproc);
    printf("     %d by %d Element Blocks\n", block_size, block_size);
    printf("\n");
  }
    
  num_rows = (int) sqrt((double) nproc);
  for (;;) {
    num_cols = nproc/num_rows;
    if (num_rows*num_cols == nproc)
      break;
    num_rows--;
  }
    
  nblocks = n/block_size;
  if (block_size * nblocks != n) {
    nblocks++;
  }
    
  edge = n%block_size;
  if (edge == 0) {
    edge = block_size;
  }
    
#ifdef DEBUG
  if(me == 2) {
    for (i=0;i<nblocks;i++) {
      for (j=0;j<nblocks;j++) 
    printf("%d ", block_owner(i, j));
      printf("\n");
    }
  }
  armci_msg_barrier();
  armci_msg_finalize();
  exit(0);
#endif
    
  for (i=0;i<nblocks;i++) {
    for (j=0;j<nblocks;j++) {
      if(block_owner(i,j) == me) {
    if ((i == nblocks-1) && (j == nblocks-1)) {
      size = edge*edge;
    }
    else if ((i == nblocks-1) || (j == nblocks-1)) {
      size = edge*block_size;
    }
    else {
      size = block_size*block_size;
    }
    proc_bytes += size*sizeof(double);
      }
    }
  }
    
  /* initialize ARMCI */
  ARMCI_Init();
  ptr = (void **)ARMCI_Malloc_local(nproc * sizeof(void *));
  ARMCI_Malloc(ptr, proc_bytes);
  
  a = (double **)ARMCI_Malloc_local(nblocks*nblocks*sizeof(double *));
  if (a == NULL) {
    fprintf(stderr, "Could not malloc memory for a\n");
    exit(-1);
  } 
  ptr_loc = (double **)ARMCI_Malloc_local(nproc*sizeof(double *));
  for(i=0; i<nproc; i++) ptr_loc[i] = (double *)ptr[i];
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
    
  /* initialize the array */
  init_array();
  
  /* barrier to ensure all initialization is done */
  armci_msg_barrier();

  /* to remove cold-start misses, all processors touch their own data */
  touch_array(block_size, me);
  armci_msg_barrier();

  if(doprint) {
    if(me == 0) {
      printf("Matrix before LU decomposition\n");
      print_array(me); 
    }
    armci_msg_barrier();
  }  

  /* Starting the timer */
  if(me == 0) start_timer();

  lu(n, block_size, me);
    
  armci_msg_barrier();

  /* Timer Stops here */
  if(me == 0) 
  printf("\nRunning time = %f milliseconds.\n\n",  elapsed_time());

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

  return 0;
}

void lu(int n, int bs, int me)
{
  int i, il, j, jl, k, kl;
  int I, J, K;
  double *A, *B, *C, *D;
  int strI, strJ, strK;
  /*unsigned int t1, t2, t3, t4, t11, t22;*/
  int diagowner;
  double *dbuf, **bufr, **bufc;
  armci_hdl_t *hdl1=NULL, *hdl2=NULL;
  int r, c, w, t1, t2, br[MAXBLOCK], bc[MAXBLOCK];  
  
  dbuf = (double *)ARMCI_Malloc_local((armci_size_t) block_size*block_size*sizeof(double));

  bufr = (double **)ARMCI_Malloc_local(nblocks * sizeof(double *));
  bufc = (double **)ARMCI_Malloc_local(nblocks * sizeof(double *));
  if (bufr == NULL || bufc == NULL)
    printf("Could not ARMCI_Malloc_local() mem\n");

  for (i = 0; i < nblocks; i++) {
    bufr[i] = (double *)ARMCI_Malloc_local(block_size*block_size*sizeof(double));
    bufc[i] = (double *)ARMCI_Malloc_local(block_size*block_size*sizeof(double));
  }
  
  for (k=0, K=0; k<n; k+=bs, K++) {
    kl = k + bs; 
    if (kl > n) {
      kl = n;
      strK = kl - k;
    } else {
      strK = bs;
    }
    
    /* factor diagonal block */
    diagowner = block_owner(K, K);
    if (diagowner == me) {
      A = a[K+K*nblocks]; 
      lu0(A, strK, strK);
    }
    armci_msg_barrier(); 
    
    /* divide column k by diagonal block */
    if(block_owner(K, K) == me)
      D = a[K+K*nblocks];
    else {
      D = dbuf;
      get_remote(D, K, K, NULL);
    }
    
    for (i=kl, I=K+1; i<n; i+=bs, I++) {
      if (block_owner(I, K) == me) {  /* parcel out blocks */
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
      if (block_owner(K, J) == me) {  /* parcel out blocks */
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
        
    armci_msg_barrier();

    /* modify subsequent block columns */
    
    t1 = t2 = 0;
    memset(br, 0, sizeof(br));
    memset(bc, 0, sizeof(bc));

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

    if (block_owner(I, J) == me) {  /* parcel out blocks */
    
      /* Pre-fetch next two blocks that will be required by me */
      /* First, identify the next IJ-th block owned by me */
      /* This caculation is for block-cyclic distribution */
     
        r = I;
        c = J + num_cols;
        if (c >= nblocks) {
          r = I + num_rows;
          w = J - (K+1);
          if (w >= num_cols)
        c = w%num_cols + (K+1);
          else
        c = J;
        }
      
      /* This processor will need the blocks [r,K] and [K, c] next */  
      /* Now, pre-fetch blocks [r,K] and [K,c] using non-blocking gets*/
      if (r <  nblocks && c < nblocks) {
        if (!br[c] && block_owner(K, c) != me) { /* if this block has not been pre-fetched yet and if I already don't own it*/
          if (hdl1 == NULL) {/* this is the first time, no previous non-blocking call */
          get_remote(bufr[c], K, c, hdl1);
        }
        else {
          if (!ARMCI_Wait(hdl1)) {/* only if previous call with hdl1 returned, then fetch next block */
        get_remote(bufr[c], K, c, hdl1);
        t1 = 1;
          }
        }
      }
      
      if (!bc[r] && block_owner(r, K) != me) {        
        if (hdl2 == NULL)
          get_remote(bufc[r], r, K, hdl2);
        else {
          if (!ARMCI_Wait(hdl2)) {
        get_remote(bufc[r], r, K, hdl2);
        t2 = 1;
          }
        }
      }
      } /* end of if (r < nblocks && c < nblocks) */

      if(block_owner(I,K) == me)
        A = a[I+K*nblocks];
      else {
        if (!t1)
          get_remote(bufc[I], I, K, NULL); /* This is the first time, so make a blocking call */
        A = bufc[I];
        bc[I] = 1;
      }
      
      if(block_owner(K,J) == me)
        B = a[K+J*nblocks];
      else {
        if (!t2)
          get_remote(bufr[J], K, J, NULL); /* This is the first time, so make a blocking call */
        B = bufr[J];
        br[J] = 1;
      }
      C = a[I+J*nblocks];
      bmod(A, B, C, strI, strJ, strK, strI, strK, strI);
    }
      }
    }
  }
  
  ARMCI_Free_local(dbuf);
  ARMCI_Free_local(bufr);
  ARMCI_Free_local(bufc);
}

void get_remote(double *buf, int I, int J, armci_hdl_t *handle)
{
  int proc_owner;
  int edge, size;
  int rc;
  
  proc_owner = block_owner(I, J);
    
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
  
  if (handle == NULL) {/* do a blocking get */
    ARMCI_Get(a[I+J*nblocks], buf, size, proc_owner);
  }
  else {
     if ((rc = ARMCI_NbGet(a[I+J*nblocks], buf, size, proc_owner, handle))) /* do a non-blocking get */
      ARMCI_Error("Error in ARMCI_NbGet", rc);
  }
  
}

void lu0(double *a, int n, int stride)
{
  int j; 
  int k; 
  /*int length;*/
  double alpha;
    
  for (k=0; k<n; k++) {
    /* modify subsequent columns */
    for (j=k+1; j<n; j++) {
      a[k+j*stride] /= a[k+k*stride];
      alpha = -a[k+j*stride];
      /*length = n-k-1;*/
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
  int j; 
  int k; 
  /*int length;*/
  double alpha;
    
  for (k=0; k<dimi; k++) {
    for (j=0; j<dimj; j++) {
      c[k+j*stride_c] /= a[k+k*stride_a];
      alpha = -c[k+j*stride_c];
      /*length = dimi - k - 1;*/
      daxpy(&c[k+1+j*stride_c], &a[k+1+k*stride_a], dimi-k-1, alpha);
    }
  }
}

void bmod(double *a, double *b, double *c, int dimi, int dimj, int dimk,
          int stridea, int strideb, int stridec)
{
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
  return((J%num_cols) + (I%num_rows)*num_cols); 
}

void init_array()
{
  int i, j;
  int ii, jj;
  int edge;
  int ibs;
  int jbs, skip;
    
  srand48((long) 1);
  edge = n%block_size;
  for (j=0; j<n; j++) {
    for (i=0; i<n; i++) {
      if(block_owner((i/block_size), (j/block_size)) == me) {
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
    /*            a[ii][jj] = ((double) lrand48())/MAXRAND; */
    a[ii][jj] = i + j*6 + 1;
    if (i == j) {
      a[ii][jj] *= 10;
    }
      }
    }
  }
}

double touch_array(int bs, int me)
{
  int i, j, I, J;
  double tot = 0.0;
  int ibs;
  int jbs;
    
  /* touch my portion of A[] */
    
  for (J=0; J<nblocks; J++) {
    for (I=0; I<nblocks; I++) {
      if (block_owner(I, J) == me) {
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
  return(tot);
}

void print_array(int myid)
{
  int i, j;
  double **buf;

  int ii, jj;
  int edge;
  int ibs, jbs, skip;

  buf = (double **)ARMCI_Malloc_local(nblocks*nblocks*sizeof(double *));

  for(i=0; i<nblocks; i++) 
    for(j=0; j<nblocks; j++) 
      if(block_owner(i, j) == myid)
    buf[i+j*nblocks] = a[i+j*nblocks];
      else {
    buf[i+j*nblocks] = (double *)ARMCI_Malloc_local(block_size*block_size*
                        sizeof(double));
    get_remote(buf[i+j*nblocks], i, j, NULL);
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
      if(block_owner(i, j) != myid) ARMCI_Free_local(buf[i+j*nblocks]);
  
  ARMCI_Free_local(buf);
}

void print_block()
{
  int i, j, k;

  for(i=0; i<nblocks; i++)
    for(j=0; j<nblocks; j++)
      if(block_owner(i,j) == me) {
    printf("Block %d (%d,%d)\t", i+j*nblocks, i, j);
    for(k=0; k<block_size*block_size; k++)
      printf("%8.1f ", a[i+j*nblocks][k]);
    printf("\t me = %d\n", me);
      }
}
