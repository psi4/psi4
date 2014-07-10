#if HAVE_CONFIG_H
#   include "config.h"
#endif

/****************************************************************************
* program: ga-mpi.c
* date:    Tue Oct  3 12:31:59 PDT 1995
* author:  Jarek Nieplocha
* purpose: This program demonstrates interface between GA and MPI. 
*          For a given square matrix, it creates a vector that contains maximum
*          elements for each matrix row. MPI group communication is used.
*
* notes:   The program can run in two modes:
*          1. Using TCGMSG calls available through the TCGMSG-MPI library
*             and MPI. In this mode initialization must be done with 
*             the TCGMSG PBEGIN call.
*          2. Using MPI calls only -- preprocessor symbol MPI must be defined.
*
****************************************************************************/

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif

#include <mpi.h>

#include "ga.h"
#include "globalp.h"
#include "macdecls.h"
#include "mp3.h"


#define N 100           /* dimension of matrices */


void do_work()
{
int ZERO=0;   /* useful constants */
int g_a, g_b;
int n=N, ndim=2,type=MT_F_DBL,dims[2]={N,N},coord[2];
int me=GA_Nodeid(), nproc=GA_Nnodes();
int row, i, j;
 int lo[2], hi[2];

/* Note: on all current platforms DoublePrecision = double */
DoublePrecision buf[N], *max_row=NULL;

MPI_Comm ROW_COMM;
int ilo,ihi, jlo,jhi, ld, prow, pcol;
int root=0, grp_me=-1;

     if(me==0)printf("Creating matrix A\n");
     dims[0]=n; dims[1]=n;
     g_a = NGA_Create(type, ndim, dims, "A", NULL);
     if(!g_a) GA_Error("create failed: A",n); 
     if(me==0)printf("OK\n");
     
     if(me==0)printf("Creating matrix B\n");
     dims[0]=n;
     g_b = NGA_Create(type, 1, dims, "B", NULL);
     if(!g_b) GA_Error("create failed: B",n); 
     if(me==0)printf("OK\n");
     
     GA_Zero(g_a);   /* zero the matrix */
     
     if(me==0)printf("Initializing matrix A\n");
     /* fill in matrix A with values: A(i,j) = (i+j) */ 
     for(row=me; row<n; row+= nproc){
    /**
     * simple load balancing: 
     * each process works on a different row in MIMD style 
     */ 
    for(i=0; i<n; i++) buf[i]=(DoublePrecision)(i+row+1); 
    lo[0]=hi[0]=row;
    lo[1]=ZERO;  hi[1]=n-1; 
    NGA_Put(g_a, lo, hi, buf, &n); 
     }
     
     /* GA_print(&g_a);*/
     NGA_Distribution(g_a, me, lo, hi);
     ilo=lo[0]; ihi=hi[0];
     jlo=lo[1]; jhi=hi[1];
     
     GA_Sync(); 
     if(ihi-ilo+1 >0){
        max_row=(DoublePrecision*)malloc(sizeof(DoublePrecision)*(ihi-ilo+1));
        if (!max_row) GA_Error("malloc 3 failed",(ihi-ilo+1));
        for (i=0; i<(ihi-ilo+1); i++) {
            max_row[i] = 0.0;
        }
     }
     NGA_Proc_topology(g_a, me, coord);  /* block coordinates */
     prow = coord[0];
     pcol = coord[1];

     if(me==0)printf("Splitting comm according to distribution of A\n");
     
     /* GA on SP1 requires synchronization before & after message-passing !!*/
     GA_Sync(); 
     
     if(me==0)printf("Computing max row elements\n");
     /* create communicator for processes that 'own' A[:,jlo:jhi] */
     MPI_Barrier(MPI_COMM_WORLD);
     if(pcol < 0 || prow <0)
    MPI_Comm_split(MPI_COMM_WORLD,MPI_UNDEFINED,MPI_UNDEFINED, &ROW_COMM);
     else
    MPI_Comm_split(MPI_COMM_WORLD, (int)pcol, (int)prow, &ROW_COMM);
     
     if(ROW_COMM != MPI_COMM_NULL){
    double *ptr;
    MPI_Comm_rank(ROW_COMM, &grp_me);
    
    /* each process computes max elements in the block it 'owns' */
    lo[0]=ilo; hi[0]=ihi;
    lo[1]=jlo; hi[1]=jhi;
    NGA_Access(g_a, lo, hi, &ptr, &ld);
    for(i=0; i<ihi-ilo+1; i++){
       for(j=0; j<jhi-jlo+1; j++)
          if(max_row[i] < ptr[i*ld + j]){
         max_row[i] = ptr[i*ld + j];
          }
    }
    MPI_Reduce(max_row, buf, ihi-ilo+1, MPI_DOUBLE, MPI_MAX,
           root, ROW_COMM);
    
     }else fprintf(stderr,"process %d not participating\n",me);
     GA_Sync(); 
     
     /* processes with rank=root in ROW_COMM put results into g_b */
     ld = 1;
     if(grp_me == root) {
    lo[0]=ilo;  hi[0]=ihi;
    NGA_Put(g_b, lo, hi, buf, &ld); 
     }
        
     GA_Sync();

     if(me==0)printf("Checking the result\n");
     if(me==0){
    lo[0]=ZERO; hi[0]=n-1;
        NGA_Get(g_b, lo, hi, buf, &n); 
        for(i=0; i< n; i++)if(buf[i] != (double)n+i){
            fprintf(stderr,"error:%d max=%f should be:%d\n",i,buf[i],n+i);
            GA_Error("terminating...",1);
        }
     }
     
     if(me==0)printf("OK\n");

     GA_Destroy(g_a);
     GA_Destroy(g_b);
}
     


int main(argc, argv)
int argc;
char **argv;
{
int heap=20000, stack=20000;
int me, nproc;

    MP_INIT(argc,argv);
    GA_INIT(argc,argv);                            /* initialize GA */
    me=GA_Nodeid();
    nproc=GA_Nnodes();
    if(me==0) printf("Using %ld processes\n",(long)nproc);

    heap /= nproc;
    stack /= nproc;
    if(! MA_init((int)MT_F_DBL, stack, heap)) 
       GA_Error("MA_init failed",stack+heap);   /* initialize memory allocator*/ 
    do_work();

    if(me==0)printf("Terminating ..\n");
    GA_Terminate();

#   ifdef MPI
      MPI_Finalize();
#   else
      tcg_pend();
#   endif

    return 0;
}

