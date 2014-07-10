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


void do_work()
{
int ONE=1 ;   /* useful constants */
int g_a, g_b;
int n=N, type=MT_F_DBL;
int me=GA_Nodeid(), nproc=GA_Nnodes();
int i, row;
int dims[2]={N,N};
int lo[2], hi[2];

/* Note: on all current platforms DoublePrecision == double */
double buf[N], err, alpha, beta;

     if(me==0)printf("Creating matrix A\n");
     g_a = NGA_Create(type, 2, dims, "A", NULL);
     if(!g_a) GA_Error("create failed: A",n); 
     if(me==0)printf("OK\n");

     if(me==0)printf("Creating matrix B\n");
     /* create matrix B  so that it has dims and distribution of A*/
     g_b = GA_Duplicate(g_a, "B");
     if(! g_b) GA_Error("duplicate failed",n); 
     if(me==0)printf("OK\n");

     GA_Zero(g_a);   /* zero the matrix */

     if(me==0)printf("Initializing matrix A\n");
     /* fill in matrix A with random values in range 0.. 1 */ 
     lo[1]=0; hi[1]=n-1;
     for(row=me; row<n; row+= nproc){
         /* each process works on a different row in MIMD style */
         lo[0]=hi[0]=row;   
         for(i=0; i<n; i++) buf[i]=sin((double)i + 0.1*(row+1));
         NGA_Put(g_a, lo, hi, buf, &n);
     }


     if(me==0)printf("Symmetrizing matrix A\n");
     GA_Symmetrize(g_a);   /* symmetrize the matrix A = 0.5*(A+A') */
   

     /* check if A is symmetric */ 
     if(me==0)printf("Checking if matrix A is symmetric\n");
     GA_Transpose(g_a, g_b); /* B=A' */
     alpha=1.; beta=-1.;
     GA_Add(&alpha, g_a, &beta, g_b, g_b);  /* B= A - B */
     err= GA_Ddot(g_b, g_b);
     
     if(me==0)printf("Error=%f\n",(double)err);
     
     if(me==0)printf("\nChecking atomic accumulate \n");

     GA_Zero(g_a);   /* zero the matrix */
     for(i=0; i<n; i++) buf[i]=(double)i;

     /* everybody accumulates to the same location/row */
     alpha = 1.0;
     row = n/2;
     lo[0]=hi[0]=row;
     lo[1]=0; hi[1]=n-1;
     NGA_Acc(g_a, lo, hi, buf, &ONE, &alpha );
     GA_Sync();

     if(me==0){ /* node 0 is checking the result */

        NGA_Get(g_a, lo, hi, buf,&ONE);
        for(i=0; i<n; i++) if(buf[i] != (double)nproc*i)
           GA_Error("failed: column=",i);
        printf("OK\n\n");

     }
     
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
    if(me==0) {
       if(GA_Uses_fapi())GA_Error("Program runs with C array API only",1);
       printf("Using %ld processes\n",(long)nproc);
       fflush(stdout);
    }

    heap /= nproc;
    stack /= nproc;
    if(! MA_init(MT_F_DBL, stack, heap)) 
       GA_Error("MA_init failed",stack+heap);  /* initialize memory allocator*/ 
    
    do_work();

    if(me==0)printf("Terminating ..\n");
    GA_Terminate();

    MP_FINALIZE();

    return 0;
}

