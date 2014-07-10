#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_MATH_H
#   include <math.h>
#endif
#if HAVE_STDLIB_H
#include <stdlib.h>
#endif
#if HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif
#if HAVE_SYS_STAT_H
#include <sys/stat.h>
#endif
#if HAVE_FCNTL_H
#include <fcntl.h>
#endif

#include "armci.h"
#include "ga.h"
#include "macdecls.h"
#include "finclude.h"
#include "mp3.h"

#define VERIFY_RESULT 1

int na,nz;
int bvec,dvec,svec,dmvec,m_dvec,amat,xvec,axvec,rvec,qvec,ridx,cidx;
int me, nproc;
int myfirstrow=0,mylastrow=0;
double epsilon=1e-4;
double time_get=0;
double *entirexvecptr,*xvecptr,*entiredvecptr,*dvecptr;
int isvectormirrored=0;
static int niter;
void read_and_create(int,char **);
void computeminverser(double *,double *, double *);
void computeminverse(double *,double *, int *,int *);
void finalize_arrays(int);
extern void matvecmul(double *,int,double *,int,int *,int *);
extern double *ga_vecptr;
void conjugate_gradient(int nit,int dopreconditioning)
{
int i,zero=0;
int lo,hi;
double d_one=1.0,d_negone=-1.0;
double delta0=0.0,deltaold=0.0,deltanew=0.0,alpha=0.0,negalpha,beta,dtransposeq;
double *axvecptr,*qvecptr,*aptr,*dmvecptr,*rvecptr,*svecptr,*bvecptr;
double time0;
int *mycp,*myrp;

 int j;
 double sum;

    NGA_Distribution(cidx,me,&lo,&hi);
    NGA_Access(cidx,&lo,&hi,&mycp,&zero);
    NGA_Access(amat,&lo,&hi,&aptr,&zero);

    NGA_Distribution(ridx,me,&lo,&hi);
    NGA_Access(ridx,&lo,&hi,&myrp,&zero);
    NGA_Access(axvec,&lo,&hi,&axvecptr,&zero);
    NGA_Access(qvec,&lo,&hi,&qvecptr,&zero);
    NGA_Access(rvec,&lo,&hi,&rvecptr,&zero);
    NGA_Access(bvec,&lo,&hi,&bvecptr,&zero);

    if(dopreconditioning){
       NGA_Access(dmvec,&lo,&hi,&dmvecptr,&zero);
       NGA_Access(svec,&lo,&hi,&svecptr,&zero);
       /* NGA_Distribution(dvec, 0, &lo, &hi);
      NGA_Access(dvec, &lo, &hi, &entiredvecptr, &zero); */
    }
    printf("\n%d:before matvecmul\n",me);fflush(stdout);
    /* compute Ax */
    f_matvecmul(aptr,entirexvecptr,axvecptr,&zero,&myfirstrow,&mylastrow,myrp,mycp);
    /* r=b-Ax */
    f_addvec(&d_one,bvecptr,&d_negone,axvecptr,rvecptr,&myfirstrow,&mylastrow); 

    if(dopreconditioning){
      f_computeminverse(dmvecptr,aptr,myrp,mycp,&myfirstrow,&mylastrow);
      f_computeminverser(dmvecptr,rvecptr,dvecptr,&myfirstrow,&mylastrow);
      NGA_Put(dvec,&lo,&hi,dvecptr,&hi);
      if (me == 0)
    printf("Doing preconditioning!\n");
    }
    else{
      if(me==0){
         na--;
         NGA_Get(rvec,&zero,&na,entiredvecptr,&na);
         NGA_Put(dvec,&zero,&na,entiredvecptr,&na);
         na++;
      }
    }

    deltanew = GA_Ddot(dvec,rvec);                /* deltanew = r.r_tranpose */

    delta0 = deltanew;                            /* delta0 = deltanew */

    if(me==0)printf("\n\tdelta0 is %f\n",delta0);
    /*if(me==0)printf("\n\titer\tbeta\tdelta");*/

    for(i=0;i<nit && deltanew>(1e-8*delta0);i++){
      na--;
      NGA_Get(dvec, &zero, &na, entiredvecptr, &na);
      na++;

       if(isvectormirrored)
         matvecmul(aptr,m_dvec,qvecptr,1,myrp,mycp);/* q = Ad */
       else{
         f_matvecmul(aptr,entiredvecptr,qvecptr,&zero,&myfirstrow,&mylastrow,myrp,mycp);

     sum = 0.0;
     for (j = 0; j < na; j++)
       if (entiredvecptr[j] != 0.0)
         sum += entiredvecptr[j];

     /* if (me == 0)
        printf("me: %d, sum: %g\n", me, sum); */
       }

       NGA_Put(dvec,&lo,&hi,dvecptr,&hi);
       dtransposeq=GA_Ddot(dvec,qvec);            /* compute d_transpose.q */

       alpha = deltanew/dtransposeq;              /* deltanew/(d_transpose.q) */

       if(i>10000 && i%25==0){
         /* compute Ax*/
         f_matvecmul(aptr,entirexvecptr,axvecptr,&zero,&myfirstrow,&mylastrow,myrp,mycp);
         /* x = x+ alpha.d*/ /* r=b-Ax*/
         f_2addvec(&d_one,xvecptr,&alpha,dvecptr,xvecptr,&d_one,bvecptr,
                         &d_negone,axvecptr,rvecptr,&myfirstrow,&mylastrow);
       }
       else{
         negalpha = 0.0-alpha;                         
         /* x = x+ alpha.d*/ /* r=r-alpha.q*/
         f_2addvec(&d_one,xvecptr,&alpha,dvecptr,xvecptr,&d_one,rvecptr,
                         &negalpha,qvecptr,rvecptr,&myfirstrow,&mylastrow);
       }

       if(dopreconditioning)
         computeminverser(dmvecptr,rvecptr,svecptr);

       deltaold = deltanew;                        /* deltaold = deltanew*/

       if(dopreconditioning)
         deltanew = GA_Ddot(svec,rvec);            /* deltanew = r_transpose.r*/
       else
         deltanew = GA_Ddot(rvec,rvec);            /* deltanew = r_transpose.r*/

       beta = deltanew/deltaold;                   /* beta = deltanew/deltaold*/

       if(dopreconditioning)
         f_addvec(&d_one,svecptr,&beta,dvecptr,dvecptr,&myfirstrow,&mylastrow);    /* d = s + beta.d */
       else
         f_addvec(&d_one,rvecptr,&beta,dvecptr,dvecptr,&myfirstrow,&mylastrow);    /* d = r + beta.d */

       if(isvectormirrored)
         GA_Copy(dvec,m_dvec);                     /*copy from distributed */

       /* if(me==0)printf("\n\t%d\t%0.4f\t%f",(i+1),beta,deltanew); */
    }
    if(i < nit && me == 0)
        printf("\n Done with CG before reaching max iter %f",sqrt(deltanew/delta0));
    niter = i;

#if VERIFY_RESULT
    GA_Zero(qvec);
    GA_Zero(rvec);
    matvecmul(aptr,xvec,qvecptr,0,myrp,mycp);
    GA_Add(&d_one,qvec,&d_negone,bvec,rvec);
    time0=GA_Ddot(rvec,rvec);
    if(me==0)printf("\n%d:error is %f",me,time0);
#endif
    
}

void initialize_arrays(int dpc)
{
double d_one=1.0;
int i;
    GA_Zero(dvec);
    GA_Fill(xvec,&d_one);
    GA_Zero(axvec);
    GA_Zero(rvec);
    GA_Zero(qvec);
    if(dpc){
       GA_Zero(dmvec);
       GA_Zero(svec);
    }
    for(i=0;i<na;i++)
            entirexvecptr[i]=1.0;
}

void **myptrarrx;
void **myptrarrd;
static void create_entire_vecs()
{
int i,lo,hi;
    myptrarrx = (void **)malloc(sizeof(void*)*nproc);
    myptrarrd = (void **)malloc(sizeof(void*)*nproc);

    i=ARMCI_Malloc(myptrarrx,na*sizeof(double));
    if(i!=0)GA_Error("malloc failed",0);
    entirexvecptr=myptrarrx[me];

    i=ARMCI_Malloc(myptrarrd,na*sizeof(double));
    if(i!=0)GA_Error("malloc failed",0);
    entiredvecptr=myptrarrd[me];
    
    NGA_Distribution(ridx,me,&lo,&hi);

    xvecptr=entirexvecptr+lo;
    dvecptr=entiredvecptr+lo;

    printf("me: %d, entiredvecptr: %p, dvecptr: %p\n",
            me, (void*)entiredvecptr, (void*)dvecptr);
}

int main(argc, argv)
int argc;
char **argv;
{
int heap=200000, stack=200000;
int dopreconditioning=1;
double time0,time1;

    MP_INIT(argc, argv);                    /* initialize message passing */
    GA_Initialize();                           /* initialize GA */

    me=GA_Nodeid(); 
    nproc=GA_Nnodes();
    if(me==0)printf("\n                          CONJUGATE GRADIENT EXAMPLE\n");
    if(argc<3){
       if(me==0){
         printf(" CORRECT USAGE IS:");
         printf("\n\n <launch commands> ga_cg.x na nz file");
         printf("\n\n where:");
         printf("\n\tna is array dimention (only square arrays supported)");
         printf("\n\tnz is number of non-zeros");
         printf("\n\tfile is either the input file or the word random");
         printf("\n\t  use the word random if you to use random input");
         printf("\n\t  input should be in row compressed format");
         printf("\n\t  file should have matrix a followed by row, col & b (Ax=b)");
         printf("\n\t  if file also has na and nz, pass them as 0's and the");
         printf("\n\t  program will read them from the file");
         printf("\n\nexample usages are:");
         printf("\n\tmpirun -np 4 ./ga_cg.x 5000 80000 /home/me/myinput.dat");
         printf("\n\tor");
         printf("\n\tmpirun -np 4 ./ga_cg.x 5000 80000 random\n\n");
         fflush(stdout);
       }
       GA_Terminate();
       MP_FINALIZE();
       return 0;
    }

    heap /= nproc;
    stack /= nproc;
    if(! MA_init(MT_F_DBL, stack, heap)) 
       GA_Error("MA_init failed",stack+heap);  /* initialize memory allocator*/ 
    
    read_and_create(argc,argv);
    create_entire_vecs();

    if(me==0)printf("\nWarmup and initialization run");
    initialize_arrays(dopreconditioning);
    conjugate_gradient(1,dopreconditioning);
    time_get =0.0;
    if(me==0)printf("\n\nStarting Conjugate Gradient ....");
    initialize_arrays(dopreconditioning);

    time0=MP_TIMER();
    conjugate_gradient(30000/*2*/,dopreconditioning);
    time1=MP_TIMER();

    /* GA_Print(xvec); */
    /* GA_Print(dvec); */
    
    if(me==0)printf("\n%d:in %d iterations time to solution=%f-%f\n",me,niter,(time1-time0),time_get);

    finalize_arrays(dopreconditioning);
    MP_BARRIER();

    if(me==0)printf("Terminating ..\n");
    GA_Terminate();
    MP_FINALIZE();
    return 0;
}


void finalize_arrays(int dpc)
{
     GA_Destroy(bvec);
     GA_Destroy(dvec);
     if(isvectormirrored)
        GA_Destroy(m_dvec);
     GA_Destroy(amat);
     GA_Destroy(xvec);
     GA_Destroy(axvec);
     GA_Destroy(rvec);
     GA_Destroy(qvec);
     GA_Destroy(ridx);
     GA_Destroy(cidx);
     if(dpc){
        GA_Destroy(svec);
        GA_Destroy(dmvec);
     }
     ARMCI_Free(myptrarrx[me]);
     ARMCI_Free(myptrarrd[me]);
}     
