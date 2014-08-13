#if HAVE_CONFIG_H
#   include "config.h"
#endif

/*$Id: computation_impact.c,v 1.1.2.1 2007-06-20 17:42:13 vinod Exp $*/

#ifndef _GNU_SOURCE
#   define _GNU_SOURCE
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_STRING_H
#   include <string.h>
#endif
#if HAVE_STRINGS_H
#   include <strings.h>
#endif
#if HAVE_UNISTD_H
#   include <unistd.h>
#endif

#include <mpi.h>
#include <sched.h>

#include "armci.h"
#include "gpc.h"

typedef struct {double real,imag;} DoubleComplex;

#define LOOP 100
int accloop;
int gpcwork_daxpy;
int gpcwork_dgemm;
int gpcwork_ddot;
int gpcwork_memcpy;

int me,nprocs;
int size=2100*1024*LOOP;
int count=2100*1024;
int stride_arr[2]={0,0};
double alpha=0.1;
double c_alpha=0.1;
char **myptrs;
double t0,t1,t2,t3;
double tmpbuf[LOOP*100];
double tmpbuf1[LOOP*100];
double tmpbuf2[LOOP*100000];
double tmpbuf3[LOOP*100000];
#define DGS 450
double dga[DGS][DGS];
double dgb[DGS][DGS];
double dgc[DGS][DGS];

int gpc_work_handler_daxpy(int to, int from, void *hdr,   int hlen,
                      void *data,  int dlen,
                      void *rhdr,  int rhlen, int *rhsize,
                      void *rdata, int rdlen, int *rdsize,
                      int rtype)
{
int *rem;
int tmp_loop;
int i,j;
int ONE=1;
int N=DGS;

     rem = (int*)ARMCI_Gpc_translate(*(void**)hdr,to,from);

     tmp_loop = *rem;

      t2=MPI_Wtime();
      for(j=0;j<tmp_loop*80;j++){
        alpha=alpha+j*rand();
        daxpy_(&N,&alpha,tmpbuf2,&ONE,tmpbuf3,&ONE);
      }
      t3=MPI_Wtime()-t2;
      printf("\n%d:Server_Compute_daxpy %d %f\n",me,tmp_loop,t3);
     return GPC_DONE;
}
int gpc_work_handler_memcpy(int to, int from, void *hdr,   int hlen,
                      void *data,  int dlen,
                      void *rhdr,  int rhlen, int *rhsize,
                      void *rdata, int rdlen, int *rdsize,
                      int rtype)
{
int *rem;
int i,j;
int tmp_loop;

     rem = (int*)ARMCI_Gpc_translate(*(void**)hdr,to,from);

     tmp_loop = *rem;
     printf("\n%d:tmp_loop=%d",me,tmp_loop);

      t2=MPI_Wtime();
      for(j=0;j<tmp_loop*1400;j++){
        memcpy(tmpbuf2,tmpbuf3,LOOP*10000);
      }
      t3=MPI_Wtime()-t2;
      printf("\n%d:Server_Compute_memcpy %d %f\n",me,tmp_loop,t3);
     return GPC_DONE;
}

int gpc_work_handler_dgemm(int to, int from, void *hdr,   int hlen,
                      void *data,  int dlen,
                      void *rhdr,  int rhlen, int *rhsize,
                      void *rdata, int rdlen, int *rdsize,
                      int rtype)
{
int i,j;
int m,n,k;
int *rem;
int tmp_loop;
char notr='n';
DoubleComplex ZERO;
     ZERO.real=0.;ZERO.imag=0.;
     rem = (int*)ARMCI_Gpc_translate(*(void**)hdr,to,from);
     m=n=k=DGS;
     tmp_loop = *rem;

      t2=MPI_Wtime();
      for(j=0;j<tmp_loop*15;j++){
         alpha=alpha+j*rand();
         dgemm_(&notr,&notr,&m,&n,&k,&alpha,dga,&m,dgb,&n,&ZERO,dgc,&k,1,1);
      }
      t3=MPI_Wtime()-t2;
      printf("\n%d:Server_Compute_dgemm %d %f\n",me,tmp_loop,t3);
     return GPC_DONE;
}

int gpc_work_handler_ddot(int to, int from, void *hdr,   int hlen,
                      void *data,  int dlen,
                      void *rhdr,  int rhlen, int *rhsize,
                      void *rdata, int rdlen, int *rdsize,
                      int rtype)
{
int *rem;
int i,j;
int tmp_loop;

    rem = (int*)ARMCI_Gpc_translate(*(void**)hdr,to,from);

    tmp_loop = *rem;

    t2=MPI_Wtime();
    for(j=0;j<tmp_loop*90;j++){
    }
    t3=MPI_Wtime()-t2;
    printf("\n%d:Server_Compute_Ddot %d %f\n",me,accloop,t3);
    return GPC_DONE;
}



double c_dga[DGS][DGS];
double c_dgb[DGS][DGS];
double c_dgc[DGS][DGS];

int main(int argc, char **argv)
{
int i,peer,j;
cpu_set_t mycpuid,new_mask;
char str[CPU_SETSIZE];
int rrr;
char cid[8];
extern char * cpuset_to_cstr(cpu_set_t *mask, char *str);
extern int cstr_to_cpuset(cpu_set_t *mask, const char* str);
gpc_hdl_t nbh;
char rheader[100];
int hlen, rhlen, rhsize;
int rdsize;
int rem;
void *header=&rem;
int locval=0;
void *loc=&locval;
int right;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&me);
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    if(nprocs<2){
      printf("\ncan run only on >=2 procs\n");
      MPI_Finalize();
      exit(1);
    }     
    right = (me+1)%nprocs; 
    hlen=sizeof(header);
    bzero(rheader,100);
    rhlen = hlen;

    ARMCI_Init();
    accloop=atoi(argv[1]);
    rem=accloop;
    myptrs = (char **)malloc(sizeof(char *)*nprocs);
    ARMCI_Malloc((void **)myptrs,size);

    MPI_Barrier(MPI_COMM_WORLD);

    gpcwork_memcpy = ARMCI_Gpc_register(gpc_work_handler_memcpy);
    gpcwork_ddot =ARMCI_Gpc_register(gpc_work_handler_ddot);
    gpcwork_daxpy = ARMCI_Gpc_register(gpc_work_handler_daxpy);
    gpcwork_dgemm = ARMCI_Gpc_register(gpc_work_handler_dgemm);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    ARMCI_Gpc_init_handle(&nbh);
    if(ARMCI_Gpc_exec(gpcwork_memcpy, right, &header, hlen, loc, sizeof(int), 
                            rheader, rhlen,loc, sizeof(int), &nbh))
       fprintf(stderr,"ARMCI_Gpc_exec failed\n");
    {
      int m,n,k;
      char notr='n';
      DoubleComplex ZERO;
      usleep(100);
      ZERO.real=0.;ZERO.imag=0.;
      m=n=k=DGS;
      t0=MPI_Wtime();
#ifdef DGEMM_WORK
      for(j=0;j<4*15;j++){
         c_alpha=c_alpha+j*rand();
         dgemm_(&notr,&notr,&m,&n,&k,&c_alpha,c_dga,&m,c_dgb,&n,&ZERO,c_dgc,&k,1,1);
      }
#elif IUNIT_WORK 
      for(j=0;j<2*LOOP*100;j++){
        for(i=0;i<LOOP*100;i++){
          tmpbuf1[i]=tmpbuf1[i]*1.1214+i/tmpbuf1[j/2];
        }
      }
#elif DAXPY_WORK
      for(j=0;j<tmp_loop*80;j++){
        alpha=alpha+j*rand();
        daxpy_(&N,&alpha,tmpbuf1,&ONE,tmpbuf2,&ONE);
      }
#endif

      t1=MPI_Wtime()-t0;
      printf("\n%d:Compute_During_memcpy %d %f\n",me,accloop,t1);
    }

    ARMCI_Gpc_wait(&nbh);
    
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    ARMCI_Gpc_init_handle(&nbh);
    if(ARMCI_Gpc_exec(gpcwork_ddot, right, &header, hlen, loc, sizeof(int), 
                            rheader, rhlen,loc, sizeof(int), &nbh))
       fprintf(stderr,"ARMCI_Gpc_exec failed\n");
    {
      int m,n,k;
      char notr='n';
      DoubleComplex ZERO;
      usleep(100);
      ZERO.real=0.;ZERO.imag=0.;
      m=n=k=DGS;
      t0=MPI_Wtime();
#ifdef DGEMM_WORK
      for(j=0;j<4*15;j++){
         c_alpha=c_alpha+j*rand();
         dgemm_(&notr,&notr,&m,&n,&k,&c_alpha,c_dga,&m,c_dgb,&n,&ZERO,c_dgc,&k,1,1);
      }
#elif IUNIT_WORK 
      for(j=0;j<2*LOOP*100;j++){
        for(i=0;i<LOOP*100;i++){
          tmpbuf1[i]=tmpbuf1[i]*1.1214+i/tmpbuf1[j/2];
        }
      }
#elif DAXPY_WORK
      for(j=0;j<tmp_loop*80;j++){
        alpha=alpha+j*rand();
        daxpy_(&N,&alpha,tmpbuf1,&ONE,tmpbuf2,&ONE);
      }
#endif
      t1=MPI_Wtime()-t0;
      printf("\n%d:Compute_During_Ddot %d %f\n",me,accloop,t1);
    }
    ARMCI_Gpc_wait(&nbh);

    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    ARMCI_Gpc_init_handle(&nbh);
    if(ARMCI_Gpc_exec(gpcwork_daxpy, right, &header, hlen, loc, sizeof(int), 
                            rheader, rhlen,loc, sizeof(int), &nbh))
       fprintf(stderr,"ARMCI_Gpc_exec failed\n");
    {
      int m,n,k;
      char notr='n';
      DoubleComplex ZERO;
      usleep(100);
      ZERO.real=0.;ZERO.imag=0.;
      m=n=k=DGS;
      t0=MPI_Wtime();
#ifdef DGEMM_WORK
      for(j=0;j<4*15;j++){
         c_alpha=c_alpha+j*rand();
         dgemm_(&notr,&notr,&m,&n,&k,&c_alpha,c_dga,&m,c_dgb,&n,&ZERO,c_dgc,&k,1,1);
      }
#elif IUNIT_WORK 
      for(j=0;j<2*LOOP*100;j++){
        for(i=0;i<LOOP*100;i++){
          tmpbuf1[i]=tmpbuf1[i]*1.1214+i/tmpbuf1[j/2];
        }
      }
#elif DAXPY_WORK
      for(j=0;j<tmp_loop*80;j++){
        alpha=alpha+j*rand();
        daxpy_(&N,&alpha,tmpbuf1,&ONE,tmpbuf2,&ONE);
      }
#endif
      t1=MPI_Wtime()-t0;
      printf("\n%d:Compute_During_Daxpy %d %f\n",me,accloop,t1);
    }
    ARMCI_Gpc_wait(&nbh);

    MPI_Barrier(MPI_COMM_WORLD);

    ARMCI_Gpc_init_handle(&nbh);
    if(ARMCI_Gpc_exec(gpcwork_dgemm, right, &header, hlen, loc, sizeof(int), 
                            rheader, rhlen,loc, sizeof(int), &nbh))
       fprintf(stderr,"ARMCI_Gpc_exec failed\n");
    {
      int m,n,k;
      char notr='n';
      DoubleComplex ZERO;
      usleep(100);
      ZERO.real=0.;ZERO.imag=0.;
      m=n=k=DGS;
      t0=MPI_Wtime();
#ifdef DGEMM_WORK
      for(j=0;j<4*15;j++){
         c_alpha=c_alpha+j*rand();
         dgemm_(&notr,&notr,&m,&n,&k,&c_alpha,c_dga,&m,c_dgb,&n,&ZERO,c_dgc,&k,1,1);
      }
#elif IUNIT_WORK 
      for(j=0;j<2*LOOP*100;j++){
        for(i=0;i<LOOP*100;i++){
          tmpbuf1[i]=tmpbuf1[i]*1.1214+i/tmpbuf1[j/2];
        }
      }
#elif DAXPY_WORK
      for(j=0;j<tmp_loop*80;j++){
        alpha=alpha+j*rand();
        daxpy_(&N,&alpha,tmpbuf1,&ONE,tmpbuf2,&ONE);
      }
#endif
      t1=MPI_Wtime()-t0;
      printf("\n%d:Compute_During_Dgemm %d %f\n",me,accloop,t1);
    }
    ARMCI_Gpc_wait(&nbh);

    MPI_Barrier(MPI_COMM_WORLD);

    ARMCI_AllFence();

    ARMCI_Finalize();
    MPI_Finalize();
}
