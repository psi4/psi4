#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Id$ */

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_UNISTD_H
#   include <unistd.h>
#endif

#include "armci.h"
#include "message.h"

int me,nprocs;
int LOOP=10;
int main(int argc, char **argv)
{
int i;
double **myptrs;
double t0,t1,tnbget=0,tnbwait=0,t2=0;
    armci_msg_init(&argc,&argv);
    nprocs = armci_msg_nproc();
    if (nprocs==1)
    {
        fprintf(stderr,"You must use more than 1 process for this test.  Exiting gently.");
        return 0;
    }
    me = armci_msg_me();
    myptrs = (double **)malloc(sizeof(double *)*nprocs);
    ARMCI_Init();
    ARMCI_Malloc((void **)myptrs,LOOP*sizeof(double)); 
    armci_msg_barrier();
    if(me==0){
       for(i=0;i<LOOP;i++){
         ARMCI_Get(myptrs[me+1]+i,myptrs[me]+i,sizeof(double),me+1);
       }
       t0 = armci_timer(); 
       for(i=0;i<LOOP;i++){
         ARMCI_Get(myptrs[me+1]+i,myptrs[me]+i,sizeof(double),me+1);
       }
       t1 = armci_timer(); 
       printf("\nGet Latency=%f\n",1e6*(t1-t0)/LOOP);fflush(stdout);
       t1=t0=0;
       for(i=0;i<LOOP;i++){
         armci_hdl_t nbh;
         ARMCI_INIT_HANDLE(&nbh);
         t0 = armci_timer(); 
         ARMCI_NbGet(myptrs[me+1]+i,myptrs[me]+i,sizeof(double),me+1,&nbh);
         t1 = armci_timer(); 
         ARMCI_Wait(&nbh);
         t2 = armci_timer();
         tnbget+=(t1-t0);
         tnbwait+=(t2-t1);
       }
       printf("\nNb Get Latency=%f Nb Wait=%f\n",1e6*tnbget/LOOP,1e6*tnbwait/LOOP);fflush(stdout);
    }
    else
      sleep(1);
    armci_msg_barrier();
    ARMCI_Finalize();
    armci_msg_finalize();
    return 0;
}
