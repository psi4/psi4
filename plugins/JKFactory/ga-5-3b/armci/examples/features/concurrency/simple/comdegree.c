#if HAVE_CONFIG_H
#   include "config.h"
#endif

/**                             Copyright (c) 2006
 *                      Pacific Northwest National Laboratory,
 *                           Battelle Memorial Institute.
 *                              All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met: - Redistributions of source code must retain the above
 * copyright notice, this list of conditions and the following disclaimer.
 * 
 * - Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the distribution.
 * - Neither the name of the Battelle nor the names of its contributors
 *   may be used to endorse or promote products derived from this software
 *   without specific prior written permission.
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * $Id: comdegree.c,v 1.1.2.1 2007-06-20 17:41:49 vinod Exp $
 *
 * This test checks the networks ability to overlap data transfers.
 * It does it both for an optimitic case (with no other communication) and
 * a more realistic case.
 * --Vinod Tipparaju
 * --Pacific Northwest National Laboratory
 * --vinod@pnl.gov
 */
#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_MATH_H
#   include <math.h>
#endif
#if HAVE_STRINGS_H
#   include <strings.h>
#endif
#if HAVE_UNISTD_H
#   include <unistd.h>
#endif
#if HAVE_ASSERT_H
#   include <assert.h>
#endif

#define DEBUG__ 

#include "armci.h"
#include "message.h"

/***************************** macros ************************/
#define COPY(src, dst, bytes) memcpy((dst),(src),(bytes))
#define ARMCI_MAX(a,b) (((a) >= (b)) ? (a) : (b))
#define ARMCI_MIN(a,b) (((a) <= (b)) ? (a) : (b))
#define ARMCI_ABS(a) (((a) <0) ? -(a) : (a))

/***************************** global data *******************/
int me, nproc;


void create_array(void *a[], int size)
{
     armci_size_t bytes=size;
     int rc;

     rc = ARMCI_Malloc(a, bytes);
     assert(rc==0);
     
#ifdef DEBUG_
     printf("%d after malloc ndim=%d b=%d ptr=%p\n",me,ndim,(int) bytes,a[me]);
     fflush(stdout);
#endif

     assert(a[me]);
     bzero(a[me],bytes);
}

void destroy_array(void *ptr[])
{
    armci_msg_barrier();

    assert(!ARMCI_Free(ptr[me]));
}



#define LCC 13
#define MAXPROC 128
    
void test_get_multidma()
{
int i,j;
void *b[MAXPROC], *a[MAXPROC];
int left = (me+nproc-1) % nproc;
int right = (me+1) % nproc;
/*int sendersright=0,sendersleft=0;*/
int loopcnt=10, itercount=5;
double tt/*, t0[LCC]*/,t1[LCC],t2[LCC],t3[LCC];
armci_hdl_t hdl1,hdl2;

    for(i=0;i<LCC;i++){
       /*t0[i]=0;*/
       t1[i]=0;
       t2[i]=0;
       t3[i]=0;
    }

    /* create shared and local arrays */
    create_array(b, 1024*1024*10);
    create_array(a, 1024*1024*10);
    /* warmup */
    ARMCI_INIT_HANDLE(&hdl1);
    ARMCI_INIT_HANDLE(&hdl2);
    ARMCI_NbGet((double*)b[left],(double*)a[me],1024,left,&hdl1);
    ARMCI_NbGet((double*)b[right]+1024,(double*)a[me]+1024,1024,
                         right,&hdl2);
    ARMCI_Wait(&hdl1);
    ARMCI_Wait(&hdl2);

    ARMCI_Barrier();

    /*start test*/
    for(j=0;j<itercount;j++){
       for(i=0;i<loopcnt;i++){
         int bytes;
        
         /*sendersright = (j+1)%nproc;*/
         /*sendersleft = (j+nproc-1)%nproc;*/

         bytes = 1024*pow(2,i);

         ARMCI_INIT_HANDLE(&hdl1);

         armci_msg_barrier();
         /*first time a regular call*/
         tt = armci_timer();
         ARMCI_NbGet((double*)b[left],(double*)a[me],bytes, left,&hdl1);
         ARMCI_Wait(&hdl1);
         t1[i] += (armci_timer()-tt);

         
         armci_msg_barrier();
         /*now time 1 left + 1 right but realize there is one xtra issue*/
         ARMCI_INIT_HANDLE(&hdl1);
         ARMCI_INIT_HANDLE(&hdl2);
         tt = armci_timer();
         ARMCI_NbGet((double*)b[left],(double*)a[me],bytes/2,left,&hdl1);
         ARMCI_NbGet((double*)b[right]+bytes/16,(double*)a[me]+bytes/16,bytes/2,
                         right,&hdl2);
         ARMCI_Wait(&hdl1);
         ARMCI_Wait(&hdl2);
         t2[i] += (armci_timer()-tt);

         ARMCI_Barrier();
         armci_msg_barrier();
         /*now time both to the left*/
         ARMCI_INIT_HANDLE(&hdl1);
         ARMCI_INIT_HANDLE(&hdl2);
         tt = armci_timer();
         ARMCI_NbGet((double*)b[left],(double*)a[me],bytes/2,left,&hdl1);
         ARMCI_NbGet((double*)b[left]+bytes/16,(double*)a[me]+bytes/16,bytes/2,
                         left,&hdl2);
         ARMCI_Wait(&hdl1);
         ARMCI_Wait(&hdl2);
         t3[i] += ( armci_timer()-tt);

         ARMCI_Barrier();
       }
    }
    armci_msg_barrier();
    if(0==me){
       for(i=0;i<loopcnt;i++){
         fprintf(stderr,"\n%.0f\t%.2e\t%.2e\t%.2e",
                         1024.0*pow(2,i),t1[i]/loopcnt,t3[i]/loopcnt,
                         t2[i]/loopcnt);
       }
    }
    fflush(stdout);
    armci_msg_barrier();
    armci_msg_barrier();
    if((nproc-3)==me){
       for(i=0;i<loopcnt;i++){
         fprintf(stderr,"\n%.0f\t%.2e\t%.2e\t%.2e",
                         1024.0*pow(2,i),t1[i]/loopcnt,t3[i]/loopcnt,
                         t2[i]/loopcnt);
       }
    }
    fflush(stdout);
    armci_msg_barrier();
#if 0
    for(j=0;j<nproc;j++) {
       if(j==me){
         for(i=0;i<loopcnt;i++){
           printf("\n%d:size=%f onesnd=%.2e twosnd=%.2e twosnddiffdir=%.2e\n",
                           me,1024.0*pow(2,i),t1[i]/loopcnt,t3[i]/loopcnt,
                           t2[i]/loopcnt);
         }
         armci_msg_barrier();
       }
       else
         armci_msg_barrier();
    }
#endif


    ARMCI_Barrier();

    destroy_array(b);
    destroy_array(a);
}


void test_put_multidma()
{
int i,j;
void *b[MAXPROC], *a[MAXPROC];
int left = (me+nproc-1) % nproc;
int right = (me+1) % nproc;
/*int sendersright=0,sendersleft=0;*/
int loopcnt=LCC, itercount=1000;
double tt/*, t0[LCC]*/,t1[LCC],t2[LCC],t3[LCC];
armci_hdl_t hdl1,hdl2;


    /* create shared and local arrays */
    create_array(b, 1024*1024*10);
    create_array(a, 1024*1024*10);
    for(i=0;i<LCC;i++){
       /*t0[i]=0;*/
       t1[i]=0;
       t2[i]=0;
       t3[i]=0;
    }

    ARMCI_Barrier();
    for(j=0;j<itercount;j++){
       for(i=0;i<loopcnt;i++){
         int wc,bytes;
         /* int lc1,rc1,wc1; */
        
         /*sendersright = (j+1)%nproc;*/
         /*sendersleft = (j+nproc-1)%nproc;*/

         bytes = 1024*pow(2,i)/8;

         ARMCI_INIT_HANDLE(&hdl1);
         ARMCI_NbPut((double*)a[me]+1024,(double*)b[left]+1024,bytes,left,&hdl1);
         ARMCI_Wait(&hdl1);
         ARMCI_INIT_HANDLE(&hdl1);
         ARMCI_NbPut((double*)a[me]+1024,(double*)b[right]+1024,bytes,right,&hdl1);
         ARMCI_Wait(&hdl1);

         ARMCI_INIT_HANDLE(&hdl1);
         armci_msg_barrier();

         tt = armci_timer();
         ARMCI_NbPut((double*)a[me],(double*)b[left],bytes, left,&hdl1);
         ARMCI_Wait(&hdl1);
         t1[i] += (armci_timer()-tt);
         /* (void)armci_notify(left); */
         /* tt = armci_timer(); */
         /* (void)armci_notify_wait(right,&wc);  */
         /* t1[i] += (armci_timer()-tt); */

         ARMCI_INIT_HANDLE(&hdl1);
         ARMCI_INIT_HANDLE(&hdl2);
         armci_msg_barrier();

         tt = armci_timer();
         ARMCI_NbPut((double*)a[me],(double*)b[left],bytes,left,&hdl1);
         ARMCI_NbPut((double*)a[me],(double*)b[right],bytes,
                         right,&hdl2);
         ARMCI_Wait(&hdl1);
         t2[i] += (armci_timer()-tt);
         /* (void)armci_notify(left); */
         /* lc1=armci_notify(right); */
         /* tt = armci_timer(); */
         /* rc1 = armci_notify_wait(left,&wc1);  */
         /* (void)armci_notify_wait(right,&wc);  */
         /* t2[i] += (armci_timer()-tt); */
         /* ARMCI_Wait(&hdl1); */
         ARMCI_Wait(&hdl2);

         ARMCI_INIT_HANDLE(&hdl1);
         ARMCI_INIT_HANDLE(&hdl2);
         armci_msg_barrier();

         tt = armci_timer();
         ARMCI_NbPut((double*)a[me],(double*)b[left],bytes/2,left,&hdl1);
         ARMCI_NbPut((double*)a[me]+bytes/16,(double*)b[left]+bytes/16,bytes/2,
                         left,&hdl2);
         /* ARMCI_Wait(&hdl1); */
         /* ARMCI_Wait(&hdl2); */
         t3[i] += ( armci_timer()-tt);
         (void)armci_notify(left);
         tt = armci_timer();
         (void)armci_notify_wait(right,&wc); 
         t3[i] += ( armci_timer()-tt);
         ARMCI_Wait(&hdl1);
         ARMCI_Wait(&hdl2);

         ARMCI_Barrier();
       }
    }
    armci_msg_barrier();
    if(0==me){
       for(i=0;i<loopcnt;i++){
         fprintf(stderr,"\n%.0f\t%.2e\t%.2e\t%.2e",
                         128.0*pow(2,i),t1[i]/itercount,t3[i]/itercount,
                         t2[i]/itercount);
       }
    }
    fflush(stdout);
    fflush(stdout);
    armci_msg_barrier();
    armci_msg_barrier();
    if((nproc-1)==me){
       for(i=0;i<loopcnt;i++){
         fprintf(stderr,"\n%.0f\t%.2e\t%.2e\t%.2e",
                         128.0*pow(2,i),t1[i]/itercount,t3[i]/itercount,
                         t2[i]/itercount);
       }
    }
    fflush(stdout);
    armci_msg_barrier();
#if 0
    for(j=0;j<nproc;j++) {
       if(j==me){
         for(i=0;i<loopcnt;i++){
           printf("\n%d:size=%f onesnd=%.2e twosnd=%.2e twosnddiffdir=%.2e\n",
                           me,1024.0*pow(2,i),t1[i]/loopcnt,t3[i]/loopcnt,
                           t2[i]/loopcnt);
         }
         armci_msg_barrier();
       }
       else
         armci_msg_barrier();
    }
#endif


    ARMCI_Barrier();

    destroy_array(b);
    destroy_array(a);
}


int main(int argc, char* argv[])
{
    armci_msg_init(&argc, &argv);
    nproc = armci_msg_nproc();
    me = armci_msg_me();

    
    ARMCI_Init();

    armci_msg_barrier();
    if(me==0){
       printf("\nTesting transfer overlap with ARMCI put calls\n");
       printf("\nsize\tone-send\ttwo-sends\ttwo-sends-diff-dir\n");
       fflush(stdout);
       sleep(1);
    }
    armci_msg_barrier();
    test_put_multidma();
    armci_msg_barrier();
    if(me==0){
       printf("\nTesting transfer overlap with ARMCI get calls\n");
       printf("\nsize\tone-send\ttwo-sends\ttwo-sends-diff-dir\n");
       fflush(stdout);
       sleep(1);
    }
    armci_msg_barrier();
    test_get_multidma();
    if(me==0)printf("\n");

    ARMCI_Finalize();
    armci_msg_finalize();
    return(0);
}
