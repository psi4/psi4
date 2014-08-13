#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Id: simple.c,v 1.1.2.1 2007-06-20 17:41:40 vinod Exp $ */

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_UNISTD_H
#   include <unistd.h>
#elif HAVE_WINDOWS_H
#   include <windows.h>
#   define sleep(x) Sleep(1000*(x))
#endif
#if HAVE_ASSERT_H
#   include <assert.h>
#endif

#include "armci.h"
#include "message.h"

#define DIM1 5
#define DIM2 3
#ifdef __sun
/* Solaris has shared memory shortages in the default system configuration */
# define DIM3 6
# define DIM4 5
# define DIM5 4
#elif defined(__alpha__)
# define DIM3 8
# define DIM4 5
# define DIM5 6
#else
# define DIM3 8
# define DIM4 9
# define DIM5 7
#endif
#define DIM6 3
#define DIM7 2


#define OFF 1
#define EDIM1 (DIM1+OFF)
#define EDIM2 (DIM2+OFF)
#define EDIM3 (DIM3+OFF)
#define EDIM4 (DIM4+OFF)
#define EDIM5 (DIM5+OFF)
#define EDIM6 (DIM6+OFF)
#define EDIM7 (DIM7+OFF)

#define DIMS 4
#define MAXDIMS 7
#define MAX_DIM_VAL 50 
#define LOOP 200

#define BASE 100.
#define MAXPROC 128
#define TIMES 100

#ifdef CRAY
# define ELEMS 800
#else
# define ELEMS 200
#endif


/***************************** macros ************************/
#define COPY(src, dst, bytes) memcpy((dst),(src),(bytes))
#define ARMCI_MAX(a,b) (((a) >= (b)) ? (a) : (b))
#define ARMCI_MIN(a,b) (((a) <= (b)) ? (a) : (b))
#define ARMCI_ABS(a) (((a) <0) ? -(a) : (a))

/***************************** global data *******************/
int me, nproc;
void* work[MAXPROC]; /* work array for propagating addresses */



#ifdef MSG_COMMS_PVM
void pvm_init(int argc, char *argv[])
{
    int mytid, mygid, ctid[MAXPROC];
    int np, i;

    mytid = pvm_mytid();
    if((argc != 2) && (argc != 1)) goto usage;
    if(argc == 1) np = 1;
    if(argc == 2)
        if((np = atoi(argv[1])) < 1) goto usage;
    if(np > MAXPROC) goto usage;

    mygid = pvm_joingroup(MPGROUP);

    if(np > 1)
        if (mygid == 0) 
            i = pvm_spawn(argv[0], argv+1, 0, "", np-1, ctid);

    while(pvm_gsize(MPGROUP) < np) sleep(1);

    /* sync */
    pvm_barrier(MPGROUP, np);
    
    printf("PVM initialization done!\n");
    
    return;

usage:
    fprintf(stderr, "usage: %s <nproc>\n", argv[0]);
    pvm_exit();
    exit(-1);
}
#endif
          
void create_array(void *a[], int elem_size, int ndim, int dims[])
{
     int bytes=elem_size, i, rc;

     assert(ndim<=MAXDIMS);
     for(i=0;i<ndim;i++)bytes*=dims[i];

     rc = ARMCI_Malloc(a, bytes);
     assert(rc==0);
     
     assert(a[me]);
     
}

void destroy_array(void *ptr[])
{
    armci_msg_barrier();

    assert(!ARMCI_Free(ptr[me]));
}

#define MAXELEMS      1000
#define MAX_REQUESTS  MAXELEMS

void test_aggregate(int dryrun) {
  
    int i, j, rc, bytes, elems[2] = {MAXPROC, MAXELEMS};
    double *ddst_put[MAXPROC];
    double *ddst_get[MAXPROC];
    double *dsrc[MAXPROC];
    armci_hdl_t aggr_hdl_put[MAXPROC];
    armci_hdl_t aggr_hdl_get[MAXPROC];
    armci_hdl_t hdl_put[MAXELEMS];
    armci_hdl_t hdl_get[MAXELEMS];
    armci_giov_t darr;
    void *src_ptr[MAX_REQUESTS], *dst_ptr[MAX_REQUESTS];
    int start = 0, end = 0;
    double start_time;
        
    create_array((void**)ddst_put, sizeof(double),2, elems);
    create_array((void**)ddst_get, sizeof(double),2, elems);
    create_array((void**)dsrc, sizeof(double),1, &elems[1]);
    
    for(i=0; i<elems[1]; i++) dsrc[me][i]=i*1.001*(me+1);
    for(i=0; i<elems[0]*elems[1]; i++) {
      ddst_put[me][i]=0.0;
      ddst_get[me][i]=0.0;
    }
    
    armci_msg_barrier();

    /* only proc 0 does the work */
    if(me == 0) {
      if(!dryrun)printf("Transferring %d doubles (Not an array of %d doubles)\n", MAXELEMS, MAXELEMS);
      
      /* initializing non-blocking handles */
      for(i=0; i<elems[1]; i++) ARMCI_INIT_HANDLE(&hdl_put[i]);
      for(i=0; i<elems[1]; i++) ARMCI_INIT_HANDLE(&hdl_get[i]);
      
      /* aggregate handles */
      for(i=0; i<nproc; i++) ARMCI_INIT_HANDLE(&aggr_hdl_put[i]);
      for(i=0; i<nproc; i++) ARMCI_INIT_HANDLE(&aggr_hdl_get[i]);
      for(i=0; i<nproc; i++) ARMCI_SET_AGGREGATE_HANDLE(&aggr_hdl_put[i]);
      for(i=0; i<nproc; i++) ARMCI_SET_AGGREGATE_HANDLE(&aggr_hdl_get[i]);    
      
      bytes = sizeof(double);
      
      /* **************** PUT **************** */    
      /* register put */
      start_time=armci_timer();
      start = 0; end = elems[1]; 
      for(i=1; i<nproc; i++) {
    for(j=start; j<end; j++) {  
      ARMCI_NbPutValueDouble(dsrc[me][j], &ddst_put[i][me*elems[1]+j], i, 
                 &hdl_put[j]);
    }
    for(j=start; j<end; j++) ARMCI_Wait(&hdl_put[j]);
      }
      if(!dryrun)printf("%d: Value Put time      = %.2es\n", me, armci_timer()-start_time);
 
      /* vector put */
      start_time=armci_timer();
      for(i=1; i<nproc; i++) {
    for(j=start; j<end; j++) {
      src_ptr[j] = (void *)&dsrc[me][j];
      dst_ptr[j] = (void *)&ddst_put[i][me*elems[1]+j];
    }
    darr.src_ptr_array = src_ptr;
    darr.dst_ptr_array = dst_ptr;
    darr.bytes = sizeof(double);
    darr.ptr_array_len = elems[1];
    if((rc=ARMCI_NbPutV(&darr, 1, i, &hdl_put[i])))
      ARMCI_Error("armci_nbputv failed\n",rc);
      }
      for(i=1; i<nproc; i++) ARMCI_Wait(&hdl_put[i]);
      if(!dryrun)printf("%d: Vector Put time     = %.2es\n", me, armci_timer()-start_time);
      
      /* regular put */
      start_time=armci_timer();    
      for(i=1; i<nproc; i++) {
    for(j=start; j<end; j++) {  
      if((rc=ARMCI_NbPut(&dsrc[me][j], &ddst_put[i][me*elems[1]+j], bytes,
                 i, &hdl_put[j])))
        ARMCI_Error("armci_nbput failed\n",rc);
    }
    for(j=start; j<end; j++) ARMCI_Wait(&hdl_put[j]);
      }
      if(!dryrun)printf("%d: Regular Put time    = %.2es\n", me, armci_timer()-start_time);
      
      /* aggregate put */
      start_time=armci_timer();
      for(i=1; i<nproc; i++) {
    for(j=start; j<end; j++) {  
      if((rc=ARMCI_NbPut(&dsrc[me][j], &ddst_put[i][me*elems[1]+j], bytes,
                 i,  &aggr_hdl_put[i])))
        ARMCI_Error("armci_nbput failed\n",rc);
    }
      }
      for(i=1; i<nproc; i++) ARMCI_Wait(&aggr_hdl_put[i]);
      if(!dryrun)printf("%d: Aggregate Put time  = %.2es\n\n", me, armci_timer()-start_time);
      
      
      /* **************** GET **************** */    
      
      /* vector get */
      start_time=armci_timer();
      for(i=1; i<nproc; i++) {
    for(j=start; j<end; j++) {
      src_ptr[j] = (void *)&dsrc[i][j];
      dst_ptr[j] = (void *)&ddst_get[me][i*elems[1]+j];
    }
    darr.src_ptr_array = src_ptr;
    darr.dst_ptr_array = dst_ptr;
    darr.bytes = sizeof(double);
    darr.ptr_array_len = elems[1];
    if((rc=ARMCI_NbGetV(&darr, 1, i, &hdl_get[i])))
      ARMCI_Error("armci_nbgetv failed\n",rc);
    ARMCI_Wait(&hdl_get[i]);
      }
      if(!dryrun)printf("%d: Vector Get time     = %.2es\n", me, armci_timer()-start_time);
      
      /* regular get */
      start_time=armci_timer();    
      for(i=1; i<nproc; i++) {
    for(j=start; j<end; j++) {  
      if((rc=ARMCI_NbGet(&dsrc[i][j], &ddst_get[me][i*elems[1]+j], bytes,
                 i, &hdl_get[j])))
        ARMCI_Error("armci_nbget failed\n",rc);
    }
    for(j=start; j<end; j++) ARMCI_Wait(&hdl_get[j]);
      }
      if(!dryrun)printf("%d: Regular Get time    = %.2es\n", me, armci_timer()-start_time);
      
      /* aggregate get */
      start_time=armci_timer();
      for(i=1; i<nproc; i++) {
    for(j=start; j<end; j++) {  
      ARMCI_NbGet(&dsrc[i][j], &ddst_get[me][i*elems[1]+j], bytes,
              i, &aggr_hdl_get[i]);
    }
      }
      for(i=1; i<nproc; i++) ARMCI_Wait(&aggr_hdl_get[i]);
      if(!dryrun)printf("%d: Aggregate Get time  = %.2es\n", me, armci_timer()-start_time);
    }

    armci_msg_barrier();
    ARMCI_AllFence();
    armci_msg_barrier();

    /* Verify */
    if(!(me==0))
      for(j=0; j<elems[1]; j++) {
    if( ARMCI_ABS(ddst_put[me][j]-j*1.001) > 0.1) {
      ARMCI_Error("aggregate put failed...1", 0);
    }
      }
    armci_msg_barrier();
    if(!dryrun)if(me==0) printf("\n  aggregate put ..O.K.\n"); fflush(stdout);

    if(me==0) {
      for(i=1; i<nproc; i++) {
    for(j=0; j<elems[1]; j++) {
      if( ARMCI_ABS(ddst_get[me][i*elems[1]+j]-j*1.001*(i+1)) > 0.1) {
        ARMCI_Error("aggregate get failed...1", 0);
      }
    }
      }
    }
    armci_msg_barrier();
    if(!dryrun)if(me==0) printf("  aggregate get ..O.K.\n"); fflush(stdout);


    ARMCI_AllFence();
    armci_msg_barrier();
    
    if(!dryrun)if(me==0){printf("O.K.\n"); fflush(stdout);}
    destroy_array((void **)ddst_put);
    destroy_array((void **)ddst_get);
    destroy_array((void **)dsrc);
}


int main(int argc, char* argv[])
{

    armci_msg_init(&argc, &argv);
    nproc = armci_msg_nproc();
    me = armci_msg_me();

/*    printf("nproc = %d, me = %d\n", nproc, me);*/
    
    if(nproc>MAXPROC && me==0)
       ARMCI_Error("Test works for up to %d processors\n",MAXPROC);

    if(me==0){
       printf("ARMCI test program (%d processes)\n",nproc); 
       fflush(stdout);
       sleep(1);
    }
    
    ARMCI_Init();

    if(me==0){
      printf("\nAggregate put/get requests\n\n");
      fflush(stdout);
    }
    test_aggregate(1); /* cold start */
    test_aggregate(0); /* warm start */
    
    ARMCI_AllFence();
    armci_msg_barrier();
    if(me==0){printf("\nSuccess!!\n"); fflush(stdout);}
    sleep(2);
    
    armci_msg_barrier();
    ARMCI_Finalize();
    armci_msg_finalize();
    return(0);
}
