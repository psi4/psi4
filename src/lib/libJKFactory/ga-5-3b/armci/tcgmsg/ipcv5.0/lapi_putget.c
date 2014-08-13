#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_SIGNAL_H
#   include <signal.h>
#endif
#if HAVE_UNISTD_H
#   include <unistd.h>
#endif

#include <lapi.h>
#include <pthread.h>

#include "tcgmsgP.h"

lapi_handle_t lapi_handle;
lapi_info_t   lapi_info;
extern ShmemBuf TCGMSG_receive_buffer[];

#define LEN 2
int nxtval_counter=0;
int *nxtval_cnt_adr = &nxtval_counter;
static lapi_cntr_t req_cnt;

#define INCR 1           /* increment for NXTVAL */
#define BUSY -1L         /* indicates somebody else updating counter*/
/*#define TRACEINFO 1*/


/**
 * initialize lapi
 */
void lapi_initialize()
{
    int myid, numtasks,rc;

    bzero(&lapi_info,sizeof(lapi_info)); /* needed under Mohonk */
    rc = LAPI_Init(&lapi_handle, &lapi_info);
    if(rc) Error("lapi_init failed",rc);

    rc=LAPI_Qenv(lapi_handle, TASK_ID, &myid);
    if(rc) Error("lapi_qenv failed",rc);
    rc=LAPI_Qenv(lapi_handle, NUM_TASKS, &numtasks);
    if(rc) Error("lapi_qenv failed 2",rc);

    TCGMSG_nodeid = (long)myid;
    TCGMSG_nnodes = (long)numtasks;

    /* disable LAPI internal error checking */
    LAPI_Senv(lapi_handle, ERROR_CHK, 0);

    rc = LAPI_Setcntr(lapi_handle, &req_cnt, 0);
    if(rc)Error("lapi_initialize: setcntr failed",rc);

#ifdef DEBUG
    printf("me=%d initialized %d processes\n", myid, numtasks);
#endif
    fflush(stdout);
}


void lapi_adr_exchg()
{
    long node, tgt;
    int rc;
    void **table;
    int i;

    table = (void **)malloc(TCGMSG_nnodes * sizeof(void *));
    if (!table) Error(" lapi_adr_exchg: malloc failed", 0);

    /* allocate and initialize send buffers */
    sendbuf_arr = (sendbuf_t*)malloc(SENDBUF_NUM*sizeof(sendbuf_t));
    if(!sendbuf_arr) Error(" lapi_adr_exchg:malloc 2 failed", 0);
    /*
       bzero(sendbuf_arr,SENDBUF_NUM*sizeof(sendbuf_t));
       */

    for(i=0; i< SENDBUF_NUM; i++){
        LAPI_Setcntr(lapi_handle,&sendbuf_arr[i].cntr, 1);
        sendbuf_arr[i].next = sendbuf_arr+i+1;
    }
    sendbuf_arr[SENDBUF_NUM-1].next = sendbuf_arr;
    localbuf = sendbuf_arr;
    if(sizeof(ShmemBuf) < sizeof(sendbuf_t))
        Error("lapi_adr_exchg: buffer size problem",0);

    /* exchange addresses */
    for(node = 0; node < TCGMSG_nnodes; node++){

        /* Lapi does not like NULL address for buffer that we have 
           for sending msg to itself - use some invalid address */
        if (node == TCGMSG_nodeid) 
            TCGMSG_proc_info[node].recvbuf = (ShmemBuf *)1; 
        else
            if(LAPI_Setcntr(lapi_handle,
                        &(TCGMSG_proc_info[node].recvbuf->cntr),0))
                Error("lapi_adr_exchg: setcntr failed",-1);

        rc = LAPI_Address_init(lapi_handle, TCGMSG_proc_info[node].recvbuf, 
                table); 
        if(rc) Error(" lapi_adr_exchg: address_init failed", node);

        if(rc) Error(" lapi_adr_exchg: cntr init failed", node);

        if(TCGMSG_nodeid == node) {
            for(tgt=0; tgt<TCGMSG_nnodes; tgt++){ 
/*
#ifdef DEBUG
                printf("%ld org=%lx final=%lx\n",
                        TCGMSG_nodeid, TCGMSG_proc_info[tgt].sendbuf,
                        table[tgt]);
                fflush(stdout);
#endif
*/
                TCGMSG_proc_info[tgt].sendbuf = (ShmemBuf *)table[tgt];
            }
        }
    }
    free(table);
#ifdef DEBUG
    printf("%ld: bufarray: (%lx, %lx) %ld bytes\n",
            TCGMSG_nodeid, TCGMSG_receive_buffer,
            TCGMSG_receive_buffer + MAX_PROC, sizeof(ShmemBuf)*MAX_PROC);
    fflush(stdout);
#endif
}


/** Error handler */
void Error(char *string, long code)
{
    char cmd[100];

    (void) fflush(stdout);

    (void) fprintf(stdout, "%3d:%s %ld(%x)\n", NODEID_(), string, code, code);
    (void) fflush(stdout);
    (void) fprintf(stderr, "%3d:%s %ld(%x)\n", NODEID_(), string, code, code);
    (void) perror("system message");

    (void) fflush(stdout);
    (void) fflush(stderr);

#ifdef TRACEINFO
    sprintf (cmd,"/u2/d3h325/bin/dbxwhere_pid %d ",getpid());
    if(-1 == system(cmd)) fprintf(stderr,"system call failed\n");
    sleep(1);
#endif
    if(-1 == kill(getpid(),SIGTERM)) fprintf(stderr,"kill failed\n");

    LAPI_Term(lapi_handle); /* JN: not sure if it should be called here */
    exit(1);
}


/**
 * Get next value of shared counter.
 *
 * mproc > 0 ... returns requested value
 * mproc < 0 ... server blocks until abs(mproc) processes are queued
 *               and returns junk
 * mproc = 0 ... indicates to server that I am about to terminate
 */
long NXTVAL_(long *mproc)
{
#define INC 1
    int local;
    long stype = INTERNAL_SYNC_TYPE;
    lapi_cntr_t req_id;
    int rc, inc = INC;

    int  server = (int)NNODES_() -1;         /* id of server process */

    if (server>0) { 
        /* parallel execution */
        if (DEBUG_) {
            (void) printf("%2ld: nxtval: mproc=%ld\n",NODEID_(), *mproc);
            (void) fflush(stdout);
        }

        if (*mproc < 0) {
            SYNCH_(&stype);
            /* reset the counter value to zero */
            if( NODEID_() == server) nxtval_counter = 0;
            SYNCH_(&stype);
        }
        if (*mproc > 0) {
            /* use atomic swap operation to increment nxtval counter */
            rc = LAPI_Setcntr(lapi_handle, &req_id, 0);
            if(rc)Error("nxtval: setcntr failed",rc);
            rc = LAPI_Rmw(lapi_handle, FETCH_AND_ADD, server, nxtval_cnt_adr,
                    &inc, &local, &req_id);
            if(rc)Error("nxtval: rmw failed",rc);
            rc = LAPI_Waitcntr(lapi_handle, &req_id, 1, NULL);
            if(rc)Error("nxtval: waitcntr failed",rc);
        }
    } else {
        /* Not running in parallel ... just do a simulation */
        static int count = 0;
        if (*mproc == 1){
            int val = count;
            count+=INCR;
            local = val;
        }else if (*mproc == -1) {
            count = 0;
            local = 0;
        }
        else
            Error("nxtval: sequential version with silly mproc ", (long) *mproc);
    }

    return (long)local;
}


/** blocking get */
void lapi_get(void* dest, void* src, long bytes, long node)
{
    int rc;

#ifdef DEBUG
    printf("%ld getting %ld bytes from addr=%lx node %ld to adr=%lx\n", 
            TCGMSG_nodeid, bytes, src, node, dest );
    fflush(stdout);
#endif

    rc = LAPI_Get(lapi_handle, (uint)node, (uint)bytes, src, dest, NULL,&req_cnt);
    if(rc)Error("lapi_get: get failed",rc);
    rc = LAPI_Waitcntr(lapi_handle, &req_cnt, 1, NULL);
    if(rc)Error("lapi_get: waitcntr failed",rc);
}


/** put with nonblocking semantics */
void lapi_put(void* dest, void* src, long bytes, long node)
{
    int rc;

    /*  LAPI_Fence(lapi_handle);*/
#ifdef DEBUG
    printf("%ld puting %ld bytes to addr=%lx node %ld\n", TCGMSG_nodeid,
            bytes, dest, node);
    fflush(stdout);
#endif
#ifdef ERR_CHECKING
    if(dest < (void*)TCGMSG_receive_buffer){
        printf("%ld: Warning: Out of range? %lx(%ld) < %lx\n",
                TCGMSG_nodeid, dest, node, TCGMSG_receive_buffer);
        fflush(stdout); 
    }
    if(dest + bytes > (void*)(TCGMSG_receive_buffer+MAX_PROC) ){
        printf("%ld: Warning: Out of range? %lx(%ld) < %lx\n",
                TCGMSG_nodeid, dest+bytes, node, TCGMSG_receive_buffer+MAX_PROC);
        fflush(stdout);
    }
#endif

    rc=LAPI_Put(lapi_handle, (uint)node, (uint)bytes, dest, src,NULL, &req_cnt,NULL);
    if(rc)Error("lapi_put: sdput failed",rc);
    rc = LAPI_Waitcntr(lapi_handle, &req_cnt, 1, NULL);
    if(rc)Error("lapi_put: waitcntr failed",rc);

}


/** put with nonblocking semantics and counter */
void lapi_put_c(void* dest, void* src, long bytes, long node, lapi_cntr_t* cntr)
{
    int rc;
    rc = LAPI_Put(lapi_handle, (uint)node, (uint)bytes, dest, src,cntr,NULL,NULL);
    if(rc)Error("lapi_put_c: put failed",rc);
}


void PBEGINF_()
{
    PBEGIN_(NULL,NULL);
}


double fred =0.;
void Busy(int n)
{
    while (n-- >= 0) fred++;
    /*   LAPI_Probe(lapi_handle); */
}


void SYNCH_(long* type)
{
    int rc;

    rc=LAPI_Gfence(lapi_handle);
    if(rc) Error("lapi_gfence failed",rc);
}


/** Interface from fortran to c error routine */
void PARERR_(long *code)
{
    Error("User detected error in FORTRAN", *code);
}
