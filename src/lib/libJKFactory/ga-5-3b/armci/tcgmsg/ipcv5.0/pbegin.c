#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STRINGS_H
#   include <strings.h>
#endif
#if HAVE_SYS_TYPES_H
#   include <sys/types.h>
#endif
#if HAVE_UNISTD_H
#   include <unistd.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif

/* extern long atol(const char *nptr); */
/* extern void exit(int status); */

/* Define PBEGIN_C so that global variables in tcgmsgP.h are defined here
   and declared extern everywhere else ... SGI linker is a whiner */
#define PBEGIN_C

#include "tcgmsgP.h"

#ifdef LAPI
ShmemBuf  TCGMSG_receive_buffer[MAX_PROC];
void lapi_initialize();
#endif

extern void TrapSigint(void);
extern void TrapSigchld(void);
extern int WaitAll(long);

static int SR_initialized=0;


long TCGREADY_()
{
    return (long)SR_initialized;
}


/* Define what was externally declared in tcgmsgP.h */
long        TCGMSG_nodeid;
long        TCGMSG_nnodes;
long        DEBUG_=0; /* debug flag ... see setdbg */
long        TCGMSG_nodeid;
long        TCGMSG_nnodes;
char*       TCGMSG_shmem;
long        TCGMSG_shmem_id;
long        TCGMSG_shmem_size;
long        TCGMSG_caught_sigint;
ProcInfo*   TCGMSG_proc_info;
SendQEntry* TCGMSG_sendq_ring;


/**
 * shared-memory version of TCGMSG
 */  
void tcgi_pbegin(int argc, char **argv)
{
    long arg, node, i, max_n_msg;

    TCGMSG_nodeid = 0;
    TCGMSG_nnodes = 1;        /* By default just sequential */

    if(SR_initialized)Error("TCGMSG initialized already???",-1);
    else SR_initialized=1;

#ifdef LAPI
    lapi_initialize();
#else /* LAPI */
    for (arg=1; arg<(argc-1); arg++)
        if (strcmp(argv[arg],"-np") == 0) {
            TCGMSG_nnodes = atol(argv[arg+1]);
            break;
        }
#endif /* LAPI */

    if (TCGMSG_nnodes > MAX_PROC){
        if(NODEID_()){
            sleep(1);
            return;
        }
        fprintf(stderr,"\nTCGMSG has been configured for up to %d processes\n",
                MAX_PROC);
        fprintf(stderr,"Please change MAX_PROC in `tcgmsgP.h` and recompile\n\n");
        sleep(1);
        Error("aborting ... ",0);
    }
    if (TCGMSG_nnodes == 1) {
        return;
    };

    /* Set up handler for SIGINT and SIGCHLD */

#ifndef LAPI
    TrapSigint();
    TrapSigchld();
#endif

    /* Allocate the process info structures */

    if (!(TCGMSG_proc_info = (ProcInfo *)
                malloc((size_t) (TCGMSG_nnodes*sizeof(ProcInfo)))))
        Error("pbegin: failed to malloc procinfo",
                (long) (TCGMSG_nnodes*sizeof(ProcInfo)));
    bzero((char *) TCGMSG_proc_info, (int) (TCGMSG_nnodes*sizeof(ProcInfo)));

    /* Allocate a ring of message q entries to avoid having a malloc/free
       pair for every message sent */

    max_n_msg = 2*TCGMSG_nnodes;
    if (max_n_msg < MAX_N_OUTSTANDING_MSG) max_n_msg = MAX_N_OUTSTANDING_MSG;

    if (!(TCGMSG_sendq_ring = (SendQEntry *)
                malloc((size_t) (max_n_msg*sizeof(SendQEntry)))))
        Error("pegin: failed to malloc entries for send q", 0L);

    for (i=0; i<max_n_msg; i++) {
        TCGMSG_sendq_ring[i].active = 0;
        TCGMSG_sendq_ring[i].next_in_ring = TCGMSG_sendq_ring + ((i+1)%max_n_msg);
    }

    /* Create the shared memory and fill with zeroes */

#ifndef LAPI
    TCGMSG_shmem_size = (long) (TCGMSG_nnodes * TCGMSG_nnodes * sizeof(ShmemBuf));
    TCGMSG_shmem_size += 64;
    TCGMSG_shmem = CreateSharedRegion(&TCGMSG_shmem_id, &TCGMSG_shmem_size);
    nxtval_shmem = (long*)(((char*)(TCGMSG_shmem))+TCGMSG_shmem_size-64); 
#else
    TCGMSG_shmem_size = (long)(TCGMSG_nnodes * sizeof(ShmemBuf));
    TCGMSG_shmem = (char *) TCGMSG_receive_buffer;
#endif

    bzero(TCGMSG_shmem, (int) TCGMSG_shmem_size);

    /* Fork the child processes */

    TCGMSG_proc_info[0].pid = getpid();

#ifndef LAPI
    for (node=1; node<TCGMSG_nnodes; node++) {
        pid_t pid = fork();

        if (pid == 0) {
            TCGMSG_nodeid = node;    /* Generate my unique id */
            break;
        }      
        else {
            TCGMSG_proc_info[node].pid = pid;
        }
    }
#endif


    /* Now everyone initializes the pointers to the shared-
       memory buffers.

       Each process has TCGMSG_nnodes buffers, one for each
       other process that it can receive from via shared
       memory (plus one extra!). */

    for (node=0; node<TCGMSG_nnodes; node++) {
        long me = TCGMSG_nodeid;
        if (me != node) {

#ifndef LAPI
            TCGMSG_proc_info[node].sendbuf = ((ShmemBuf *) TCGMSG_shmem) +
                (node*TCGMSG_nnodes + me);
            TCGMSG_proc_info[node].recvbuf = ((ShmemBuf *) TCGMSG_shmem) +
                (me*TCGMSG_nnodes + node);
#else
            TCGMSG_proc_info[node].sendbuf = ((ShmemBuf *) TCGMSG_shmem) + me;
            TCGMSG_proc_info[node].recvbuf = ((ShmemBuf *) TCGMSG_shmem) + node;
#endif
        }
    }

#ifdef LAPI
    lapi_adr_exchg();
#endif

    /* At this point communication is possible. 

       Synchronize and continue. */

    {
        long type = 1;

        SYNCH_(&type);
    }
}


void tcgi_alt_pbegin(int *argc, char **argv[])
{ 
    tcgi_pbegin(*argc, *argv);
} 


void PEND_(void)
{
    long type = 999;
#ifndef LAPI
    (void) signal(SIGCHLD, SIG_DFL); /* Death of children now OK */
#endif

    SYNCH_(&type);
    SR_initialized = 0;

#ifndef LAPI
    if (TCGMSG_nodeid == 0 && TCGMSG_nnodes > 1) {
        int status;
        int rc;
        status = WaitAll(TCGMSG_nnodes-1);       /* Wait for demise of children */
        rc=DeleteSharedRegion(TCGMSG_shmem_id);
        if(rc)printf("DeleteSharedMem returned %d\n",rc);
        if (status) exit(1);
    }
#endif
}
