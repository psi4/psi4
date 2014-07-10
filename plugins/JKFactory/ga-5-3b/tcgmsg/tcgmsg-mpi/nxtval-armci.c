#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <mpi.h>

#if HAVE_MALLOC_H
#   include <malloc.h>
#endif

#include "armci.h"
#include "message.h"
#include "tcgmsgP.h"

#define LEN 2
static long pnxtval_counter_val;
static long *pnxtval_counter=&pnxtval_counter_val;
int nxtval_installed=0;
extern int     *tcgi_argc;
extern char  ***tcgi_argv;
#define INCR 1   /**< increment for NXTVAL */
#define BUSY -1L /**< indicates somebody else updating counter*/
#define NXTV_SERVER ((int)NNODES_() -1)
static int ARMCIinitialized = 0;

extern void make_tcgmsg_comm(void);

/**
 *  Get next value of shared counter.
 *  mproc > 0 ... returns requested value
 *  mproc < 0 ... server blocks until abs(mproc) processes are queued
 *                and returns junk
 *  mproc = 0 ... indicates to server that I am about to terminate
 */
long NXTVAL_(long *mproc)
{
    long local;
    int rc;
    int server = NXTV_SERVER;         /* id of server process */

    install_nxtval(tcgi_argc, tcgi_argv);

    if (SR_parallel) {
        if (DEBUG_) {
            (void) printf(FMT_INT ": nxtval: mproc=" FMT_INT "\n",
                          NODEID_(), *mproc);
            (void) fflush(stdout);
        }

        if (*mproc < 0) {
            rc=MPI_Barrier(TCGMSG_Comm); 
            if(rc!=MPI_SUCCESS) {
                Error("nxtval: barrier failed",0);
            }

            /* reset the counter value to zero */
            if( NODEID_() == server) {
                *pnxtval_counter = 0;
            }

            rc=MPI_Barrier(TCGMSG_Comm); 
            if(rc!=MPI_SUCCESS) {
                Error("nxtval: barrier failed",0);
            }
        }
        if (*mproc > 0) {
#if   SIZEOF_F77_INTEGER == SIZEOF_INT
            int op = ARMCI_FETCH_AND_ADD;
#elif SIZEOF_F77_INTEGER == SIZEOF_LONG
            int op = ARMCI_FETCH_AND_ADD_LONG;
#else
#   error
#endif
            rc = ARMCI_Rmw(op,(void*)&local,(void*)pnxtval_counter,1,server);
        }
    } else {
        /* Not running in parallel ... just do a simulation */
        static int count = 0;
        if (*mproc == 1) {
            local = count++;
        } else if (*mproc == -1) {
            count = 0;
            local = 0;
        } else {
            Error("nxtval: sequential version with silly mproc ", (long) *mproc);
        }
    }

    return local;
}


/**
 * initialization for nxtval
 */
void install_nxtval(int *argc, char **argv[])
{
    int rc;
    int me = (int)NODEID_(), bytes, server;
    void **ptr_ar;

    if (nxtval_installed) {
        return;
    }
    nxtval_installed = 1;

#if HAVE_ARMCI_INITIALIZED_FUNCTION
    if (!ARMCI_Initialized())
#else
    if (!ARMCIinitialized)
#endif
    {
        ARMCI_Init_args(argc, argv);
        ARMCIinitialized = 1;
        make_tcgmsg_comm();
    }

    ptr_ar = (void **)malloc(sizeof(void *)*(int)NNODES_());
    if(!ptr_ar) {
        Error("malloc failed in install_nxtval", (long)NNODES_());  
    }

    server = NXTV_SERVER;

    if(me== server) {
        bytes = sizeof(long);
    } else {
        bytes =0;
    }

    rc = ARMCI_Malloc(ptr_ar,bytes);
    if(rc) {
        Error("nxtv: armci_malloc failed",rc);
    }

    pnxtval_counter = (long*) ptr_ar[server];

    if(me==server) {
        *pnxtval_counter = (long)0;
    }

    free(ptr_ar);
    rc=MPI_Barrier(TCGMSG_Comm); 
    if(rc!=MPI_SUCCESS) {
        Error("init_nxtval: barrier failed",0);
    }
}


void finalize_nxtval()
{
/*
 * Cannot call ARMCI functions here as ARMCI might have been terminated
 * by now. NOTE: finalize_nxtval is called in pend(), which is called after
 * GA_Terminate/ARMCI_Finalize.
 */    
#if 0
    if(NODEID_() == NXTV_SERVER)ARMCI_Free(pnxtval_counter);
#endif
    ARMCI_Finalize();
}
