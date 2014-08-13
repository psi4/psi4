#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <mpi.h>

extern void exit(int status);

#include "tcgmsgP.h"
#include "srftoc.h"
#include "armci.h"

char     tcgmsg_err_string[ERR_STR_LEN];
MPI_Comm TCGMSG_Comm=MPI_COMM_WORLD;
int      _tcg_initialized=0;
long  DEBUG_;
int      SR_parallel; 
int      SR_single_cluster =1;
int      tcgi_argc=0;
char   **tcgi_argv=NULL;

static int SR_initialized=0;
extern int nxtval_installed;

extern MPI_Comm armci_group_comm(ARMCI_Group *group);

long TCGREADY_()
{
    return (long)SR_initialized;
}


/**
 * number of processes
 */
long NNODES_()
{
    int numprocs;

    MPI_Comm_size(TCGMSG_Comm, &numprocs);
#ifdef NXTVAL_SERVER
    if(SR_parallel) {
        return((long)numprocs-1);
    }
#endif
    return((long)numprocs);
}


/**
 * Get calling process id
 */
long NODEID_()
{
    int myid;

    MPI_Comm_rank(TCGMSG_Comm,&myid);
    return((long)myid);
}


void Error(char *string, long code)
{
    fprintf(stdout, FMT_INT ": %s " FMT_INT " (%#lx).\n",
            NODEID_(), string, code, (long unsigned int)code);
    fflush(stdout);
    fprintf(stderr, FMT_INT ": %s " FMT_INT " (%#lx).\n",
            NODEID_(), string, code, (long unsigned int)code);

    finalize_nxtval(); /* clean nxtval resources */
    MPI_Abort(TCGMSG_Comm,(int)code);
}


/**
 * this is based on the MPI Forum decision that MPI_COMM_WORLD is a C constant 
 */
void make_tcgmsg_comm()
{
    extern int single_cluster();

#ifdef NXTVAL_SERVER
    if( SR_parallel ){   
        /* data server for a single process */
        int server;
        MPI_Group MPI_GROUP_WORLD, tcgmsg_grp;

        MPI_Comm_size(MPI_COMM_WORLD, &server);
        server --; /* the highest numbered process will be excluded */
        MPI_Comm_group(MPI_COMM_WORLD, &MPI_GROUP_WORLD);
        MPI_Group_excl(MPI_GROUP_WORLD, 1, &server, &tcgmsg_grp); 
        MPI_Comm_create(MPI_COMM_WORLD, tcgmsg_grp, &TCGMSG_Comm); 
    }else
#endif
    {
#if HAVE_ARMCI_INITIALIZED
        if (ARMCI_Initialized())
#else
        if (nxtval_installed)
#endif
        {
            ARMCI_Group group;
            ARMCI_Group_get_world(&group);
#   if HAVE_ARMCI_GROUP_COMM_MEMBER
            TCGMSG_Comm = group.comm;
#   else
            TCGMSG_Comm = armci_group_comm(&group);
#   endif
        }
        else
        {
            TCGMSG_Comm = MPI_COMM_WORLD;
        }
    }
}


/**
 * Alternative initialization for C programs
 * used to address argv/argc manipulation in MPI
 */
void tcgi_alt_pbegin(int *argc, char **argv[])
{
    int numprocs, myid;
    int init=0;

    if(SR_initialized) {
        Error("TCGMSG initialized already???",-1);
    } else {
        SR_initialized=1;
        _tcg_initialized=1;
    }

    /* check if another library initialized MPI already */
    MPI_Initialized(&init);

    if(!init){ 
        /* nope */
#if defined(DCMF) || defined(MPI_MT)
        int provided;
        MPI_Init_thread(argc, argv, MPI_THREAD_MULTIPLE, &provided);
#else
        MPI_Init(argc, argv);
#endif
        MPI_Errhandler_set(TCGMSG_Comm, MPI_ERRORS_RETURN);
    }

    MPI_Comm_size(TCGMSG_Comm, &numprocs);
    MPI_Comm_rank(TCGMSG_Comm, &myid);
    SR_parallel = numprocs > 1 ? 1 : 0;

#if NEED_DELAY_TCGMSG_MPI_STARTUP
    /* printf("%d:ready to go\n",NODEID_()); */
    /* wait until the last possible moment to call install_nxtval
     * it could be called by ARMCI_Init
     * or is called the first time nxtval is invoked (yuck) */
    /*install_nxtval(argc, argv);*/
#else
    install_nxtval(argc, argv); /* which calls ARMCI_Init(), if needed  */
#endif

    make_tcgmsg_comm();
    MPI_Barrier(TCGMSG_Comm);
}


/**
 * Initialization for C programs
 */
void tcgi_pbegin(int argc, char* argv[])
{
    tcgi_argc = argc;
    tcgi_argv = argv;
    tcgi_alt_pbegin(&argc, &argv);
}


/**
 * shut down message-passing library
 */ 
void PEND_()
{
#ifdef NXTVAL_SERVER
    long zero=0;
    if( SR_parallel ) {
        (void) NXTVAL_(&zero);
    }
    MPI_Barrier(TCGMSG_Comm);
#endif
    finalize_nxtval();
    MPI_Finalize();
    exit(0);
}


double TCGTIME_()
{
    static int first_call = 1;
    static double first_time, last_time, cur_time;
    double diff;

    if (first_call) {
        first_time = MPI_Wtime();
        first_call = 0;
        last_time  = -1e-9; 
    }

    cur_time = MPI_Wtime();
    diff = cur_time - first_time;

    /* address crappy MPI_Wtime: consectutive calls must be at least 1ns apart  */
    if(diff - last_time < 1e-9) {
        diff +=1e-9;
    }
    last_time = diff;

    return diff;                  /* Add logic here for clock wrap */
}


long MTIME_()
{
    return (long) (TCGTIME_()*100.0); /* time in centiseconds */
}



/**
 * longerface from Fortran to C error routine
 */
void PARERR_(long *code)
{
    Error("User detected error in FORTRAN", *code);
}


void SETDBG_(long *onoff)
{
    DEBUG_ = *onoff;
}

void STATS_()
{
    printf("STATS not implemented\n");
} 
