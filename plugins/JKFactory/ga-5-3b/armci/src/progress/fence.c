#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include "armcip.h"
#include "armci.h"
#include "copy.h"
#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if defined(PVM)
#   include <pvm3.h>
#elif defined(TCGMSG)
#   include <sndrcv.h>
static void tcg_synch(long type)
{
    long atype = type;

    SYNCH_(&atype);
}
#elif defined(BGML)
#   include "bgml.h"
#else
#   include <mpi.h>
#endif

char *_armci_fence_arr;

void armci_init_fence()
{
#if defined (DATA_SERVER) || defined(PORTALS)
#if defined(THREAD_SAFE)
     _armci_fence_arr = calloc(armci_nproc*armci_user_threads.max,1);
#else
     _armci_fence_arr=calloc(armci_nproc,1);
#endif
     if(!_armci_fence_arr)
         armci_die("armci_init_fence: calloc failed",0);
#endif
}

void armci_finalize_fence()
{
#if defined (DATA_SERVER) || defined(PORTALS)
     free(_armci_fence_arr);
     _armci_fence_arr = NULL;
#endif
}

#ifdef PORTALS
void armci_update_fence_array(int proc, int inc)
{
    if (inc)
        FENCE_ARR(proc)++;
    else
        FENCE_ARR(proc)--;
}
#endif


void PARMCI_Fence(int proc)
{
#if defined(DATA_SERVER) && !(defined(GM) && defined(ACK_FENCE))
     if(FENCE_ARR(proc) && (armci_nclus >1)){

           int cluster = armci_clus_id(proc);
           int master = armci_clus_info[cluster].master;

           armci_rem_ack(cluster);

           bzero(&FENCE_ARR(master),
                   armci_clus_info[cluster].nslave);
     }
#elif defined(ARMCIX)
     ARMCIX_Fence (proc);
#elif defined(BGML)
     BGML_WaitProc(proc);
     MEM_FENCE;
#else
     FENCE_NODE(proc);
     MEM_FENCE;
#endif
}


void PARMCI_AllFence()
{
#if defined(ARMCIX)
    ARMCIX_AllFence ();
#elif defined(BGML)
    BGML_WaitAll();
#elif defined(LAPI) || defined(CLUSTER)
    int p;

    for(p = 0;p < armci_nproc; p++) {
        PARMCI_Fence(p); 
    }
#endif
    MEM_FENCE;
}

void PARMCI_Barrier()
{
    if (armci_nproc==1)
        return;
#if defined(BGML)
    BGML_WaitAll();
    bgml_barrier(3);
#else
    PARMCI_AllFence();
#  ifdef MPI
    MPI_Barrier(ARMCI_COMM_WORLD);
#  else
    {
       long type=ARMCI_TAG;
       tcg_synch(type);
    }
#  endif
#endif
    MEM_FENCE;
}
