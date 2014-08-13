#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Id: fence.c,v 1.25.4.6 2007-08-30 19:17:02 manoj Exp $ */
#include "armcip.h"
#include "armci.h"
#include "copy.h"
#include <stdio.h>
#if defined(PVM)
#   include <pvm3.h>
#elif defined(TCGMSG)
#   include <tcgmsg.h>
#elif defined(BGML)
#   include "bgml.h"
#else
#   include <mpi.h>
#endif

char *_armci_fence_arr;

void armci_init_fence()
{
#if defined (DATA_SERVER)
     _armci_fence_arr=calloc(armci_nproc,1);
     if(!_armci_fence_arr)armci_die("armci_init_fence: calloc failed",0);
#endif
}

void ARMCI_DoFence(int proc)
{
int i;
  if(!SAMECLUSNODE(proc) && (armci_nclus >1)){
  int cluster = armci_clus_id(proc);
    armci_rem_ack(cluster);
  }
}

void PARMCI_Fence(int proc)
{
int i;

#if defined(DATA_SERVER) && !(defined(GM) && defined(ACK_FENCE))
//   printf("%d [cp] fence_arr(%d)=%d\n",armci_me,proc,FENCE_ARR(proc));
     if(FENCE_ARR(proc) && (armci_nclus >1)){

           int cluster = armci_clus_id(proc);
           int master=armci_clus_info[cluster].master;

           armci_rem_ack(cluster);

           /* one ack per cluster node suffices */
           /* note, in multi-threaded case it will only clear for current thread */
           bzero(&FENCE_ARR(master),armci_clus_info[cluster].nslave);
     }
#elif defined(BGML)
     BGML_WaitProc(proc);
     MEM_FENCE;
#else
     FENCE_NODE(proc);
     MEM_FENCE;
#endif
}


/*
    portals developers' note:
    armci fence is not guaranteed to be correct unless PUT_START events are captured
    PUT_ENDs do NOT guarantee order; only PUT_STARTs
*/
void PARMCI_AllFence()
{
#if defined(CLUSTER)
     { int p; for(p=0;p<armci_nproc;p++)PARMCI_Fence(p); }
#endif
     MEM_FENCE;
}

void PARMCI_Barrier()
{
    if(armci_nproc==1)return;
    PARMCI_AllFence();
#  ifdef MPI
    MPI_Barrier(ARMCI_COMM_WORLD);
#  else
    {
       long type=ARMCI_TAG;
       tcg_synch(type);
    }
#  endif
    MEM_FENCE;
}
