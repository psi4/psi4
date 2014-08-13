#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Id: rmw.c,v 1.24.2.5 2007-08-29 17:32:47 manoj Exp $ */
#include "armcip.h"
#include "locks.h"
#include "copy.h"
#include <stdio.h>
#if (defined(__i386__) || defined(__x86_64__)) && !defined(_CRAYC)
#  include "atomics-i386.h"
#endif

#ifdef LIBELAN_ATOMICS 

ELAN_ATOMIC *a;

int elan_int_fadd(int *target, int inc, int vp)
{
    int result;

    elan_wait(elan_atomic32(a, ELAN_ATOMIC_ADD, target, inc, 0, vp, &result), elan_base->waitType);
    return(result);
}

int elan_long_fadd(long *target, long inc, int vp)
{
    long result;
    
#ifdef _LP64
    elan_wait(elan_atomic64(a, ELAN_ATOMIC_ADD, target, inc, 0, vp, &result), elan_base->waitType);
#else
    elan_wait(elan_atomic32(a, ELAN_ATOMIC_ADD, target, inc, 0, vp, &result), elan_base->waitType);
#endif

    return(result);
}

int elan_int_swap(int *target, int value, int vp)
{
    int result;

    elan_wait(elan_atomic32(a, ELAN_ATOMIC_SWAP, target, value, 0, vp, &result), elan_base->waitType);
    return(result);
}

int elan_long_swap(long *target, long value, int vp)
{
    long result;
    
#ifdef _LP64
    elan_wait(elan_atomic64(a, ELAN_ATOMIC_SWAP, target, value, 0, vp, &result), elan_base->waitType);
#else
    elan_wait(elan_atomic32(a, ELAN_ATOMIC_SWAP, target, value, 0, vp, &result), elan_base->waitType);
#endif

    return(result);
}
#endif /* LIBELAN_ATOMICS */

/* enable use of newer interfaces in SHMEM */
#ifndef CRAY
#ifndef LIBELAN_ATOMICS
/* manpages for shmem_fadd exist on the T3E but library code does not */
#define SHMEM_FADD 
#endif
#endif


/* global scope to prevent compiler optimization of volatile code */
int  _a_temp;
long _a_ltemp;

void armci_generic_rmw(int op, void *ploc, void *prem, int extra, int proc)
{
#if defined(CLUSTER) && !defined(SGIALTIX)
    int lock = (proc-armci_clus_info[armci_clus_id(proc)].master)%NUM_LOCKS;
#else
    int lock = 0;
#endif

    ARMCI_PR_DBG("enter",0);
    NATIVE_LOCK(lock,proc);
    switch (op) {
      case ARMCI_FETCH_AND_ADD:
                armci_get(prem,ploc,sizeof(int),proc);
                _a_temp = *(int*)ploc + extra;
                armci_put(&_a_temp,prem,sizeof(int),proc);
           break;
      case ARMCI_FETCH_AND_ADD_LONG:
                armci_get(prem,ploc,sizeof(long),proc);
                _a_ltemp = *(long*)ploc + extra;
                armci_put(&_a_ltemp,prem,sizeof(long),proc);
           break;
      case ARMCI_SWAP:
#if (defined(__i386__) || defined(__x86_64__)) && !defined(_CRAYC)
        if(SERVER_CONTEXT || armci_nclus==1){
	  atomic_exchange(ploc, prem, sizeof(int));
        }
        else 
#endif
        {
	  armci_get(prem,&_a_temp,sizeof(int),proc);
	  armci_put(ploc,prem,sizeof(int),proc);
	  *(int*)ploc = _a_temp; 
        }
	break;
      case ARMCI_SWAP_LONG:
                armci_get(prem,&_a_ltemp,sizeof(long),proc);
                armci_put(ploc,prem,sizeof(long),proc);
                *(long*)ploc = _a_ltemp;
           break;
      default: armci_die("rmw: operation not supported",op);
    }
    /*TODO memfence here*/
    NATIVE_UNLOCK(lock,proc);
    ARMCI_PR_DBG("exit",0);
}


int PARMCI_Rmw(int op, void *ploc, void *prem, int extra, int proc)
{
    if(!SAMECLUSNODE(proc)){
    # if defined CRAY_REGISTER_ARMCI_MALLOC && HAVE_ONESIDED_FADD
      if(op == ARMCI_FETCH_AND_ADD_LONG) {
         armci_onesided_fadd(ploc, prem, extra, proc);
      } else {
    # endif
         armci_rem_rmw(op, ploc, prem,  extra, proc);
    # if defined CRAY_REGISTER_ARMCI_MALLOC && HAVE_ONESIDED_FADD
      }
    # endif
      return 0;
    }

    switch (op) {
      case ARMCI_FETCH_AND_ADD_LONG:
      # if defined CRAY_REGISTER_ARMCI_MALLOC && HAVE_ONESIDED_FADD
        armci_onesided_fadd(ploc, prem, extra, proc);
        break;
      # endif
      case ARMCI_FETCH_AND_ADD:
      case ARMCI_SWAP:
      case ARMCI_SWAP_LONG:
           armci_generic_rmw(op, ploc, prem,  extra, proc);
        break;
      default: armci_die("rmw: operation not supported",op);
    }

    return 0;
}

