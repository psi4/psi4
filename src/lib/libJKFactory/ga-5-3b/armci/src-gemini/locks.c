#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Id: locks.c,v 1.15.6.1 2006-12-14 13:24:36 manoj Exp $ */
#define _LOCKS_C_
#include "armcip.h"
#include "locks.h"
#ifndef WIN32
#   include <unistd.h>
#endif
#include <stdio.h>


extern void armci_die(char*,int);

#if defined(SPINLOCK) || defined(PMUTEXES)

void **ptr_arr;

void CreateInitLocks(int num_locks, lockset_t *plockid)
{
int locks_per_proc, size;
extern void armci_set_serv_mutex_arr(void *);
    ARMCI_PR_DBG("enter",0);
  ptr_arr = (void**)malloc(armci_nproc*sizeof(void*));
  locks_per_proc = (num_locks*armci_nclus)/armci_nproc + 1;
  size=locks_per_proc*sizeof(PAD_LOCK_T);
  PARMCI_Malloc(ptr_arr, size);
  _armci_int_mutexes = (PAD_LOCK_T*) ptr_arr[armci_master];
# ifdef PORTALS_SPECIFIC_QUESTION
  if(armci_me==armci_master)armci_set_serv_mutex_arr(_armci_int_mutexes);
# endif

  if(!_armci_int_mutexes) armci_die("Failed to create spinlocks",size);

#ifdef PMUTEXES
  if(armci_me == armci_master) {
       int i;
       pthread_mutexattr_t pshared;
       if(pthread_mutexattr_init(&pshared))
            armci_die("armci_allocate_locks: could not init mutex attr",0);
#      ifndef LINUX
         if(pthread_mutexattr_setpshared(&pshared,PTHREAD_PROCESS_SHARED))
            armci_die("armci_allocate_locks: could not set PROCESS_SHARED",0);
#      endif

       for(i=0; i< locks_per_proc*armci_clus_info[armci_clus_me].nslave; i++){
             if(pthread_mutex_init(_armci_int_mutexes+i,&pshared))
                armci_die("armci_allocate_locks: could not init mutex",i);
       }
  }
#else

  bzero((char*)ptr_arr[armci_me],size);
    ARMCI_PR_DBG("exit",0);
#endif
} 

void InitLocks(int num_locks, lockset_t lockid)
{
    /* what are you doing here ? 
       All processes should've called CreateInitLocks().
       Check preprocessor directtives and see lock allocation in armci_init */
    armci_die("InitLocks(): what are you doing here ?",armci_me);
}


void DeleteLocks(lockset_t lockid)
{
  _armci_int_mutexes = (PAD_LOCK_T*)0;
}

#else
/*********************** every thing else *************************/

void CreateInitLocks(int num_locks, lockset_t  *lockid)
{}

void InitLocks(int num_locks, lockset_t lockid)
{
}


void DeleteLocks(lockset_t lockid)
{
}

#endif

