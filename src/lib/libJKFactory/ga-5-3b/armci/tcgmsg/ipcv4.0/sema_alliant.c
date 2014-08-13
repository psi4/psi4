#if HAVE_CONFIG_H
#   include "config.h"
#endif

/*
  These routines simplify the interface to semaphores for use in mutual
  exclusion and queuing. Hopefully I can also make this portable.

  An external routine Error is assumed which is called upon an error
  and tidies up by calling SemSetDestroyAll.

  In most cases errors cause an internal hard failure (by calling Error).

  1) make an array of n_sem semaphores, returning the id associated
     with the entire set. All the semaphore values are initialized to value
     which should be a positve integer (queuing) or 0 (synchronization).
     The semaphores in the set are indexed from 0 to n_sem-1.

     long SemSetCreate(long n_sem, long value)

  2) Decrement and test the value associated with the semaphore specified by 
     (sem_set_id, sem_num). In effect this:

     if (value >= 0) {
        continue execution
     }
     else {
        wait in queue for the semaphore
     }
     decrement value

     void SemWait(long sem_set_id, long sem_num)

  3) Increment the value associated with the semaphore specified by
     (sem_set_id, sem_num). If value <= 0 (i.e. there are processes
     in the queue) this releases the next process.

     void SemPost(long sem_set_id, long sem_num)
     
  4) Return the current value associated with the semaphore sepcified by
     (sem_set_id, sem_num).

     long SemValue(long sem_set_id, long sem_num)

  5) Destroy the set of semaphores. Any other processes that are accessing
     or try to access the semaphore set should get an error.
     On the SUN (all system V machines?) the semaphore sets should
     be destroyed explicitly before the final process exits.
     0 is returned if OK. -1 implies an error.

     long SemSetDestroy(long sem_set_id)

  6) Destroy all the semaphore sets that are known about. This is really
     meant for an error routine to call to try and tidy up. Though all
     applications could call it before the last process exits.
     0 is returned if OK. -1 implies an error.

     long SemSetDestroyAll()
*/

extern void Error();

/*************************************************************
    Alliant Concentrix 5.0 and Concentrix FX/2800
 *************************************************************/

/* This is very specific to the Alliant. */

#include <sys/rtsem.h>
#include <sys/errno.h>

extern int errno;

/* On the alliant semaphores are handed out one at a time rather than
   in sets, so have to maintain sets manually */

#define MAX_SEM_SETS 20
#define MAX_N_SEM 128

static struct sem_set_list_struct {
  int id[MAX_N_SEM];                       /* alliant semaphore id */
  int n_sem;                               /* no. of semaphores in set */
} sem_set_list[MAX_SEM_SETS];

static int num_sem_set = 0;


void InitSemSetList()
/* Initialise sem_set_list */
{
  int i, j;
  
  for (i=0; i<MAX_SEM_SETS; i++) {
    sem_set_list[i].n_sem = 0;
    for (j=0; j<MAX_N_SEM; j++)
      sem_set_list[i].id[j] = -1;
  }
}

long SemSetCreate(n_sem, value)
     long n_sem;
     long value;
{
  int semid, i, j;

  /* Check for errors and initialise data if first entry */

  if ( (n_sem <= 0) || (n_sem >= MAX_N_SEM) )
    Error("SemSetCreate: n_sem has invalid value", (long) n_sem);

  if (num_sem_set == 0)
    InitSemSetList();
  else if (num_sem_set >= MAX_SEM_SETS)
    Error("SemSetCreate: Exceeded man no. of semaphore sets",
          (long) num_sem_set);

  /* Find first empty slot in sem_set_list */

  for (i=0; i < MAX_SEM_SETS; i++) 
    if (sem_set_list[i].n_sem == 0)
      break;

  if (i == MAX_SEM_SETS)
    Error("SemSetCreate: internal error puting semid in list", (long) i);

  /* Actually make the semaphore set */

  for (j=0; j<n_sem; j++) {
    if ( (semid = sem_create(value, value, SEMQUE_FIFO, 0)) < 0)
      Error("SemSetCreate: failed to create semaphore", (long) j);
    sem_set_list[i].id[j] = semid;
  }

  num_sem_set++;
  
  return i;
}

void SemWait(sem_set_id, sem_num)
     long sem_set_id;
     long sem_num;
{
interrupted:
  if (sem_wait(sem_set_list[sem_set_id].id[sem_num]) < 0) {
    if (errno == EINTR)
      goto interrupted;   /* got zapped by a signal ... try again */
    else
      Error("SemWait: error from sem_wait", (long) -1);
  }
}

void SemPost(sem_set_id, sem_num)
     long sem_set_id;
     long sem_num;
{
  if (sem_post(sem_set_list[sem_set_id].id[sem_num]) < 0)
    Error("SemPost: error from sem_post", (long) -1);
}

long SemValue(sem_set_id, sem_num)
     long sem_set_id;
     long sem_num;
{
  SEM_INFO info;

  if (sem_info(sem_set_list[sem_set_id].id[sem_num], &info, sizeof info) < 0)
    Error("SemValue: error from sem_info", (long) -1);

  return info.curval;
}

long SemSetDestroy(sem_set_id)
     long sem_set_id;
{
  int status=0, i;

  /* Close each semaphore in the set */

  for (i=0; i<sem_set_list[sem_set_id].n_sem; i++) {
    status += sem_destroy(sem_set_list[sem_set_id].id[i]);
    sem_set_list[sem_set_id].id[i] = -1;
  }

  sem_set_list[sem_set_id].n_sem = 0;

  num_sem_set--;

  if (status)
    status = -1;

  return (long) status;
}
  
long SemSetDestroyAll()
{
  int i, status=0;

  for (i=0; i<MAX_SEM_SETS; i++)
    if (sem_set_list[i].n_sem)
      status += SemSetDestroy(i);

  if (status)
    status = -1;

  return (long) status;
}
