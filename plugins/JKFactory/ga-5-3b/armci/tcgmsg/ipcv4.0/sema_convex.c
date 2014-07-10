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

#include <stdio.h>
#include <sys/param.h>
#include <sys/file.h>
#include <sys/mman.h>
#include <sys/types.h>

#define MAX_SEM_SETS 20
#define MAX_N_SEM 100

/* On the convex a semaphore is a structure but on the apollo
   it is an array which does not need dereferencing.  Use ADDR
   to generate the address of a semaphore */
#ifdef APOLLO
#define ADDR(x) x
#else
#define ADDR(x) &x
#endif

extern char *mktemp();

struct sem_set_struct {
  int n_sem;                    /* no. of semaphores in set */
  semaphore lock[MAX_N_SEM];    /* locks for changing value */
  semaphore wait[MAX_N_SEM];    /* locks for queing */
  int value[MAX_N_SEM];         /* values */
};

static int num_sem_set = 0;
static struct sem_set_struct *sem_sets;
static int fd = -1;
static char template[] = "/tmp/SEMA.XXXXXX";
static char *filename = (char *) NULL;

void InitSemSets()
/* Initialise sem_sets and allocate associated shmem region */
{
  int i, j;
  unsigned size = sizeof(struct sem_set_struct) * MAX_SEM_SETS;

#ifndef APOLLO
  /* Generate scratch file to identify region ... mustn't do this
     on the APOLLO */

  filename = mktemp(template);
  if ( (fd = open(filename, O_RDWR|O_CREAT, 0666)) < 0 )
    Error("InitSemSets: failed to open temporary file",0);
#endif

  sem_sets = (struct sem_set_struct *) mmap((caddr_t) 0, &size,
                     PROT_READ|PROT_WRITE,
                     MAP_ANON|MAP_HASSEMAPHORE|MAP_SHARED, fd, 0);

#ifdef APOLLO
  if (sem_sets == (struct sem_set_struct *) 0)
    Error("InitSemSets: mmap failed", (long) -1);
#else
  if (sem_sets == (struct sem_set_struct *) -1)
    Error("InitSemSets: mmap failed", (long) -1);
#endif

  for (i=0; i<MAX_SEM_SETS; i++) {
    sem_sets[i].n_sem = 0;
    for (j=0; j<MAX_N_SEM; j++) {
      mclear(ADDR(sem_sets[i].lock[j]));
      mclear(ADDR(sem_sets[i].wait[j]));
      sem_sets[i].value[j] = 0;
    }
  }
}

long SemSetCreate(n_sem, value)
     long n_sem;
     long value;
{
  int i;

  /* Check for errors and initialise data if first entry */

  if ( (n_sem <= 0) || (n_sem >= MAX_N_SEM) )
    Error("SemSetCreate: n_sem has invalid value",n_sem);

  if (num_sem_set == 0)
    InitSemSets();
  else if (num_sem_set >= MAX_SEM_SETS)
    Error("SemSetCreate: Exceeded man no. of semaphore sets",
          num_sem_set);

  /* Initialize the values */

  for (i=0; i<n_sem; i++)
    sem_sets[num_sem_set].value[i] = value;

  sem_sets[num_sem_set].n_sem = n_sem;

  num_sem_set++;

  return (long) (num_sem_set - 1);
}

void SemWait(sem_set_id, sem_num)
     long sem_set_id;
     long sem_num;
{
  if ( (sem_set_id < 0) || (sem_set_id >= num_sem_set) )
    Error("SemWait: invalid sem_set_id",sem_set_id);
  if ( (sem_num < 0) || (sem_num >= sem_sets[sem_set_id].n_sem) )
    Error("SemWait: invalid semaphore number in set",sem_num);

  while (1) {

    /* Get the lock around the whole semaphore */

    (void) mset(ADDR(sem_sets[sem_set_id].lock[sem_num]), 1);

    /* If the value is positive fall thru, else wait */

    if (sem_sets[sem_set_id].value[sem_num] > 0)
      break;
    else {
      (void) mclear(ADDR(sem_sets[sem_set_id].lock[sem_num]));
      (void) mset(ADDR(sem_sets[sem_set_id].wait[sem_num]), 1);
    }
  }

  /* Are ready to go ... decrement the value and release lock */

  sem_sets[sem_set_id].value[sem_num]--;
  (void) mclear(ADDR(sem_sets[sem_set_id].lock[sem_num]));

}

void SemPost(sem_set_id, sem_num)
     long sem_set_id;
     long sem_num;
{
  int i;

  if ( (sem_set_id < 0) || (sem_set_id >= num_sem_set) )
    Error("SemPost: invalid sem_set_id",sem_set_id);
  if ( (sem_num < 0) || (sem_num >= sem_sets[sem_set_id].n_sem) )
    Error("SemPost: invalid semaphore number in set",sem_num);

  /* Get the lock around the whole semaphore */

  (void) mset(ADDR(sem_sets[sem_set_id].lock[sem_num]), 1);
  
  /* Read and increment the value. If is now zero wake up
     up the queue */

  sem_sets[sem_set_id].value[sem_num]++;
  i = sem_sets[sem_set_id].value[sem_num];

  (void) mclear(ADDR(sem_sets[sem_set_id].lock[sem_num]));
  if (i >= 0)
    (void) mclear(ADDR(sem_sets[sem_set_id].wait[sem_num]));
}

long SemValue(sem_set_id, sem_num)
     long sem_set_id;
     long sem_num;
{
  int i;

  if ( (sem_set_id < 0) || (sem_set_id >= num_sem_set) )
    Error("SemValue: invalid sem_set_id",sem_set_id);
  if ( (sem_num < 0) || (sem_num >= sem_sets[sem_set_id].n_sem) )
    Error("SemValue: invalid semaphore number in set",sem_num);

  /* There seems no point in getting the lock just to read
     the value and it seems more useful not to (e.g. debugging) */

  i = sem_sets[sem_set_id].value[sem_num];

  return (long) (i-1);
}

long SemSetDestroy(sem_set_id)
     long sem_set_id;
{

  if ( (sem_set_id < 0) || (sem_set_id >= num_sem_set) )
    return -1;

  sem_sets[sem_set_id].n_sem = 0;

  return (long) 0;
}
  
long SemSetDestroyAll()
{
  long i, status=0;

  for (i=0; i<num_sem_set; i++)
    if (sem_sets[i].n_sem)
      status += SemSetDestroy(i);

  if (fd >= 0) {
    (void) close(fd);
    fd = -1;
    (void) unlink(filename);
  }

  status += munmap((char *) sem_sets, 0);

  if (status)
    status = -1;

  return status;
}
