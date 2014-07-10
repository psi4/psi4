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

#include <unistd.h>
#include <sys/fcntl.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/semaphore.h>

#define MAX_SEMA 32
static int fd = -1;
static char template[] = "/tmp/SEMA.XXXXXX";
static char *filename = (char *) NULL;
static unsigned shmem_size;

#if defined(NAMED_SEMAPHORES_SUPPORTED)

static sem_t *sem_arr;

long SemSetCreate(long n_sem, long value)
{
  int i;
  shmem_size = MAX_SEMA*sizeof(sem_t);

  if ( (n_sem <= 0) || (n_sem >= MAX_SEMA) )
    Error("SemSetCreate: n_sem has invalid value",n_sem);

  /* allocate shared memory for locks and semaphore val */
  filename = mktemp(template);
  if ( (fd = shm_open(filename, O_CREAT|O_RDWR, 0666)) < 0 )
    Error("SemSetCreate: failed to open temporary shm file",0);
  sem_arr = (sem_t*) mmap((caddr_t)0, shmem_size, PROT_READ|PROT_WRITE,
                     MAP_ANON|MAP_HASSEMAPHORE|MAP_SHARED, fd, (off_t)0);
  if(!sem_arr)Error("SemSetCreate: failed to mmap",0);

  /* initialize locks and semaphore values */
  for (i=0; i<n_sem; i++) {
      if(sem_init(sem_arr+i,1,1)<0)
             Error("SemSetCreate: sem_init failed",(long)i);
  }
  return 1L;
}

long SemSetDestroyAll()
{
  long status=0;
  if((int)unlink(filename)==-1)Error("SemSetDestroyAll: unlink failed",0);
  status = munmap((char *) shmem_size, 0);
  if(status)status = -1;
  return status;
}


void SemWait(long sem_set_id, long sem_num)
{
  if ( (sem_num < 0) || (sem_num >= MAX_SEMA) )
    Error("SemWait: invalid sem_num",sem_num);
  if(sem_wait(sem_arr+sem_num)<0)
    Error("SemWait: failed",sem_num);
}

void SemPost(long sem_set_id, long sem_num)
{ 
  if ( (sem_num < 0) || (sem_num >= MAX_SEMA) )
    Error("SemPost: invalid sem_num",sem_num);
  if(sem_post(sem_arr+sem_num)<0)
    Error("SemPost: failed",sem_num);
 
}
  
long SemValue(long sem_set_id, long sem_num)  
{ 
  Error("SemValue: not implemented", sem_num);
  return 1L;
}

#else


typedef struct{
        int state;
        int pad[15];
} lock_t;
static lock_t *locks;


static char template1[] = "/tmp/SEMA1.XXXXXX";
static char *filename1 = (char *) NULL;
static sem_t *sem;
static lock_t *locks;

#include <stdio.h>

long SemSetCreate(long n_sem, long value)
{
  int i;
  shmem_size = MAX_SEMA*sizeof(lock_t);

  if ( (n_sem <= 0) || (n_sem >= MAX_SEMA) )
    Error("SemSetCreate: n_sem has invalid value",n_sem);

  /* allocate shared memory for locks and semaphore val */
  locks = (lock_t*) mmap((caddr_t)0, shmem_size, PROT_READ|PROT_WRITE,
                     MAP_ANON|MAP_SHARED, -1, (off_t)0);
  if(locks == (lock_t*)-1)Error("SemSetCreate: failed to mmap",shmem_size);

  filename1 = mktemp(template1);
  sem = sem_open(filename1, O_CREAT|O_EXCL, 0666, 1); 
  if(!sem)Error("SemSetCreate: failed to sem_open",0);
  
  /* initialize locks and semaphore values */
  bzero(locks,shmem_size);
  return 1L;
}

long SemSetDestroyAll()
{
  long status=0;
  status = munmap((char *) locks, shmem_size);
  if(status)status = -1;
  sem_unlink(filename1);
  return status;
}

double __tcgmsg_fred__=0.0;
void Dummy()
{
  int n = 200;                  /* This seems optimal */  
  while(n--) 
    __tcgmsg_fred__++;
}   
  
void SemWait(long sem_set_id, long sem_num)
{ 
  int value = 0, count=0;
  
  if ( (sem_num < 0) || (sem_num >= MAX_SEMA) )
    Error("SemWait: invalid sem_num",sem_num);
   
  while (value<=0) {
    if(sem_wait(sem)<0)Error("SemWait: sem_op error",sem_num);;
    value = locks[sem_num].state;
    if (value>0)
      locks[sem_num].state--;
    if(sem_post(sem)<0)Error("SemWait: sem_op error",sem_num);;
    if (value<=0) Dummy();
    count++;
    if(count%1000 == 999)usleep(1);
  }
} 
  
void SemPost(long sem_set_id, long sem_num)
{ 
  if ( (sem_num < 0) || (sem_num >= MAX_SEMA) )
    Error("SemPost: invalid sem_num",sem_num); 
    
  if(sem_wait(sem)<0)Error("SemPost: sem_op error",sem_num);;
      locks[sem_num].state++;
  if(sem_post(sem)<0)Error("SemWait: sem_op error",sem_num);;
} 

long SemValue(long sem_set_id, long sem_num)
{ 
  if ( (sem_num < 0) || (sem_num >= MAX_SEMA) )
    Error("SemVal: invalid sem_num",sem_num); 
  return (long)locks[sem_num].state;
}   
#endif  
