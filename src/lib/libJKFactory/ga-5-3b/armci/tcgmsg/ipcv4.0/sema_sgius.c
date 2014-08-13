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

/*
  SGI fast US library semaphores ... aren't any faster
  than system V semaphores ... implement using spin locks
*/

#include <stdio.h>
#include <unistd.h>  
#define MAX_SEMA 512
static volatile int *val;
#define NAME_LEN 200

#ifdef SGI
# include <ulocks.h>
  static usptr_t *arena_ptr;
  static ulock_t *locks[MAX_SEMA];
  static char arena_name[NAME_LEN];
# define EIGHT 8
# define LOCK ussetlock
# define UNLOCK usunsetlock
#define  JUMP EIGHT

#include "sndrcvP.h"

extern char *getenv(const char *);

long SemSetCreate(long n_sem, long value)
{
  int i;
  char *tmp;
  if (!(tmp = getenv("ARENA_DIR"))) tmp = "/tmp";

   sprintf(arena_name,"%s/tcgmsg.arena.%ld",tmp, (long)getpid()); 
#ifdef PRIVATE_ARENA
   (void) usconfig(CONF_ARENATYPE, US_SHAREDONLY);
#endif
   (void) usconfig(CONF_INITUSERS, (unsigned int)SR_clus_info[SR_clus_id].nslave );
#ifdef SGI
    (void) usconfig(CONF_INITSIZE, 1024*1024);
#endif

  if (!(arena_ptr = usinit(arena_name)))
    Error("SemSetCreate: failed to create arena", 0L);

  /* Magic factors of EIGHT here to ensure that values are
     in different cache lines to avoid aliasing -- good on SGI and Convex */

  if (!(val = (int *) usmalloc(EIGHT*MAX_SEMA*sizeof(int), arena_ptr)))
    Error("SemSetCreate: failed to get shmem", (long) (MAX_SEMA*sizeof(int)));

  for (i=0; i<n_sem; i++) {
    if (!(locks[i] = usnewlock(arena_ptr)))
      Error("SemSetCreate: failed to create lock", (long) i);
    val[i*EIGHT] = (int) value;
  }
  return 1L;
}

long SemSetDestroyAll()
{
/*  usdetach (arena_ptr);*/
  arena_ptr = 0;
  unlink(arena_name);
  return 0;
}

#endif


#ifdef SPPLOCKS
#include <sys/param.h>
#include <sys/file.h>
#include <sys/cnx_mman.h>
#include <sys/mman.h>
#include <sys/types.h>

#define SIXTEEN 16
#define JUMP SIXTEEN
typedef struct{
        int state;
        int pad[15]; 
} lock_t;

static lock_t *locks;

#  define    LOCK(x) set_lock(&x.state)
#  define  UNLOCK(x) unset_lock(&x.state)
#  define INILOCK(x) init_lock(&x.state)


void init_lock(int * volatile ip)
{
    *ip = 1;
}

void set_lock(int * volatile ip)
{
    while (1) {
        while (!(*ip));
        if (__ldcws32(ip))
                        break;
    }
}

void unset_lock(int *ip)
{
  *ip = 1;
  asm("sync");
}

static int fd = -1;
static char template[] = "/tmp/SEMA.XXXXXX";
static char *filename = (char *) NULL;
static unsigned shmem_size;

long SemSetCreate(long n_sem, long value)
{
  int i;
  shmem_size = SIXTEEN*MAX_SEMA*sizeof(int)+MAX_SEMA*sizeof(lock_t);

  if ( (n_sem <= 0) || (n_sem >= MAX_SEMA) )
    Error("SemSetCreate: n_sem has invalid value",n_sem);

  /* allocate shared memory for locks and semaphore val */
  filename = mktemp(template);
  if ( (fd = open(filename, O_RDWR|O_CREAT, 0666)) < 0 )
    Error("SemSetCreate: failed to open temporary file",0);
  val = (int *) mmap((caddr_t) 0, shmem_size,
                     PROT_READ|PROT_WRITE,
                     MAP_ANONYMOUS|CNX_MAP_SEMAPHORE|MAP_SHARED, fd, 0);
  locks = (lock_t*)( val + SIXTEEN*MAX_SEMA);

  /* initialize locks and semaphore values */
  for (i=0; i<n_sem; i++) {
    INILOCK(locks[i]);
    val[i*SIXTEEN] = (int) value;
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

#endif


double __tcgmsg_fred__=0.0;

void Dummy()
{
  int n = 200;            /* This seems optimal */
  while(n--)
    __tcgmsg_fred__++;
}

void SemWait(long sem_set_id, long sem_num)
{
  int value = 0;
  int off = sem_num*JUMP;

  if ( (sem_num < 0) || (sem_num >= MAX_SEMA) )
    Error("SemWait: invalid sem_num",sem_num);

  while (value<=0) {
    LOCK(locks[sem_num]);
    value = val[off];
    if (value>0)
      val[off]--;
    UNLOCK(locks[sem_num]);
    if (value<=0) 
      Dummy();
  }
}

void SemPost(long sem_set_id, long sem_num)
{
  int off = sem_num*JUMP;
  if ( (sem_num < 0) || (sem_num >= MAX_SEMA) )
    Error("SemPost: invalid sem_num",sem_num);

  LOCK(locks[sem_num]);
  val[off]++;
  UNLOCK(locks[sem_num]);
}

long SemValue(long sem_set_id, long sem_num)
{
  Error("SemValue: not implemented", sem_num);
  return 1;
}

long SemSetDestroy(long sem_set_id)
{
  return(SemSetDestroyAll());
}
