/* $Header: /tmp/hpctools/ga/tcgmsg/ipcv4.0/sema.h,v 1.4 1995-02-24 02:17:41 d3h325 Exp $ */

/* Header file declaring stubs for semaphore routines. */

/*
  These routines simplify the interface to semaphores for use in mutual
  exclusion and queuing. Hopefully I can also make this portable.

  Interruption by signals is not tested for.

  An external routine Error is assumed which is called upon an error
  and tidies up by calling SemSetDestroyAll.

  In most cases errors cause an internal hard failure (by calling Error).
*/

/*
  1) make an array of n_sem semaphores, returning the id associated
     with the entire set. All the semaphore values are initialized to value
     which should be a positve integer (queuing) or 0 (synchronization).
     The semaphores in the set are indexed from 0 to n_sem-1.

     long SemSetCreate(long n_sem, long value)
*/
extern long SemSetCreate();


/*
  2) Decrement and test the value associated with the semaphore specified by 
     (sem_set_id, sem_num). In effect this:

     decrement value

     if (value >= 0) {
        continue execution
	}
     else {
        wait in queue for the semaphore
	}

     void SemWait(long sem_set_id, long sem_num)
*/
extern void SemWait();


/*
  3) Increment the value associated with the semaphore specified by
     (sem_set_id, sem_num). If value <= 0 (i.e. there are processes
     in the queue) this releases the next process.

     void SemPost(long sem_set_id, long sem_num)
*/
extern void SemPost();

     
/*
  4) Return the current value associated with the semaphore sepcified by
     (sem_set_id, sem_num).

     long SemValue(long sem_set_id, long sem_num)
*/
extern long SemValue();


/*
  5) Destroy the set of semaphores. Any other processes that are accessing
     or try to access the semaphore set should get an error.
     On the SUN (all system V machines?) the semaphore sets should
     be destroyed explicitly before the final process exits.
     0 is returned if OK. -1 implies an error.

     long SemSetDestroy(long sem_set_id)
*/
extern long SemSetDestroy();


/*
  6) Destroy all the semaphore sets that are known about. This is really
     meant for an error routine to call to try and tidy up. Though all
     applications could call it before the last process exits.
     0 is returned if OK. -1 implies an error.

     long SemSetDestroyAll()
*/
extern long SemSetDestroyAll();
