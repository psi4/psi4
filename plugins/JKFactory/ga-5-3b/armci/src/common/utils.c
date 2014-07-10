#if HAVE_CONFIG_H
#   include "config.h"
#endif

/*
 * A barrier causes threads to wait until a set of threads has
 * all "reached" the barrier. The number of threads required is
 * set when the barrier is initialized, and cannot be changed
 * except by reinitializing.
 *
 * The barrier_init() and barrier_destroy() functions,
 * respectively, allow you to initialize and destroy the
 * barrier.
 *
 * The barrier_wait() function allows a thread to wait for a
 * barrier to be completed. One thread (the one that happens to
 * arrive last) will return from barrier_wait() with the status
 * -1 on success -- others will return with 0. The special
 * status makes it easy for the calling code to cause one thread
 * to do something in a serial region before entering another
 * parallel section of code.
 */

#if HAVE_SYS_TYPES_H
#   include <sys/types.h>
#endif
#if HAVE_SYS_ERRNO_H
#   include <sys/errno.h>
#endif
#if HAVE_SYS_TIME_H
#   include <sys/time.h>
#endif
#if HAVE_UNISTD_H
#   include <unistd.h>
#endif
#if HAVE_ERRNO_H
#   include <errno.h>
#endif
#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_STRING_H
#   include <string.h>
#endif
#include "utils.h"

#define DEBUG_

int mt_size;   /* number of processes: needed for collective mt ops */
int mt_tpp; /* number of threads used for collective ops */
thread_barrier_t mt_barrier; /* static barrier used for multi-threaded MT_BARRIER */

int armci_malloc_mt(void *ptr[], int bytes)
{
    int rc, th_size, i, j;

    th_size = mt_size * mt_tpp;
    if (thread_barrier_wait(&mt_barrier)==-1) {
        rc = PARMCI_Malloc(ptr, bytes * mt_tpp);
#ifdef DEBUG
        printf("bytes=%d\n", bytes);
        for (i = 0; i < mt_size; i++) printf("ptr[%d]=%p\n",i,ptr[i]);
#endif
        /* at this point proc ptrs are at beggining of the list */
        for (i = mt_size - 1; i >= 0; i--) for (j = mt_tpp - 1; j >= 0; j--) {
#ifdef DEBUG
            printf("mt_size=%d,mt_tpp=%d,i=%d,j=%d,ptr[%d]=%p+%d\n",
                    mt_size,mt_tpp,i,j,i*mt_tpp+j,ptr[i],j*bytes);
            fflush(stdout);
#endif
            ptr[i * mt_tpp + j] = ((char*)ptr[i]) + j * bytes;
        }
    }
    thread_barrier_wait(&mt_barrier);

    return rc;
}

int armci_free_mt(void *ptr, int th_idx)
{
}

#ifdef POSIX_THREADS
/*
 * Initialize a barrier for use.
 */
int thread_barrier_init (thread_barrier_t *barrier, int count)
{
    int status;

    barrier->threshold = barrier->counter = count;
    barrier->cycle = 0;
    status = pthread_mutex_init (&barrier->mutex, NULL);
    if (status != 0)
        return status;
    status = pthread_cond_init (&barrier->cv, NULL);
    if (status != 0) {
        pthread_mutex_destroy (&barrier->mutex);
        return status;
    }
    barrier->valid = BARRIER_VALID;
    return 0;
}

/*
 * Destroy a barrier when done using it.
 */
int thread_barrier_destroy (thread_barrier_t *barrier)
{
    int status, status2;

    if (barrier->valid != BARRIER_VALID)
        return EINVAL;

    status = pthread_mutex_lock (&barrier->mutex);
    if (status != 0)
        return status;

    /*
     * Check whether any threads are known to be waiting; report
     * "BUSY" if so.
     */
    if (barrier->counter != barrier->threshold) {
        pthread_mutex_unlock (&barrier->mutex);
        return EBUSY;
    }

    barrier->valid = 0;
    status = pthread_mutex_unlock (&barrier->mutex);
    if (status != 0)
        return status;

    /*
     * If unable to destroy either 1003.1c synchronization
     * object, return the error status.
     */
    status = pthread_mutex_destroy (&barrier->mutex);
    status2 = pthread_cond_destroy (&barrier->cv);
    return (status == 0 ? status : status2);
}

/*
 * Wait for all members of a barrier to reach the barrier. When
 * the count (of remaining members) reaches 0, broadcast to wake
 * all threads waiting.
 */
int thread_barrier_wait (thread_barrier_t *barrier)
{
    int status, cancel, tmp, cycle;

    if (barrier->valid != BARRIER_VALID)
        return EINVAL;

    status = pthread_mutex_lock (&barrier->mutex);
    if (status != 0)
        return status;

    cycle = barrier->cycle;   /* Remember which cycle we're on */

    if (--barrier->counter == 0) {
        barrier->cycle = !barrier->cycle;
        barrier->counter = barrier->threshold;
        status = pthread_cond_broadcast (&barrier->cv);
        /*
         * The last thread into the barrier will return status
         * -1 rather than 0, so that it can be used to perform
         * some special serial code following the barrier.
         */
        if (status == 0)
            status = -1;
    } else {
        /*
         * Wait with cancellation disabled, because barrier_wait
         * should not be a cancellation point.
         */
        pthread_setcancelstate (PTHREAD_CANCEL_DISABLE, &cancel);

        /*
         * Wait until the barrier's cycle changes, which means
         * that it has been broadcast, and we don't want to wait
         * anymore.
         */
        while (cycle == barrier->cycle) {
            status = pthread_cond_wait (
                    &barrier->cv, &barrier->mutex);
            if (status != 0) break;
        }

        pthread_setcancelstate (cancel, &tmp);
    }
    /*
     * Ignore an error in unlocking. It shouldn't happen, and
     * reporting it here would be misleading -- the barrier wait
     * completed, after all, whereas returning, for example,
     * EINVAL would imply the wait had failed. The next attempt
     * to use the barrier *will* return an error, or hang, due
     * to whatever happened to the mutex.
     */
    pthread_mutex_unlock (&barrier->mutex);
    return status;          /* error, -1 for waker, or 0 */
}
#endif

#if 0

/***
   NAME
     timing.c
   PURPOSE
     Timing routines for calculating the execution time:
       void start_timer(void);  Set the timer.
       double elapsed_time(void);  Return the timing elapsed since
                                   the timer has been set.
   NOTES
     Jialin Ju - Oct 16, 1995 Created.
***/

/* Timing routines that use standard Unix gettingofday() */
static struct timezone tz;
static struct timeval start_time, finish_time;

/* Start measuring a time delay */
void start_timer(void)
{
    gettimeofday( &start_time, &tz);
}

/* Retunrn elapsed time in milliseconds */
double elapsed_time(void)
{
    gettimeofday( &finish_time, &tz);
    return(1000.0*(finish_time.tv_sec - start_time.tv_sec) +
           (finish_time.tv_usec - start_time.tv_usec)/1000.0 );
}

/* Return the stopping time in milliseconds */
double stop_time(void)
{
    gettimeofday( &finish_time, &tz);
    return(1000.0*finish_time.tv_sec + finish_time.tv_usec/1000.0);
}

#endif
