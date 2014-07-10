/* $Id: utils.h,v 1.1.2.3 2007-07-02 05:35:31 d3p687 Exp $
 *
 * primitives for transparent handling of multi-threading
 */

#ifndef UTILS_H
#define UTILS_H

/*
 * This header file describes the "barrier" synchronization
 * construct. The type barrier_t describes the full state of the
 * barrier including the POSIX 1003.1c synchronization objects
 * necessary.
 *
 * A barrier causes threads to wait until a set of threads has
 * all "reached" the barrier. The number of threads required is
 * set when the barrier is initialized, and cannot be changed
 * except by reinitializing.
 */


#ifdef THREAD_SAFE
#   ifdef POSIX_THREADS

#       include <pthread.h>

#if 1
        typedef pthread_mutex_t thread_lock_t;
#       define THREAD_LOCK_INIT(x)    pthread_mutex_init(&x,NULL)
#       define THREAD_LOCK_DESTROY(x) pthread_mutex_destroy(&x)
#       define THREAD_LOCK(x)         pthread_mutex_lock(&x)
#       define THREAD_UNLOCK(x)       pthread_mutex_unlock(&x)
#else

#ifndef INLINE
#       define INLINE
#       include "spinlock.h"
#       undef INLINE
#else
#       include "spinlock.h"
#endif

        typedef LOCK_T thread_lock_t;
#       define THREAD_LOCK_INIT(x)    armci_init_spinlock(&x)
#       define THREAD_LOCK_DESTROY(x) 0
#       define THREAD_LOCK(x)         armci_acquire_spinlock(&x)
#       define THREAD_UNLOCK(x)       armci_release_spinlock(&x)
#endif
        typedef pthread_t thread_t;
#       define THREAD_CREATE(th_,func_,arg_) pthread_create(th_,NULL,func_,arg_)
#       define THREAD_JOIN(th_,ret_) pthread_join(th_,ret_)

        /* structure describing a barrier */
        typedef struct thread_barrier_tag {
            pthread_mutex_t     mutex;          /* Control access to barrier */
            pthread_cond_t      cv;             /* wait for barrier */
            int                 valid;          /* set when valid */
            int                 threshold;      /* number of threads required */
            int                 counter;        /* current number of threads */
            int                 cycle;          /* alternate wait cycles (0 or 1) */
        } thread_barrier_t;

#       define BARRIER_VALID   0xdbcafe

        /* support static initialization of barriers */
#       define BARRIER_INITIALIZER(cnt) {\
            PTHREAD_MUTEX_INITIALIZER, PTHREAD_COND_INITIALIZER,\
            BARRIER_VALID, cnt, cnt, 0}

#   else
#       error ONLY PTHREADS SUPPORT HAS BEEN IMPLEMENTED
#   endif

#   define TH2PROC(th_) (th_/mt_tpp) /* computes processor from thread id */

    /* barrier functions */
    int thread_barrier_init (thread_barrier_t *barrier, int count);
    int thread_barrier_destroy (thread_barrier_t *barrier);
    int thread_barrier_wait (thread_barrier_t *barrier);

    /* multi-threaded memory functions */
    int armci_malloc_mt(void *ptr[], int bytes);
    int armci_free_mt(void *ptr, int th_idx);
#   define ARMCI_MALLOC_MT armci_malloc_mt
#   define ARMCI_FREE_MT armci_free_mt


#   define TH_INIT(p_,t_)   mt_size=p_;mt_tpp=t_;\
                            thread_barrier_init(&mt_barrier,mt_tpp)
#   define TH_FINALIZE()    thread_barrier_destroy(&mt_barrier)
#   define MT_BARRIER()     if (thread_barrier_wait(&mt_barrier)==-1) armci_msg_barrier();\
                                thread_barrier_wait(&mt_barrier)

    extern int mt_size;
    extern int mt_tpp;
    extern thread_barrier_t mt_barrier;
#else
#   define THREAD_LOCK_INIT(x)
#   define THREAD_LOCK_DESTROY(x)
#   define THREAD_LOCK(x)
#   define THREAD_UNLOCK(x)
#   define TH_INIT(p_,t_)
#   define TH_FINALIZE()
#   define MT_BARRIER armci_msg_barrier
#   define ARMCI_MALLOC_MT PARMCI_Malloc
#   define ARMCI_FREE_MT(p_,th_) PARMCI_Free(p_)
#endif






#endif/*UTILS_H*/


