#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Id: threads.c,v 1.1.2.5 2007-08-28 21:29:46 manoj Exp $ */

#if 0
#   define PRNDBG3(m,a1,a2,a3) \
           fprintf(stderr,"DBG %d: " m,armci_me,a1,a2,a3);fflush(stderr)
#   define PRNDBG(m) PRNDBG3(m,0,0,0)
#   define PRNDBG1(m,a1) PRNDBG3(m,a1,0,0)
#   define PRNDBG2(m,a1,a2) PRNDBG3(m,a1,a2,0)
#else
#   define PRNDBG(m)
#   define PRNDBG1(m,a1)
#   define PRNDBG2(m,a1,a2)
#   define PRNDBG3(m,a1,a2,a3)
#endif


#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#include "armcip.h"

armci_user_threads_t armci_user_threads;

void armci_init_threads()
{
    int i, bytes;
    char *uval = getenv("ARMCI_MAX_THREADS");
    
    armci_user_threads.max   = 1;
    armci_user_threads.avail = 0;

    if (uval != NULL) sscanf(uval, "%d", &armci_user_threads.max);
    
    if (armci_user_threads.max < 1 ||
        armci_user_threads.max > ARMCI_THREADS_LIMIT)
    {
       printf("Error: Only 1-%d threads are supported. ",ARMCI_THREADS_LIMIT);
       printf("Set ARMCI_MAX_THREADS appropriately\n"); fflush(stdout);
       armci_die("armci_init_threads: failed", 0);
    }

    bytes = sizeof(thread_id_t) * armci_user_threads.max;
    if ( !(armci_user_threads.ids = (thread_id_t*) malloc(bytes)) )
    {
       armci_die("armci_init_threads: armci_user_threads.ids malloc failed",
                 armci_user_threads.max);
    }
    memset(armci_user_threads.ids, 0, bytes);

#if 0 /* spinlock has void return value */
    if (THREAD_LOCK_INIT(armci_user_threads.lock)     ||
        THREAD_LOCK_INIT(armci_user_threads.buf_lock) ||
        THREAD_LOCK_INIT(armci_user_threads.net_lock))
        armci_die("armci_init_threads:locks initialization failed", 0);
#else
    THREAD_LOCK_INIT(armci_user_threads.lock);
    THREAD_LOCK_INIT(armci_user_threads.buf_lock);
    THREAD_LOCK_INIT(armci_user_threads.net_lock);
#endif
    
#if 0
    /* using one lock per socket for now, it might be feasible (and usefull)
     * to use two (one for sending and one for receiving) */
    armci_user_threads.sock_locks = malloc(armci_nclus *sizeof(thread_lock_t));
    for (i = 0; i < armci_nclus; i++)
       if (THREAD_LOCK_INIT(armci_user_threads.sock_locks[i]))
          armci_die("armci_init_threads:sock locks initialization failed", i);
#endif
}

void armci_finalize_threads()
{
    THREAD_LOCK_DESTROY(armci_user_threads.lock);
    THREAD_LOCK_DESTROY(armci_user_threads.net_lock);
    THREAD_LOCK_DESTROY(armci_user_threads.buf_lock);
    free(armci_user_threads.ids);
}

/* calling armci_thread_idx for every function that accesses thread-private data
 * might be expensive -- needs optiomization */
INLINE int armci_thread_idx()
{
    int i, n = ARMCI_MIN(armci_user_threads.avail, armci_user_threads.max);
    thread_id_t id = THREAD_ID_SELF();

    for (i = 0; i < n; i++) if (id == armci_user_threads.ids[i]) {
        /*PRNDBG2("thread id=%ld already registered, idx=%d\n", id, i);*/
        return i;
    }

    /* see this thread for the first time */
    return armci_register_thread(id);
}

INLINE int armci_register_thread(thread_id_t id)
{
    int i;

    THREAD_LOCK(armci_user_threads.lock);

    i = armci_user_threads.avail;
    armci_user_threads.avail++;

    THREAD_UNLOCK(armci_user_threads.lock);

    if (i < armci_user_threads.max)
        armci_user_threads.ids[i] = id;
    else
        armci_die("armci_thread_idx: too many threads, adjust ARMCI_MAX_THREADS",
                  armci_user_threads.avail);

    PRNDBG2("registered a new thread: idx=%d, id=%ld\n", i, id);
    return i;
}

