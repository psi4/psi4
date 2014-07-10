/**
 * @file spinlock.h
 *
 * This file attempts to implement spin locks for various platforms and/or CPU
 * instruction sets. 
 */
#ifndef SPINLOCK_H
#define SPINLOCK_H

#define DEBUG_SPINLOCK 0

#define OPENPA 0

#if OPENPA
#   if DEBUG_SPINLOCK
#       warning SPINLOCK: openpa
#   endif
#   define SPINLOCK
#   include "opa_primitives.h"
#   define LOCK_T OPA_int_t
#   define TESTANDSET(x) OPA_swap_int((x), 1)
#   define MEMORY_BARRIER OPA_read_write_barrier

#elif (defined(PPC) || defined(__PPC__) || defined(__PPC))
#   if DEBUG_SPINLOCK
#       warning SPINLOCK: PPC
#   endif
#   define SPINLOCK
#   include "asm-ppc.h"
//#   define TESTANDSET testandset
//#   define TESTANDSET acquireLock
#   define armci_acquire_spinlock acquire_spinlock
#   define armci_release_spinlock release_spinlock
#   define MEMORY_BARRIER memory_barrier
static int testandset(void *spinlock) { 
    int v=1;
    atomic_exchange(&v,spinlock,sizeof(int));
    return v;
}
static void memory_barrier() {
    __asm__ __volatile__ ("sync" : : : "memory");
}

#elif defined(__i386__) || defined(__x86_64__)
#   if DEBUG_SPINLOCK
#       warning SPINLOCK: x86_64
#   endif
#   define SPINLOCK
#   include "atomics-i386.h"
static int testandset(void *spinlock) { 
    int v=1;
    atomic_exchange(&v,spinlock,sizeof(int));
    return v;
}
#   define TESTANDSET testandset

#elif defined(__ia64)
#   if DEBUG_SPINLOCK
#       warning SPINLOCK: ia64
#   endif
#   define SPINLOCK
#   include "atomic_ops_ia64.h"
static int testandset(void *spinlock) { 
    int val=1;
    int res;
    atomic_swap_int(spinlock, val, &res);
    return res;
}
#   define TESTANDSET testandset

#elif defined(DECOSF)
#   if DEBUG_SPINLOCK
#       warning SPINLOCK: DECOSF
#   endif
#   error "no implementation"

#elif defined(SGI)
#   if DEBUG_SPINLOCK
#       warning SPINLOCK: SGI
#   endif
#   include <mutex.h>
#   define SPINLOCK 
#   define TESTANDSET(x) __lock_test_and_set((x), 1)
#   define RELEASE_SPINLOCK __lock_release

/*#elif defined(AIX)*/
#elif HAVE_SYS_ATOMIC_OP_H
#   if DEBUG_SPINLOCK
#       warning SPINLOCK: sys/atomic_op.h (AIX)
#   endif
#   include <sys/atomic_op.h>
#   define SPINLOCK 
#   define TESTANDSET(x) (_check_lock((x), 0, 1)==TRUE) 
#   define RELEASE_SPINLOCK(x) _clear_lock((x),0) 

#elif defined(SOLARIS)
#   if DEBUG_SPINLOCK
#       warning SPINLOCK: SOLARIS
#   endif
#   include <sys/atomic.h>
#   include <sys/machlock.h>
#   define SPINLOCK  
#   define TESTANDSET(x) (!_lock_try((x))) 
#   define RELEASE_SPINLOCK _lock_clear 

#elif defined(MACX)

#elif defined(HPUX__)
#   if DEBUG_SPINLOCK
#       warning SPINLOCK: HPUX__
#   endif
extern int _acquire_lock();
extern void _release_lock();
#   define SPINLOCK  
#   define TESTANDSET(x) (!_acquire_lock((x))) 
#   define RELEASE_SPINLOCK _release_lock 

#elif defined(HPUX) && defined(__ia64) /* HPUX on IA64, non gcc */
#   if DEBUG_SPINLOCK
#       warning SPINLOCK: HPUX ia64
#   endif
#   define SPINLOCK
typedef unsigned int slock_t;
#   include <ia64/sys/inline.h>
#   define TESTANDSET(lock) _Asm_xchg(_SZ_W, lock, 1, _LDHINT_NONE)
#   define RELEASE_SPINLOCK(lock) (*((volatile LOCK_T *) (lock)) = 0)

#elif defined(NEC)
#   if DEBUG_SPINLOCK
#       warning SPINLOCK: NEC
#   endif
extern ullong ts1am_2me();
#   define LOCK_T ullong
#   define _LKWD (1ULL << 63)
#   define SPINLOCK  
#   define TESTANDSET(x) ((_LKWD & ts1am_2me(_LKWD, 0xffULL, (ullong)(x))))
#   define MEMORY_BARRIER mpisx_clear_cache 
extern void mpisx_clear_cache();
#   define RELEASE_SPINLOCK(x) ts1am_2me(0ULL, 0xffULL, (ullong)x); 

#endif

#ifdef SPINLOCK

#if DEBUG_
#   if HAVE_STDIO_H
#       include <stdio.h>
#   endif
#endif

#if HAVE_UNISTD_H
#   include <unistd.h>
#endif

#ifndef DBL_PAD
#   define DBL_PAD 16
#endif

/* make sure that locks are not sharing the same cache line */
typedef struct{
    double  lock[DBL_PAD];
}pad_lock_t;

#ifndef LOCK_T
#   define LOCK_T int
#endif
#define PAD_LOCK_T pad_lock_t

static inline void armci_init_spinlock(LOCK_T *mutex)
{
#if OPENPA
    OPA_store_int(mutex, 0);
#else
    *mutex =0;
#endif
}

#ifdef TESTANDSET

static inline void armci_acquire_spinlock(LOCK_T *mutex)
{
#if defined(BGML) || defined(DCMF)
    return;
#else
    int loop=0, maxloop =10;

    while (TESTANDSET(mutex)){
        loop++;
        if(loop==maxloop){ 
#   if DEBUG_
            extern int armci_me;
            printf("%d:spinlock sleeping\n",armci_me); fflush(stdout);
#   endif
            usleep(1);
            loop=0;
        }
    }
#endif
}

#ifdef RELEASE_SPINLOCK
#   ifdef MEMORY_BARRIER
#       define armci_release_spinlock(x) MEMORY_BARRIER(); RELEASE_SPINLOCK(x)
#   else
#       define armci_release_spinlock(x) RELEASE_SPINLOCK(x)
#   endif
#else
static inline void armci_release_spinlock(LOCK_T *mutex)
{
#if defined(BGML) || defined(DCMF)
    return;
#else
#   ifdef MEMORY_BARRIER
    MEMORY_BARRIER ();
#   endif
#if OPENPA
    OPA_store_int(mutex, 0);
#else
    *mutex =0;
#endif
#   ifdef MEMORY_BARRIER
    MEMORY_BARRIER ();
#   endif
#   if (defined(MACX)||defined(LINUX)) && defined(__GNUC__) && defined(__ppc__)  
    __asm__ __volatile__ ("isync" : : : "memory");
#   endif
#endif
}
#endif /* RELEASE_SPINLOCK */

#endif /* TESTANDSET */

#endif /* SPINLOCK */

#endif /* SPINLOCK_H */
