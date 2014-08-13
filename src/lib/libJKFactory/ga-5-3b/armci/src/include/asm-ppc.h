/** Atomic instructions for ppc. To be populated as need arises.
 * @author Sriram Krishnamoorthy
 */
#ifndef __ATOMICS_PPC__
#define __ATOMICS_PPC__

#if HAVE_ASSERT_H
#   include <assert.h>
#endif

#define v4b (volatile unsigned int *)

/* sriram's original asm (didn't work) */
static inline void atomic_exchange(void *val, void *ptr, int size) {
    int ret;
    assert(size==4);
    assert(sizeof(unsigned int)==4);
    assert((((unsigned)ptr)&3)==0); /*make sure it is 4-byte aligned*/

    asm volatile(
        "loop: lwarx   %[ret],0,%[ptr] \n\t" /*Load & reserve*/
        "      stwcx.  %[val],0,%[ptr] \n\t" /*Store if still reserved*/
        "      bne-    loop            \n\t" /*Loop otherwise*/
        : [ret]"=&r"(ret)
        : [val]"r"(*v4b(val)), [ptr]"r"(v4b(ptr))
        : "memory", "cc"
        );
    *v4b(val) = ret;
}

/* adapted from "Synchronising C/C++ and POWER" by Sarkar et al, appearing in
 * PLDI'12 */
static inline void acquire_spinlock(void *lock) {
    int temp;
    int free = 0;
    int held = 1;

    asm volatile(
        "0: lwarx  %0,0,%3   \n" /*load-reserve lock into temp*/
        "   cmpw   %1,%0     \n" /*lock is free?*/
        "   bne-   0b        \n" /*loop if lock not free*/
        "   stwcx. %2,0,%3   \n" /*store if still reserved*/
        "   bne-   0b        \n" /*loop if lost reservation*/
        "   isync            \n" /*import barrier*/
        : "+r"(temp)
        : "r"(free), "r"(held), "r"(lock)
        : "cr0", "memory" );
}

/* adapted from "Synchronising C/C++ and POWER" by Sarkar et al, appearing in
 * PLDI'12 */
static inline void release_spinlock(void *lock) {
    int free = 0;
    assert(sizeof(unsigned int)==4);
    assert((((unsigned)lock)&3)==0); /*make sure it is 4-byte aligned*/

    /* TODO: lwsync might not be strong enough, consider 'sync' */
    asm volatile(
        "lwsync       \n" /*export barrier*/
        "stw %0,0(%1) \n" /*normal store to release lock*/
        :
        : "r"(free), "r"(lock)
        : "memory" );
}

/* adapted from "Synchronising C/C++ and POWER" by Sarkar et al, appearing in
 * PLDI'12 */
static inline int fetch_and_add(void *addr, int inc) {
    int ret;
    int tmp;

    asm volatile(
        "0: lwarx  %0,0,%2   \n" /*load-reserve addr into tmp*/
        "   mr     %1,%0     \n" /*copy tmp into ret*/
        "   add    %0,%3,%0  \n" /*add inc to addr, store in addr*/
        "   stwcx. %0,0,%2   \n" /*store if still reserved*/
        "   bne-   0b        \n" /*loop if lost reservation*/
        : "+r"(tmp), "+r"(ret)
        : "r"(addr), "r"(inc) );
}

/* from http://www.ibm.com/developerworks/rational/library/inline-assembly-C-Cpp-guide/index.html */
/* this did NOT work */
static inline int acquireLock(int *lock){
    int returnvalue = 0;
    int lockval;
    asm volatile (
    "0: lwarx %0,0,%2  \n" //load lock and reserve
    "   cmpwi 0,%0,0   \n" //compare the lock value to 0
    "   bne 1f         \n" //not 0 then exit function
    "   ori %0,%0,1    \n" //set the lock to 1
    "   stwcx. %0,0,%2 \n" //try to acquire the lock
    "   bne 0b         \n" //reservation lost, try again
    "   sync           \n" //import barrier
    "   ori  %1,%1,1   \n" //set the return value to true
    "1:                \n" //didn't get lock, return false
    : "+r" (lockval), "+r" (returnvalue)
    : "r"(lock)            //parameter lock is an address
    : "cr0" );             //cmpwi, stwcx both clobber cr0
   return returnvalue;
}

#undef v4b

#endif /* __ATOMICS_PPC__ */
