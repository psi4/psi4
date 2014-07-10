/* Excerpt from:
 *
 * Implementing Spinlocks
 * on the Intel Itanium Architecture
 * and PA-RISC
 * 
 * Tor Ekqvist and David Graves
 * 
 * Version 1.0
 * 
 * June 30, 2003
 * 
 * Copyright Hewlett-Packard Company 2003
 *
 * Notice
 * 
 * The software examples and their documentation have not been subjected
 * quality control and are not a Hewlett-Packard Company product. This
 * information is provided as a guide to allow you to create your own
 * spinlock implementation, which you would maintain and own.
 * 
 * The software and documentation is provided "as is". HP does not
 * warrant that the software or documentation is error free. HP
 * disclaims all warranties, express and implied, with regard to the
 * software and the documentation. HP specifically disclaims all
 * warranties of merchantability and fitness for a particular purpose.
 * 
 * Hewlett-Packard Company will not in any event be liable for any
 * direct, indirect, special, incidental or consequential damages
 * (including lost profits) related to any use, reproduction,
 * modification, or distribution of the software or documentation.
 */
#include <ia64/sys/inline.h>

#define LOCKED      1 /* spinlock is locked */
#define UNLOCKED    0 /* spinlock is free */
#define LOCK_ERROR  2 /* spinlock is in a bad state */

/* Memory fence to be used, when consistency needs to be
 * guaranteed.
 */
#define MEMORY_FENCE _Asm_mf()

/* Compile fence that can be used to make sure that the
 * compiler does not reorder memory operations around the
 * fence
 */
#define COMP_FENCE _Asm_fence(_UP_MEM_FENCE | _DOWN_MEM_FENCE)

/* get_lock attempts to quickly get a spoinlock by doing a
 * atomic 'xchg' on the word containing the lock
 * If the lock attempt is successfull, the value
 * 'LOCKED' will be returened, indicating that thie lock
 * is now loked to me, the calling process
 *
 * Note, that the 'xchg' instruction will blindly set the
 * memory location to the 'LOCKED' state, even if the
 * lock word contains an illegal value!
 *
 * The xchg instruction always has the 'acquire' semantics,
 * so the necessary memory ordering is enforced by this
 * instruction.
 *
 * this function returns '1', if the lock was free,
 * '0' if the lock was not free.
 */
#define get_lock_1(lock)                                \
(                                                       \
    _Asm_xchg((_Asm_sz)_SZ_W, (int *)lock,              \
        (int)LOCKED, (_Asm_ldhint)_LDHINT_NONE) > 0 ?   \
        (int)UNLOCKED : (int)LOCKED                     \
)

/* free_lock_1 frees a spinlock by executing a ordered store
 * with the release semantics on the word containing the
 * lock.
 *
 * The only possible value of the memory location (int *)lock
 * is really 'LOCKED' when this macro is called, as the calling
 * process/thread should have obtained teh lock by a call to
 * 'get_lock', setting the lock to 'LOCKED'.
 *
 * Just to make sure, we check for the correct value, to make
 * sure that no other process or thread has messed up the
 * shared memory location containing the lock.
 */
#define free_lock(lock)                 \
(                                       \
    ( *lock == LOCKED ) ?               \
    (_Asm_st_volatile((_Asm_sz)_SZ_W,   \
    (_Asm_sthint)_STHINT_NONE,          \
    (int *)lock,                        \
    UNLOCKED ),UNLOCKED) : LOCK_ERROR   \
)

/* A lock free macro using the 'fetchadd' instruction
 * to release the lock instead of a volatile store
 * This assumes, that a free lock is indicated by the
 * integer value '0'
 */
#define free_lock_fetchadd(lock)                            \
(                                                           \
    ( *lock == LOCKED ) ?                                   \
    (_Asm_fetchadd((_Asm_fasz)_SZ_W, (_Asm_sem)_SEM_REL,    \
    (int *)lock, -1, (_Asm_ldhint)_LDHINT_NONE ), UNLOCKED) \
    : LOCK_ERROR                                            \
)

/* Compile fence to use with the store to application
 * register, in order to make it possible for the
 * compiler to more efficiently schedule instructions
 */
#define FENCE                                       \
(_Asm_fence)(_UP_CALL_FENCE | _UP_SYS_FENCE |       \
            _DOWN_CALL_FENCE | _DOWN_SYS_FENCE )

/* The cas_acq macro performs an atomic compare and
 * exchange instruction on a 64-bit memory location
 * If the compare succeedes, the value '1' is returned
 * otherwise a '0' is returned
 */
#define cas_acq(location, currvalue, newvalue)          \
(                                                       \
    _Asm_mov_to_ar((_Asm_app_reg)_AREG_CCV,             \
        (long)currvalue, FENCE),                        \
    _Asm_cmpxchg((_Asm_sz)_SZ_D, (_Asm_sem)_SEM_ACQ,    \
        (long *)location, (long)newvalue,               \
        (_Asm_ldhint)_LDHINT_NONE)                      \
)

/* The cas_rel macro performs an atomic comapare and
 * exchange instruction on a 64-bit memory location.
 * If the compare succeedes, '1' is returned, otherwise
 * a '0' is returned.
 * The 'release' semantincs is used for the cmpxchg
 * instruction.
 */
#define cas_rel(location, currvalue, newvalue)          \
(                                                       \
    _Asm_mov_to_ar((_Asm_app_reg)_AREG_CCV,             \
        (long)currvalue, FENCE),                        \
    _Asm_cmpxchg((_Asm_sz)_SZ_D, (_Asm_sem)_SEM_REL,    \
        (long *)location, (long)newvalue,               \
        (_Asm_ldhint)_LDHINT_NONE)                      \
)

/* The cas_mf macro performs an atomic compare and exchange
 * on a memory location. The insttruction is fenced by
 * a 'mf' instruction, making it exepsive to execute,
 * but safe to execute in situations, where it is
 * impossible to make a choice between the acquire and
 * releas semantics.
 */
#define cas_mf(location, currvalue, newvalue)       \
(                                                   \
    _Asm_mov_to_ar((_Asm_app_reg)_AREG_CCV,         \
        (long)currvalue, FENCE),                    \
    _Asm_mf(),                                      \
    _Asm_cmpxchg((_Asm_sz)SZ_D, (_Asm_sem)_SEM_ACQ, \
        (long *)location, (long)newvalue,           \
        (_Asm_ldhint)_LDHINT_NONE)                  \
)

/* The count_acq macro adds a value to the specified
 * memory location with acquire semantics. Only the
 * values 0,1,2,4,8 and 16 are supported (positive and
 * negative */
#define count_acq(location, value)                          \
( _Asm_fetchadd( (_Asm_fasz)_FASZ_D, (_Asm_sem)_SEM_ACQ,    \
    (long *)location, -1,                                   \
    (_Asm_ldhint)_LDHINT_NONE )                             \
)

/* The count_rel macro adds a value to the specified
 * memory location with release semantics. Only the
 * values 0,1,2,4,8 and 16 are supported (positive and
 * negative */
#define count_rel(location, value)                          \
( _Asm_fetchadd( (_Asm_fasz)_FASZ_D, (_Asm_sem)_SEM_REL,    \
    (long *)location, (int)value,                           \
    (_Asm_ldhint)_LDHINT_NONE )                             \
)

/* The count_mf macro adds a value to the specified
 * memory location with a memory fence. Only the
 * values 0,1,2,4,8 and 16 are supported (positive and
 * negative */
#define count_mf(location, value)                           \
( _Asm_mf(),                                                \
    _Asm_fetchadd( (_Asm_fasz)_FASZ_D, (_Asm_sem)_SEM_REL,  \
    (long *)location, (int)value,                           \
    (_Asm_ldhint)_LDHINT_NONE                               \
)

/* spin_get tries to get the spinlock in the location defined
 * by the argument 'lock'. The number of spins to spin on the
 * lock before giving up is passed as the second argument
 * THe value of the spin count is dependent on your environment,
 * on a single CPY system, it is proabaly not worth spinning
 * at all
 */
__inline int spin_get(int *lock, int spin)
{
int i;
for ( i = spin; i > 0; i-- ) {
    /* Test the lock before actually try to exchange it's value,
     * to reduce memory traffic for highly contended locks
     */
    if (*lock == UNLOCKED && get_lock_1(lock) == LOCKED )
        return LOCKED;
    }
return UNLOCKED;
}
