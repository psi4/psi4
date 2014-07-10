/** Atomic instructions for i386. To be populated as need arises. 
 * @author Sriram Krishnamoorthy
 */
#ifndef __ATOMICS_I386__
#define __ATOMICS_I386__

#include <assert.h>

#define v4b (volatile unsigned int *)

static inline void atomic_exchange(void *val, void *ptr, int size) {
  assert(size == 4);
  __asm__ __volatile__ ("xchgl %0, %1"
			: "=r"(*v4b(val)), "+m"(*v4b(ptr))
			: "0"(*v4b(val))
			: "memory");
}

#undef v4b

#endif /*__ATOMICS_I386__*/

