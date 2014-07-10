/** Atomic instructions for i386. To be populated as need arises. 
 * @author Sriram Krishnamoorthy
 */
#ifndef __ATOMICS_I386__
#define __ATOMICS_I386__

#if HAVE_ASSERT_H
#   include <assert.h>
#endif

#define v4b (volatile unsigned int *)

static inline void atomic_exchange(void *val, void *ptr, int size) {
  assert(size == 4);
  __asm__ __volatile__ ("xchgl %0, %1"
			: "=r"(*v4b(val)), "+m"(*v4b(ptr))
			: "0"(*v4b(val))
			: "memory");
}

/*SK: fixme. only available in 486+. This breaks i386 compatibility.
  atomic: *(type*)ploc = *(type*)prem; *(type*)prem += extra*/ 
static inline void atomic_fetch_and_add(void *prem, void *ploc, int extra, int size) {
  int _a_temp;
  assert(size == 4);
#if 0
  *(int*)ploc = __sync_fetch_and_add((int*)prem, extra);
#elif 0
  __asm__ __volatile__ ("movq %2, %%rax; \
                         lock xaddl %0, (%%rax);"
			: "=r"(_a_temp)
			: "0"(extra), "m"((int*)prem)
			: "memory", "rax");
  *(int *)ploc = _a_temp;
#else
  __asm__ __volatile__ ("lock xaddl %0, (%2)"
                        : "=r"(_a_temp)
                        : "0"(extra), "r"((int*)prem)
                        : "memory");
  *(int *)ploc = _a_temp;
#endif
}

#undef v4b

#endif /*__ATOMICS_I386__*/

