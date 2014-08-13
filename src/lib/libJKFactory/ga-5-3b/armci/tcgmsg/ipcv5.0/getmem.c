#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_MALLOC_H
#   include <malloc.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif

#define getmem_ F77_FUNC(getmem,GETMEM)

/**
 * getmem gets n real*8 storage locations and returns its
 * address (iaddr) and offset (ioff) within the real*8 array work
 * so that the usable memory is (work(i+ioff),i=1,n).
 * e.g. 
 *      call getmem(n,work,iaddr,ioff)
 *      if (iaddr.eq.0) call error
 *
 * Mods are needed to release this later.
 */
void getmem_(
        unsigned long *pn,
        double *pwork,
        unsigned long *paddr,
        unsigned long *pioff)
{
    double *ptemp;
    unsigned int size = 8;

#if HAVE_MEMALIGN
    ptemp = (double *) memalign(size, (unsigned) size* *pn);
#else
    ptemp = (double *) malloc((unsigned) size* *pn);
#endif
    *paddr = (unsigned long) ptemp;
    *pioff = ptemp - pwork;
}
