#if HAVE_CONFIG_H
#   include "config.h"
#endif

/*$Id: getmem.c,v 1.2 1995-02-02 23:24:10 d3g681 Exp $*/
extern char * memalign();

#if (defined(AIX) || defined(NEXT) || defined(HPUX)) && !defined(EXTNAME)
#define getmem_ getmem
#endif

#if defined(CRAY) || defined(ARDENT)
#define getmem_ GETMEM
#endif

/* getmem gets n real*8 storage locations and returns its
   address (iaddr) and offset (ioff) within the real*8 array work
   so that the usable memory is (work(i+ioff),i=1,n).
   e.g. 
        call getmem(n,work,iaddr,ioff)
        if (iaddr.eq.0) call error

   Mods are needed to release this later. */

void getmem_(pn,pwork,paddr,pioff)
     unsigned long *pn,*paddr,*pioff;
     double *pwork;
{
  double *ptemp;
  unsigned int size = 8;

  ptemp = (double *) memalign(size, (unsigned) size* *pn);
  *paddr = (unsigned long) ptemp;
  *pioff = ptemp - pwork;
}
