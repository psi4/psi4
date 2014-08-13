#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Header: /tmp/hpctools/ga/tcgmsg/ipcv4.0/srmover.c,v 1.4 1995-02-24 02:17:53 d3h325 Exp $ */

#if defined(SEQUENT) || defined(CONVEX)
#define memcpy(a ,b ,c) bcopy((b), (a), (c))
#else
#include <memory.h>
#endif

void SRmover(a, b, n)
     char *a, *b;
     long n;
/*
  Move n bytes from b to a
*/
{
#if defined(ALLIANT) || defined(IBM) || defined(IBMNOEXT) || \
    defined(CRAY)    || defined(CONVEX) || defined(APOLLO)
  /* memcpy is fast, Cray is not actually used but
     alignment crap below won't work in anycase */
  (void) memcpy(a, b, (int) n);
#else
#define UNALIGNED(a) (((unsigned long) (a)) % sizeof(int))

  if (UNALIGNED(a) || UNALIGNED(b))
    (void) memcpy(a, b, (int) n);      /* abdicate responsibility */
  else {
    /* Data is integer aligned ... move first n/sizeof(int) bytes
       as integers and the remainder as bytes */

    int ni = n/sizeof(int);
    int *ai = (int *) a;
    int *bi = (int *) b;
    int i;

#ifdef ARDENT
#pragma ivdep
#endif
    for (i=0; i<ni; i++)
      ai[i] = bi[i];

    /* Handle the remainder */

    a += ni*sizeof(int);
    b += ni*sizeof(int);
    n -= ni*sizeof(int);

    for (i=0; i<n; i++)
      a[i] = b[i];
  }
#endif
}
