#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Header: /tmp/hpctools/ga/tcgmsg/ipcv4.0/mdtoi.c,v 1.4 1995-02-24 02:17:24 d3h325 Exp $ */

#include "sndrcv.h"

/*
  These routines use C's knowledge of the sizes of data types
  to generate a portable mechanism for FORTRAN to translate
  between bytes, integers and doubles. Note that we assume that
  FORTRAN integers are the same size as C longs.
*/

long MDTOI_(n)
     long *n;
/*
  Return the minimum no. of integers which will hold n doubles.
*/
{
  if (*n < 0)
    Error("MDTOI_: negative argument",*n);

   return (long) ( (MDTOB_(n) + sizeof(long) - 1) / sizeof(long) );
}
