#if HAVE_CONFIG_H
#   include "config.h"
#endif

/** @file
 *
 * $Header: /tmp/hpctools/ga/tcgmsg/ipcv5.0/mitod.c,v 1.1 1997-03-05 18:42:31 d3e129 Exp $
 *
 * These routines use C's knowledge of the sizes of data types
 * to generate a portable mechanism for FORTRAN to translate
 * between bytes, integers and doubles. Note that we assume that
 * FORTRAN integers are the same size as C longs.
 */
#include "sndrcv.h"


/**
 * Return the minimum no. of doubles in which we can store n longs
 */
long MITOD_(long *n)
{
    if (*n < 0)
        Error("MITOD_: negative argument",*n);

    return (long) ( (MITOB_(n) + sizeof(double) - 1) / sizeof(double) );
}
