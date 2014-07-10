#if HAVE_CONFIG_H
#   include "config.h"
#endif

/** @file
 * These routines use C's knowledge of the sizes of data types
 * to generate a portable mechanism for FORTRAN to translate
 * between bytes, integers and doubles. 
 */

#include "sndrcv.h"

/**
 * Return the no. of bytes that n doubles occupy
 */
long MDTOB_(long *n)
{
    if (*n < 0) {
        Error("MDTOB_: negative argument",*n);
    }

    return (long) (*n * sizeof(double));
}


/**
 * Return the minimum no. of integers which will hold n doubles.
 */
long MDTOI_(long *n)
{
    if (*n < 0) {
        Error("MDTOI_: negative argument",*n);
    }

    return (long) ( (MDTOB_(n) + sizeof(long) - 1) / sizeof(long) );
}


/**
 * Return the no. of bytes that n ints=longs occupy
 */
long MITOB_(long *n)
{
    if (*n < 0) {
        Error("MITOB_: negative argument",*n);
    }

    return (long) (*n * sizeof(long));
}


/**
 * Return the minimum no. of doubles in which we can store n longs
 */
long MITOD_(long *n)
{
    if (*n < 0) {
        Error("MITOD_: negative argument",*n);
    }

    return (long) ( (MITOB_(n) + sizeof(double) - 1) / sizeof(double) );
}
