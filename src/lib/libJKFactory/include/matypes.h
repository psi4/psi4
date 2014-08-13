/** @file
 * Private header file containing C type definitions.
 *
 * This file should only be included directly by internal C
 * header files (e.g., macdecls.h).  It may be included indirectly
 * by external C files that include the appropriate header
 * file (e.g., macdecls.h).
 */
#ifndef _MATYPES_H
#define _MATYPES_H

/**
 ** types
 **/

#include "typesf2c.h"

typedef Integer Boolean; /* MA_TRUE or MA_FALSE */
typedef char * Pointer;  /* generic pointer */

/* not all C compilers support long double */
typedef long double MA_LongDouble;

/* no C compilers support complex types */
typedef struct {float dummy[2];} MA_SingleComplex;
typedef struct {double dummy[2];} MA_DoubleComplex;
typedef struct {double dummy[4];} MA_LongDoubleComplex;

typedef long MA_AccessIndex;

#endif /* _matypes_h */
