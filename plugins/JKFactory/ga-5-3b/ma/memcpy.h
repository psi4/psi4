/** @file
 * Private header file containing symbolic constants, type declarations,
 * and macro definitions for OS memory routines, to provide a level of
 * abstraction between them and routines that use them.
 *
 * This file should only be included by internal C files.
 */
#ifndef _memcpy_h
#define _memcpy_h

/**
 ** constants
 **/

/* ensure that NULL is defined */
#ifndef NULL
#define NULL 0
#endif

/**
 ** macros
 **/

/* allocate bytes */
#define bytealloc(nbytes)    malloc((unsigned long)(nbytes))

/* deallocate bytes */
#define bytefree(pointer)    (void)free((char *)(pointer))

#if HAVE_STRING_H
#   include <string.h>
#else
#   ifndef HAVE_STRCHR
#       define strchr index
#       define strrchr rindex
#   endif
    char *strchr (), *strrchr ();
#   ifndef HAVE_MEMCPY
#       define memcpy(d, s, n) bcopy ((s), (d), (n))
#       define memmove(d, s, n) bcopy ((s), (d), (n))
#   endif
#endif

#define bytecopy(from,to,nbytes) \
    ((void)memcpy((char *)(to), (char *)(from), (int)(nbytes)))

#endif /* _memcpy_h */
