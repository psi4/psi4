#if HAVE_CONFIG_H
#   include "config.h"
#endif

/*$Id: string-util.c,v 1.2 1995-02-02 23:18:32 d3g681 Exp $*/
/*
 * string-util.c
 */

/*
 * <Global Comment 1> :
 *
 * The standard string functions (strcat, strcmp, strcpy, strlen, etc.)
 * dereference their pointer arguments without first checking to see if
 * the arguments are NULL (0).  On some machines (e.g., the VAX) this is
 * harmless, but on others (e.g., the DECstation 3100) it causes a
 * segmentation fault.
 *
 * The string utilities defined here provide functionality similar to that
 * of their standard counterparts, but they are careful not to dereference
 * arguments which are NULL.
 *
 * <Global Comment 2> :
 *
 * A string is of type "pointer to char" (char *); a string is terminated
 * by a byte whose value is 0 (a null character).  The term "Nstring" will
 * be used to mean "NULL string" (pointer value of 0); "Zstring" will be
 * used to mean "zero-length string" (nonzero pointer value, but no nonzero
 * bytes pointed to).  <Global Comment 1> is concerned with Nstrings.
 */

#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif
#if HAVE_STRING_H
#   include <string.h>
#endif
#include "string-util.h"

/**
 ** public routines
 **/

/* ------------------------------------------------------------------------- */
/*
 * Return the length of (i.e., number of nonzero bytes in) the given string.
 */
/* ------------------------------------------------------------------------- */

unsigned int str_len(s)
    char        *s;        /* string */
{
    int            length = 0;

    /* see <Global Comment 1> */
    if (!s)
    {
        /* length of Nstring is 0 */
        return(0);
    }

    while (*s++)
    {
        length++;
    }

    return(length);
}

/* ------------------------------------------------------------------------- */
/*
 * Return the index of the string in slist for which s is a match; s "matches"
 * a string t in slist if (a) s and t are equal, or (b) s is a prefix of t but
 * is not a prefix of any other string in slist.
 *
 * If s matches no string in slist, return SM_NONE.  If s matches more than
 * one string in slist, return SM_MANY.
 *
 * Restriction: slist can contain no Nstrings.
 */
/* ------------------------------------------------------------------------- */

int str_match(s, slist, n)
    char        *s;        /* string to match */
    char        *slist[];    /* list of strings to search */
    unsigned int    n;        /* # of strings in slist */
{
    size_t     i;        /* loop index */
    size_t     length;        /* of s */
    int        match;        /* index of string in slist matched by s */

    /* see <Global Comment 1> */
    if (!s)
    {
        return(SM_NONE);
    }

    /* we know s is not an Nstring */
    length = strlen(s);
    match = -1;

    for (i = 0; i < n; i++)
    {
        /* s, slist[i] are not Nstrings */
        if (!strncmp(s, slist[i], length))
        {
            /* s is at least a prefix */
            if (length == strlen(slist[i]))
            {
                /* exact match */
                return(i);
            }

            /* now we know s is a proper prefix */
            if (match < 0)
            {
                /* first match we've seen */
                match = i;
            }
            else
            {
                /* s is a proper prefix of more than one string in slist */
                return(SM_MANY);
            }
        }
    }

    if (match < 0)
    {
        /* s is not a prefix of any string in slist */
        return(SM_NONE);
    }
    else
    {
        /* s is a proper prefix of exactly one string in slist */
        return(match);
    }
}
