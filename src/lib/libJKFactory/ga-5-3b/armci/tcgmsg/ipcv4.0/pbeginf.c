#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STRING_H
#   include <string.h>
#endif

#include "srftoc.h"
#include "sndrcv.h"


/**
 * Hewlett Packard Risc box and new SparcWorks F77 2.* compilers.
 * Have to construct the argument list by calling FORTRAN.
 */
void PBEGINF_()
{
}


/**
 * Alternative entry for those senstive to FORTRAN making reference
 * to 7 character external names
 */
void PBGINF_()
{
    PBEGINF_();
}
