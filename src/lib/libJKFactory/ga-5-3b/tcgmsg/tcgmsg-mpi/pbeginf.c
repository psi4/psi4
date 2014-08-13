#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STRING_H
#   include <string.h>
#endif

#include "farg.h"
#include "sndrcv.h"
#include "srftoc.h"
#define LEN 255

extern int    tcgi_argc;
extern char **tcgi_argv;

/**
 * Hewlett Packard Risc box, SparcWorks F77 2.* and Paragon compilers.
 * Have to construct the argument list by calling FORTRAN.
 */
void PBEGINF_()
{
    ga_f2c_get_cmd_args(&tcgi_argc, &tcgi_argv);
    tcgi_pbegin(tcgi_argc, tcgi_argv);
}


/**
 * Alternative entry for those senstive to FORTRAN making reference
 * to 7 character external names
 */
void PBGINF_()
{
    PBEGINF_();
}
