#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif

/** @file
 * This checks the functioning of the include file gaf2cp.h.
 */

#include "farg.h"
#include "typesf2c.h"

#define ARGLIMIT 256

#define PARG_ F77_FUNC(parg,PARG)

void PARG_()
{
    Integer i;
    Integer limit = F2C_IARGC();
    printf("argc=%d\n", (int)limit);

    for (i=0; i<limit; i++) {
        char fstring[ARGLIMIT];
        char cstring[ARGLIMIT];
        F2C_GETARG(&i, fstring, ARGLIMIT);
        ga_f2cstring(fstring, ARGLIMIT, cstring, ARGLIMIT);
        (void) printf("argv(%ld)=%s\n", (long)i, cstring);
    }
}
