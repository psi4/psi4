#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif

#include "srftoc.h"
#include "tcgmsgP.h"


/**
 * Define value of debug flag
 */
void SETDBG_(long *onoff)
{
    DEBUG_ = *onoff;
}


/**
 * Print out statistics for communications ... not yet implemented
 */
void STATS_()
{
    (void) fprintf(stderr,"STATS_ not yet supported\n");
    (void) fflush(stderr);
}
