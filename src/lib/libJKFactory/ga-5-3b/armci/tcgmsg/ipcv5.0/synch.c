#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include "tcgmsgP.h"
#include "sndrcv.h"

void SYNCH_(long *ptype)
{
    long junk = 0, n = 1;

    IGOP_(ptype, &junk, &n, "+", 1);
}
