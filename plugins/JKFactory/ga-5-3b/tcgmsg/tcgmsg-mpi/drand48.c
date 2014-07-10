#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>

#include "sndrcv.h"


double DRAND48_()
{
    double val=((double) random() ) * 4.6566128752458e-10;
    return val;
}


void SRAND48_(long *seed)
{
    (void) srandom(*seed);
}
