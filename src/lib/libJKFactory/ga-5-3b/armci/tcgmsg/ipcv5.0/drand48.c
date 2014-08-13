#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Header: /tmp/hpctools/ga/tcgmsg/ipcv5.0/drand48.c,v 1.2 1994-12-30 20:55:40 d3h325 Exp $ */

#include "srftoc.h"

extern long random();
extern int srandom();

double DRAND48_()
{
    return ( (double) random() ) * 4.6566128752458e-10;
}

void SRAND48_(seed)
    unsigned *seed;
{
    (void) srandom(*seed);
}
