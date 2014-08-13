#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Header: /tmp/hpctools/ga/tcgmsg/ipcv4.0/drand48.c,v 1.4 1995-02-24 02:17:14 d3h325 Exp $ */

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
