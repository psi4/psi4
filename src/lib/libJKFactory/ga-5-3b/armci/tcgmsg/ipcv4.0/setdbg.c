#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Header: /tmp/hpctools/ga/tcgmsg/ipcv4.0/setdbg.c,v 1.4 1995-02-24 02:17:42 d3h325 Exp $ */

#include "sndrcv.h"
#include "sndrcvP.h"

void SETDBG_(value)
    long *value;
/*
  set global debug flag for this process to value
*/
{
  SR_debug = *value;
}

