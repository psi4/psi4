#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Header: /tmp/hpctools/ga/tcgmsg/ipcv4.0/waitcom.c,v 1.3 1995-02-24 02:18:06 d3h325 Exp $ */

#include "sndrcv.h"

/*ARGSUSED*/
void WAITCOM_(node)
  long *node;
/*
  Wait for async communications to complete ... null operation in
  the UNIX environment
*/
{
}
