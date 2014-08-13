#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Header: /tmp/hpctools/ga/tcgmsg/ipcv4.0/brdcst.c,v 1.6 2002-07-17 17:20:11 vinod Exp $ */

#include "sndrcv.h"

#include <stdio.h>
#include "sndrcvP.h"

void BRDCST_(type, buf, lenbuf, originator)
     long *type;
     void *buf;
     long *lenbuf;
     long *originator;
/*
  broadcast buffer to all other processes from process originator
  ... all processes call this routine specifying the same
  orginating process.

  Optimized for communicating clusters of processes ... broadcast
  amoung cluster masters, and then amoung slaves in a cluster.
*/
{
  long me = NODEID_();
  long master = SR_clus_info[SR_clus_id].masterid;
  long nslave = SR_clus_info[SR_clus_id].nslave;
  long slaveid = me - master;
  long synch = 1;
  long lenmes, from, up, left, right;

  /* Process zero is at the top of the broadcast tree */

  if ((me == *originator) && (me != 0)) {
    long zero = 0;
    SND_(type, buf, lenbuf, &zero, &synch);
  }
  else if ((*originator != 0) && (me == 0)) {
    RCV_(type, buf, lenbuf, &lenmes, originator, &from, &synch);
  }

  if ((*originator != 0) && (SR_n_proc == 2)) return;	/* Special case */

  /* Broadcast amoung cluster masters */

  if (me == master) {
    up    = (SR_clus_id-1)/2;
    left  = 2*SR_clus_id + 1;
    right = 2*SR_clus_id + 2;
    up = SR_clus_info[up].masterid;
    left = (left < SR_n_clus) ? SR_clus_info[left].masterid : -1;
    right = (right < SR_n_clus) ? SR_clus_info[right].masterid : -1;
  
    if (me != 0)
      RCV_(type, buf, lenbuf, &lenmes, &up, &from, &synch);
    if (left > 0)
      SND_(type, buf, lenbuf, &left, &synch);
    if (right > 0)
      SND_(type, buf, lenbuf, &right, &synch);
  }

  /* Broadcast amoung local slaves */

  up    = master + (slaveid-1)/2;
  left  = master + 2*slaveid + 1;
  right = master + 2*slaveid + 2;

  if (me != master)
    RCV_(type, buf, lenbuf, &lenmes, &up, &from, &synch);
  if (left < (master+nslave))
    SND_(type, buf, lenbuf, &left, &synch);
  if (right < (master+nslave))
    SND_(type, buf, lenbuf, &right, &synch);
}  
