#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Header: /tmp/hpctools/ga/tcgmsg/ipcv4.0/nxtsrv.c,v 1.4 1995-02-24 02:17:33 d3h325 Exp $ */

#include <stdio.h>
#include "sndrcv.h"

#define DEBUG_ 0
#define TYPE_NXTVAL 32767

int main()
/*
  This runs as highest no. process(es)
*/	   
{
  long cnt     = 0;            /* actual counter */
  long lencnt  = sizeof cnt;   /* length of cnt */
  long ndone   = 0;            /* no. finished for this loop */
  long ntermin = 0;            /* no. terminated so far (pend) */
  long node    = -1;           /* select any node */
  long type    = TYPE_NXTVAL;  /* message type */
  long buf[2];                 /* buffer to get values */
  long lenbuf  = sizeof buf;   /* length of buffer */
  long mproc;                  /* no. of processes running loop */
  long nval;                   /* no. of values requested */
  long done_list[16384];       /* list of processes finished with this loop */
  long sync = 1;               /* all info goes synchronously */
  long on=0;
  long lenmes, nodefrom;


  PBEGIN_();
  SETDBG_(&on);

  while (1) {

    /* Wait for input from any node */
    
    RCV_(&type, (char *) buf, &lenbuf, &lenmes, &node, &nodefrom, &sync);

    if (lenmes != lenbuf) 
      Error("NextValueServer: lenmes != lenbuf", lenmes);


    mproc = buf[0];
    nval = buf[1];
    if (DEBUG_)
      (void) printf("NVS: from=%d, mproc=%d, ndone=%d, ntermin=%d\n",
		    nodefrom, mproc, ndone, ntermin);

    if (mproc == 0) {

      /* Sending process is about to terminate. Send reply and disable
	 sending to him. If all processes have finished return. */

      SND_(&type, (char *) &cnt, &lencnt, &nodefrom, &sync);

      if (++ntermin == NNODES_())
	return 0;
    }
    else if (mproc > 0) {
      
      /* This is what we are here for */

      SND_(&type, (char *) &cnt, &lencnt, &nodefrom, &sync);
      cnt += nval;
    }
    else if (mproc < 0) {

      /* This process has finished the loop. Wait until all mproc
	 processes have finished before releasing it */

      done_list[ndone++] = nodefrom;

      if (ndone == -mproc) {
	while (ndone--) {
	  nodefrom = done_list[ndone];
	  SND_(&type, (char *) &cnt, &lencnt, &nodefrom, &sync);
	}
	cnt = 0;
	ndone = 0;
      }
    }
  }
  return 0;
}
