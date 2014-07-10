#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Header: /tmp/hpctools/ga/tcgmsg/ipcv4.0/nxtval.c,v 1.6 2004-04-01 02:04:57 manoj Exp $ */

#include <stdio.h>
#include <setjmp.h>
#include <signal.h>
#include <unistd.h>
#include "sndrcvP.h"
#include "sndrcv.h"

jmp_buf SR_jmp_buf;   /* Jumped to on soft error */

void NextValueServer()
/*
  This runs as process SR_n_proc and provides load balancing service.

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
  long done_list[MAX_PROCESS]; /* list of processes finished with this loop */
  long sync = 1;               /* all info goes synchronously */
  long on=0;
  long lenmes, nodefrom;

  SR_exit_on_error = FALSE; /* Want to return no matter what */

  if (setjmp(SR_jmp_buf)) {  /* Error should long jump to here */
/*    (void) printf("Error long jumped to NXTVAL ... returning.\n"); */
    SR_exit_on_error = TRUE;
    return;
  }

  SETDBG_(&on);

  while (1) {

    /* Wait for input from any node */
    
    RCV_(&type, (char *) buf, &lenbuf, &lenmes, &node, &nodefrom, &sync);

    if (lenmes != lenbuf) {
      Error("NextValueServer: lenmes != lenbuf", lenmes);
      return;   /* Never actually gets here as does long jump */
    }

    mproc = buf[0];
    nval = buf[1];
    if (DEBUG_)
      (void) printf("NVS: from=%ld, mproc=%ld, ndone=%ld, ntermin=%ld\n",
		    nodefrom, mproc, ndone, ntermin);

    if (mproc == 0) {

      /* Sending process is about to terminate. Send reply and disable
	 sending to him. If all processes have finished return. 

         Modified so that all processes block on waiting for message
         from nxtval server before terminating. nxtval only lets
         everyone go when all have registered termination. 
	 This is so that processes do not close their sockets
	 while another process is doing a RCV from any node (which
	 results in an unavoidable error condition). */

      if (++ntermin == NNODES_()) {
	(void) signal(SIGCHLD, SIG_DFL); /* Will be dying naturally */
	for (node=0; node<NNODES_(); node++) {
	  SND_(&type, (char *) &cnt, &lencnt, &node, &sync);
	  (void) close(SR_proc_info[node].sock);
	  SR_proc_info[node].sock = -1;
	}
	return;
      }
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
}

long NXTVAL_(mproc)
     long *mproc;
/*
  Get next value of shared counter.

  mproc > 0 ... returns requested value
  mproc < 0 ... server blocks until abs(mproc) processes are queued
                and returns junk
  mproc = 0 ... indicates to server that I am about to terminate

  this needs to be extended so that clusters of processes with
  shared memory collectively get a bunch of values from the server
  thus reducing the overhead of calling nextvalue.
*/
{
  long server = NNODES_();         /* id of server process */
  long buf[2];
  long lenbuf = sizeof buf;
  long type = TYPE_NXTVAL;
  long lenmes, nodefrom;
  long sync = 1;
  long result=0;

  if (SR_parallel) {
    buf[0] = *mproc;
    buf[1] = 1;

    if (DEBUG_) {
      (void) printf("%2ld: nxtval: mproc=%ld\n",NODEID_(), *mproc);
      (void) fflush(stdout);
    }

    SND_(&type, (char *) buf, &lenbuf, &server, &sync);
    RCV_(&type, (char *) buf, &lenbuf, &lenmes, &server, &nodefrom, &sync);
    result = buf[0];
  }
  else {
    /* Not running in parallel ... just do a simulation */
    static int count = 0;
    if (*mproc == 1)
      result = count++;
    else if (*mproc == -1) {
      count = 0;
      result = 0;
    }
    else
      Error("nxtval: sequential version with silly mproc ", (long) *mproc);
  }

  return result;
}
