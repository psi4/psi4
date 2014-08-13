#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Header: /tmp/hpctools/ga/tcgmsg/ipcv4.0/waitall.c,v 1.3 1995-02-24 02:18:05 d3h325 Exp $ */

#include <stdio.h>
#if defined(SUN) || defined(ALLIANT) || defined(ENCORE) || defined(SEQUENT) || \
    defined(AIX) || defined(NEXT)    || defined(DECOSF) || defined(LINUX)
#include <sys/wait.h>
#endif

int WaitAll(nchild)
     long nchild;
/*
  Wait for all children to finish and return appropriate status
  0 = OK
  1 = bad news
*/
{
  int status, pid, child, stat=0, lo, hi;
  
#if defined(ALLIANT) || defined(ENCORE) || defined(SEQUENT) || defined(NEXT)
  union wait ustatus;
#endif

  for (child=0; child<nchild; child++) {
#if defined(ALLIANT) || defined(ENCORE) || defined(SEQUENT) || defined(NEXT)
      pid = wait(&ustatus);
      status = ustatus.w_status;
#else
      pid = wait(&status);
#endif
/*
(void) printf("Child finished pid=%d, status=0x%x\n",pid, status);
(void) fflush(stdout);
*/
      if ( pid < 0 ) {
        (void) fprintf(stderr,"WaitAll: No children or error in wait?\n");
		return 1;
      }

      if (status != 0) {
	(void) fprintf(stderr, "WaitAll: Child (%d) finished, status=0x%x ",
		       pid, status);

	lo = status & 0xff;
	hi = (status >> 8) & 0xff;
	
	if ( lo == 0177 )
	  (void) fprintf(stderr, "(stopped by signal %d).\n", hi);
	else if ( (lo != 0) && (lo & 0200) )
	  (void) fprintf(stderr, "(killed by signal %d, dumped core).\n", 
			 lo & 0100);
	else if ( lo != 0 )
	  (void) fprintf(stderr, "(killed by signal %d).\n",lo);
	else
	  (void) fprintf(stderr, "(exited with code %d).\n",hi);

	(void) fflush(stderr);
	stat = 1;
      }

    }

  return stat;
}
