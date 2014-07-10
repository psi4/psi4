#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Header: /tmp/hpctools/ga/tcgmsg/ipcv4.0/error.c,v 1.9 2003-05-27 22:06:51 edo Exp $ */

#include <stdio.h>
#include <setjmp.h>
#include <signal.h>

#include "sndrcvP.h"
#include "sndrcv.h"
#include "signals.h"
#include "tcgsockets.h"

#if defined(SHMEM) || defined(SYSV)
#include "sema.h"
#include "tcgshmem.h"
#endif

extern jmp_buf SR_jmp_buf;   /* Jumped to on soft error */ 

#include <errno.h>

extern void exit();
extern int SR_caught_sigint;

void Error(string, integer)
     char *string;
     long integer;
{
  (void) signal(SIGCHLD, SIG_DFL); /* Death of children to be expected */
  (void) signal(SIGINT, SIG_IGN);

  (void) fflush(stdout);
  if (SR_caught_sigint) {
    (void) fprintf(stderr,"%3ld: interrupt(%d)\n",NODEID_(), SR_caught_sigint);
    (void) fflush(stderr);
  }
  else {
    (void) fprintf(stdout,"%3ld: %s %ld (%#lx).\n", NODEID_(), string,
		   integer,integer);
    (void) fflush(stdout);
    (void) fprintf(stderr,"%3ld: %s %ld (%#lx).\n", NODEID_(), string,
		   integer,integer);
    if (errno != 0)
      perror("system error message");
    if (DEBUG_)
      PrintProcInfo();
  }
  (void) fflush(stdout);
  (void) fflush(stderr);

  /* Shut down the sockets and remove shared memory and semaphores to
     propagate an error condition to anyone that is trying to communicate
     with me */

   ZapChildren();  /* send interrupt to children which should trap it
		     and call Error in the handler */

#if defined(SHMEM) || defined(SYSV)
#   if (defined(SGI_N32) || defined(SGITFP))
#       define PARTIALSPIN
#   else
#       define NOSPIN
#   endif
#endif

#if defined(SHMEM) || defined(SYSV)
#if defined(NOSPIN) || defined(PARTIALSPIN)
  (void) SemSetDestroyAll();
#endif
  (void) DeleteSharedRegion(SR_proc_info[NODEID_()].shmem_id);
#endif
  ShutdownAll();         /* Close sockets for machines with static kernel */

/*  abort(); */

  if (SR_exit_on_error)
    exit(1);
  else {
    SR_error = 1;
    (void) longjmp(SR_jmp_buf, 1); /* For NXTVAL server */
  }
}

void PARERR_(code)
   long *code;
/*
  Interface from fortran to c error routine
*/
{
  Error("User detected error in FORTRAN", *code);
}
