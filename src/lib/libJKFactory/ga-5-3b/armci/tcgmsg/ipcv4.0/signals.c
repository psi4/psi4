#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Header: /tmp/hpctools/ga/tcgmsg/ipcv4.0/signals.c,v 1.11 2004-04-01 02:04:57 manoj Exp $ */

#include <signal.h>
#include "sndrcvP.h"
#if defined(SUN) || defined(ALLIANT) || defined(ENCORE) || defined(SEQUENT) || \
    defined(AIX) || defined(NEXT)
#include <sys/wait.h>
#else
#include <sys/types.h>
#include <sys/wait.h>
#endif

extern void Error();
int SR_caught_sigint = 0;

#if (defined(ENCORE) || defined(SEQUENT) || defined(ARDENT))
#   define SigType  int
#else
#   define SigType  void
#endif

#ifndef SIG_ERR
#   define SIG_ERR         (SigType (*)())-1
#endif


#if (defined(ENCORE) || defined(SEQUENT) || defined(ARDENT)) || (defined(SUN) && !defined(SOLARIS))
SigType SigintHandler(sig, code, scp, addr)
     int code;
     struct sigcontext *scp;
     char *addr;
#else
SigType SigintHandler(sig)
#endif
     int sig;
{
  SR_caught_sigint = 1;
  Error("SigintHandler: signal was caught",(long) sig);
}


void TrapSigint()
/*
  Trap the signal SIGINT so that we can propagate error
  conditions and also tidy up shared system resources in a
  manner not possible just by killing everyone
*/
{
  if ( signal(SIGINT, SigintHandler) == SIG_ERR)
       Error("TrapSigint: error from signal setting SIGINT",(long) SIGINT);
}


void ZapChildren()
/*
  kill -SIGINT all of my beloved children
*/
{
  while (SR_numchild--)
    (void) kill((int) SR_pids[SR_numchild], SIGINT);
}

#if (defined(ENCORE) || defined(SEQUENT) || defined(ARDENT)) || (defined(SUN) && !defined(SOLARIS))
SigType SigchldHandler(sig, code, scp, addr)
     int code;
     struct sigcontext *scp;
     char *addr;
#else
void SigchldHandler(sig)
#endif
     int sig;
{
  int status;
  
#if defined(ALLIANT) || defined(ENCORE) || defined(SEQUENT) || defined(NEXT)
  union wait ustatus;
#endif

#if defined(ALLIANT) || defined(ENCORE) || defined(SEQUENT) || defined(NEXT)
  (void) wait(&ustatus);
  status = ustatus.w_status;
#else
  (void) wait(&status);
#endif
  SR_caught_sigint = 1;
  Error("Child process terminated prematurely, status=",(long) status);
}

void TrapSigchld()
/*
  Trap SIGCHLD so that can tell if children die unexpectedly.
*/
{
  if ( signal(SIGCHLD, SigchldHandler) == SIG_ERR)
    Error("TrapSigchld: error from signal setting SIGCHLD", (long) SIGCHLD);
}

#if (defined(ENCORE) || defined(SEQUENT) || defined(ARDENT)) || (defined(SUN) && !defined(SOLARIS))
SigType SigsegvHandler(sig, code, scp, addr)
     int code;
     struct sigcontext *scp;
     char *addr;
#else
SigType SigsegvHandler(sig)
#endif
     int sig;
{
  SR_caught_sigint = 1;
  Error("SigsegvHandler: signal was caught",(long) sig);
}


void TrapSigsegv()
/*
  parallel needs to trap the signal SIGSEGV under Solaris 
  that is generated when interrupted in NxtVal  
*/
{
  if ( signal(SIGSEGV, SigsegvHandler) == SIG_ERR)
       Error("TrapSigsegv: error from signal setting SIGSEGV", (long) SIGSEGV);
}

#if (defined(ENCORE) || defined(SEQUENT) || defined(ARDENT)) || (defined(SUN) && !defined(SOLARIS))
SigType SigtermHandler(sig, code, scp, addr)
     int code;
     struct sigcontext *scp;
     char *addr;
#else
SigType SigtermHandler(sig)
#endif
     int sig;
{
  SR_caught_sigint = 1;
  Error("SigtermHandler: signal was caught",(long) sig);
}

void TrapSigterm()
/*
  parallel needs to trap the SIGTERM for batch jobs
*/
{
  if ( signal(SIGTERM, SigtermHandler) == SIG_ERR)
       Error("TrapSigterm: error from signal setting SIGTERM", (long) SIGTERM);
}


