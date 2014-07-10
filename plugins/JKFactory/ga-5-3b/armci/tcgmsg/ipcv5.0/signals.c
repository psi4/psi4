#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Id: signals.c,v 1.3 2000-11-14 20:43:56 d3h325 Exp $ */

#if HAVE_SIGNAL_H
#   include <signal.h>
#endif
#if HAVE_SYS_WAIT_H
#   include <sys/wait.h>
#endif
#if HAVE_SYS_TYPES_H
#   include <sys/types.h>
#endif

#include "tcgmsgP.h"

#ifndef SIG_ERR
#   define SIG_ERR (RETSIGTYPE (*)(int))-1
#endif


RETSIGTYPE SigintHandler(int sig)
{
    TCGMSG_caught_sigint = 1L;
    Error("SigintHandler: signal was caught",0L);
}


/**
 * Trap the signal SIGINT so that we can propagate error
 * conditions and also tidy up shared system resources in a
 * manner not possible just by killing everyone
 */
void TrapSigint()
{
    if ( signal(SIGINT, SigintHandler) == SIG_ERR)
        Error("TrapSigint: error from signal setting SIGINT",(long) SIGINT);
}


/**
 * kill -SIGINT all of my beloved children
 */
void ZapChildren()
{
    long node;

    for (node=0; node<TCGMSG_nnodes; node++)
        if (node != TCGMSG_nodeid)
            (void) kill((int) TCGMSG_proc_info[node].pid, SIGINT);
}


RETSIGTYPE SigchldHandler(int sig)
{
    int status;

    (void) wait(&status);
    TCGMSG_caught_sigint = 1;
    Error("Child process terminated prematurely, status=",(long) status);
}


/**
 * Trap SIGCHLD so that can tell if children die unexpectedly.
 */
void TrapSigchld()
{
    if ( signal(SIGCHLD, SigchldHandler) == SIG_ERR)
        Error("TrapSigchld: error from signal setting SIGCHLD", (long) SIGCHLD);
}
