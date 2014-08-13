#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_ERRNO_H
#   include <errno.h>
#endif
#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_SIGNAL_H
#   include <signal.h>
#endif

#include "sndrcv.h"
#include "tcgmsgP.h"

extern void perror(const char *);
extern void exit(int);
extern void ZapChildren(void);

#define DEV stderr


void Error(char *string, long integer)
{
    (void) signal(SIGINT, SIG_IGN);
    (void) signal(SIGCHLD, SIG_DFL); /* Death of children to be expected */

    (void) fflush(stdout);
    if (TCGMSG_caught_sigint) {
        (void) fprintf(DEV,"%2ld: interrupt\n",(long)NODEID_());
    }
    else {
        (void) fprintf(DEV,"%3ld: %s %ld (%#lx).\n", (long)NODEID_(), string,
                       (long)integer,(long)integer);
        if (errno != 0)
            perror("system error message");
    }
    (void) fflush(DEV);

    /* Shut down the sockets and remove shared memory and semaphores to
       propagate an error condition to anyone that is trying to communicate
       with me */

#ifndef LAPI
    ZapChildren();  /* send interrupt to children which should trap it
                       and call Error in the handler */
    DeleteSharedRegion(TCGMSG_shmem_id);
#endif

    exit(1);
}


/**
 * Interface from fortran to c error routine
 */
void PARERR_(long *code)
{
    long lcode = (long)(*code);
    Error("User detected error in FORTRAN", lcode);
}
