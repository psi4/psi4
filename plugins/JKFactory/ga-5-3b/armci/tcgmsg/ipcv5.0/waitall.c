#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $$ */

#if HAVE_SYS_TYPES_H
#   include <sys/types.h>
#endif
#if HAVE_SYS_WAIT_H
#   include <sys/wait.h>
#endif
#if HAVE_STDIO_H
#   include <stdio.h>
#endif

#include "sndrcv.h"


/**
 * Wait for all children to finish and return appropriate status
 * 0 = OK
 * 1 = bad news
 */
int WaitAll(long nchild)
{
    int status, pid, child, stat=0, lo, hi;

    for (child=0; child<nchild; child++) {
        pid = wait(&status);
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
