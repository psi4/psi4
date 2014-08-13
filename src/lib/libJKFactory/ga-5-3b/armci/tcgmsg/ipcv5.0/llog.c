#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_SYS_TYPES_H
#   include <sys/types.h>
#endif
#if HAVE_TIME_H
#   include <time.h>
#endif
#if HAVE_SYS_TIME_H
#   include <sys/time.h>
#endif

#include "sndrcv.h"

extern void Error();

/**
 * close and open stdin and stdout to append to a local logfile
 * with the name log.<process#> in the current directory
*/
void LLOG_()
{
    char name[12];
    time_t t;

    (void) sprintf(name, "log.%03ld",(long)NODEID_());

    (void) fflush(stdout);
    (void) fflush(stderr);

    if (freopen(name, "a", stdout) == (FILE *) NULL) {
        Error("LLOG_: error re-opening stdout", (long) -1);
    }

    if (freopen(name, "a", stderr) == (FILE *) NULL) {
        Error("LLOG_: error re-opening stderr", (long) -1);
    }

    (void) time(&t);
    (void) printf("\n\nLog file opened : %s\n\n",ctime(&t));
    (void) fflush(stdout);
}
