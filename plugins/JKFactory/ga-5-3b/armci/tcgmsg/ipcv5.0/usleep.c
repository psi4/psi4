#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Header: /tmp/hpctools/ga/tcgmsg/ipcv5.0/usleep.c,v 1.4 1997-11-07 23:44:20 d3h325 Exp $ */

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_SYS_SELECT_H
#   include <sys/select.h>
#endif
#if HAVE_SYS_TYPES_H
#   include <sys/types.h>
#endif
#if HAVE_SYS_TIME_H
#   include <sys/time.h>
#endif

#include "tcgmsgP.h"

#ifdef STUPIDUSLEEP
void USleep(long us)
{
    int s = us/1000000;
    if (s == 0)
        s = 1;
    (void) sleep(s);
}
#else /* STUPIDUSLEEP */
/**
 * Sleep for the specified no. of micro-seconds ... uses the timeout
 * on select ... it seems to be accurate to about a few centiseconds
 * on a sun.  I don't know how much system resources it eats.
 */
void USleep(long us)
{
    int width=0;
    struct timeval timelimit;

    /*  printf("%2ld: sleeping for %ldus\n", TCGMSG_nodeid, us);
        fflush(stdout);*/

    timelimit.tv_sec = (int) (us/1000000);
    timelimit.tv_usec = (int) (us - timelimit.tv_sec*1000000);

    (void) select(width, (fd_set *) 0, (fd_set *) 0, (fd_set *) 0,
                  &timelimit);
}
#endif /* STUPIDUSLEEP */
