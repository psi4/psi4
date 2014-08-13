#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Header: /tmp/hpctools/ga/tcgmsg/ipcv5.0/mtime.c,v 1.5 2002-03-12 18:59:31 d3h325 Exp $ */

#if HAVE_STDIO_H
#   include <stdio.h>
#endif

#include "srftoc.h"


/**
 * return wall clock time in centiseconds
 */
long MTIME_()
{
    double TCGTIME_();
    return (long) (TCGTIME_()*100.0);
}

#if defined(LAPI) && defined(AIX)
#   define LAPI_AIX
#endif

#ifndef LAPI_AIX

#if HAVE_SYS_TYPES_H
#   include <sys/types.h>
#endif
#if HAVE_SYS_TIME_H
#   include <sys/time.h>
#endif

static unsigned firstsec=0;     /* Reference for timer */
static unsigned firstusec=0;    /* Reference for timer */

void MtimeReset()               /* Sets timer reference */
{
    struct timeval tp;
    struct timezone tzp;

    (void) gettimeofday(&tp,&tzp);

    firstsec = tp.tv_sec;
    firstusec = tp.tv_usec;
}


/**
 * Return wall clock time in seconds as accurately as possible
 */
double TCGTIME_()
{
    static int firstcall=1;
    double low, high;

    struct timeval tp;
    struct timezone tzp;

    if (firstcall) {
        MtimeReset();
        firstcall = 0;
    }

    (void) gettimeofday(&tp,&tzp);

    low = (double) (tp.tv_usec>>1) - (double) (firstusec>>1);
    high = (double) (tp.tv_sec - firstsec);

    return high + 1.0e-6*(low+low);
}

#endif

#ifdef LAPI_AIX

#if HAVE_SYS_TIME_H
#   include <sys/time.h>
#endif
#if HAVE_SYS_SYSTEMCFG_H
#   include <sys/systemcfg.h>
#endif

static int firstsec=0;    /* Reference for timer */
static int firstnsec=0;    

void MtimeReset()               /* Sets timer reference */
{
    timebasestruct_t t;
    read_real_time(&t, TIMEBASE_SZ);
    time_base_to_time(&t, TIMEBASE_SZ);

    firstsec = t.tb_high;
    firstnsec = t.tb_low;
}


/**
 * Return wall clock time in seconds as accurately as possible
 */
double TCGTIME_()
{
    static int firstcall=1;
    timebasestruct_t t;
    int low, high;
    int secs, nsecs;

    if (firstcall) {
        MtimeReset();
        firstcall = 0;
    }


    read_real_time(&t, TIMEBASE_SZ);
    time_base_to_time(&t, TIMEBASE_SZ);

    secs = t.tb_high - firstsec;
    nsecs = t.tb_low - firstnsec;

    /* If there was a carry from low-order to high-order during
       the measurement, we have to undo it */
    if(nsecs < 0){
        secs--;
        nsecs+= 1000000000;
    }
    return (double)(secs + 1.0e-9*nsecs);
}

#endif /* LAPI_AIX */
