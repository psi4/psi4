#if HAVE_CONFIG_H
#   include "config.h"
#endif

/** @file
 * Crude FORTRAN interface to C event logging routines.
 *  See evlog.c for more details.
 *
 *  FORTRAN character variables are so unportable that guaranteeing
 *  that U can parse a variable length argument list is next to impossible.
 *
 *  This provides very basic event logging functionality.
 *
 *  CALL EVON()
 *     enable logging.
 *
 *  CALL EVOFF()
 *     disable logging.
 *
 *  CALL EVBGIN("event description")
 *     push event onto state stack
 *
 *  CALL EVEND("event description")
 *     pop event off state stack
 *
 *  CALL EVENT("event description")
 *     log occurence of event that doesn't change state stack
 */

#include "evlog.h"

/* These to get portable FORTRAN interface ... these routines
 * will not be called from C which has the superior evlog interface */

#define evon_   F77_FUNC(evon,EVON)
#define evoff_  F77_FUNC(evoff,EVOFF)
#define evbgin_ F77_FUNC(evbgin,EVBGIN)
#define evend_  F77_FUNC(evend,EVEND)
#define event_  F77_FUNC(event,EVENT)

void evon_()
{
#ifdef EVENTLOG
    evlog(EVKEY_ENABLE, EVKEY_LAST_ARG);
#endif
}

void evoff_()
{
#ifdef EVENTLOG
    evlog(EVKEY_DISABLE, EVKEY_LAST_ARG);
#endif
}

void evbgin_(char *string, long len)
{
#ifdef EVENTLOG
    char *value = malloc( (unsigned) (len+1) );

    if (value) {
        (void) bcopy(string, value, len);
        value[len] = '\0';
        evlog(EVKEY_BEGIN, value, EVKEY_LAST_ARG);
        (void) free(value);
    }
#endif
}

void evend_(char *string, long len)
{
#ifdef EVENTLOG
    char *value = malloc( (unsigned) (len+1) );

    if (value) {
        (void) bcopy(string, value, len);
        value[len] = '\0';
        evlog(EVKEY_END, value, EVKEY_LAST_ARG);
        (void) free(value);
    }
#endif
}

void event_(char *string, long len)
{
#ifdef EVENTLOG
    char *value = malloc( (unsigned) (len+1) );

    if (value) {
        (void) bcopy(string, value, len);
        value[len] = '\0';
        evlog(EVKEY_EVENT, value, EVKEY_LAST_ARG);
        (void) free(value);
    }
#endif
}
