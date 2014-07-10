#if HAVE_CONFIG_H
#   include "config.h"
#endif

/* $Header: /tmp/hpctools/ga/tcgmsg/ipcv4.0/evon.c,v 1.4 1995-02-24 02:17:17 d3h325 Exp $ */

/* Crude FORTRAN interface to C event logging routines.
   See evlog.c for more details.

   FORTRAN character variables are so unportable that guaranteeing
   that U can parse a variable length argument list is next to impossible.

   This provides very basic event logging functionality.

   CALL EVON()

      enable logging.

   CALL EVOFF()
 
      disable logging.

   CALL EVBGIN("event description")

      push event onto state stack

   CALL EVEND("event description")

      pop event off state stack

   CALL EVENT("event description")

      log occurence of event that doesn't change state stack
*/

#include <stdlib.h>

#ifdef IPSC
#define bcopy(a, b, n) memcpy((b), (a), (n))
#endif

#if 0
#if defined(ULTRIX) || defined(SGI) || defined(NEXT) || defined(HPUX) || \
    defined(KSR)    || defined(DECOSF)
extern void *malloc();
#else
extern char *malloc();
#endif
#endif

#include "evlog.h"

/* These to get portable FORTRAN interface ... these routines
   will not be called from C which has the superior evlog interface */

#if (defined(AIX) || defined(NEXT) || defined(HPUX)) && !defined(EXTNAME)
#define evon_     evon
#define evoff_    evoff
#define evbgin_   evbgin
#define evend_    evend
#define event_    event
#endif

#if (defined(CRAY) || defined(ARDENT))
#define evon_     EVON
#define evoff_    EVOFF
#define evbgin_   EVBGIN
#define evend_    EVEND
#define event_    EVENT
#endif

/* Define crap for handling FORTRAN character arguments */

#ifdef CRAY
#include <fortran.h>
#endif
#ifdef ARDENT
struct char_desc {
  char *string;
  int len;
};
#endif

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

#ifdef ARDENT
void evbgin_(arg)
     struct char_desc *arg;
{
  char *string = arg->string;
  int   len = arg->len;
#endif
#ifdef CRAY
void evbgin_(arg)
     _fcd arg;
{
  char *string = _fcdtocp(arg);
  int len = _fcdlen(arg);
#endif
#if !defined(ARDENT) && !defined(CRAY)
void evbgin_(string, len)
  char *string;
  int   len;
{
#endif
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

#ifdef ARDENT
void evend_(arg)
     struct char_desc *arg;
{
  char *string = arg->string;
  int   len = arg->len;
#endif
#ifdef CRAY
void evend_(arg)
     _fcd arg;
{
  char *string = _fcdtocp(arg);
  int len = _fcdlen(arg);
#endif
#if !defined(CRAY) && !defined(ARDENT)
void evend_(string, len)
  char *string;
  int   len;
{
#endif
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
  
#ifdef ARDENT
void event_(arg)
     struct char_desc *arg;
{
  char *string = arg->string;
  int   len = arg->len;
#endif
#ifdef CRAY
void event_(arg)
     _fcd arg;
{
  char *string = _fcdtocp(arg);
  int len = _fcdlen(arg);
#endif
#if !defined(ARDENT) && !defined(CRAY)
void event_(string, len)
  char *string;
  int   len;
{
#endif
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
