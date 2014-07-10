#if HAVE_CONFIG_H
#   include "config.h"
#endif

/** $Header: /tmp/hpctools/ga/tcgmsg/ipcv5.0/evlog.c,v 1.3 2003-06-27 13:53:12 manoj Exp $ */

/** Event logging routine with key driven varargs interface */

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDARG_H
#   include <stdarg.h>
#endif
#if HAVE_STRINGS_H
#   include <strings.h>
#endif
#if HAVE_STRING_H
#   include <string.h>
#endif

extern long nodeid_();

#include "evlog.h"
#include "sndrcv.h"

#define ERROR_RETURN() do { \
    error = 1; \
        return; \
} while (0)

#define DUMPBUF() do { \
    (void) fputs(buffer, file); \
        (void) fflush(file); \
        if (ferror(file)) { \
            ERROR_RETURN; \
        } \
    bufpt = buffer; \
        left = BUFLEN; \
} while (0)

#define RECORD(A) do { \
    A; \
        nchars = strlen(bufpt); \
        bufpt += nchars; \
        left -= nchars; \
} while (0)

static double walltime();


/**
 * The format of the argument list is as follows:
 *
 * evlog([(int) key, [values, ...]], ..., EVKEY_LAST_ARG)
 *
 * Arguments are read as keys with corresponding values. Recognised keys
 * are defined in evlog.h and are described in detail below.
 *
 * Logging is enabled/disabled by calling evlog with one of EVKEY_ENABLE
 * or EVKEY_DISABLE specified. Note that EVKEY_ENABLE must be the
 * first key specified for it to be recognized and that all keys
 * in the argument list after EVKEY_DISABLE are ignored. By default
 * events are logged in the file events. This can be overridden with
 * the key EVKEY_FILENAME, which takes the filename as its value.
 *
 * The model for logging events assumed by the post-analysis routines
 * assumes that upon logging an event:
 *
 *   a) no state chage occurs (EVKEY_EVENT). The event is just recorded.
 *
 *   b) the process changes state by pushing the event onto the state stack
 *      (EVKEY_BEGIN).
 *
 *   c) the process changes state by popping an event off the state stack
 *      (EVKEY_END). If the event or state popped off the stack does not
 *      match the specified event then the post-analysis may get confused
 *      but this does not interfere with the actual logging.
 *
 * EVKEY_EVENT, EVKEY_BEGIN or EVKEY_END must be the first key specified other
 * than a possible EVKEY_ENABLE.
 *
 * Internally an event is stored as a large character string to simplify
 * post-analysis. Users specify data for storage in addition to
 * that which is automatically stored (only the time and process) with
 * key, value combinations (EVKEY_STR_INT, EVKEY_STR_DBL, EVKEY_STR).
 * Many such key-value combinations as required may be specified.
 * Since the internal data format uses colons ':', double quotation
 * marks '"' and carriage returns users should avoid these in their
 * string data.
 *
 * ----------------------------
 * Sample calling sequence:
 *
 * evlog(EVKEY_ENABLE, EVKEY_FILENAME, "events.log", EVKEY_LAST_ARG);
 *
 * evlog(EVKEY_EVENT, "Finished startup code",
 *       EVKEY_STR, "Now do some real work",
 *       EVKEY_LAST_ARG);
 *
 * evlog(EVKEY_BEGIN, "Get Matrix", EVKEY_LAST_ARG);
 *
 * evlog(EVKEY_END, "Get matrix",
 *       EVKEY_STR_INT, "Size of matrix", (int) N,
 *       EVKEY_STR_DBL, "Norm of matrix", (double) matrix_norm,
 *       EVKEY_LAST_ARG);
 *
 * evlog(EVKEY_BEGIN, "Transform matrix",
 *       EVKEY_STR_DBL, "Recomputed norm", (double) matrix_norm,
 *       EVKEY_LAST_ARG);
 *
 * evlog(EVKEY_END, "Transform matrix",
 *       EVKEY_STR_INT, "No. of iterations", (int) niters,
 *       EVKEY_LAST_ARG);
 *
 * evlog(EVKEY_DUMP, EVKEY_DISABLE, EVKEY_LAST_ARG);
 *
 * evlog(EVKEY_EVENT, "Logging is disabled ... this should not print",
 *       EVKEY_DUMP, EVKEY_LAST_ARG);
 *
 *  ----------------------------
 *
 * EVKEY_LAST_ARG
*      Terminates list ... takes no value and must be present
*
* EVKEY_EVENT, (char *) event
*      Simply log occurence of the event
*
* EVKEY_BEGIN, (char *) event
*      Push event onto process state stack
*
* EVKEY_END, (char *) event
*      Pop event off process state stack
*
* EVKEY_MSG_LEN, (int) length
*      Value is (int) mesage length SND/RCV only
*
* EVKEY_MSG_TO, (int) to
*      Value is (int) to process id SND/RCV only
*
* EVKEY_MSG_FROM, (int) from
*      Value is (int) from process  SND/RCV only
*
* EVKEY_MSG_TYPE, (int) type
*      Value is (int) message type  SND/RCV only
*
* EVKEY_STR_INT, (char *) string, (int) data
*      User data value pair
*
    * EVKEY_STR_DBL, (char *) string, (double) data
*      User data value pair (char *), (double)
    *
    * EVKEY_STR, (char *) string
*      User data value (char *)
    *
    * EVKEY_ENABLE
    *      Enable logging
    *
    * EVKEY_DISABLE
    *      Disable logging
    *
    * EVKEY_DUMP
    *      Dump out the current buffer to disk
    *
    * EVKEY_FILE, (char *) filename
    *      Use specified file to capture events. Default is "events".
    */
void evlog(int farg_key, ...)
{
    static int   logging=0;  /* Boolean flag for login enabled/disabled */
    static int   error=0;    /* Boolean flag for error detected         */
    static int   ncall=0;    /* Need to do stuff on first entry         */
    static char *buffer;     /* Logging buffer ... null terminated      */
    static char *bufpt;      /* Pointer to next free char in buffer     */
    static int   left;       /* Amount of free space in buffer          */
#define  BUFLEN 262144     /* Size allocated for buffer ... biggish */
#define  MAX_EV_LEN 1000   /* Assumed maximum size of single event record */
    static FILE *file;       /* File where events will be dumped */
    static char *filename = "events";   /* Default name of events file */

    va_list  ap;             /* For variable argument list */
    int      key;            /* Hold key being processed   */
    int      nchars;         /* No. of chars printed by sprintf call */
    char    *temp;           /* Temporary copy of bufpt */
    char    *string;         /* Temporary */
    int      integer;        /* Temporary */
    double   dbl;            /* Temporary */
    int      valid;          /* Temporary */

    /* If an error was detected on a previous call don't even try to
       do anything */

    if (error) {
        ERROR_RETURN();
    }

    /* First time in need to allocate the buffer, open the file etc */

    if (ncall == 0) {
        ncall = 1;
        if (!(bufpt = buffer = malloc((unsigned) BUFLEN))) {
            ERROR_RETURN();
        }
        left = BUFLEN;

        if (!(file = fopen(filename, "w"))) {
            ERROR_RETURN();
        }
    }

    /* Parse the arguments */

    temp = bufpt; /* Save to check if anything has been logged */
    valid = 0;    /* One of BEGIN, END or EVENT must preceed most keys */

    va_start(ap, farg_key);
    key = farg_key;
    while (key != EVKEY_LAST_ARG) {

        if ( (!logging) && (key != EVKEY_ENABLE) )
            return;

        switch (key) {

            case EVKEY_ENABLE:
                logging = 1;
                break;

            case EVKEY_DISABLE:
                logging = 0;
                goto done;
                /*      break; */

            case EVKEY_FILENAME:
                if (!(filename = strdup(va_arg(ap, char *))))
                {ERROR_RETURN();}
                if (!(file = freopen(filename, "w", file))) {ERROR_RETURN();}
                break;

            case EVKEY_BEGIN:
                valid = 1;
                RECORD(sprintf(bufpt, ":BEGIN:%s", va_arg(ap, char *)));
                RECORD(sprintf(bufpt, ":TIME:%.2f", walltime()));
                break;

            case EVKEY_END:
                valid = 1;
                RECORD(sprintf(bufpt, ":END:%s", va_arg(ap, char *)));
                RECORD(sprintf(bufpt, ":TIME:%.2f", walltime()));
                break;

            case EVKEY_EVENT:
                valid = 1;
                RECORD(sprintf(bufpt, ":EVENT:%s", va_arg(ap, char *)));
                RECORD(sprintf(bufpt, ":TIME:%.2f", walltime()));
                break;

            case EVKEY_MSG_LEN:
                if (!valid) {ERROR_RETURN();}
                RECORD(sprintf(bufpt, ":MSG_LEN:%d", va_arg(ap, int)));
                break;

            case EVKEY_MSG_TO:
                if (!valid) {ERROR_RETURN();}
                RECORD(sprintf(bufpt, ":MSG_TO:%d", va_arg(ap, int)));
                break;

            case EVKEY_MSG_FROM:
                if (!valid) {ERROR_RETURN();}
                RECORD(sprintf(bufpt, ":MSG_FROM:%d", va_arg(ap, int)));
                break;

            case EVKEY_MSG_TYPE:
                if (!valid) {ERROR_RETURN();}
                RECORD(sprintf(bufpt, ":MSG_TYPE:%d", va_arg(ap, int)));
                break;

            case EVKEY_MSG_SYNC:
                if (!valid) {ERROR_RETURN();}
                RECORD(sprintf(bufpt, ":MSG_SYNC:%d", va_arg(ap, int)));
                break;

            case EVKEY_STR_INT:
                if (!valid) {ERROR_RETURN();}
                string  = va_arg(ap, char *);
                integer =    va_arg(ap, int);
                RECORD(sprintf(bufpt, ":STR_INT:%s:%d", string, integer));
                break;

            case EVKEY_STR_DBL:
                if (!valid) {ERROR_RETURN();}
                string  = va_arg(ap, char *);
                dbl = va_arg(ap, double);
                RECORD(sprintf(bufpt, ":STR_DBL:%s:%g", string, dbl));
                break;

            case EVKEY_STR:
                if (!valid) {ERROR_RETURN();}
                RECORD(sprintf(bufpt, ":STR:%s", va_arg(ap, char *)));
                break;

            case EVKEY_DUMP:
                {DUMPBUF();}
                if (temp != bufpt) {
                    RECORD(sprintf(bufpt, "\n"));
                    temp = bufpt;
                }
                break;

            default:
                {DUMPBUF();}
                {ERROR_RETURN();}
        }
        key = va_arg(ap, int);
    }

done:
    va_end(ap);

    /* Put a linefeed on the end of the record if something is written */

    if (temp != bufpt) {
        RECORD(sprintf(bufpt, "\n"));
        temp = bufpt;
    }

    /* Should really check on every access to the buffer that there is
       enough space ... however just assume a very large maximum size
       for a single event log entry and check here */

    if (left <= 0) {
        ERROR_RETURN();
    }

    if (left < MAX_EV_LEN) {
        DUMPBUF();
    }
}


/**
 * return the wall time in seconds as a double
 */
static double walltime()
{
    return ((double) MTIME_()) * 0.01;
}

/*
int main(int argc, char **argv)
{
   int N = 19;
   double matrix_norm = 99.1;
   int niters = 5;

   evlog(EVKEY_ENABLE, EVKEY_FILENAME, "events.log", EVKEY_LAST_ARG);

   evlog(EVKEY_EVENT, "Finished startup code",
   EVKEY_STR, "Now do some real work",
   EVKEY_LAST_ARG);

   evlog(EVKEY_BEGIN, "Get Matrix", EVKEY_LAST_ARG);

   evlog(EVKEY_END, "Get matrix",
   EVKEY_STR_INT, "Size of matrix", (int) N,
   EVKEY_STR_DBL, "Norm of matrix", (double) matrix_norm,
   EVKEY_LAST_ARG);

   evlog(EVKEY_BEGIN, "Transform matrix",
   EVKEY_STR_DBL, "Recomputed norm", (double) matrix_norm,
   EVKEY_LAST_ARG);

   evlog(EVKEY_END, "Transform matrix",
   EVKEY_STR_INT, "No. of iterations", (int) niters,
   EVKEY_LAST_ARG);

   evlog(EVKEY_DUMP, EVKEY_LAST_ARG);

   evlog(EVKEY_EVENT, "Logging is disabled ... this should not print",
   EVKEY_DUMP, EVKEY_LAST_ARG);

   return 0;
}
*/
