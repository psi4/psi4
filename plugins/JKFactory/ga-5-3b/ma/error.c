#if HAVE_CONFIG_H
#   include "config.h"
#endif

/** @file
 * Error handling module.
 */

#if HAVE_STDIO_H
#   include <stdio.h>
#endif
#if HAVE_STDLIB_H
#   include <stdlib.h>
#endif

#include "error.h"
#include "scope.h"

/** default # of initial table entries */
#define MA_EBUF_SIZE 1024

/** buffer for error messages */
public char ma_ebuf[MA_EBUF_SIZE];

/** print error messages for nonfatal errors? */
public Boolean ma_error_print = MA_TRUE;

/** terminate execution upon any error? */
public Boolean ma_hard_fail = MA_FALSE;

void (*ma_func_terminate)() = 0;


void MA_set_error_callback(void (*func)())
{
    ma_func_terminate = func;
}


/**
 * Depending on the given arguments and certain global parameters,
 * possibly print a message to stderr and/or terminate the program.
 *
 * @param elevel severity of error
 * @param etype category of error
 * @param func name of routine in which error was found
 * @param emsg msg describing error
 */
public void ma_error(ErrorLevel elevel, ErrorType etype, char *func, char *emsg)
{
    /* print a message? */
    if ((elevel == EL_Fatal) || ma_hard_fail || ma_error_print)
    {
        char *s1; /* internal or not */
        char *s2; /* class of error */

        /* set s1 */
        if (etype == ET_Internal) {
            s1 = "internal ";
        } else {
            s1 = "";
        }

        /* set s2 */
        if (elevel == EL_Fatal) {
            s2 = "fatal error";
        } else if (ma_hard_fail) {
            s2 = "hard failure";
        } else {
            s2 = "error";
        }

        /* print the message */
        (void)fflush(stdout);
        (void)fflush(stderr);
        (void)fprintf(stderr, "MA %s%s: %s: %s\n", s1, s2, func, emsg);
        (void)fflush(stderr);
    }

    /* terminate execution? */
    if ((elevel == EL_Fatal) || ma_hard_fail) {
        if(ma_func_terminate) {
            ma_func_terminate("MA aborting",0);
        } else {
            exit(1);
        }
    }
}
