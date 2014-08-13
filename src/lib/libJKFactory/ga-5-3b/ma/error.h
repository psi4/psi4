/** @file
 * Private header file containing symbolic constants and type declarations
 * for the error handling module.
 *
 * This file should only be included by internal C files.
 */
#ifndef _error_h
#define _error_h

#include "macdecls.h"

/**
 ** types
 **/

/* severity levels of errors */
typedef enum
{
    EL_Fatal,        /* unrecoverable */
    EL_Nonfatal        /* recoverable */
} ErrorLevel;

/* categories of errors */
typedef enum
{
    ET_External,    /* anticipated, client-level problem */
    ET_Internal        /* unanticipated problem internal to MA */
} ErrorType;

/**
 ** function types
 **/

extern void ma_error(ErrorLevel elevel, ErrorType etype, char *func, char *emsg);

/**
 ** variables
 **/

/* buffer for error messages */
extern char ma_ebuf[];

/* print error messages for nonfatal errors? */
extern Boolean ma_error_print;

/* terminate execution upon any error? */
extern Boolean ma_hard_fail;

#endif /* _error_h */
