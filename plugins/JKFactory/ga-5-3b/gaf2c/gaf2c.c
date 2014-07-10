#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STRING_H
#   include <string.h>
#endif
#if HAVE_STDIO_H
#   include <stdio.h>
#endif

#include "farg.h"
#include "message.h"


/**
 * Converts C strings to Fortran strings.
 *
 * @param cstring C buffer
 * @param fstring Fortran string
 * @param flength length of fstring
 */
void ga_c2fstring(char *cstring, char *fstring, int flength)
{
    int clength = strlen(cstring);

    /* remove terminal \n character if any */
    if(cstring[clength] == '\n') {
        clength--;
    }

    /* Truncate C string into Fortran string */
    if (clength > flength) {
        clength = (int)flength;
    }

    /* Copy characters over */
    flength -= clength;
    while (clength--) {
        *fstring++ = *cstring++;
    }

    /* Now terminate with blanks */
    while (flength--) {
        *fstring++ = ' ';
    }
}


/**
 * Converts Fortran strings to C strings.
 *
 * Strip trailing blanks from fstring and copy it to cstring,
 * truncating if necessary to fit in cstring, and ensuring that
 * cstring is NUL-terminated.
 *
 * @param fstring Fortran string
 * @param flength length of fstring
 * @param cstring C buffer
 * @param clength max length (including NUL) of cstring
 */
void ga_f2cstring(char *fstring, int flength, char *cstring, int clength)
{
    /* remove trailing blanks from fstring */
    while (flength-- && fstring[flength] == ' ') ;

    /* the postdecrement above went one too far */
    flength++;

    /* truncate fstring to cstring size */
    if (flength >= clength) {
        flength = clength - 1;
    }

    /* ensure that cstring is NUL-terminated */
    cstring[flength] = '\0';

    /* copy fstring to cstring */
    while (flength--) {
        cstring[flength] = fstring[flength];
    }
}

void ga_f2c_get_cmd_args(int *argc, char ***argv)
{
    Integer i=0;
    int iargc=F2C_IARGC();
    char **iargv=NULL;

    if (iargc > F2C_GETARG_ARGV_MAX) {
        printf("ga_f2c_get_cmd_args: too many cmd line args");
        armci_msg_abort(1);
    }
    iargv = malloc(sizeof(char*)*F2C_GETARG_ARGV_MAX);
    if (!iargv) {
        printf("ga_f2c_get_cmd_args: malloc iargv failed");
        armci_msg_abort(1);
    }
    for (i=0; i<iargc; i++) {
        char fstring[F2C_GETARG_ARGLEN_MAX];
        char cstring[F2C_GETARG_ARGLEN_MAX];
        F2C_GETARG(&i, fstring, F2C_GETARG_ARGLEN_MAX);
        ga_f2cstring(fstring, F2C_GETARG_ARGLEN_MAX,
                cstring, F2C_GETARG_ARGLEN_MAX);
        iargv[i] = strdup(cstring);
    }
    *argc = iargc;
    *argv = iargv;
}


#if NOFORT

/* To avoid missing symbols even though these should never be called. */

/**
 */
void F2C_GETARG(Integer *a, char *b, int c)
{
    printf("GA was built without support for Fortran.  You have attempted "
            "to retreive command line arguments from a Fortran program. "
            "Please recompile GA and avoid using --disable-f77");
    armci_msg_abort(1);
}


/**
 */
Integer F2C_IARGC()
{
    printf("GA was built without support for Fortran.  You have attempted "
            "to retreive command line arguments from a Fortran program. "
            "Please recompile GA and avoid using --disable-f77");
    armci_msg_abort(1);
    return 0;
}

#endif /* NOFORT */
