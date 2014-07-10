#if HAVE_CONFIG_H
#   include "config.h"
#endif

/** @file
 * Primitive version of strtok for alliant etc who don't have it.
 * I think it works .... ?
 */

#undef NULL
#define NULL 0

/**
 * Return 1 if in set, 0 otherwise.
 */
static int InSet(char *a, char *set)
{
    register char test;
    register char b = (*a);

    while ( (test = *set++) ) {
        if (test == b) {
            return 1;
        }
    }

    return 0;
}

/*
 * Return pointer to next character in string not in set or
 * return NULL pointer if no such character.
 */
static char *NextNotInSet(char *string, char *set)
{
    /* Return NULL if given NULL */
    if (string == (char *) NULL) {
        return (char *) NULL;
    }

    while (*string) {
        if (InSet(string, set)) {
            string++;
        } else {
            break;
        }
    }

    if (*string) {
        return string;
    } else {
        return (char *) NULL;
    }
}

/**
 * Return pointer to next character in string in set or
 * return NULL pointer if no such character.
 */
static char *NextInSet(char *string, char *set)
{
    /* Return NULL if given NULL */
    if (string == (char *) NULL) {
        return (char *) NULL;
    }

    while (*string) {
        if (InSet(string, set)) {
            break;
        } else {
            string++;
        }
    }

    if (*string) {
        return string;
    } else {
        return (char *) NULL;  
    }
}

/**
 * Naive version of strtok for the alliant.
 */
char *strtok(char *s1, char *s2)
{
    static char *string = (char *) NULL;
    char *start, *end;

    /* Initialize on first call */
    if (s1 != (char *) NULL) {
        string = s1;
    }

    start = NextNotInSet(string, s2); /* Find start of next token */

    end = NextInSet(start, s2);       /* Find end of this token */

    if (end == (char *) NULL) {
        string = (char *) NULL;
    } else {
        string = end + 1;
        *end = '\0';
    }

    return start;
}
