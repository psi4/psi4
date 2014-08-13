#if HAVE_CONFIG_H
#   include "config.h"
#endif

#include <stdlib.h>
extern char *strcpy();
extern size_t strlen();

char *strdup(char *s)
{
    char *new;

    if ((new = malloc((size_t) (strlen(s)+1))))
        (void) strcpy(new,s);

    return new;
}
