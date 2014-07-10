#if HAVE_CONFIG_H
#   include "config.h"
#endif

#if HAVE_STDIO_H
#   include <stdio.h>
#endif

#include "tcgmsg.h"

/**
 * Traditional first parallel program.
 */
int main(int argc, char **argv)
{
    tcg_pbegin(argc, argv);

    (void) printf("Hello from node %ld\n",tcg_nodeid());

    tcg_pend();

    return 0;
}
