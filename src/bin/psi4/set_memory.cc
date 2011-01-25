/*!
** \file
** \brief Get the amount of core memory available from input
** \ingroup CIOMR
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libparallel/parallel.h>
#include <psi4-dec.h>

namespace psi {

#define DEF_MAXCRR (256000000)  // default maxcor 256 M bytes

void set_memory(FILE *outfile)
{
    long int maxcrr;

    maxcrr = DEF_MAXCRR; // set to default

    if (maxcrr < 1e9) {
        if(Communicator::world->me() == 0) fprintf(outfile,"    Memory level set to %.3lf MB\n", maxcrr / 1e6 );
    }
    else {
        if(Communicator::world->me() == 0) fprintf(outfile,"    Memory level set to %.3lf GB\n", maxcrr / 1e9 );
    }

    Process::environment.set_memory(maxcrr);

    return;
}

}
