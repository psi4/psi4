#include <stdio.h>
#include <stdlib.h>
#include "psi4-dec.h"
#include "libqt/qt.h"
#include "libchkpt/chkpt.hpp"
#include "libpsio/psio.h"
#include "libchkpt/chkpt.h"
#include "defines.h"
#include "dcft.h"

using namespace psi;


namespace psi{ namespace dcft{

PsiReturnType
dcft(Options &options, int argc, char *argv[])
{
    // Start the timers
    tstart();
    timer_init();

    psio_init();
    psio_ipv1_config();
    chkpt_init(PSIO_OPEN_OLD);

    fprintf(outfile,"\n\n\t\t*******************************************\n");
    fprintf(outfile,    "\t\t*  DCFT - A Density Cumulant Functional   *\n");
    fprintf(outfile,    "\t\t*  Theory code written by Andy Simmonett  *\n");
    fprintf(outfile,    "\t\t*******************************************\n\n");

    // The solver object doesn't really need to be an object at this point
    // but one fine day, main might do something a bit more elaborate...
    DCFTSolver dcft(options);
    dcft.compute();

    // Shut down the timers
    timer_done();
    tstop();

    return Success;   
}

}} // End Namespaces
