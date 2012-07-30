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
using namespace boost;

namespace psi{ namespace dcft{

PsiReturnType
dcft(Options &options)
{

    // Start the timers
    tstart();

    fprintf(outfile,"\n\n\t\t*******************************************\n");
    fprintf(outfile,    "\t\t*  DCFT - A Density Cumulant Functional   *\n");
    fprintf(outfile,    "\t\t*  Theory code written by Andy Simmonett  *\n");
    fprintf(outfile,    "\t\t*******************************************\n\n");

    boost::shared_ptr<Wavefunction> dcft = boost::shared_ptr<Wavefunction>(new DCFTSolver(Process::environment.wavefunction(), options));
    Process::environment.set_wavefunction(dcft);

    dcft->compute_energy();

    // Shut down the timers
    tstop();

    return Success;
}

}} // End Namespaces
