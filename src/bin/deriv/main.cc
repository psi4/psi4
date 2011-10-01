#include <stdio.h>
#include <stdlib.h>

#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libmints/mints.h>

#include <psi4-dec.h>

using namespace boost;

namespace psi { namespace deriv {

PsiReturnType deriv(Options &)
{
    tstart();

    fprintf(outfile, " DERIV: Derivative code.\n   by Justin Turney\n\n");

    Deriv test(Process::environment.reference_wavefunction(),
               0x1,
               false,
               false);
    test.compute();

    // Shut down psi
    tstop();

    return Success;
}

}} // namespaces
