#include <libplugin/plugin.h>
#include <liboptions/liboptions.h>
#include <libciomr/libciomr.h>

#include "omp2wave.h"

using namespace boost;

namespace psi{ namespace omp2wave {

PsiReturnType
omp2wave(Options &options)
{
    // Start the timers
    tstart();
   
    boost::shared_ptr<Wavefunction> omp2 = boost::shared_ptr<Wavefunction>(new OMP2Wave(Process::environment.reference_wavefunction(), options));
    Process::environment.set_reference_wavefunction(omp2);
    omp2->compute_energy();
    
    
    // Shut down the timers
    tstop();

    return Success;
}

}} // End Namespaces


