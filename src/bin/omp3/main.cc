#include <libplugin/plugin.h>
#include <liboptions/liboptions.h>
#include <libciomr/libciomr.h>

#include "omp3wave.h"

using namespace boost;

namespace psi{ namespace omp3wave {

PsiReturnType
omp3wave(Options &options)
{
    // Start the timers
    tstart();
   
    boost::shared_ptr<Wavefunction> omp3 = boost::shared_ptr<Wavefunction>(new OMP3Wave(Process::environment.reference_wavefunction(), options));
    Process::environment.set_reference_wavefunction(omp3);
    omp3->compute_energy();
    
    
    // Shut down the timers
    tstop();

    return Success;
}

}} // End Namespaces


