#include <libplugin/plugin.h>
#include <liboptions/liboptions.h>
#include <libciomr/libciomr.h>

#include "ocepawave.h"

using namespace boost;

namespace psi{ namespace ocepawave {

PsiReturnType
ocepawave(Options &options)
{
    // Start the timers
    tstart();
   
    boost::shared_ptr<Wavefunction> ocepa = boost::shared_ptr<Wavefunction>(new OCEPAWave(Process::environment.wavefunction(), options));
    Process::environment.set_wavefunction(ocepa);
    ocepa->compute_energy();
    
    // Shut down the timers
    tstop();

    return Success;
}

}} // End Namespaces


