#include <libplugin/plugin.h>
#include <liboptions/liboptions.h>
#include <libciomr/libciomr.h>

#include "occwave.h"

using namespace boost;

namespace psi{ namespace occwave {

PsiReturnType
occwave(Options &options)
{
    // Start the timers
    tstart();
   
    boost::shared_ptr<Wavefunction> occ = boost::shared_ptr<Wavefunction>(new OCCWave(Process::environment.wavefunction(), options));
    Process::environment.set_wavefunction(occ);
    occ->compute_energy();
    
    // Shut down the timers
    tstop();

    return Success;
}
}} // End Namespaces


