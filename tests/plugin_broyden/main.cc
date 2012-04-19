#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include "broyden.h"

INIT_PLUGIN

namespace psi{ namespace scf {

extern "C" 
int read_options(std::string name, Options& options)
{
    if (name == "BROYDEN"|| options.read_globals()) {

        /*- Step Type -*/
        options.add_str("SCF_STEP", "BROYDEN", "BROYDEN BFGS SR1 DFP ROOTHAAN"); 
        /*- Pade Expansion Size in expm -*/
        options.add_int("BROYDEN_PADE_N", 2);
        /*- Broyden start iteration -*/
        options.add_int("BROYDEN_START", 1);
        /*- Broyden maximum iterations before restart -*/
        options.add_int("BROYDEN_MAXITER", 10);

    }

    return true;
}


extern "C" 
PsiReturnType plugin_broyden(Options& options)
{
    boost::shared_ptr<Wavefunction> hf(new BroydenRHF(options));
    hf->compute_energy();
    return Success;
}


}} // End namespaces

