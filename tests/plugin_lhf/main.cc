#include <libplugin/plugin.h>
#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.hpp>
#include "lhf.h"

INIT_PLUGIN

namespace psi{ namespace scf {

extern "C" 
int read_options(std::string name, Options& options)
{
    if (name == "LHF"|| options.read_globals()) {

    }

    return true;
}


extern "C" 
PsiReturnType plugin_lhf(Options& options)
{
    boost::shared_ptr<Wavefunction> hf(new LHF(options));
    hf->compute_energy();
    return Success;
}


}} // End namespaces

