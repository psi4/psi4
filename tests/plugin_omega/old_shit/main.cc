#include <libplugin/plugin.h>
#include "psi4-dec.h"
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.hpp>
#include <libqt/qt.h>
#include "omega.h"

INIT_PLUGIN

using namespace psi::scf;

namespace psi{ namespace omega_wrapper{

extern "C" int
read_options(std::string name, Options &options)
{
    if (name == "PLUGIN_OMEGA"|| options.read_globals()) {
        /*- The amount of information printed
            to the output file -*/
        options.add_int("PRINT", 1);
        /*- The amount of information printed
            to the output file -*/
        options.add_int("DEBUG", 1);
    }
    return true;
}

extern "C" PsiReturnType
plugin_omega(Options &options)
{
    tstart();

    // Initialize the psi3 timer library.
    timer_init();

    boost::shared_ptr<PSIO> psio = PSIO::shared_object();
    boost::shared_ptr<Molecule> molecule = Process::environment.molecule();

    // Build the OmegaKS object
    boost::shared_ptr<OmegaKS> scf(new OmegaKS(options, psio));

    // Run the OmegaKS procedure
    double energy = scf->compute_energy();

    // Set some environment variables
    Process::environment.globals["SCF ENERGY"] = energy;
    Process::environment.globals["CURRENT ENERGY"] = energy;

    // Shut down psi.
    timer_done();

    tstop();

    return Success;
}
}} // End Namespaces
