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
        options.add_int("DEBUG", 0);
        /*- The Omega optimization procedure
            -*/
        options.add_str("OMEGA_PROCEDURE", "IP", "IP");
        /*- The Omega optimization root-finding algorithm
            -*/
        options.add_str("OMEGA_ROOT_ALGORITHM", "REGULA_FALSI3", "BISECTION REGULA_FALSI REGULA_FALSI2 REGULA_FALSI3 BRENT");
        /*- Maximum number of omega iterations to perform
            -*/
        options.add_int("OMEGA_MAXITER", 30);
        /*- Procedure used to guess initial omega
            -*/
        options.add_str("OMEGA_GUESS", "HOMO_SIZE" , "HOMO_SIZE DEFAULT");
        /*- a in w_0^-1 = a <R>_HOMO + b
            -*/
        options.add_double("OMEGA_GUESS_A", 2.0);
        /*- b in w_0^-1 = a <R>_HOMO + b
            -*/
        options.add_double("OMEGA_GUESS_B", 0.0);
        /*- Multiplier to use to backet omega (>1)
            -*/
        options.add_double("OMEGA_BRACKET_ALPHA", 2.0);
        /*- Convergence threshold for omega, 10^-thresh
            -*/
        options.add_int("OMEGA_CONVERGE", 3);
        /*- Interpolate Fock matrices at omega steps?
            -*/
        options.add_bool("OMEGA_GUESS_INTERPOLATE", true);
        /*- The Omega functional selected for this procedure
            -*/
        options.add_str("DFT_FUNCTIONAL", "");

    }
    return true;
}

extern "C" PsiReturnType
plugin_omega(Options &options)
{
    // Initialize the psi3 timer library.
//    timer_init();

    boost::shared_ptr<PSIO> psio = PSIO::shared_object();
    boost::shared_ptr<Molecule> molecule = Process::environment.molecule();

    // Build the OmegaKS object
    boost::shared_ptr<OmegaKS> scf;

    if (options.get_str("OMEGA_PROCEDURE") == "IP") {
        scf = boost::shared_ptr<OmegaKS>(new OmegaIPKS(options, psio));
    }

    // Run the OmegaKS procedure
    scf->run_procedure();

    // Shut down psi.
//    timer_done();

    return Success;
}

}} // End Namespaces
