#include <libplugin/plugin.h>
#include "psi4-dec.h"
#include <libparallel/parallel.h>
#include <liboptions/liboptions.h>
#include <libmints/mints.h>
#include <libpsio/psio.h>
#include <libchkpt/chkpt.hpp>
#include <libqt/qt.h>
#include "mad_mp2.h"

INIT_PLUGIN

namespace psi{ 

extern "C" int
read_options(std::string name, Options &options)
{
    if (name == "MAD_MP2"|| options.read_globals()) {
        /*- The amount of information printed
            to the output file -*/
        options.add_int("PRINT", 1);
        /*- The amount of information printed
            to the output file -*/
        options.add_int("MADMP2_DEBUG", 1);
        /*- How long to sleep before running -*/
        options.add_int("MADMP2_SLEEP", 0);
        /*- The schwarz cutoff -*/
        options.add_double("SCHWARZ_CUTOFF", 0.0);
        /*- The same-spin scale factor -*/
        options.add_double("SCALE_SS", 1.0/3.0);
        /*- The opposite-spin scale factor -*/
        options.add_double("SCALE_OS", 6.0/5.0);
        /*- Denominator algorithm for PT methods -*/
        options.add_str("DENOMINATOR_ALGORITHM", "LAPLACE", "LAPLACE CHOLESKY");
        /*- Maximum denominator error allowed (Max error norm in Delta tensor) -*/
        options.add_double("DENOMINATOR_DELTA", 1.0E-6);
        /*- MP2 Algorithm -*/
        options.add_str("MP2_ALGORITHM", "DFMP2", "DFMP2 DFMP2J");
    }
    return true;
}

extern "C" PsiReturnType
plugin_mad_mp2(Options &options)
{
    tstart();

    boost::shared_ptr<PSIO> psio = PSIO::shared_object();

    boost::shared_ptr<psi::mad_mp2::MAD_MP2> mp2(new psi::mad_mp2::MAD_MP2(options, psio));
   
    mp2->compute_energy();

    tstop();

    return Success;
}

} // End Namespaces

