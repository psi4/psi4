#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <libmints/wavefunction.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include "plugin.h"

using namespace boost;
using namespace psi;

extern "C" void init_plugin()
{
    // Prepares the plugin's environment to a state that it can expect when it is
    // fully incorporated into the psi4 source tree.
}
