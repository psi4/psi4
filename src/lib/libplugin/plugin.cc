#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <libmints/wavefunction.h>
#include "plugin.h"

using namespace psi;

extern "C" void init_plugin(const shared_ptr<Communicator>& comm, const Process::Environment& env)
{
    // Prepares the plugin's environment to a state that it can expect when it is
    // fully incorporated into the psi4 source tree.
    Wavefunction::initialize_singletons();
    Communicator::world = comm;
    Process::environment = env;
}
