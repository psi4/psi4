#include <psi4-dec.h>
#include <libparallel/parallel.h>
#include <libmints/wavefunction.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp>
#include <libyeti/env.h>
#include "plugin.h"

using namespace boost;
using namespace psi;

extern "C" void init_plugin(const shared_ptr<Communicator>& comm, const Process::Environment& env,
                            const shared_ptr<Chkpt> &chkpt, const shared_ptr<PSIO> &psio,
                            const yeti::Env& yetiEnv)
{
    // Prepares the plugin's environment to a state that it can expect when it is
    // fully incorporated into the psi4 source tree.
//    Wavefunction::initialize_singletons();
//    Communicator::world = comm;
//    Process::environment = env;
//    _default_psio_lib_ = psio;
//    _default_chkpt_lib_ = chkpt;
//    psi::yetiEnv = yetiEnv;
}
