#include "cc.h"

using namespace std;
using namespace psi;

namespace psi { namespace dfcc {

CC::CC(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt)
    : Wavefunction(options, psio, chkpt)
{
}

CC::~CC()
{
}

}}
