#include "mp2.h"

using namespace std;
using namespace psi;

namespace psi { namespace dfcc {

MP2::MP2(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt)
  : CC(options, psio, chkpt)
{
}

MP2::~MP2()
{
}

double MP2::compute_energy()
{
  return(0.0);
}

void MP2::print_header()
{
}

}}
