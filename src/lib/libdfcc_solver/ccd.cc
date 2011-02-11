#include "ccd.h"

using namespace std;
using namespace psi;

namespace psi { namespace dfcc {

CCD::CCD(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt)
    : CC(options, psio, chkpt)
{
}

CCD::~CCD()
{
}

double CCD::compute_energy()
{
    return(0.0);
}

}}
