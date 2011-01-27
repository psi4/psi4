#include <libciomr/libciomr.h>
#include <libmints/mints.h>

#include <psi4-dec.h>

using namespace boost;

namespace psi { namespace mints {

PsiReturnType mints(Options & options)
{
    tstart();

    shared_ptr<MintsHelper> mints(new MintsHelper(options));
    mints->integrals();

    tstop();

    return Success;
}

}}
