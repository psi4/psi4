#include<libciomr/libciomr.h>
#include <libefp_solver/efp_solver.h>
#include "efp.h"

namespace psi { namespace efp {

PsiReturnType efp_init(Options & options)
{
    // new efp object
    boost::shared_ptr<EFP> myefp(new EFP(options));

    // set efp object in process environment
    Process::environment.set_efp(myefp);

    return Success;
}

}}
