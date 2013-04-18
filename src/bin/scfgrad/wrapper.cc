#include "psi4-dec.h"
#include <libmints/mints.h>
#include <liboptions/liboptions.h>
#include <libchkpt/chkpt.h>
#include <libciomr/libciomr.h>
#include "scf_grad.h"

namespace psi{ 
namespace scfgrad {

PsiReturnType scfgrad(Options &options)
{
    tstart();

    boost::shared_ptr<SCFGrad> grad(new SCFGrad());
    SharedMatrix G = grad->compute_gradient();

    Process::environment.set_gradient(G); 
    Process::environment.wavefunction()->set_gradient(G);

    tstop();

    return Success;
}

}} // End Namespaces
