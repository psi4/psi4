#include<psi4-dec.h>
#include<libmints/wavefunction.h>
#include"coupledpair.h"

using namespace boost;

namespace psi{ namespace cepa{

PsiReturnType cepa(Options&options){

  boost::shared_ptr<psi::Wavefunction> wfn = Process::environment.wavefunction();
  boost::shared_ptr<Wavefunction> cepa = boost::shared_ptr<Wavefunction>(new CoupledPair(wfn,options));

  //Process::environment.set_wavefunction(cepa);

  cepa->compute_energy();

  return Success;
}

}}


