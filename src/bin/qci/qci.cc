#include"ccsd.h"
#include"frozen_natural_orbitals.h"

using namespace boost;

namespace psi{ namespace qci{

PsiReturnType qci(Options &options) {

  boost::shared_ptr<Wavefunction> wfn = Process::environment.wavefunction();

  // frozen natural orbital ccsd(t)
  if (options.get_bool("NAT_ORBS")) {
      boost::shared_ptr<FrozenNO> fno(new FrozenNO(options));
      wfn = (boost::shared_ptr<Wavefunction>)fno;
  }

  boost::shared_ptr<CoupledCluster> ccsd(new CoupledCluster(wfn,options));
  ccsd->compute_energy();

  return  Success;
} // end qci

}} // end namespaces
