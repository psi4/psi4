#include"ccsd.h"
#include"frozen_natural_orbitals.h"
#include<libciomr/libciomr.h>

using namespace boost;

namespace psi{ namespace fnocc{

PsiReturnType fnocc(Options &options) {

  boost::shared_ptr<Wavefunction> wfn;

  if ( !options.get_bool("DFCC") ){

      // frozen natural orbital ccsd(t)
      if (options.get_bool("NAT_ORBS")) {
          boost::shared_ptr<FrozenNO> fno(new FrozenNO(Process::environment.wavefunction(),options));
          fno->ComputeNaturalOrbitals();
          wfn = (boost::shared_ptr<Wavefunction>)fno;
      }else {
          wfn = Process::environment.wavefunction();
      }

      if ( !options.get_bool("RUN_CEPA") ) {
          boost::shared_ptr<CoupledCluster> ccsd(new CoupledCluster(wfn,options));
          Process::environment.set_wavefunction(ccsd);
          ccsd->compute_energy();
      } else {
          boost::shared_ptr<CoupledPair> cepa (new CoupledPair(wfn,options));
          Process::environment.set_wavefunction(cepa);
          cepa->compute_energy();
      }

  }else {

      // generate frozen natural orbitals?
      if (options.get_bool("NAT_ORBS")){
          tstart();
          boost::shared_ptr<DFFrozenNO> fno(new DFFrozenNO(Process::environment.wavefunction(),options));
          fno->ComputeNaturalOrbitals();
          wfn = (boost::shared_ptr<Wavefunction>)fno;
          tstop();
      }else {
          wfn = Process::environment.wavefunction();
      }

      // ccsd(t)!
      boost::shared_ptr<DFCoupledCluster> ccsd (new DFCoupledCluster(wfn,options));
      ccsd->compute_energy();
  }

  return  Success;
} // end fnocc

}} // end namespaces
