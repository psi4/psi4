#include"ccsd.h"
#include"frozen_natural_orbitals.h"
#include<libciomr/libciomr.h>
#include<libtrans/integraltransform.h>
#include<libtrans/mospace.h>

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

          // transform integrals
          tstart();
          fprintf(outfile,"        ==> Transform all two-electron integrals <==\n");
          fprintf(outfile,"\n");

          std::vector<shared_ptr<MOSpace> > spaces;
          spaces.push_back(MOSpace::all);
          boost::shared_ptr<IntegralTransform> ints(new IntegralTransform(wfn, spaces, IntegralTransform::Restricted,
                     IntegralTransform::IWLOnly, IntegralTransform::QTOrder, IntegralTransform::OccAndVir, false));
          ints->set_dpd_id(0);
          ints->set_keep_iwl_so_ints(true);
          ints->set_keep_dpd_so_ints(true);
          ints->initialize();
          ints->transform_tei(MOSpace::all, MOSpace::all, MOSpace::all, MOSpace::all);
          tstop();

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
