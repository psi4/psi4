#include"psi4-dec.h"
#include<libciomr/libciomr.h>
#include<libmints/wavefunction.h>

#include"coupledpair.h"

using namespace psi;
using namespace boost;

namespace psi{ namespace cepa{
  void OPDM(boost::shared_ptr<psi::cepa::CoupledPair>cepa,Options&options);
  void RunCoupledPair(Options &options,boost::shared_ptr<psi::Wavefunction> wfn);
  void TransformIntegrals(boost::shared_ptr<Wavefunction>wfn,Options&options);
  void SortIntegrals(int nfzc,int nfzv,int norbs,int ndoccact,int nvirt,bool isdirect,bool islocal);
}}

namespace psi{ namespace cepa{

PsiReturnType cepa(Options&options){
  boost::shared_ptr<psi::Wavefunction> ref = Process::environment.reference_wavefunction();
  RunCoupledPair(options,ref);
  return Success;
}

void RunCoupledPair(Options &options,boost::shared_ptr<psi::Wavefunction> wfn){

  boost::shared_ptr<CoupledPair> cepa(new CoupledPair(wfn,options));
  PsiReturnType status;

  // integral transformation.  only needed if not cim or integral direct
  if ( !wfn->isCIM() && options.get_bool("CEPA_VABCD_DIRECT")) {
     tstart();
     TransformIntegrals(wfn,options);
     tstop();
  }

  // integral sort
  tstart();
  SortIntegrals(cepa->nfzc,cepa->nfzv,cepa->nmo+cepa->nfzc+cepa->nfzv,cepa->ndoccact,cepa->nvirt,options.get_bool("CEPA_VABCD_DIRECT"),wfn->isCIM());
  tstop();

  // solve cepa equations
  tstart();
  cepa->WriteBanner(options);
  cepa->AllocateMemory(options);
  status = cepa->CEPAIterations(options);
  tstop();

  // mp2 energy
  Process::environment.globals["MP2 CORRELATION ENERGY"] = cepa->emp2;
  Process::environment.globals["MP2 TOTAL ENERGY"] = cepa->emp2 + cepa->escf;

  // cepa energy
  char*cepatype = (char*)malloc(100*sizeof(char));
  if (cepa->cepa_level == 0){
     Process::environment.globals["CEPA(0) CORRELATION ENERGY"] = cepa->ecepa;
     Process::environment.globals["CEPA(0) TOTAL ENERGY"] = cepa->ecepa + cepa->escf;
  }
  if (cepa->cepa_level == 1){
     Process::environment.globals["CEPA(1) CORRELATION ENERGY"] = cepa->ecepa;
     Process::environment.globals["CEPA(1) TOTAL ENERGY"] = cepa->ecepa + cepa->escf;
  }
  if (cepa->cepa_level == 2){
     Process::environment.globals["CEPA(2) CORRELATION ENERGY"] = cepa->ecepa;
     Process::environment.globals["CEPA(2) TOTAL ENERGY"] = cepa->ecepa + cepa->escf;
  }
  if (cepa->cepa_level == 3){
     Process::environment.globals["CEPA(3) CORRELATION ENERGY"] = cepa->ecepa;
     Process::environment.globals["CEPA(3) TOTAL ENERGY"] = cepa->ecepa + cepa->escf;
  }
  if (cepa->cepa_level == -1){
     Process::environment.globals["CISD CORRELATION ENERGY"] = cepa->ecepa;
     Process::environment.globals["CISD TOTAL ENERGY"] = cepa->ecepa + cepa->escf;
  }
  if (cepa->cepa_level == -2){
     Process::environment.globals["ACPF CORRELATION ENERGY"] = cepa->ecepa;
     Process::environment.globals["ACPF TOTAL ENERGY"] = cepa->ecepa + cepa->escf;
  }
  if (cepa->cepa_level == -3){
     Process::environment.globals["AQCC CORRELATION ENERGY"] = cepa->ecepa;
     Process::environment.globals["AQCC TOTAL ENERGY"] = cepa->ecepa + cepa->escf;
  }
  Process::environment.globals["CURRENT ENERGY"] = cepa->ecepa + cepa->escf;
  Process::environment.globals["CURRENT CORRELATION ENERGY"] = cepa->ecepa;


  // dipole moments for CISD, ACPF, or AQCC
  if (options.get_bool("DIPMOM")){
     if (cepa->cepa_level>0) 
        throw PsiException("coupled-pair dipole moments available only for CEPA(0), CISD, ACFP, and AQCC",__FILE__,__LINE__);
     OPDM(cepa,options);
  }

  // free memory
  free(cepa->integrals);
  free(cepa->tempt);
  free(cepa->tempv);
  if (!cepa->t2_on_disk){
     free(cepa->tb);
  }
  free(cepa->w1);
  free(cepa->t1);
  free(cepa->I1);
  free(cepa->I1p);
  free(cepa->diisvec);

  free(cepatype);
  cepa.reset();
}

}}


