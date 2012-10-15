#include<libmints/wavefunction.h>
#include"coupledpair.h"

using namespace boost;

namespace psi{ namespace cepa{

PsiReturnType cepa(Options&options){

  boost::shared_ptr<Wavefunction> wfn  = Process::environment.wavefunction();
  boost::shared_ptr<Wavefunction> cepa = boost::shared_ptr<Wavefunction>(new CoupledPair(wfn, options));
  Process::environment.set_wavefunction(cepa);
  cepa->compute_energy();

  return Success;
}

}}
//
//  Comments so that autodoc utility will find these PSI variables
//
//  Process::environment.globals["CEPA(0) DIPOLE X"] =
//  Process::environment.globals["CISD DIPOLE X"] =
//  Process::environment.globals["ACPF DIPOLE X"] =
//  Process::environment.globals["AQCC DIPOLE X"] =
//
//  Process::environment.globals["CEPA(0) DIPOLE Y"] =
//  Process::environment.globals["CISD DIPOLE Y"] =
//  Process::environment.globals["ACPF DIPOLE Y"] =
//  Process::environment.globals["AQCC DIPOLE Y"] =
//
//  Process::environment.globals["CEPA(0) DIPOLE Z"] =
//  Process::environment.globals["CISD DIPOLE Z"] =
//  Process::environment.globals["ACPF DIPOLE Z"] =
//  Process::environment.globals["AQCC DIPOLE Z"] =
//
//  Process::environment.globals["CEPA(0) QUADRUPOLE XX"] =
//  Process::environment.globals["CISD QUADRUPOLE XX"] =
//  Process::environment.globals["ACPF QUADRUPOLE XX"] =
//  Process::environment.globals["AQCC QUADRUPOLE XX"] =
//
//  Process::environment.globals["CEPA(0) QUADRUPOLE YY"] =
//  Process::environment.globals["CISD QUADRUPOLE YY"] =
//  Process::environment.globals["ACPF QUADRUPOLE YY"] =
//  Process::environment.globals["AQCC QUADRUPOLE YY"] =
//
//  Process::environment.globals["CEPA(0) QUADRUPOLE ZZ"] =
//  Process::environment.globals["CISD QUADRUPOLE ZZ"] =
//  Process::environment.globals["ACPF QUADRUPOLE ZZ"] =
//  Process::environment.globals["AQCC QUADRUPOLE ZZ"] =
//
//  Process::environment.globals["CEPA(0) QUADRUPOLE XY"] =
//  Process::environment.globals["CISD QUADRUPOLE XY"] =
//  Process::environment.globals["ACPF QUADRUPOLE XY"] =
//  Process::environment.globals["AQCC QUADRUPOLE XY"] =
//
//  Process::environment.globals["CEPA(0) QUADRUPOLE XZ"] =
//  Process::environment.globals["CISD QUADRUPOLE XZ"] =
//  Process::environment.globals["ACPF QUADRUPOLE XZ"] =
//  Process::environment.globals["AQCC QUADRUPOLE XZ"] =
//
//  Process::environment.globals["CEPA(0) QUADRUPOLE YZ"] =
//  Process::environment.globals["CISD QUADRUPOLE YZ"] =
//  Process::environment.globals["ACPF QUADRUPOLE YZ"] =
//  Process::environment.globals["AQCC QUADRUPOLE YZ"] =
//
