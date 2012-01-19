#include <libmints/mints.h>
#include <libqt/qt.h>

#include <libsapt_solver/sapt0.h>
#include <libsapt_solver/sapt2.h>
#include <libsapt_solver/sapt2p.h>
#include <libsapt_solver/sapt2p3.h>
//#include <libsapt_solver/sapt_dft.h>
#include "wrapper.h"

namespace psi { namespace sapt {

std::string to_string(const int val);   // In matrix.cpp

PsiReturnType sapt(Options & options)
{
  tstart();

  boost::shared_ptr<PSIO> psio(new PSIO);
  boost::shared_ptr<Chkpt> chkpt(new Chkpt(psio, PSIO_OPEN_OLD));

  if (options.get_str("SAPT_LEVEL") == "SAPT0") {
    SAPT0 sapt(options, psio, chkpt);
    sapt.compute_energy();
  } else if (options.get_str("SAPT_LEVEL") == "SAPT2") {
    SAPT2 sapt(options, psio, chkpt);
    sapt.compute_energy();
  } else if (options.get_str("SAPT_LEVEL") == "SAPT2+") {
    SAPT2p sapt(options, psio, chkpt);
    sapt.compute_energy();
  } else if (options.get_str("SAPT_LEVEL") == "SAPT2+3") {
    SAPT2p3 sapt(options, psio, chkpt);
    sapt.compute_energy();
//  } else if (options.get_str("SAPT_LEVEL") == "MP2C") {
//    MP2C sapt(options, psio, chkpt);
//    sapt.compute_energy();
  } else {
    throw PSIEXCEPTION("Unrecognized SAPT type");
  }

  // Shut down psi.

  tstop();

  return Success;
}

}}

