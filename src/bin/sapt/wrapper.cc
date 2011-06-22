#include <libmints/mints.h>
#include <libqt/qt.h>

#include <libsapt_solver/sapt0.h>
#include <libsapt_solver/sapt_dft.h>
#include "wrapper.h"

using namespace boost;

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
  } else if (options.get_str("SAPT_LEVEL") == "MP2C") {
    MP2C sapt(options, psio, chkpt);
    sapt.compute_energy();
  } else {
    throw PSIEXCEPTION("Unrecognized SAPT type");
  }


  // Shut down psi.

  tstop();

  return Success;
}

}}

