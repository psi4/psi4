#include <libmints/mints.h>
#include <libqt/qt.h>

#include <libsapt_solver/sapt0.h>
#include "wrapper.h"

using namespace boost;

namespace psi { namespace sapt {

std::string to_string(const int val);   // In matrix.cpp

PsiReturnType sapt(Options & options)
{
  tstart();

  shared_ptr<PSIO> psio(new PSIO);

  shared_ptr<Chkpt> chkpt(new Chkpt(psio, PSIO_OPEN_OLD));

  // Initialize the psi3 timer library.
  timer_init();

  if (options.get_str("SAPT_LEVEL") == "SAPT0") {
    SAPT0 sapt(options, psio, chkpt);
    sapt.compute_energy();
  }

  // Shut down psi.
  timer_done();

  tstop();

  return Success;
}

}}

