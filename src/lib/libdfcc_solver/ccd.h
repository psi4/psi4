#ifndef CCD_H
#define CCD_H

#include <libpsio/psio.hpp>
#include <libmints/wavefunction.h>
#include <libmints/basisset.h>
#include <psi4-dec.h>

#include "cc.h"

#define DFCC_INT_FILE 56  // temporary

using namespace psi;

namespace psi { namespace dfcc {

class CCD : public CC {
private:
  void df_integrals();

protected:

public:
  CCD(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
  virtual ~CCD();

  virtual double compute_energy();

};

}}

#endif
