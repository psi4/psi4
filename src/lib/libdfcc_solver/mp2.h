#ifndef MP2_H
#define MP2_H

#include <libpsio/psio.hpp>
#include <libmints/wavefunction.h>
#include <libmints/basisset.h>
#include <psi4-dec.h>

#include "cc.h"

using namespace psi;

namespace psi { namespace dfcc {

class MP2 : public CC {
private:
  virtual void print_header();

protected:

public:
  MP2(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
  virtual ~MP2();

  virtual double compute_energy();

};

}}

#endif
