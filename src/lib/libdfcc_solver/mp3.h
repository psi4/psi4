#ifndef MP3_H
#define MP3_H

#include <libpsio/psio.hpp>
#include <libmints/wavefunction.h>
#include <libmints/basisset.h>
#include <lib3index/3index.h>
#include <psi4-dec.h>

#include "cc.h"

using namespace psi;

namespace psi { namespace dfcc {

class MP3 : public CC {
private:
  void print_header();

protected:
  shared_ptr<DFTensor> dfints_;

public:
  MP3(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
  virtual ~MP3();

  virtual double compute_energy();

};

}}

#endif
