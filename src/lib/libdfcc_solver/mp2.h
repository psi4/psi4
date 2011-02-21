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

protected:
  void print_header();
  std::string mp2_algorithm_;
  double E_MP2J_;
  double E_MP2K_;  
 
  void compute_DF_MP2();
  void compute_OS_MP2();
  void compute_PS_MP2();
  void compute_PS2_MP2();
  void compute_PS3_MP2();
  void common_init();

public:
  MP2(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
  virtual ~MP2();

  virtual double compute_energy();

};

}}

#endif
