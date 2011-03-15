#ifndef MP2_H
#define MP2_H

#include "cc.h"

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class PSIO;
class Chkpt;

namespace dfcc {

class MP2 : public CC {
private:

protected:
  void print_header();
  void print_energies();
  std::string mp2_algorithm_;

  void compute_DF_MP2();
  void compute_OS_MP2();
  void compute_PS_MP2();
  void compute_PS2_MP2();
  void compute_PS3_MP2();
  void common_init();

public:
  MP2(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt);
  virtual ~MP2();

  virtual double compute_energy();

};

}}

#endif
