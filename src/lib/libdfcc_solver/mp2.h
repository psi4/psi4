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

  void check_integrals(SharedMatrix E); 
  void test_denominators();
  void test_df();
  void test_ps();
  void test_ps_omega();
  void test_dps_omega();
  void compute_MP2();
  void compute_DF_MP2();
  void compute_PS_MP2();
  void compute_DF_MP2J();
  void compute_PS_MP2J();
  void compute_PS_MP2K();
  void compute_PS_DF_MP2K();
  void common_init();

public:
  MP2(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt);
  virtual ~MP2();

  virtual double compute_energy();

};

}}

#endif
