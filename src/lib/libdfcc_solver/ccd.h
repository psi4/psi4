#ifndef CCD_H
#define CCD_H

#include "cc.h"

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class PSIO;
class Chkpt;

namespace dfcc {

class CCD : public CC {
private:
  void print_header();

protected:
  boost::shared_ptr<DFCCDIIS> diis_;

  double *tIAJB_;
  double *t2IAJB_;
  double *vIAJB_;
  double *xIAJB_;

  void apply_denom();
  void symmetrize();
  void term_1();
  void term_2();
  void term_3();
  void term_4();
  void term_5();
  void term_6();
  void term_7();
  void term_8();
  void term_9();
  void term_10();
  void term_11();
  void term_12();

  double energy();
  double store_error_vecs();

public:
  CCD(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt);
  virtual ~CCD();

  virtual double compute_energy();

};

}}

#endif
