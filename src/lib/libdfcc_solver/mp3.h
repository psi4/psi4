#ifndef MP3_H
#define MP3_H

#include "cc.h"

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {

class PSIO;
class Chkpt;
class DFTensor;

namespace dfcc {

class MP3 : public CC {
private:
  void print_header();

protected:
  double *tIAJB_;
  double *t2IAJB_;
  double *vIAJB_;

  void apply_denom(double *);
  void symmetrize();
  void term_1();
  void term_2();
  void term_3();
  void term_4();
  void term_5();

  double energy(double *);

public:
  MP3(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt);
  virtual ~MP3();

  virtual double compute_energy();

};

#if 0
class PSMP3 : public CC {
private:
  void print_header();

protected:
  boost::shared_ptr<DFTensor> dfints_;

  double *tIAJB_;
  double *t2IAJB_;
  double *vIAJB_;

  void apply_denom(double *);
  void symmetrize();
  void term_1();
  void term_2();
  void term_3();
  void term_4();
  void term_5();
  void term_6();

  double energy(double *);

  long int mem_; // Memory (in doubles)

public:
  PSMP3(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt);
  virtual ~PSMP3();

  virtual double compute_energy();

};

#endif

}}

#endif
