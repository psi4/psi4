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
  MP3(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
  virtual ~MP3();

  virtual double compute_energy();

};

class PSMP3 : public CC {
private:
  void print_header();

protected:
  shared_ptr<DFTensor> dfints_;

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
  PSMP3(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
  virtual ~PSMP3();

  virtual double compute_energy();

};

}}

#endif
