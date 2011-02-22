#ifndef CCD_H
#define CCD_H

#include <libpsio/psio.hpp>
#include <libmints/wavefunction.h>
#include <libmints/basisset.h>
#include <psi4-dec.h>

#include "cc.h"

#define DFCC_INT_FILE 56  // temporary
#define DFCC_DIIS_FILE 42  // temporary

using namespace psi;

namespace psi { namespace dfcc {

class CCD : public CC {
private:
  void print_header();

  void df_integrals();
  void mo_integrals();

protected:
  shared_ptr<DFCCDIIS> diis_;

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
  CCD(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
  virtual ~CCD();

  virtual double compute_energy();

};

}}

#endif
