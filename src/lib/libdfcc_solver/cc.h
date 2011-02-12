#ifndef CC_H
#define CC_H

#include <libpsio/psio.hpp>
#include <libmints/wavefunction.h>
#include <libmints/basisset.h>
#include <psi4-dec.h>

using namespace psi;

namespace psi { namespace dfcc {

class CC : public Wavefunction {
private:
  void get_options();
  void get_params();
  void get_ribasis();

  void print_header();

protected:
  shared_ptr<BasisSet> ribasis_;
  shared_ptr<BasisSet> zero_;

  // MO basis parameters
  int nso_;
  int nmo_;
  int nocc_; // Total occupied orbitals
  int nvir_; // Total virtual orbitals
  int nfocc_; // Frozen occupieds
  int nfvir_; // Frozen virtuals
  int naocc_; // Active occupieds
  int navir_; // Active virtuals
  int nirrep_; // Number of irreps (must be 1 for now)

  // DF basis parameters
  int ndf_; // Number of DF auxiliary functions

  // Reference wavefunction
  double Eref_;
  double *evals_;
  double **C_;

public:
  CC(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
  virtual ~CC();

  virtual double compute_energy()=0;
};

}}

#endif
