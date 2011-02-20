#ifndef CC_H
#define CC_H

#include <libmints/wavefunction.h>
#include <psi4-dec.h>

using namespace psi;

namespace psi { 

class BasisSet; 
class Matrix; 
class Vector; 
class PSIO;
class Options;

namespace dfcc {

class CC : public Wavefunction {
private:
  void get_options();
  void get_params();
  void get_ribasis();
  void get_dealiasbasis();

protected:
  shared_ptr<BasisSet> ribasis_;
  shared_ptr<BasisSet> dealias_;
  shared_ptr<BasisSet> zero_;

  virtual void print_header(); 

  // MO basis parameters
  int nso_; // Total basis functions (AO)
  int nmo_; // Total basis functions (MO)
  int namo_; // Total active MO's
  int nocc_; // Total occupied orbitals
  int nvir_; // Total virtual orbitals
  int nfocc_; // Frozen occupieds
  int nfvir_; // Frozen virtuals
  int naocc_; // Active occupieds
  int navir_; // Active virtuals
  int nirrep_; // Number of irreps (must be 1 for now)

  // DF basis parameters
  int ndf_; // Number of DF auxiliary functions
  int ndealias_; // Number of dealias functions

  // Reference wavefunction
  double Eref_; // HF Reference energy

  // New libmints objects for HF information 
  shared_ptr<Vector> evals_; // HF eigenvalues (all)
  shared_ptr<Vector> evals_aocc_; // HF active occupied eigenvalues
  shared_ptr<Vector> evals_avir_; // HF active virtual eigenvalues
  shared_ptr<Matrix> C_; // HF eigenvectors (all)
  shared_ptr<Matrix> C_aocc_; // HF eigenvectors (active occupieds) 
  shared_ptr<Matrix> C_avir_; // HF eigenvectors (active virtuals) 

  // Pointers to the new libmints objects so's Ed doesn't get bitchy
  double *evalsp_;
  double *evals_aoccp_;
  double *evals_avirp_;
  double **Cp_;
  double **C_aoccp_;
  double **C_avirp_;

  // Functions to permute array indices
  void iajb_ibja(double *);
  void iajb_ijab(double *);

public:
  CC(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
  virtual ~CC();

  virtual double compute_energy()=0;
} ;

class DFCCDIIS {

private:
    int filenum_;
    char *vec_label_;
    char *err_label_;
    int max_diis_vecs_;

    int diis_file_;
    int vec_length_;

    int curr_vec_;
    int num_vecs_;

    char *get_err_label(int);
    char *get_vec_label(int);

protected:
    shared_ptr<PSIO> psio_;

public:
    DFCCDIIS(int, int, char *, char *, int, int, shared_ptr<PSIO>);
    ~DFCCDIIS();

    void store_vectors();
    void get_new_vector();
};

}}

#endif
