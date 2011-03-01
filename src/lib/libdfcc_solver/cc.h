#ifndef CC_H
#define CC_H

#include <libmints/wavefunction.h>
#include <psi4-dec.h>

#define DFCC_INT_FILE 56  // temporary
#define DFCC_DIIS_FILE 42  // temporary

// This is a dirty hack...bite me
#define INDEX(i,j) ((i>=j) ? (ioff[i] + j) : (ioff[j] + i))

using namespace psi;

namespace psi { 

class BasisSet; 
class Matrix; 
class Vector; 
class PSIO;
class Options;
class PseudoGrid;

namespace dfcc {

class CC : public Wavefunction {
private:
  void get_options();
  void get_params();
  void get_ribasis();
  void get_dealiasbasis();
  void get_pseudogrid();

protected:
  shared_ptr<BasisSet> ribasis_;
  shared_ptr<BasisSet> dealias_;
  shared_ptr<BasisSet> zero_;
  shared_ptr<PseudoGrid> grid_;

  virtual void print_header(); 

  // Available memory in doubles
  unsigned long int doubles_;

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

  // Energies
  double Eref_; // HF Reference energy
  std::map<std::string, double> energies_; // Table of energies

  // Spin Treatment Parameters
  double sss_; // Same-spin scaling parameter
  double oss_; // Opposite-spin scaling parameter

  // Angles and dangles
  std::string denominator_algorithm_;
  double denominator_delta_;

  // Schwarz cutoff 
  double schwarz_cutoff_;

  // Maximum allowed fitting metric condition
  std::string fitting_algorithm_;
  double fitting_condition_;

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

  // Simple function to compute DF integrals
  // To be replaced with lib3index calls
  void df_integrals();
  void mo_integrals();

  // Functions to permute array indices
  void iajb_ibja(double *);
  void iajb_ijab(double *);
  void ijab_iajb(double *);

  // Fun with disk I/O
  void zero_disk(int, const char *, char *, int, int);

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
    DFCCDIIS(int, int, int, shared_ptr<PSIO>);
    ~DFCCDIIS();

    void store_vectors(double *, double *);
    void get_new_vector(double *, double *);
};

}}

#endif
