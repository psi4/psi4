#ifndef CC_H
#define CC_H

#include <map>
#include <libmints/wavefunction.h>

#define DFCC_INT_FILE 56  // temporary
#define DFCC_DIIS_FILE 42  // temporary

// This is a dirty hack...bite me
#define INDEX(i,j) ((i>=j) ? (ioff[i] + j) : (ioff[j] + i))

namespace boost {
template<class T> class shared_ptr;
}

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

protected:

  // Print level
  int print_;
  // Debug level
  int debug_;

  boost::shared_ptr<BasisSet> ribasis_;
  boost::shared_ptr<BasisSet> zero_;

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
  boost::shared_ptr<Vector> evals_; // HF eigenvalues (all)
  boost::shared_ptr<Vector> evals_aocc_; // HF active occupied eigenvalues
  boost::shared_ptr<Vector> evals_avir_; // HF active virtual eigenvalues
  SharedMatrix C_; // HF eigenvectors (all)
  SharedMatrix C_aocc_; // HF eigenvectors (active occupieds)
  SharedMatrix C_avir_; // HF eigenvectors (active virtuals)

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
  CC(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt);
  virtual ~CC();
  virtual bool same_a_b_orbs() const { return reference_wavefunction_->same_a_b_orbs(); }
  virtual bool same_a_b_dens() const { return reference_wavefunction_->same_a_b_dens(); }

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
    boost::shared_ptr<PSIO> psio_;

public:
    DFCCDIIS(int, int, int, boost::shared_ptr<PSIO>);
    ~DFCCDIIS();

    /* Stores new vectors and error vectors in one call */
    void store_vectors(double *, double *);

    /* Stores new vectors and error vectors separately
       these should be called as a pair, followed by increment_vectors()
       or the DIIS will not work properly. */
    void store_current_vector(char *t_vec);
    void store_error_vector(char *err_vec);
    void increment_vectors();

    /* Performs the DIIS extrapolation with two arrays in core */
    void get_new_vector(double *, double *);

    /* Performs the DIIS extrapolation in chunks with one array */
    void get_new_vector(double **vec_i, int cols);

    char *get_last_vec_label();

};

}}

#endif
