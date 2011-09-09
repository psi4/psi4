#ifndef _psi_src_bin_mcscf_scf_h_
#define _psi_src_bin_mcscf_scf_h_

#include "sblock_vector.h"
#include "sblock_matrix.h"
#include <libmoinfo/libmoinfo.h>
#include <libmints/wavefunction.h>
#include <libutil/memory_manager.h>

#define STORE_TEI 0

namespace boost {
template<class T> class shared_ptr;
}

namespace psi{ namespace mcscf{

enum ReferenceType {rhf, rohf, uhf, tcscf};

class SCF  : public Wavefunction
{
public:
  explicit SCF(Options& options_, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt_);
  ~SCF();
  double compute_energy();

private:
  ReferenceType reference;
  static const int maxci   = 8;
  static const int maxdiis = 10;
  static const int maxbatches = 50;

  int  nirreps;
  int  nso;
  size_t* ioff;
  intvec sopi;
  intvec docc;
  intvec actv;
  int  turn_on_actv;
  int  root;
  int  ndiis;
  int  current_diis;
  bool use_diis;
  double total_energy;

  // Out-of-core algorithm
  bool   out_of_core;
  int    nbatch;
  size_t batch_pq_min[maxbatches];
  size_t batch_pq_max[maxbatches];
  size_t batch_index_min[maxbatches];
  size_t batch_index_max[maxbatches];
  size_t batch_size[maxbatches];
  size_t total_symmetric_block_size;
  size_t nin_core;                    // Number of matrix elements of the PK and K matrices held in of core

  // Addressing routines
  int* block_offset;
  int  npairs;              // Total number of pairs
  int* pairpi;              // Number of pairs per irrep
  int* pairs;               // The pairs stored as [ p1 q1 p2 q2 ... ]
  int** pair;               // Maps absolute p and q to the pair index
  int** pair_sym;           // Maps absolute p and q to the pair symmetry
  int*  pair_offset;        // Offset for first pair of a given irrep

  double*  PK;              // PK(pq|rs)
  double*  K;               // K(pq|rs)

  SBlockVector epsilon;     // Effective Fock matrix eigenvalues


  SBlockMatrix C;           // MO coefficients matrix
  SBlockMatrix C_t;         // transformed C matrix
  SBlockMatrix C_T;         // MO coefficients matrix transposed
  SBlockMatrix Dc;          // The density matrix (closed)
  SBlockMatrix Do;          // The density matrix (open)
  SBlockMatrix Dtc[maxci];  // The density matrix (tcscf)
  SBlockMatrix Dsum[maxci]; // The density matrix (closed + tcscf)
  SBlockMatrix Fc;          // The Fock matrix (closed)
  SBlockMatrix Fc_t;        // The transformed Fock matrix (closed)
  SBlockMatrix Fo;          // The Fock matrix (open)
  SBlockMatrix Fo_t;        // The transformed Fock matrix (open)
  SBlockMatrix Favg;        // The Fock matrix (average)
  SBlockMatrix Favg_t;      // The transformed Fock matrix (average)
  SBlockMatrix Ftc[maxci];  // The Fock matrix (tcscf)
  SBlockMatrix Ftc_t[maxci];// The transformed Fock matrix (tcscf)
  SBlockMatrix Feff_t;      // The transformed effective Fock matrix
  SBlockMatrix Feff_t_old;  // The transformed effective Fock matrix
  SBlockMatrix Feff_oAO;    // The effective Fock matrix in the orthogonal AO basis
  SBlockMatrix G;           // The G matrix
  SBlockMatrix T;           // a temp matrix
  SBlockMatrix H;           // one electron integrals
  SBlockMatrix O;           // The occupation matrix
  SBlockMatrix S;           // overlap integrals
  SBlockMatrix S_sqrt_inv;  // S^-1/2
  SBlockMatrix S_sqrt;      // S^1/2

  SBlockMatrix e;                         // MO coefficients error matrix
  SBlockMatrix diis_F[maxdiis];           // The Fock matrices saved for DIIS
  SBlockMatrix diis_e[maxdiis];           // The error matrices saved for DIIS
  double       diis_ci[maxci][maxdiis];   // The ci vector saved for DIIS

  // TWOCON Specific
  int          nci;         // Number of references
  double   norm_ci_grad;    // Norm of the CI gradient (sum_I |C_I|)
  double*       ci;         // TWOCON CI coefficients
  double*       ci_grad;    // TWOCON CI coefficients gradient
  double** H_tcscf;         // TWOCON Hamiltonian
  int tcscf_mos[maxci];     // Number of the TWOCON mos (relative to the irrep)
  int tcscf_sym[maxci];     // Symmetry of the TWOCON mos

  // Private functions
  void startup();
  void cleanup();

  void generate_pairs();
  void read_so_oei();
  void read_so_tei();
  void read_so_tei_form_PK();
  void read_so_tei_form_PK_and_K();
  void read_Raffanetti(const char* integral_type, double* integrals, int batch);
  void write_Raffanetti(const char* integral_type, double* integrals, int batch);

  void construct_S_inverse_sqrt();
  void iterate_scf_equations();
  void guess_occupation();
  void canonicalize_MO();
  void save_info();
  void check_orthonormality();
  void print_eigenvectors_and_MO();

  void initial_guess();
  void density_matrix();
  void construct_G(SBlockMatrix& density,SBlockMatrix& G,double* integrals,int batch, double factor);
  void construct_G(SBlockMatrix& density,SBlockMatrix& G,double* integrals,int batch);
  void construct_F();
  void construct_Favg();
  void construct_Feff(int cycle);
  double energy(int cycle,double old_energy);

  void diis(int cycle);

  // Auxiliary functions
  void transform(SBlockMatrix& Initial, SBlockMatrix& Final, SBlockMatrix& Transformation);
};

}} /* End Namespaces */

#endif // _psi_src_bin_mcscf_scf_h_

/*
//   void read_C();



  int* so_sym;

  // TWOCON



  // DIIS
  int maxdiis;
  int current_diis;

  double*      epsilon;     // Eigenvalues
  SBlockMatrix* C;           // MO coefficients matrix
  SBlockMatrix* C_T;         // Transpose coefficients matrix
  SBlockMatrix* C_t;         // transformed C matrix
  SBlockMatrix* D;           // density matrix
  SBlockMatrix* Dc;          // closed-shell density matrix
  SBlockMatrix* Do;          // open-shell density matrix
  SBlockMatrix* Dtc[2];      // tcscf density matrices
  SBlockMatrix* Ex;          // tcscf exchange matrix
  SBlockMatrix* F;           // Fock matrix
  SBlockMatrix* F_t;         // transformed Fock matrix
  SBlockMatrix* Fc;          // closed-shell Fock matrix
  SBlockMatrix* Fo;          // open-shell   Fock matrix

  SBlockMatrix* Fc_t;        // transformed closed-shell Fock matrix
  SBlockMatrix* Fo_t;        // transformed open-shell   Fock matrix
  SBlockMatrix* Feff_t;      // transformed effective Fock matrix
  SBlockMatrix* G;           // two electron contribution to F

  SBlockMatrix* SDF;         // SDF: used to compute the diis error
  SBlockMatrix* FDS;         // FDS: used to compute the diis error
  SBlockMatrix* CSC;         // CSC = 1 for orthornormal MOs


  SBlockMatrix* Temp;        // A temporary matrix



*/
