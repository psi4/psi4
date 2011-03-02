/*
 *  hf.h
 *  matrix
 *
 *  Created by Justin Turney on 4/9/08.
 *  Copyright 2008 by Justin M. Turney, Ph.D.. All rights reserved.
 *
 */
//#define _MKL

#ifndef HF_H
#define HF_H

#include <libpsio/psio.hpp>
#include <libmints/wavefunction.h>
#include <libmints/basisset.h>
#include <libdiis/diismanager.h>
#include <psi4-dec.h>

//#ifndef MEMORY_SAFETY_FACTOR
#define MEMORY_SAFETY_FACTOR 0.20

using namespace psi;

namespace psi {
    class TwoBodySOInt;
    namespace scf {

class PseudospectralHF;
class DFHF; 

class HF : public Wavefunction {
protected:
    SharedMatrix H_;
    SharedMatrix S_;
    SharedMatrix Shalf_;
    SharedMatrix X_;
    SharedMatrix Sphalf_;

    /// Previous iteration's energy and current energy
    double Eold_;
    double E_;

    /// The RMS error in the density
    double Drms_;

    /// Max number of iterations for HF
    int maxiter_;

    /// Current Iteration
    int iteration_;

    /// Nuclear repulsion energy
    double nuclearrep_;

    /// DOCC vector from input (if found)
    bool input_docc_;

    /// SOCC vector from input (if found)
    bool input_socc_;

    /// Number of alpha and beta electrons
    int nalpha_, nbeta_;
    /// Number of alpha and beta electrons per irrep
    int nalphapi_[8], nbetapi_[8];

    //Initial SAD doubly occupied may be more than ndocc
    int sad_nocc_[8];

    //Canonical or Symmetric orthogonalization?
    bool canonical_X_;

    /// Mapping arrays
    int *so2symblk_;
    int *so2index_;

    /// Pairs needed for PK supermatrix
    size_t pk_pairs_;
    size_t pk_size_;
    int *pk_symoffset_;

    /// SCF algorithm type
    std::string scf_type_;

    /// Perturb the Hamiltonian?
    int perturb_h_;
    /// How big of a perturbation
    double lambda_;
    /// With what...
    enum perturb { nothing, dipole_x, dipole_y, dipole_z };
    perturb perturb_;

    /// The value below which integrals are neglected
    double integral_threshold_;


    /// The SO integral generator.  Only ever constructed if needed
    boost::shared_ptr<TwoBodySOInt> eri_;

    /// Pseudospectral stuff 
    shared_ptr<PseudospectralHF> pseudospectral_;
    /// DF stuff 
    shared_ptr<DFHF> df_;

    /// DF Storage Scheme
    enum df_storage { double_core, core, flip_B_core, flip_B_disk, k_incore, disk};
    df_storage df_storage_;
    unsigned long int df_memory_;

    //Density Fitting?
    int naux_raw_; //Number of functions in the raw auxiliary basis
    int naux_fin_; //Number of functions in the finished auxiliary basis
    int *ri_back_map_;
    int *ri_pair_nu_;
    int *ri_pair_mu_;
    int sig_fun_pairs_;
    int sig_shell_pairs_;

    int *schwarz_shell_pairs_;
    int *schwarz_fun_pairs_;

    //Local K? (only with DF)
    bool local_K_;

    /// do we need Coulomb?
    bool J_is_required_;
    /// do we need Exchange?
    bool K_is_required_;

    /// Three Index tensor for DF-SCF
    double **B_ia_P_;
    double **A_ia_P_;

    // RI Basis
    shared_ptr<BasisSet> ribasis_;
    /// Poisson Basis
    shared_ptr<BasisSet> poissonbasis_;

    /// Fitting metric (J-matrix) for DF-SCF
    double **W_;
    /// Fitting metric decomposition (varies by fitting algorithm)
    double **Winv_;

    double schwarz_; //Current Schwarz magnitude (static for now)
    int ntri_; //Number of function pairs after schwarz sieve and subsequent sieves
    int ntri_naive_; //Number of function pairs after schwarz sieve


    /// save SCF Cartesian Grid
    bool save_grid_;

    /// DIIS manager intiialized?
    bool initialized_diis_manager_;
    /// DIIS manager for all SCF wavefunctions
    boost::shared_ptr<DIISManager> diis_manager_;

    /// How many min vectors for DIIS
    int min_diis_vectors_;
    /// How many max vectors for DIIS
    int max_diis_vectors_;
    /// When do we start collecting vectors for DIIS
    int diis_start_;
    /// Are we even using DIIS?
    int diis_enabled_;

public:
    /// Exactly what their name says
    std::vector<SharedSimpleMatrix> Dipole_;
    std::vector<SharedSimpleMatrix> Quadrupole_;

    /// Nuclear contributions
    SimpleVector nuclear_dipole_contribution_;
    SimpleVector nuclear_quadrupole_contribution_;

    /// Formation of H is the same regardless of RHF, ROHF, UHF
    void form_H();

    /// Formation of S^+1/2 and S^-1/2 are the same
    void form_Shalf();

    /// Perform casting of basis set if desired.
    shared_ptr<Matrix> dualBasisProjection(SharedMatrix _C, int *_noccpi, shared_ptr<BasisSet> _old, shared_ptr<BasisSet> _new);

    /// UHF Atomic Density Matrix for SAD
    /// returns atomic_basis->nbf() x atomic_basis_->nbf() double array of approximate atomic density (summed over spin)
    void getUHFAtomicDensity(shared_ptr<BasisSet> atomic_basis, int n_electrons, int multiplicity, double** D);
    // Computes the C and D matrix in place for SAD Atomic UHF
    void atomicUHFHelperFormCandD(int nelec, int norbs,double** Shalf, double**F, double** C, double** D);

    /// Set the amount of information to print
    void set_print(const int n) {print_ = n;}

    /// The number of iterations needed to reach convergence
    int iterations_needed() {return iterations_needed_;}

    /// The RMS error in the density
    double rms_density_error() {return Drms_;}
protected:
    /// Common initializer
    void common_init();

    /// Clears memory and closes files (Should they be open) prior to correlated code execution
    virtual void finalize();

    /// Compute multipole integrals
    void form_multipole_integrals();

    /// Figure out how to occupy the orbitals in the absence of DOCC and SOCC
    void find_occupation();

    /// Determine how many core and virtual orbitals to freeze
    int *compute_fcpi(int nfzc, SharedVector &eigvalues);
    int *compute_fvpi(int nfvc, SharedVector &eigvalues);

    /// Forms the _so2* mapping arrays and determines _pk_pairs
    void form_indexing();

    /// Prints some opening information
    void print_header();

    /// The amout of information to print
    int print_;

    /// The number of electrons
    int nelectron_;

    /// The charge of the system
    int charge_;

    /// The multiplicity of the system (specified as 2 Ms + 1)
    int multiplicity_;

    /// The number of iterations need to reach convergence
    int iterations_needed_;

    /// Compute energy for the iteration.
    virtual double compute_E() = 0;

    /** Read in C from checkpoint. Default implementation works for RHF and ROHF. UHF needs to read in additional C.
     *  If unable to load C from checkpoint, will call form_C to compute the value.
     *  If unable to load call compute_initial_E(), else loads SCF energy from checkpoint.
     * @returns true if successfully loaded from checkpoint, else false.
     */
    virtual bool load_or_compute_initial_C();

    /** Computes the density matrix (D_) */
    virtual void form_D() {}

    /** Compute the MO coefficients (C_) */
    virtual void form_C() {}

    /** Computes the initial MO coefficients (default is to call form_C) */
    virtual void form_initial_C() { form_C(); }

    /** Computes the initial energy. */
    virtual double compute_initial_E() { return 0.0; }

    void sort_cholesky(double*, int*, int);
    /** Form canonical three-index Cholesky tensor */
    void form_CD();
    /** Form canonical three-index DF tensor */
    void form_B();
    /** Form three-index DF tensor with Poisson fitting */
    void form_B_Poisson();
    /** Form B without metric transform for local K (makes J go crazy fast)*/
    void form_A();
    /** Computes the J and/or K matrices according to the scf_type keyword and the functor passed in*/
    template <class JKFunctor> void process_tei(JKFunctor & functor);

    /** Write tensor from memory to disk */
    void write_B();
    /** Free all memory associated with DF */
    void free_B();
    /** Free all memory associated with DF */
    void free_A();

    /** Form DF fitting tensor **/
    void form_W();
    /** Form DF fitting tensor w/ Poisson fitting **/
    void form_W_Poisson();
    /** Form DF fitting tensor inverse square root without conditioning **/
    void form_Wm12_raw();
    /** Form DF fitting tensor inverse square root with preconditioning via change of basis **/
    void form_Wm12_fin();
    /** Form DF fitting tensor square root cholesky decomposition **/
    void form_Wp12_chol();


    inline int integral_type(int i, int j, int k, int l)
    {
        int type;

        if (i == j && i == k && i == l)     // (ij|kl)  (11|11)
            type = 1;
        else if (i == j && k == l && i > k) // (ij|kl)  (22|11)
            type = 2;
        else if (i == j && i == k && i > l) // (ij|kl)  (22|21)
            type = 3;
        else if (j == k && j == l && i > j) // (ij|kl)  (21|11)
            type = 4;
        else if (i == k && j == l)          // (ij|kl)  (21|21)
            type = 5;
        else if (i == j)                    // (ij|kl)  (33|21)
            type = 6;
        else if (j >  k && k == l)          // (ij|kl)  (32|11)
            type = 7;
        else if (k == l)                    // (ij|kl)  (31|22)
            type = 8;
        else if (i == k)                    // (ij|kl)  (32|31)
            type = 9;
        else if (j == k)                    // (ij|kl)  (32|21)
            type = 10;
        else if (j == l)                    // (ij|kl)  (31|21)
            type = 11;
        else if (j >  k)                    // (ij|kl)  (43|21)
            type = 12;
        else if (j >  l)                    // (ij|kl)  (42|31)
            type = 13;
        else                                // (ij|kl)  (41|32)
            type = 14;

        return type;
    }
public:
    HF(Options& options, shared_ptr<PSIO> psio, shared_ptr<Chkpt> chkpt);
    HF(Options& options, shared_ptr<PSIO> psio);

    virtual ~HF();
};

}} // Namespaces

#endif
