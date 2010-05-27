/*
 *  hf.h
 *  matrix
 *
 *  Created by Justin Turney on 4/9/08.
 *  Copyright 2008 by Justin M. Turney, Ph.D.. All rights reserved.
 *
 */

#ifndef HF_H
#define HF_H

#include <libpsio/psio.hpp>
#include <libmints/wavefunction.h>
#include <libdiis/diismanager.h>
#include <psi4-dec.h>

//#ifndef MEMORY_SAFETY_FACTOR
#define MEMORY_SAFETY_FACTOR 0.20

using namespace psi;

namespace psi { namespace scf {

class HF : public Wavefunction {
protected:
    SharedMatrix H_;
    SharedMatrix S_;
    SharedMatrix Shalf_;
    SharedMatrix Sphalf_;
    SharedMatrix C_;

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

    /// Number of atoms and their Z value (used in find_occupation)
    int natom_;
    double *zvals_;

    /// DOCC vector from input (if found)
    int doccpi_[8];
    bool input_docc_;

    /// SOCC vector from input (if found)
    int soccpi_[8];
    bool input_socc_;

    /// Number of alpha and beta electrons
    int nalpha_, nbeta_;
    /// Number of alpha and beta electrons per irrep
    int nalphapi_[8], nbetapi_[8];

    /// Mapping arrays
    int *so2symblk_;
    int *so2index_;

    /// Pairs needed for PK supermatrix
    size_t pk_pairs_;
    size_t pk_size_;
    int *pk_symoffset_;

    /// Perturb the Hamiltonian?
    int perturb_h_;
    /// How big of a perturbation
    double lambda_;
    /// With what...
    enum perturb { nothing, dipole_x, dipole_y, dipole_z };
    perturb perturb_;

    /// Using direct integrals?
    int direct_integrals_;

    /// DF Storage Scheme
    enum df_storage { double_full, full, flip_B_core, flip_B_disk, k_incore, disk};
    df_storage df_storage_;

    //Density Fitting?
    bool ri_integrals_;
    int ri_nbf_; //Number of functions in the auxiliary basis
    int *ri_pair_nu_;
    int *ri_pair_mu_;

    //Local K? (only with DF)
    bool local_K_;
    
    /// do we need Coulomb?
    bool J_is_required_;
    /// do we need Exchange?
    bool K_is_required_;

    /// Three Index tensor for DF-SCF
    double **B_ia_P_;
    double **A_ia_P_;

    /// Fitting Tensor Inverse (J-matrix) for DF-SCF
    double **Jinv_;

    double schwarz_; //Current Schwarz magnitude (static for now)
    int ntri_; //Number of function pairs after schwarz sieve and subsequent sieves
    int ntri_naive_; //Number of function pairs after schwarz sieve
   
    
    /// save SCF Cartesian Grid
    bool save_grid_;

    /// DIIS manager for all SCF wavefunctions
    boost::shared_ptr<DIISManager> diis_manager_;

    /// How many vectors for DIIS
    int num_diis_vectors_;
    /// Are we even using DIIS?
    int diis_enabled_;

public:
    /// Exactly what their name says
    vector<SharedSimpleMatrix> Dipole_;
    vector<SharedSimpleMatrix> Quadrupole_;

    /// Nuclear contributions
    SimpleVector nuclear_dipole_contribution_;
    SimpleVector nuclear_quadrupole_contribution_;

    /// Formation of H is the same regardless of RHF, ROHF, UHF
    void form_H();

    /// Formation of S^+1/2 and S^-1/2 are the same
    void form_Shalf();


    /// Set the amount of information to print
    void set_print(const int n) {print_ = n;}

    /// The number of iterations needed to reach convergence
    int iterations_needed() {return iterationsNeeded_;}

    /// The RMS error in the density
    double rms_density_error() {return Drms_;}
protected:
    /// Common initializer
    void common_init();

    /// Compute multipole integrals
    void form_multipole_integrals();

    /// Figure out how to occupy the orbitals in the absence of DOCC and SOCC
    void find_occupation(Vector& evals);

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
    int nElec_;

    /// The charge of the system
    int charge_;

    /// The multiplicity of the system (specified as 2 Ms + 1)
    int multiplicity_;

    /// The number of iterations need to reach convergence
    int iterationsNeeded_;

    /// Whether to add in an external potential to the fock matrix
    bool addExternalPotential_;

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

    /** Perform full localization procedure */
    virtual void fully_localize_mos() {}

    /** Propagate previous localization guess a. la. R. Polly et. al. */
    virtual void propagate_local_mos(){}

    /** Form canonical three-index DF tensor */
    void form_B();
    /** Form B without metric transform for local K (makes J go crazy fast)*/
    void form_A();
        
    /** Write tensor from memory to disk */
    void write_B();
    /** Free all memory associated with DF */
    void free_B();
    /** Free all memory associated with DF */
    void free_A();

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

    virtual ~HF();
};

}}

#endif
