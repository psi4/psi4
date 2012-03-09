#ifndef HF_H 
#define HF_H
/*
 *  hf.h
 *  matrix
 *
 *  Created by Justin Turney on 4/9/08.
 *  Copyright 2008 by Justin M. Turney, Ph.D.. All rights reserved.
 *
 */


#include <libpsio/psio.hpp>
#include <libmints/wavefunction.h>
#include <libmints/basisset.h>
#include <libmints/vector.h>
#include <libdiis/diismanager.h>
#include <libdiis/diisentry.h>
#include <psi4-dec.h>

namespace boost {
template<class T> class shared_ptr;
}

namespace psi {
class Matrix;
class Vector;
class SimpleVector;
class TwoBodySOInt;
class JK;
namespace scf {

class PKIntegrals;

class HF : public Wavefunction {
protected:

    /// The kinetic energy matrix
    SharedMatrix T_;
    /// The 1e potential energy matrix
    SharedMatrix V_;
    /// The core hamiltonian
    SharedMatrix H_;
    /// The orthogonalization matrix (symmetric or canonical)
    SharedMatrix X_;
    /// Temporary matrix for diagonalize_F
    SharedMatrix diag_temp_;
    /// Temporary matrix for diagonalize_F
    SharedMatrix diag_F_temp_;
    /// Temporary matrix for diagonalize_F
    SharedMatrix diag_C_temp_;

    /// Old C Alpha matrix (if needed for MOM)
    SharedMatrix Ca_old_;
    /// Old C Beta matrix (if needed for MOM)
    SharedMatrix Cb_old_;

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

    /// Whether DIIS was performed this iteration, or not
    bool diis_performed_;

    /// DOCC vector from input (if found)
    bool input_docc_;

    /// SOCC vector from input (if found)
    bool input_socc_;

    //Initial SAD doubly occupied may be more than ndocc
    int sad_nocc_[8];

    /// Mapping arrays
    int *so2symblk_;
    int *so2index_;

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

    /// The soon to be ubiquitous JK object
    boost::shared_ptr<JK> jk_;

    /// The SO integral generator.  Only ever constructed if needed
    boost::shared_ptr<TwoBodySOInt> eri_;
    /// PK Matrix approach
    boost::shared_ptr<PKIntegrals> pk_integrals_;

    /// Are we to do MOM?
    bool MOM_enabled_;
    /// Are we to do excited-state MOM?
    bool MOM_excited_;
    /// MOM started?
    bool MOM_started_;
    /// MOM performed?
    bool MOM_performed_;

    /// Are we to fractionally occupy?
    bool frac_enabled_;
    /// Frac started? (Same thing as frac_performed_)
    bool frac_performed_;

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

    /// The amount (%) of the previous orbitals to mix in during SCF damping
    double damping_percentage_;
    /// The energy convergence at which SCF damping is disabled
    double damping_convergence_;
    /// Whether to use SCF damping
    bool damping_enabled_;
    /// Whether damping was actually performed this iteration
    bool damping_performed_;
public:
    /// Nuclear contributions
    Vector nuclear_dipole_contribution_;
    Vector nuclear_quadrupole_contribution_;

    /// The number of iterations needed to reach convergence
    int iterations_needed() {return iterations_needed_;}

    /// The JK object (or null if it has been deleted)
    boost::shared_ptr<JK> jk() const { return jk_; }

    /// The RMS error in the density
    double rms_density_error() {return Drms_;}
protected:

    /// Formation of H is the same regardless of RHF, ROHF, UHF
    void form_H();

    /// Formation of S^+1/2 and S^-1/2 are the same
    void form_Shalf();

    /// Prints the orbital occupation
    void print_occupation();

    /// Perform casting of basis set if desired.
    SharedMatrix dualBasisProjection(SharedMatrix Cold, int* napi, boost::shared_ptr<BasisSet> old_basis, boost::shared_ptr<BasisSet> new_basis);

    /// Common initializer
    void common_init();

    /// Clears memory and closes files (Should they be open) prior to correlated code execution
    virtual void finalize();

    /// Figure out how to occupy the orbitals in the absence of DOCC and SOCC
    void find_occupation();

    /// Maximum overlap method for prevention of oscillation/excited state SCF
    void MOM();
    /// Start the MOM algorithm (requires one iteration worth of setup)
    void MOM_start();

    /// Fractional occupation UHF/UKS
    void frac();
    /// Renormalize orbitals to 1.0 before saving to chkpt
    void frac_renormalize();

    /// Determine how many core and virtual orbitals to freeze
    void compute_fcpi();
    void compute_fvpi();

    /// Prints the orbitals energies and symmetries (helper method)
    void print_orbitals(const char* header, std::vector<std::pair<double,
                        std::pair<const char*, int> > > orbs);

    /// Prints the orbitals in arbitrary order (works with MOM)
    void print_orbitals();

    /// Prints some opening information
    void print_header();

    /// Prints some details about nsopi/nmopi, and initial occupations
    void print_preiterations();

    /// Do any needed integral setup
    virtual void integrals();

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

    /// Save the current density and energy.
    virtual void save_density_and_energy() = 0;

    /// Check MO phases
    void check_phases();

    /// SAD Guess and propagation
    void compute_SAD_guess();

    /// Reset to regular occupation from the fractional occupation
    void reset_SAD_occupation();

    /// Form the guess (gaurantees C, D, and E)
    virtual void guess();

    /** Computes the density matrix (D_) */
    virtual void form_D() =0;

    /** Applies damping to the density update */
    virtual void damp_update();

    /** Compute the MO coefficients (C_) */
    virtual void form_C() =0;

    /** Transformation, diagonalization, and backtransform of Fock matrix */
    virtual void diagonalize_F(const SharedMatrix& F, SharedMatrix& C, boost::shared_ptr<Vector>& eps);

    /** Computes the Fock matrix */
    virtual void form_F() =0;

    /** Computes the initial MO coefficients (default is to call form_C) */
    virtual void form_initial_C() { form_C(); }

    /** Forms the G matrix */
    virtual void form_G() =0;

    /** Computes the initial energy. */
    virtual double compute_initial_E() { return 0.0; }

    /** Test convergence of the wavefunction */
    virtual bool test_convergency() { return false; }

    /** Compute/print spin contamination information (if unrestricted) **/
    virtual void compute_spin_contamination();

    /** Saves information to the checkpoint file */
    virtual void save_information() {}

    /** Handles forming the PK matrix */
    virtual void form_PK() {}

    /** Save the Fock matrix to the DIIS object */
    virtual void save_fock() {}

    /** Performs DIIS extrapolation */
    virtual bool diis() { return false; }

    /** Form Fia (for DIIS) **/
    virtual SharedMatrix form_Fia(SharedMatrix Fso, SharedMatrix Cso, int* noccpi);

    /** Form X'(FDS - SDF)X (for DIIS) **/
    virtual SharedMatrix form_FDSmSDF(SharedMatrix Fso, SharedMatrix Dso);

    /** Save orbitals to File 100 **/
    virtual void save_orbitals();

    /** Load orbitals from File 100, projecting if needed **/
    virtual void load_orbitals();

    /** Save SAPT info (TODO: Move to Python driver **/
    virtual void save_sapt_info() {}

    /** Saves all wavefunction information to the checkpoint file*/
    void dump_to_checkpoint();

    /** Computes the J and/or K matrices according to the scf_type keyword and the functor passed in*/
    template <class JKFunctor> void process_tei(JKFunctor & functor);

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
    HF(Options& options, boost::shared_ptr<PSIO> psio, boost::shared_ptr<Chkpt> chkpt);
    HF(Options& options, boost::shared_ptr<PSIO> psio);

    virtual ~HF();

    virtual double compute_energy();
};

}} // Namespaces

#endif
