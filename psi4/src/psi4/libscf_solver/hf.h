/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef HF_H
#define HF_H

#include <vector>
#include "psi4/libpsio/psio.hpp"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/vector.h"
#include "psi4/libdiis/diismanager.h"
#include "psi4/libdiis/diisentry.h"
#include "psi4/psi4-dec.h"
#include "psi4/libqt/qt.h"

namespace psi {
class Matrix;
class Vector;
class SimpleVector;
class TwoBodySOInt;
class JK;
class MinimalInterface;
class SOSCF;
class PCM;
class SuperFunctional;
class VBase;
namespace scf {

class HF : public Wavefunction {
protected:

    /// The kinetic energy matrix
    SharedMatrix T_;
    /// The 1e potential energy matrix
    SharedMatrix V_;
    /// The DFT potential matrices (nice naming scheme)
    SharedMatrix Va_;
    SharedMatrix Vb_;
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

    /// User defined orbitals
    SharedMatrix guess_Ca_;
    SharedMatrix guess_Cb_;

    /// Energy convergence threshold
    double energy_threshold_;

    /// Density convergence threshold
    double density_threshold_;

    /// Previous iteration's energy and current energy
    double Eold_;
    double E_;

    /// Table of energy components
    std::map<std::string, double> energies_;

    /// Basis list for SAD
    std::vector<std::shared_ptr<BasisSet>> sad_basissets_;
    std::vector<std::shared_ptr<BasisSet>> sad_fitting_basissets_;

    /// The RMS error in the density
    double Drms_;

    /// Max number of iterations for HF
    int maxiter_;

    /// Fail if we don't converge by maxiter?
    bool ref_C_;

    /// Fail if we don't converge by maxiter?
    bool fail_on_maxiter_;

    /// Current Iteration
    int iteration_;

    /// Did the SCF converge?

    bool converged_;

    /// Nuclear repulsion energy
    double nuclearrep_;

    /// Whether DIIS was performed this iteration, or not
    bool diis_performed_;

    /// DOCC vector from input (if found)
    bool input_docc_;

    /// SOCC vector from input (if found)
    bool input_socc_;

    /// Whether its broken symmetry solution or not
    bool broken_symmetry_;

    //Initial SAD doubly occupied may be more than ndocc
    // int sad_nocc_[8];
    Dimension original_doccpi_;
    Dimension original_soccpi_;
    int original_nalpha_;
    int original_nbeta_;
    bool reset_occ_;

    /// Mapping arrays
    int *so2symblk_;
    int *so2index_;

    /// SCF algorithm type
    std::string scf_type_;

    /// Old SCF type for DF guess trick
    /// TODO We should really get rid of that and put it in the driver
    std::string old_scf_type_;

    /// Perturb the Hamiltonian?
    int perturb_h_;
    /// How big of a perturbation
    Vector3 perturb_dipoles_;
    /// With what...
    enum perturb { nothing, dipole_x, dipole_y, dipole_z, dipole, embpot, dx, sphere };
    perturb perturb_;

    /// The value below which integrals are neglected
    double integral_threshold_;

    /// The soon to be ubiquitous JK object
    std::shared_ptr<JK> jk_;

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
    std::shared_ptr<DIISManager> diis_manager_;

    /// How many min vectors for DIIS
    int min_diis_vectors_;
    /// How many max vectors for DIIS
    int max_diis_vectors_;
    /// When do we start collecting vectors for DIIS
    int diis_start_;
    /// Are we even using DIIS?
    int diis_enabled_;

    /// Are we doing second-order convergence acceleration?
    bool soscf_enabled_;
    /// What is the gradient threshold that we should start?
    double soscf_r_start_;
    /// Maximum number of iterations
    int soscf_min_iter_;
    /// Minimum number of iterations
    int soscf_max_iter_;
    /// Break if the residual RMS is less than this
    double soscf_conv_;
    /// Do we print the microiterations?
    double soscf_print_;

    /// The amount (%) of the previous orbitals to mix in during SCF damping
    double damping_percentage_;
    /// The energy convergence at which SCF damping is disabled
    double damping_convergence_;
    /// Whether to use SCF damping
    bool damping_enabled_;
    /// Whether damping was actually performed this iteration
    bool damping_performed_;

    // parameters for hard-sphere potentials
    double radius_; // radius of spherical potential
    double thickness_; // thickness of spherical barrier
    int r_points_; // number of radial integration points
    int theta_points_; // number of colatitude integration points
    int phi_points_; // number of azimuthal integration points

    /// DFT variables
    std::shared_ptr<SuperFunctional> functional_;
    std::shared_ptr<VBase> potential_;

public:
    /// Nuclear contributions
    Vector nuclear_dipole_contribution_;
    Vector nuclear_quadrupole_contribution_;

    /// The number of iterations needed to reach convergence
    int iterations_needed() {return iterations_needed_;}

    /// The JK object (or null if it has been deleted)
    std::shared_ptr<JK> jk() const { return jk_; }

    /// The DFT Functional object (or null if it has been deleted)
    std::shared_ptr<SuperFunctional> functional() const { return functional_; }

    /// The DFT Potential object (or null if it has been deleted)
    std::shared_ptr<VBase> V_potential() const { return potential_; }

    /// The RMS error in the density
    double rms_density_error() {return Drms_;}

    /// Returns the occupation vectors
    std::shared_ptr<Vector> occupation_a() const;
    std::shared_ptr<Vector> occupation_b() const;

    // PCM interface
    bool pcm_enabled_;
    std::shared_ptr<PCM> hf_pcm_;

protected:

    /// Formation of H is the same regardless of RHF, ROHF, UHF
    // Temporarily converting to virtual function for testing embedding
    // potentials.  TDC, 5/23/12.
    virtual void form_H();

    /// Formation of S^+1/2 and S^-1/2 are the same
    void form_Shalf();

    /// Edit matrices if we are doing canonical orthogonalization
    virtual void prepare_canonical_orthogonalization() { return; }

    /// Prints the orbital occupation
    void print_occupation();

    /// Common initializer
    void common_init();

    /// Figure out how to occupy the orbitals in the absence of DOCC and SOCC
    void find_occupation();

    /// Maximum overlap method for prevention of oscillation/excited state SCF
    void MOM();
    /// Start the MOM algorithm (requires one iteration worth of setup)
    void MOM_start();

    /// Fractional occupation UHF/UKS
    void frac();
    /// Renormalize orbitals to 1.0 before saving
    void frac_renormalize();

    /// Check the stability of the wavefunction, and correct (if requested)
    virtual bool stability_analysis();
    void print_stability_analysis(std::vector<std::pair<double, int> > &vec);


    /// Determine how many core and virtual orbitals to freeze
    void compute_fcpi();
    void compute_fvpi();

    /// Prints the orbitals energies and symmetries (helper method)
    void print_orbitals(const char* header, std::vector<std::pair<double,
                        std::pair<const char*, int> > > orbs);

    /// Prints the orbitals in arbitrary order (works with MOM)
    void print_orbitals();

    /// Prints the energy breakdown from this SCF
    void print_energies();

    /// Prints some opening information
    void print_header();

    /// Prints some details about nsopi/nmopi, and initial occupations
    void print_preiterations();

    /// Do any needed integral setup
    virtual void integrals();

    /// Which set of iterations we're on in this computation, e.g., for stability
    /// analysis, where we want to retry SCF without going through all of the setup
    int attempt_number_;
    /// Maximum number of macroiterations to take in e.g. a stability analysis
    int max_attempts_;

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
    void reset_occupation();

    /// Form the guess (gaurantees C, D, and E)
    virtual void guess();

    /** Applies damping to the density update */
    virtual void damp_update();

    /** Applies second-order convergence acceleration */
    virtual int soscf_update();

    /** Rotates orbitals inplace C' = exp(U) C, U = antisymmetric matrix from x */
    void rotate_orbitals(SharedMatrix C, const SharedMatrix x);

    /** Transformation, diagonalization, and backtransform of Fock matrix */
    virtual void diagonalize_F(const SharedMatrix& F, SharedMatrix& C, std::shared_ptr<Vector>& eps);

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

    /** Compute the orbital gradient */
    virtual void compute_orbital_gradient(bool) {}

    /** Performs DIIS extrapolation */
    virtual bool diis() { return false; }

    /** Form Fia (for DIIS) **/
    virtual SharedMatrix form_Fia(SharedMatrix Fso, SharedMatrix Cso, int* noccpi);

    /** Form X'(FDS - SDF)X (for DIIS) **/
    virtual SharedMatrix form_FDSmSDF(SharedMatrix Fso, SharedMatrix Dso);

    /** Performs any operations required for a incoming guess **/
    virtual void format_guess();

    /** Save orbitals to use later as a guess **/
    // virtual void save_orbitals();

    /** Tells whether or not to read Fock matrix as a guess **/
    // bool do_use_fock_guess();

    /** Load fock matrix from previous computation to form guess MO coefficients **/
    // virtual void load_fock();

    /** Load orbitals from previous computation, projecting if needed **/
    // virtual void load_orbitals();

public:
    HF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> funct,
       Options& options, std::shared_ptr<PSIO> psio);

    virtual ~HF();

    /// Specialized initialization, compute integrals and does everything to prepare for iterations
    virtual void initialize();

    /// Performs the actual SCF iterations
    virtual void iterations();

    /// Performs stability analysis and calls back SCF with new guess if needed,
    /// Returns the SCF energy
    ///  This function should be called once orbitals are ready for energy/property computations,
    /// usually after iterations() is called.
    virtual double finalize_E();

    /// Base class Wavefunction requires this function. Here it is simply a wrapper around
    /// initialize(), iterations(), finalize_E(). It returns the SCF energy computed by
    /// finalize_E()
    virtual double compute_energy();

    /// Clears memory and closes files (Should they be open) prior to correlated code execution
    /// Derived classes override it for additional operations and then call HF::finalize()
    virtual void finalize();

    /// Semicanonicalizes ROHF/CUHF orbitals, breaking the alpha-beta degeneracy
    /// On entrance, there's only one set of orbitals and orbital energies.  On
    /// exit, the alpha and beta Fock matrices correspond to those in the semicanonical
    /// basis, and there are distinct alpha and beta C and epsilons, also in the
    /// semicanonical basis.
    virtual void semicanonicalize();

    /// Compute the MO coefficients (C_)
    virtual void form_C();

    /// Computes the density matrix (D_)
    virtual void form_D();

    /// Computes the density matrix (V_)
    virtual void form_V();

    // Return the DFT potenitals
    SharedMatrix Va() { return Va_; }
    SharedMatrix Vb() { return Vb_; }

    // Set guess occupied orbitals, nalpha and nbeta will be taken from the number of passed in eigenvectors
    void guess_Ca(SharedMatrix Ca) { guess_Ca_ = Ca; }
    void guess_Cb(SharedMatrix Cb) { guess_Cb_ = Cb; }

    // Expert option to reset the occuption or not at iteration zero
    void reset_occ(bool reset) { reset_occ_ = reset; }

    void set_sad_basissets(std::vector<std::shared_ptr<BasisSet>> basis_vec) { sad_basissets_ = basis_vec; }
    void set_sad_fitting_basissets(std::vector<std::shared_ptr<BasisSet>> basis_vec) { sad_fitting_basissets_ = basis_vec; }
};

}} // Namespaces

#endif
