/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#ifndef HF_H
#define HF_H

#include <vector>
#include <functional>
#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/vector3.h"
#include "psi4/psi4-dec.h"

#include "psi4/pybind11.h"

namespace psi {
using PerturbedPotentialFunction = std::function<SharedMatrix(SharedMatrix)>;
using PerturbedPotentials = std::map<std::string, PerturbedPotentialFunction>;
class Vector;
class JK;
class PCM;
class SuperFunctional;
class VBase;
class BasisSet;
class DIISManager;
class PSIO;
namespace scf {

class HF : public Wavefunction {
   protected:
    /// The kinetic energy matrix
    SharedMatrix T_;
    /// The 1e potential energy matrix
    SharedMatrix V_;
    /// A temporary spot for the H matrix
    SharedMatrix Horig_;
    /// The DFT potential matrices (nice naming scheme)
    SharedMatrix Va_;
    SharedMatrix Vb_;
    /// The orthogonalization matrix (symmetric or canonical)
    SharedMatrix X_;
    /// List of external potentials to add to Fock matrix and updated at every iteration
    /// e.g. PCM potential
    std::vector<SharedMatrix> external_potentials_;

    /// Map of external potentials/perturbations to add to the CPSCF two-electron contribution
    /// e.g. PCM or PE potential
    PerturbedPotentials external_cpscf_perturbations_;

    /// Old C Alpha matrix (if needed for MOM)
    SharedMatrix Ca_old_;
    /// Old C Beta matrix (if needed for MOM)
    SharedMatrix Cb_old_;

    /// User defined orbitals
    SharedMatrix guess_Ca_;
    SharedMatrix guess_Cb_;

    // Q: right now, thresholds are removed from Wfn since only appear once, py-side.
    //    should we instead store here the E & D to which SCF was converged?

    /// Table of energy components
    std::map<std::string, double> energies_;

    /// Basis list for SAD
    std::vector<std::shared_ptr<BasisSet>> sad_basissets_;
    std::vector<std::shared_ptr<BasisSet>> sad_fitting_basissets_;

    /// Current Iteration
    int iteration_;

    /// Did the SCF converge?
    bool converged_;

    /// Nuclear repulsion energy
    double nuclearrep_;

    /// DOCC vector from input (if found)
    bool input_docc_;

    /// SOCC vector from input (if found)
    bool input_socc_;

    /// Whether its broken symmetry solution or not
    bool broken_symmetry_;

    // Initial SAD doubly occupied may be more than ndocc
    Dimension original_nalphapi_;
    Dimension original_nbetapi_;
    int original_nalpha_;
    int original_nbeta_;
    // Reset occupations in SCF iteration?
    bool reset_occ_;
    // SAD guess, non-idempotent guess density?
    bool sad_;

    /// Mapping arrays
    int* so2symblk_;
    int* so2index_;

    /// SCF algorithm type
    std::string scf_type_;

    /// The value below which integrals are neglected
    double integral_threshold_;

    /// The soon to be ubiquitous JK object
    std::shared_ptr<JK> jk_;

    /// Are we to do MOM?
    bool MOM_enabled_;
    /// Are we to do excited-state MOM?
    bool MOM_excited_;
    /// MOM performed?
    bool MOM_performed_;

    /// Frac started? (Same thing as frac_performed_)
    bool frac_performed_;
    /// The orbitals _before_ scaling needed for Frac
    SharedMatrix unscaled_Ca_;
    SharedMatrix unscaled_Cb_;

    /// DIIS manager intiialized?
    bool initialized_diis_manager_;
    /// DIIS manager for all SCF wavefunctions
    py::object diis_manager_;

    /// When do we start collecting vectors for DIIS
    int diis_start_;
    /// Are we even using DIIS?
    int diis_enabled_;

    // parameters for hard-sphere potentials
    double radius_;     // radius of spherical potential
    double thickness_;  // thickness of spherical barrier
    int r_points_;      // number of radial integration points
    int theta_points_;  // number of colatitude integration points
    int phi_points_;    // number of azimuthal integration points

    /// DFT variables
    std::shared_ptr<SuperFunctional> functional_;

    // CPHF info
    int cphf_nfock_builds_;
    bool cphf_converged_;

    /// Edit matrices if we are doing canonical orthogonalization
    virtual void prepare_canonical_orthogonalization() { return; }

    /// Prints the orbital occupation
    void print_occupation();

    /// Common initializer
    void common_init();
    /// Part of the common initializer that runs after subclass specific tasks
    void subclass_init();
    /// Construct the DFT potential.
    virtual void setup_potential() { throw PSIEXCEPTION("setup_potential virtual"); };

    /// Maximum overlap method for prevention of oscillation/excited state SCF
    void MOM();
    /// Start the MOM algorithm (requires one iteration worth of setup)
    void MOM_start();
    /// Perform MOM operations for a single spincase
    void MOM_spincase(const Dimension& npi, Vector& orb_energies, Matrix& old_C, Matrix& new_C);

    /// Fractional occupation UHF/UKS
    void frac();

    /// Determine how many core and virtual orbitals to freeze
    void compute_fcpi();

    /// Prints the orbitals energies and symmetries (helper method)
    void print_orbital_pairs(const char* header, std::vector<std::pair<double, std::pair<std::string, int>>> orbs);

    /// Which set of iterations we're on in this computation, e.g., for stability
    /// analysis, where we want to retry SCF without going through all of the setup
    int attempt_number_;

    /// The number of electrons
    int nelectron_;

    /// The charge of the system
    int charge_;

    /// The multiplicity of the system (specified as 2 Ms + 1)
    int multiplicity_;

    /// SAD Guess and propagation
    virtual void compute_SAD_guess(bool natorb);
    /// Huckel guess
    virtual void compute_huckel_guess(bool updated_rule);
    /// Forms the SAPGAU guess
    virtual void compute_sapgau_guess();

    /** Transformation, diagonalization, and backtransform of Fock matrix */
    virtual void diagonalize_F(const SharedMatrix& F, SharedMatrix& C, std::shared_ptr<Vector>& eps);

    /** Form Fia (for DIIS) **/
    virtual SharedMatrix form_Fia(SharedMatrix Fso, SharedMatrix Cso, int* noccpi);

    /** Performs any operations required for a incoming guess **/
    virtual void format_guess();

   public:
    HF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> funct, Options& options,
       std::shared_ptr<PSIO> psio);

    ~HF() override;

    /// Get and set current iteration
    int iteration() const { return iteration_; }
    void set_iteration(int iter) { iteration_ = iter; }

    /// Are we even using DIIS?
    bool diis_enabled() const { return bool(diis_enabled_); }
    void set_diis_enabled(bool tf) { diis_enabled_ = int(tf); }

    /// When do we start collecting vectors for DIIS
    int diis_start() const { return diis_start_; }
    void set_diis_start(int iter) { diis_start_ = iter; }

    /// Frac performed current iteration?
    bool frac_performed() const { return frac_performed_; }
    void set_frac_performed(bool tf) { frac_performed_ = tf; }

    /// Runs the SCF using OpenOrbitalOptimizer
    virtual void openorbital_scf() { throw PSIEXCEPTION("openorbital_scf is virtual; it has not been implemented for your class"); };

    /// Are we to do excited-state MOM?
    bool MOM_excited() const { return MOM_excited_; }
    void set_MOM_excited(bool tf) { MOM_excited_ = tf; }

    /// MOM performed?
    bool MOM_performed() const { return MOM_performed_; }
    void set_MOM_performed(bool tf) { MOM_performed_ = tf; }

    // Q: MOM_started_ was ditched b/c same info as MOM_performed_

    /// Which set of iterations we're on in this computation, e.g., for stability
    /// analysis, where we want to retry SCF without going through all of the setup
    int attempt_number() const { return attempt_number_; }
    void set_attempt_number(int an) { attempt_number_ = an; }

    /// Check the stability of the wavefunction, and correct (if requested)
    /// For UHF, this is defined Python-side. The other methods should be joining it.
    virtual bool stability_analysis();

    /** Computes the initial energy. */
    virtual double compute_initial_E() { return 0.0; }

    const std::string& scf_type() const { return scf_type_; }

    /// Check MO phases
    void check_phases();

    /// Prints the orbitals in arbitrary order (works with MOM)
    void print_orbitals();

    /// Prints some opening information
    void print_header();

    /** Compute/print spin contamination information (if unrestricted) **/
    void print_stability_analysis(std::vector<std::pair<double, int>>& vec) const;

    virtual void compute_spin_contamination();

    /// The DIIS object
    // std::shared_ptr<py::object> is probably saner, but that hits a compile error.
    // Quite probably https://github.com/pybind/pybind11/issues/787
    py::object& diis_manager() { return diis_manager_; }
    void set_diis_manager(py::object& manager) { diis_manager_ = manager; }
    bool initialized_diis_manager() const { return initialized_diis_manager_; }
    void set_initialized_diis_manager(bool tf) { initialized_diis_manager_ = tf; }

    /// The JK object (or null if it has been deleted)
    std::shared_ptr<JK> jk() const { return jk_; }

    /// Sets the internal JK object (expert)
    void set_jk(std::shared_ptr<JK> jk);

    /// The DFT Functional object (or null if it has been deleted)
    std::shared_ptr<SuperFunctional> functional() const { return functional_; }

    /// The DFT Potential object (or null if it has been deleted)
    /// This needs to be virtual so that subclasses can enforce their
    /// particular potential's derived class.
    virtual std::shared_ptr<VBase> V_potential() const = 0;

    /// Returns the occupation vectors
    std::shared_ptr<Vector> occupation_a() const;
    std::shared_ptr<Vector> occupation_b() const;

    /// Save the current density and energy.
    virtual void save_density_and_energy();

    /// Reset to the user-specified DOCC/SOCC if any, and zero's otherwise.
    /// Fractional occupation requires this.
    void reset_occupation();

    /// Compute energy for the iteration.
    virtual double compute_E();

    /** Applies second-order convergence acceleration */
    virtual int soscf_update(double soscf_conv, int soscf_min_iter, int soscf_max_iter, int soscf_print);

    /// Figure out how to occupy the orbitals in the absence of DOCC and SOCC
    void find_occupation();

    /** Performs DIIS extrapolation */
    virtual bool diis(double dnorm) { return false; }

    /** Compute the orbital gradient */
    virtual double compute_orbital_gradient(bool save_diis, int max_diis_vectors) { return 0.0; }

    /** Applies damping to the density update */
    virtual void damping_update(double);

    /// Clears memory and closes files (Should they be open) prior to correlated code execution
    /// Derived classes override it for additional operations and then call HF::finalize()
    virtual void finalize();

    /// Semicanonicalizes ROHF/CUHF orbitals, breaking the alpha-beta degeneracy
    /// On entrance, there's only one set of orbitals and orbital energies.  On
    /// exit, the alpha and beta Fock matrices correspond to those in the semicanonical
    /// basis, and there are distinct alpha and beta C and epsilons, also in the
    /// semicanonical basis.
    virtual void semicanonicalize();

    /// Renormalize orbitals to 1.0 before saving
    void frac_renormalize();
    void frac_helper();

    /// Formation of H is the same regardless of RHF, ROHF, UHF
    // Temporarily converting to virtual function for testing embedding
    // potentials.  TDC, 5/23/12.
    virtual void form_H();

    /// Do any needed integral JK setup
    virtual void initialize_gtfock_jk();

    /// Formation of S^+1/2 and S^-1/2 are the same
    void form_Shalf();

    /// Form the guess (guarantees C, D, and E)
    virtual void guess();

    /// Compute the MO coefficients (C_) using level shift
    virtual void form_C(double shift = 0.0);
    /** Computes the initial MO coefficients (default is to call form_C) */
    virtual void form_initial_C() { form_C(); }

    /// Computes the density matrix (D_)
    virtual void form_D();

    /// Computes the density matrix (V_)
    virtual void form_V();

    /** Computes the Fock matrix */
    virtual void form_F();
    /** Computes the initial Fock matrix (default is to call form_F) */
    virtual void form_initial_F() { form_F(); }

    /** Forms the G matrix */
    virtual void form_G();

    /** Form X'(FDS - SDF)X (for DIIS) **/
    virtual SharedMatrix form_FDSmSDF(SharedMatrix Fso, SharedMatrix Dso);

    /** Rotates orbitals inplace C' = C exp(U), U = antisymmetric matrix from x */
    void rotate_orbitals(SharedMatrix C, const SharedMatrix x);

    /// Hessian-vector computers and solvers
    virtual std::vector<SharedMatrix> onel_Hx(std::vector<SharedMatrix> x);
    virtual std::vector<SharedMatrix> twoel_Hx(std::vector<SharedMatrix> x, bool combine = true,
                                               std::string return_basis = "MO");
    virtual std::vector<SharedMatrix> cphf_Hx(std::vector<SharedMatrix> x);
    virtual std::vector<SharedMatrix> cphf_solve(std::vector<SharedMatrix> x_vec, double conv_tol = 1.e-4,
                                                 int max_iter = 10, int print_lvl = 1);

    // CPHF data
    bool cphf_converged() { return cphf_converged_; }
    int cphf_nfock_builds() { return cphf_nfock_builds_; }

    // Return the DFT potenitals
    SharedMatrix Va() { return Va_; }
    SharedMatrix Vb() { return Vb_; }

    // Set guess occupied orbitals, nalpha and nbeta will be taken from the number of passed in eigenvectors
    void guess_Ca(SharedMatrix Ca) { guess_Ca_ = Ca; }
    void guess_Cb(SharedMatrix Cb) { guess_Cb_ = Cb; }

    // Expert option to reset the occuption or not at iteration zero
    bool reset_occ() const { return reset_occ_; }
    void set_reset_occ(bool reset) { reset_occ_ = reset; }
    // Expert option to toggle non-idempotent density matrix or not at iteration zero
    bool sad() const { return sad_; }
    void set_sad(bool sad) { sad_ = sad; }

    // SAD information
    void set_sad_basissets(std::vector<std::shared_ptr<BasisSet>> basis_vec) { sad_basissets_ = basis_vec; }
    void set_sad_fitting_basissets(std::vector<std::shared_ptr<BasisSet>> basis_vec) {
        sad_fitting_basissets_ = basis_vec;
    }

    // Energies data
    void set_energies(std::string key, double value) { energies_[key] = value; }
    double get_energies(std::string key) { return energies_[key]; }

    // External potentials
    void clear_external_potentials() { external_potentials_.clear(); }
    void push_back_external_potential(const SharedMatrix& Vext) { external_potentials_.push_back(Vext); }
    void set_external_cpscf_perturbation(const std::string name, PerturbedPotentialFunction fun) {
        external_cpscf_perturbations_[name] = fun;
    }
    void clear_external_cpscf_perturbations() { external_cpscf_perturbations_.clear(); }
    void compute_fvpi();
};
}  // namespace scf
}  // namespace psi

#endif
