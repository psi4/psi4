/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

#include "psi4/libscf_solver/basehf.h"

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

class HF : public BaseHF, public Wavefunction {
   protected:
    // Prefer BaseHF's copies over BaseWavefunction's for these names.
    using BaseHF::nelectron_;
    using BaseHF::multiplicity_;

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

    /// Basis list for SAD
    std::vector<std::shared_ptr<BasisSet>> sad_basissets_;
    std::vector<std::shared_ptr<BasisSet>> sad_fitting_basissets_;

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

    /// Mapping arrays
    int* so2symblk_;
    int* so2index_;

    /// The soon to be ubiquitous JK object
    std::shared_ptr<JK> jk_;

    /// Frac started? (Same thing as frac_performed_)
    bool frac_performed_;
    /// The orbitals _before_ scaling needed for Frac
    SharedMatrix unscaled_Ca_;
    SharedMatrix unscaled_Cb_;

    // parameters for hard-sphere potentials
    double radius_;     // radius of spherical potential
    double thickness_;  // thickness of spherical barrier
    int r_points_;      // number of radial integration points
    int theta_points_;  // number of colatitude integration points
    int phi_points_;    // number of azimuthal integration points

    // CPHF info
    int cphf_nfock_builds_;
    bool cphf_converged_;

    /// Edit matrices if we are doing canonical orthogonalization
    virtual void prepare_canonical_orthogonalization() { return; }

    /// Prints the orbital occupation
    void print_occupation();

    /// Common initializer
    void common_init() override;
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

    /// SAD Guess and propagation
    virtual void compute_SAD_guess(bool natorb);
    /// Huckel guess
    virtual void compute_huckel_guess(bool updated_rule);
    /// Forms the SAPGAU guess
    virtual void compute_sapgau_guess();

    /** Transformation, diagonalization, and backtransform of Fock matrix */
    virtual void diagonalize_F(const SharedMatrix& F, SharedMatrix& C, std::shared_ptr<Vector>& eps);

    /** Form Fia (for DIIS) **/
    virtual SharedMatrix form_Fia(SharedMatrix Fso, SharedMatrix Cso, const Dimension& noccpi);

    /** Performs any operations required for a incoming guess **/
    virtual void format_guess();

   public:
    HF(SharedWavefunction ref_wfn, std::shared_ptr<SuperFunctional> funct, Options& options,
       std::shared_ptr<PSIO> psio);

    ~HF() override;

    /// Frac performed current iteration?
    bool frac_performed() const { return frac_performed_; }
    void set_frac_performed(bool tf) { frac_performed_ = tf; }

    /// Runs the SCF using OpenOrbitalOptimizer
    virtual void openorbital_scf() { throw PSIEXCEPTION("openorbital_scf is virtual; it has not been implemented for your class"); };

    // Q: MOM_started_ was ditched b/c same info as MOM_performed_

    /// Check the stability of the wavefunction, and correct (if requested)
    /// For UHF, this is defined Python-side. The other methods should be joining it.
    virtual bool stability_analysis();

    /// Check MO phases
    void check_phases();

    /// Prints the orbitals in arbitrary order (works with MOM)
    void print_orbitals() override;

    /// Prints some opening information
    void print_header();

    /** Compute/print spin contamination information (if unrestricted) **/
    void print_stability_analysis(std::vector<std::pair<double, int>>& vec) const;

    virtual void compute_spin_contamination();

    /// The JK object (or null if it has been deleted)
    std::shared_ptr<JK> jk() const { return jk_; }

    /// Sets the internal JK object (expert)
    void set_jk(std::shared_ptr<JK> jk);

    /* Builds the correct JK object. This is a new function defined so that the
    * python side no longer constructs the JK object separately from wfn.
    * ComplexJKs should be built for ComplexWavefunctions. However, any mechanism
    * to determine if ComplexJK is requested within JK::build_JK is smelly.
    *   1. We could infer a ComplexJK with `REFERENCE CGHF`, but every new complex
    *      reference would require an addition to this list.
    *   2. We could add another option for `IS_COMPLEX`, but being a global option
    *      anyone could change it. (not obvious what the source of truth is)
    *   3. We could pass in the wfn directly to Build_JK and do RTTI. This works
    *      but is awfuly smelly too *insert foghorn sound*
    * Here, we build the correct JK via an overriden pure virtual function.
    * Notice the const qualifier: This doesn't change `HF::jk_`. A downside
    * to this approach is that shared ptrs are not covariant, meaning
    * a bit more work is required to upcast JK/ComplexJK --> BaseJK.
    */
    std::shared_ptr<JK> build_jk(size_t memory) const;

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
    double compute_E() override;

    /** Applies second-order convergence acceleration */
    virtual int soscf_update(double soscf_conv, int soscf_min_iter, int soscf_max_iter, bool soscf_print);

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

    /// Form core Hamiltonian
    void form_H() override;

    /// Form the S^{1/2} orthogonalization matrix
    void form_Shalf() override;

    /// Compute MO coefficients from the current Fock matrix
    void form_C(double shift = 0.0) override;

    /// Compute initial MO coefficients (default calls form_C)
    void form_initial_C() override { form_C(); }

    /// Form the density matrix from the current orbitals
    void form_D() override;

    /// Form the Kohn-Sham potential matrix from the current density
    void form_V() override;

    /// Form the Fock matrix
    void form_F() override;

    /// Form the initial Fock matrix (default calls form_F)
    void form_initial_F() override { form_F(); }

    /// Form the G (J-K) matrix
    void form_G() override;

    /// Do any needed integral JK setup
    virtual void initialize_gtfock_jk();

    /// Form the guess (guarantees C, D, and E)
    virtual void guess();

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

    // SAD information
    void set_sad_basissets(std::vector<std::shared_ptr<BasisSet>> basis_vec) { sad_basissets_ = basis_vec; }
    void set_sad_fitting_basissets(std::vector<std::shared_ptr<BasisSet>> basis_vec) {
        sad_fitting_basissets_ = basis_vec;
    }

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
