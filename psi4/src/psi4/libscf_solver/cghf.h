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

#ifndef psi4_libscf_solver_cghf_h
#define psi4_libscf_solver_cghf_h

#include "psi4/libmints/complexwavefunction.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libscf_solver/basehf.h"
#include "psi4/libpsi4util/exception.h"

#include <tuple>
#include <vector>

namespace psi {
class BasisSet;
class SuperFunctional;
class VBase;
class ComplexJK;

namespace scf {

#ifndef USING_Einsums

/*! Stub CGHF: constructors throw unless Psi4 is built with Einsums. */
class CGHF : public BaseHF, public ComplexWavefunction {
   public:
    CGHF(std::shared_ptr<ComplexWavefunction> /*ref_wfn*/, std::shared_ptr<SuperFunctional> functional)
        : ComplexWavefunction(), BaseHF(functional) {
        throw PSIEXCEPTION("Psi4 not built with Einsums. CGHF not available. Recompile with -DENABLE_Einsums=ON.");
    }
    CGHF(std::shared_ptr<ComplexWavefunction> /*ref_wfn*/, std::shared_ptr<SuperFunctional> functional, Options& options,
         std::shared_ptr<PSIO> /*psio*/)
        : ComplexWavefunction(options), BaseHF(functional) {
        throw PSIEXCEPTION("Psi4 not built with Einsums. CGHF not available. Recompile with -DENABLE_Einsums=ON.");
    }
    ~CGHF() override = default;

    // Pure virtual stubs so the type is concrete for pybind.
    void form_H() override {}
    void form_Shalf() override {}
    void form_C(double /*shift*/ = 0.0) override {}
    void form_D() override {}
    void form_V() override {}
    void form_F() override {}
    void form_G() override {}
    void print_orbitals() override {}
};

#else

class CGHF : public BaseHF, public ComplexWavefunction {
   protected:
    // Prefer BaseHF's copies over BaseWavefunction's for these names.
    using BaseHF::nelectron_;
    using BaseHF::multiplicity_;

    /// Nothing to do with this object as we do not support DFT yet
    std::shared_ptr<VBase> potential_;

    /// Basis list for SAD
    std::vector<std::shared_ptr<BasisSet>> sad_basissets_;
    std::vector<std::shared_ptr<BasisSet>> sad_fitting_basissets_;

    /// JK object
    std::shared_ptr<ComplexJK> jk_;

    /// The kinetic energy matrix
    SharedComplexMatrix T_;
    /// The 1e potential energy matrix
    SharedComplexMatrix V_;
    /// The orthogonalization matrix (not the Lagrangian)
    SharedComplexMatrix X_;

    /// Temporary matrices
    SharedComplexMatrix G_;
    SharedComplexMatrix J_;
    SharedComplexMatrix K_;

    /// Occupied (or full) spinor MO guess from a previous computation / READ
    SharedComplexMatrix guess_C_;

    void common_init() override;

    /// Fills real SADGuess into D_/C_
    void compute_SAD_guess();

    // Dimension of irreps. Each irrep (h) size will be 2*nsopi_[h].
    Dimension irrep_sizes_;
   public:
    CGHF(std::shared_ptr<ComplexWavefunction> ref_wfn, std::shared_ptr<SuperFunctional> functional);
    CGHF(std::shared_ptr<ComplexWavefunction> ref_wfn, std::shared_ptr<SuperFunctional> functional, Options& options,
         std::shared_ptr<PSIO> psio);
    ~CGHF() override;

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
    void form_V() override { throw PSIEXCEPTION("CGHF::form_V() DFT for CGHF not available."); }

    /// Form the Fock matrix
    void form_F() override;

    /// Form the initial Fock matrix (default calls form_F)
    void form_initial_F() override { form_F(); }

    /// Form the G (J-K) matrix
    void form_G() override;

    /// Runs the SCF using OpenOrbitalOptimizer
    void openorbital_scf();

    /// Form the guess (guarantees C, D, and E)
    void guess();

    /// Supply occupied (or full) spinor MOs for READ / cast-up style guesses
    void guess_C(SharedComplexMatrix C) { guess_C_ = C; }

    /// Sets nelecpi_ from orbital energies (Aufbau)
    void find_occupation();

    /// SAD in Python-side requires this method to be defined
    /// nelectron_ is not mutated by SAD (unlike HF nalpha_/nbeta_)
    void reset_occupation() {}

    /// Compute 1e energy + nuc
    double compute_initial_E() override;

    /// Compute energy for the iteration.
    double compute_E() override;

    /// SAD information
    void set_sad_basissets(std::vector<std::shared_ptr<BasisSet>> basis_vec) { sad_basissets_ = basis_vec; }
    void set_sad_fitting_basissets(std::vector<std::shared_ptr<BasisSet>> basis_vec) {
        sad_fitting_basissets_ = basis_vec;
    }

    /// Clean up things after SCF
    void finalize();

    void print_orbitals() override;

    /// The JK object (or null if it has been deleted)
    std::shared_ptr<ComplexJK> jk() const { return jk_; }

    /// Sets the internal JK object (expert)
    void set_jk(std::shared_ptr<ComplexJK> jk);

    /// Builds the correct JK object (placeholder for future ComplexJK)
    std::shared_ptr<ComplexJK> build_jk(size_t memory) const;

    /// Construct the DFT potential.
    void setup_potential();

    /// The DFT Potential object (or nullptr if it has been deleted)
    std::shared_ptr<VBase> V_potential() const { return potential_; };

    /// Accessors for CGHF-only matrices (H/S/C/D/F live on ComplexWavefunction).
    SharedComplexMatrix get_T() const { return T_; }
    SharedComplexMatrix get_V() const { return V_; }
    SharedComplexMatrix get_X() const { return X_; }
    SharedComplexMatrix get_G() const { return G_; }
    SharedComplexMatrix get_J() const { return J_; }
    SharedComplexMatrix get_K() const { return K_; }

    /// Compute ⟨S²⟩ and multiplicity (2S+1) for the complex GHF determinant.
    /// Returns {spin_square, multiplicity}.
    std::tuple<double, double> spin_square() const;

    /// Fix MO column phases in place: first |C_μp| > 1e-3 becomes real and positive.
    void check_phases();

    /// AO DIIS residual FDS − SDF, transformed to the X-orthogonal basis (X^H e X).
    SharedComplexMatrix form_FDSmSDF(SharedComplexMatrix Fso, SharedComplexMatrix Dso);
};

#endif  // USING_Einsums

}  // namespace scf
}  // namespace psi

#endif
