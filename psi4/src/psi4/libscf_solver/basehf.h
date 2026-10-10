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

#ifndef psi4_libscf_solver_basehf_h
#define psi4_libscf_solver_basehf_h

#include "psi4/pybind11.h"

#include <array>
#include <map>
#include <memory>
#include <string>

namespace psi {

class Molecule;
class Options;
class SuperFunctional;
class BaseJK;

namespace scf {

/*! \brief Mixin holding SCF-loop state shared by real (HF) and complex (CGHF) solvers.
 *
 *  HF / CGHF inherit this alongside Wavefunction / ComplexWavefunction so common
 *  iteration, DIIS, MOM, and DFT-functional bookkeeping is not duplicated.
 *
 *  JK storage lives on the derived class (typed JK vs BaseJK/ComplexJK); BaseHF
 *  only exposes a polymorphic accessor interface.
 */
class BaseHF {
   protected:
    /// Current Iteration
    int iteration_;

    /// Did the SCF converge?
    bool converged_;

    /// Nuclear repulsion energy
    double nuclearrep_;

    /// Reset occupations in SCF iteration?
    bool reset_occ_;
    /// SAD guess, non-idempotent guess density?
    bool sad_;

    /// SCF algorithm type
    std::string scf_type_;

    /// Table of energy components
    std::map<std::string, double> energies_;

    /// Are we to do MOM?
    bool MOM_enabled_;
    /// Are we to do excited-state MOM?
    bool MOM_excited_;
    /// MOM performed?
    bool MOM_performed_;

    /// DIIS manager intiialized?
    bool initialized_diis_manager_;
    /// DIIS manager for all SCF wavefunctions
    py::object diis_manager_;

    /// When do we start collecting vectors for DIIS
    int diis_start_;
    /// Are we even using DIIS?
    int diis_enabled_;

    /// Which set of iterations we're on in this computation, e.g., for stability
    /// analysis, where we want to retry SCF without going through all of the setup
    int attempt_number_;

    /// The number of electrons
    int nelectron_;

    /// The charge of the system
    int charge_;

    /// The multiplicity of the system (specified as 2 Ms + 1)
    int multiplicity_;

    /// The value below which integrals are neglected
    double integral_threshold_;

    /// DFT variables
    std::shared_ptr<SuperFunctional> functional_;

    BaseHF() = default;
    explicit BaseHF(std::shared_ptr<SuperFunctional> functional) : functional_(std::move(functional)) {}

    /// Shared SCF bookkeeping
    /// `module` is BaseWavefunction::module_
    void common_init(Options& options, std::string& module, const std::shared_ptr<Molecule>& molecule,
                     const std::array<double, 3>& dipole_field_strength);

   public:
    virtual ~BaseHF() = default;

    /// Get and set current iteration
    int iteration() const { return iteration_; }
    void set_iteration(int iter) { iteration_ = iter; }

    /// Are we even using DIIS?
    bool diis_enabled() const { return bool(diis_enabled_); }
    void set_diis_enabled(bool tf) { diis_enabled_ = int(tf); }

    /// When do we start collecting vectors for DIIS
    int diis_start() const { return diis_start_; }
    void set_diis_start(int iter) { diis_start_ = iter; }

    /// Are we to do excited-state MOM?
    bool MOM_excited() const { return MOM_excited_; }
    void set_MOM_excited(bool tf) { MOM_excited_ = tf; }

    /// MOM performed?
    bool MOM_performed() const { return MOM_performed_; }
    void set_MOM_performed(bool tf) { MOM_performed_ = tf; }

    /// Which set of iterations we're on in this computation, e.g., for stability
    /// analysis, where we want to retry SCF without going through all of the setup
    int attempt_number() const { return attempt_number_; }
    void set_attempt_number(int an) { attempt_number_ = an; }

    const std::string& scf_type() const { return scf_type_; }

    /// The DIIS object
    // std::shared_ptr<py::object> is probably saner, but that hits a compile error.
    // Quite probably https://github.com/pybind/pybind11/issues/787
    py::object& diis_manager() { return diis_manager_; }
    void set_diis_manager(py::object& manager) { diis_manager_ = manager; }
    bool initialized_diis_manager() const { return initialized_diis_manager_; }
    void set_initialized_diis_manager(bool tf) { initialized_diis_manager_ = tf; }

    /// The DFT Functional object (or null if it has been deleted)
    std::shared_ptr<SuperFunctional> functional() const { return functional_; }

    // Expert option to reset the occupation or not at iteration zero
    bool reset_occ() const { return reset_occ_; }
    void set_reset_occ(bool reset) { reset_occ_ = reset; }
    // Expert option to toggle non-idempotent density matrix or not at iteration zero
    bool sad() const { return sad_; }
    void set_sad(bool sad) { sad_ = sad; }

    // Energies data
    void set_energies(std::string key, double value) { energies_[key] = value; }
    double get_energies(std::string key) { return energies_[key]; }

    /// Form core Hamiltonian
    virtual void form_H() = 0;

    /// Form the S^{1/2} orthogonalization matrix
    virtual void form_Shalf() = 0;

    /// Compute MO coefficients from the current Fock matrix
    virtual void form_C(double shift = 0.0) = 0;

    /// Compute initial MO coefficients (default calls form_C)
    virtual void form_initial_C() { form_C(); };

    /** Computes the initial energy. */
    virtual double compute_initial_E() { return 0.0; }

    /// Compute energy for the iteration.
    virtual double compute_E() { return 0.0; }

    /// Form the density matrix from the current orbitals
    virtual void form_D() = 0;

    /// Form the Kohn-Sham potential matrix from the current density
    virtual void form_V() = 0;

    /// Form the Fock matrix
    virtual void form_F() = 0;

    /// Form the initial Fock matrix (default calls form_F)
    virtual void form_initial_F() { form_F(); };

    /// Form the G (J-K) matrix
    virtual void form_G() = 0;

    virtual void print_orbitals() = 0;
};

}  // namespace scf
}  // namespace psi

#endif
