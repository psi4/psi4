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

#ifndef _psi_src_lib_libmints_complexwavefunction_h
#define _psi_src_lib_libmints_complexwavefunction_h

#include "psi4/libmints/basewavefunction.h"
#include "psi4/libmints/typedefs.h"
#include "psi4/libpsi4util/exception.h"

#include <complex>
#include <map>
#include <memory>
#include <string>
#include <vector>

#ifdef USING_Einsums
#include <Einsums/Config.hpp>
#include <Einsums/Tensor.hpp>

#include <Einsums/LinearAlgebra.hpp>
#include <Einsums/TensorAlgebra.hpp>
#include <Einsums/Runtime.hpp>
#endif

namespace psi {

class Molecule;
class BasisSet;
class Options;

#ifdef USING_Einsums

using ComplexMatrix = einsums::TiledTensor<std::complex<double>, 2>;
using SharedComplexMatrix = std::shared_ptr<ComplexMatrix>;

/*! \ingroup MINTS
 *  \class ComplexWavefunction
 *  \brief Home for wavefunctions built on complex molecular orbitals.
 */
class PSI_API ComplexWavefunction : public BaseWavefunction {
   protected:
    /// Total number of electrons
    int nelec_;

    /// Number of electrons per irrep
    std::vector<size_t> nelecpi_;

    /// Overlap matrix
    SharedComplexMatrix S_;

    /// Core Hamiltonian matrix
    SharedComplexMatrix H_;

    /// MO coefficients
    SharedComplexMatrix C_;

    /// Density matrix
    SharedComplexMatrix D_;
    /// Lagrangian matrix
    SharedComplexMatrix Lagrangian_;

    /// Fock matrix
    SharedComplexMatrix F_;

    /// Orbital energies
    SharedVector epsilon_;

    /// gradient, if available, as natom_ x 3 SharedComplexMatrix
    SharedComplexMatrix gradient_;

    /// Hessian, if available, as natom_*3 x natom_*3 SharedComplexMatrix (NOT mass-weighted!)
    SharedComplexMatrix hessian_;

    /// Should molecular orbital extents be available, they will be here
    std::vector<SharedVector> mo_extents_;

    // Collection of Matrix variables
    // * any '<mtd> GRADIENT' is an energy derivative w.r.t. nuclear perturbations (a.u.) as a (nat, 3) Matrix
    // * any '<mtd> DIPOLE GRADIENT' is a dipole derivative w.r.t. nuclear perturbations (a.u.) as a degree-of-freedom
    //   by dipole component (3 * nat, 3) Matrix
    std::map<std::string, SharedComplexMatrix> arrays_;

    void common_init() override;

   public:
    ComplexWavefunction();

    /// Constructor for an entirely new wavefunction with an existing basis
    ComplexWavefunction(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basis, Options& options);

    /// Constructor for an entirely new wavefunction with an existing basis and global options
    ComplexWavefunction(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basis);

    /// Blank constructor for derived classes
    ComplexWavefunction(Options& options);

    ~ComplexWavefunction() override;

    /// Compute energy. Subclasses override this function to compute its energy.
    virtual double compute_energy() {
        throw PSIEXCEPTION("Compute energy has not been implemented for this wavefunction.");
    }

    /// Compute gradient.  Subclasses override this function to compute the gradient.
    virtual SharedComplexMatrix compute_gradient() {
        throw PSIEXCEPTION("Analytic gradients are not available for this wavefunction.");
    }

    /// Compute Hessian.  Subclasses override this function to compute the Hessian.
    virtual SharedComplexMatrix compute_hessian() {
        throw PSIEXCEPTION("Analytic Hessians are not available for this wavefunction.");
    }

    /// Returns the number of electrons per irrep array.
    const std::vector<size_t>& nelecpi() const { return nelecpi_; }
    /// Return the number of electrons
    int nelec() const { return nelec_; }

    /// Returns the MO coefficients
    SharedComplexMatrix C() const;
    /// Set MO coefficients
    void set_C(SharedComplexMatrix C) { C_ = C; }
    /// Returns the (SO basis) Fock matrix
    SharedComplexMatrix F() const { return F_; }
    /// Set Fock matrix
    void set_F(SharedComplexMatrix F) { F_ = F; }
    /// Returns the orbital energies
    SharedVector epsilon() const { return epsilon_; }
    /// Set orbital energies
    void set_epsilon(SharedVector epsilon) { epsilon_ = epsilon; }

    /// Returns the OPDM for the wavefunction
    const SharedComplexMatrix D() const { return D_; }
    /// Set density matrix
    void set_D(SharedComplexMatrix D) { D_ = D; }

    /// Returns the SO basis Lagrangian
    SharedComplexMatrix lagrangian() const { return Lagrangian_; }
    /// Set Lagrangian matrix in SO basis
    void set_lagrangian(SharedComplexMatrix X) { Lagrangian_ = X; }

    /// Returns the overlap matrix
    SharedComplexMatrix S() const { return S_; }
    /// Set overlap matrix
    void set_S(SharedComplexMatrix S) { S_ = S; }

    /// Returns the core Hamiltonian matrix
    SharedComplexMatrix H() const { return H_; }
    /// Set core Hamiltonian matrix
    void set_H(SharedComplexMatrix H) { H_ = H; }

    /// Returns the gradient
    SharedComplexMatrix gradient() const { return gradient_; }
    /// Set the gradient for the wavefunction
    void set_gradient(SharedComplexMatrix grad);

    /// Returns the Hessian
    SharedComplexMatrix hessian() const { return hessian_; }
    /// Set the Hessian for the wavefunction
    void set_hessian(SharedComplexMatrix hess);

    /// Returns Molecular orbital extents
    std::vector<SharedVector> mo_extents() const { return mo_extents_; }
    /// Sets molecular orbital extents
    void set_mo_extents(const std::vector<SharedVector> mo_es) { mo_extents_ = mo_es; }

    /// Save the wavefunction to checkpoint
    virtual void save() const {}

    /// Get and set array variable dictionaries
    bool has_array_variable(const std::string& key);
    SharedComplexMatrix array_variable(const std::string& key);
    void set_array_variable(const std::string& key, SharedComplexMatrix value);
    int del_array_variable(const std::string& key);
    std::map<std::string, SharedComplexMatrix> array_variables();

    /// Vector of density matrices
    std::map<std::string, SharedComplexMatrix> density_map_;
};

#else  // !USING_Einsums

/// Stub type so pybind can expose ComplexMatrix without Einsums.
class PSI_API ComplexMatrix {};
using SharedComplexMatrix = std::shared_ptr<ComplexMatrix>;

/*! Stub ComplexWavefunction: public mol/basis constructors throw unless built with Einsums.
 *  Blank / Options-only constructors remain for derived stubs (e.g. CGHF). */
class PSI_API ComplexWavefunction : public BaseWavefunction {
   public:
    ComplexWavefunction() : BaseWavefunction() {}
    ComplexWavefunction(Options& options) : BaseWavefunction(options) {}
    ComplexWavefunction(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basis, Options& options)
        : BaseWavefunction(molecule, basis, options) {
        throw PSIEXCEPTION(
            "Psi4 not built with Einsums. ComplexWavefunction not available. Recompile with -DENABLE_Einsums=ON.");
    }
    ComplexWavefunction(std::shared_ptr<Molecule> molecule, std::shared_ptr<BasisSet> basis)
        : BaseWavefunction(molecule, basis) {
        throw PSIEXCEPTION(
            "Psi4 not built with Einsums. ComplexWavefunction not available. Recompile with -DENABLE_Einsums=ON.");
    }
    ~ComplexWavefunction() override = default;
};

#endif  // USING_Einsums

}  // namespace psi

#endif
