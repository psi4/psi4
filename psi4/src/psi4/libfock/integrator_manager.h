/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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

#ifndef LIBFOCK_INTEGRATOR_MANAGER_H
#define LIBFOCK_INTEGRATOR_MANAGER_H

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"

#include <map>

namespace psi {
class BasisSet;
class Options;
class SuperFunctional;

class IntegratorManager {
   protected:
    /// Options object, used to build grid
    Options& options_;
    /// Basis set used in the integration
    std::shared_ptr<BasisSet> primary_;
    /// Number of basis functions. Populated by set_D.
    size_t nbf_;
    /// Quadrature values obtained during integration
    std::map<std::string, double> quad_values_;
    /// Vector of C1 D matrices (built by USO2AO)
    std::vector<SharedMatrix> D_AO_;
    /// AO2USO matrix (if not C1)
    SharedMatrix AO2USO_;
    SharedMatrix USO2AO_;
    bool can_compute_gradient_ = false;
    bool can_compute_Vx_ = false;
    bool can_compute_Vx_triplet_ = false;
    bool can_compute_hessian_ = false;
    bool can_compute_grac_ = false;

   public:
    IntegratorManager(std::shared_ptr<BasisSet> primary, Options& options) : primary_(primary), options_(options) {};
    void set_D(std::vector<SharedMatrix> Dvec);
    virtual void initialize() = 0;
    virtual void print_header() const = 0;
    /// Must be implemented for any manager.
    virtual std::map<std::string, double> compute_V(std::vector<SharedMatrix> ret) = 0;
    /// Need to distinguish between not being defined and throwing, so proceed to next
    virtual SharedMatrix compute_gradient() { return nullptr; };
    /// Either works or it doesn't.
    virtual void compute_Vx(const std::vector<SharedMatrix> Dx, std::vector<SharedMatrix> ret) {};
    /// Either works or it doesn't.
    virtual void compute_Vx_triplet(const std::vector<SharedMatrix> Dx, std::vector<SharedMatrix> ret) {};
    /// Need to distinguish between not being defined and throwing, so proceed to next
    virtual SharedMatrix compute_hessian() { return nullptr; };
    virtual std::vector<SharedMatrix> compute_fock_derivatives() { return {}; };
    const std::shared_ptr<Molecule> molecule();
    virtual bool can_compute_gradient() { return can_compute_gradient_; } ;
    virtual bool can_compute_Vx() { return can_compute_Vx_; } ;
    virtual bool can_compute_Vx_triplet() { return can_compute_Vx_triplet_; } ;
    virtual bool can_compute_hessian() { return can_compute_hessian_; } ;
    virtual bool can_compute_grac() { return can_compute_grac_; } ;
    std::map<std::string, double>& quadrature_values() { return quad_values_; } ;

    /// Size of collocation matrices. Default to 0, which means no information.
    virtual size_t collocation_size() const { return 0; } ;

    const std::vector<SharedMatrix>& Dao() const { return D_AO_; } ;
    std::shared_ptr<BasisSet> basis() const { return primary_; } ;

    // Return the functional object, if any.
    virtual std::shared_ptr<SuperFunctional> functional() const { return nullptr; } ;
};

}

#endif

