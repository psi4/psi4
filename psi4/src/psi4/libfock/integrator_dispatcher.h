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

#ifndef LIBFOCK_INTEGRATOR_DISPATCHER_H
#define LIBFOCK_INTEGRATOR_DISPATCHER_H

#include "psi4/libmints/matrix.h"

#include <map>

namespace psi {

class IntegratorManager;
class SuperFunctional;
class VBase;

class IntegratorDispatcher {
  private:
    std::vector<std::shared_ptr<IntegratorManager>> managers_; 
    std::map<std::string, double> quad_values_;
  public:
    IntegratorDispatcher(std::vector<std::shared_ptr<IntegratorManager>> managers) : managers_(managers) {};
    IntegratorDispatcher(std::shared_ptr<IntegratorManager> manager) : managers_({manager}) {};
    void initialize();
    void print_header() const;
    std::map<std::string, double> compute_V(std::vector<SharedMatrix> ret);
    SharedMatrix compute_gradient();
    void compute_Vx(const std::vector<SharedMatrix> Dx, std::vector<SharedMatrix> ret);
    void compute_Vx_triplet(const std::vector<SharedMatrix> Dx, std::vector<SharedMatrix> ret);
    SharedMatrix compute_hessian();
    std::vector<SharedMatrix> compute_fock_derivatives(); 
    void set_D(std::vector<SharedMatrix> Dvec);
    void set_grac_shift(double);
    std::map<std::string, double>& quadrature_values() { return quad_values_; };
    // Right now, only VBase has a cache. If that changes, rethink the design.
    size_t collocation_size() const;
    void build_collocation_cache(size_t doubles_per_integrator);
    void clear_collocation_cache();
    // The functional, if any.
    std::shared_ptr<SuperFunctional> functional();
    std::shared_ptr<VBase> psi_manager() const;
};

}

#endif

