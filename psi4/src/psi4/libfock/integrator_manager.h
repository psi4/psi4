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
    /// Desired superfunctional kernel
    std::shared_ptr<SuperFunctional> functional_;
    /// Vector of C1 D matrices (built by USO2AO)
    std::vector<SharedMatrix> D_AO_;
    /// AO2USO matrix (if not C1)
    SharedMatrix AO2USO_;
    SharedMatrix USO2AO_;

   public:
    IntegratorManager(std::shared_ptr<SuperFunctional> functional, std::shared_ptr<BasisSet> primary, Options& options) : functional_(functional), primary_(primary), options_(options) {};
    void set_D(std::vector<SharedMatrix> Dvec);
    virtual void initialize() = 0;
    /// Throws by default
    virtual std::map<std::string, double> compute_V(std::vector<SharedMatrix> ret) = 0;
};

}

#endif

