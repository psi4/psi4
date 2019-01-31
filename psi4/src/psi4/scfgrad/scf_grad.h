/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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

#ifndef SCF_GRAD_H
#define SCF_GRAD_H

#include "psi4/libmints/wavefunction.h"
#include "psi4/libmints/typedefs.h"

namespace psi {
class SuperFunctional;
class VBase;

namespace scfgrad {

class SCFGrad : public Wavefunction {

protected:

    /// Common initialization
    void common_init();
    std::shared_ptr<SuperFunctional> functional_;
    std::shared_ptr<VBase> potential_;
    std::map<std::string, SharedMatrix> gradients_;
    std::map<std::string, SharedMatrix> hessians_;
    SharedMatrix dipole_gradient_;

public:
    SCFGrad(SharedWavefunction ref_wfn, Options& options);
    ~SCFGrad() override;

    double compute_energy() override { throw PSIEXCEPTION("SCFGrad needs a rehash, call Rob."); }

    SharedMatrix compute_gradient() override;

    SharedMatrix compute_hessian() override;

    SharedMatrix rhf_hessian_response();

    SharedMatrix dipole_gradient() const { return dipole_gradient_; }
};

}} // Namespaces

#endif
