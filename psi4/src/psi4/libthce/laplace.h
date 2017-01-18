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

#ifndef LAPLACE_H
#define LAPLACE_H

#include "psi4/libmints/typedefs.h"

namespace psi {

class Tensor;
class CoreTensor;

/**
 * LaplaceDenom is the second-generation Laplace Denominator class,
 * specifically developed for THC applications in arbitrary-order PT.
 *
 * The decomposition has the topology
 *
 * [(e_a + e_b + ...)_o - (e_i + e_j + ...)_r + \omega]^{-1} \approx
 * \tau_w^i \tau_w^j ... \tau_w^a \tau_w^b ...
 *
 * e_i - occupied eigenvalues
 * e_a - virtual eigenvalues
 * \omega - bias, e.g., in EOM-CC2 applications
 * r - rank, e.g., 2 in MP2
 *
 * delta is the maximum infinity norm in the rescaled problem 1/r in [1, R]
 **/
class LaplaceDenom {

protected:

    // => Input Spec <= //

    /// Active occupied orbital eigenvalues
    std::shared_ptr<Vector> eps_occ_;
    /// Active virtual orbital eigenvalues
    std::shared_ptr<Vector> eps_vir_;
    /// Allowed maximum error
    double delta_;
    /// Bias (careful if negative)
    double omega_;
    /// Number of times the occupied or virtual indices occur (e.g., 2 for MP2)
    int rank_;

    // => Output Spec <= //

    /// Number of Laplace points selected
    int npoints_;
    /// Laplace factor for occupied space, N_w x N_i
    std::shared_ptr<Tensor> tau_occ_;
    /// Laplace factor for virtual space, N_w x N_a
    std::shared_ptr<Tensor> tau_vir_;

public:

    // => Constructors <= //

    /// Master constructor
    LaplaceDenom(std::shared_ptr<Vector> eps_occ,
                 std::shared_ptr<Vector> eps_vir,
                 double delta = 1.0E-6,
                 double omega = 0.0,
                 int rank = 2);

    /// Master destructor
    virtual ~LaplaceDenom();

    // => Computers <= //

    /// Compute the Laplace factors, assigning tensor names by occ_name and vir_name
    void compute(const std::string& occ_name = "U", const std::string& vir_name = "V", double power = 1.0);

    // => Accessors <= //

    /// Number of Laplace points
    int npoints() const { return npoints_; }
    /// Laplace factor for occupied space, N_w x N_i
    std::shared_ptr<Tensor> tau_occ() const { return tau_occ_; }
    /// Laplace factor for virtual space, N_w x N_a
    std::shared_ptr<Tensor> tau_vir() const { return tau_vir_; }

};


} // End namespace

#endif
