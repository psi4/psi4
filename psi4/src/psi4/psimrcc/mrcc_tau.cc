/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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

/***************************************************************************
 *  PSIMRCC : Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/
#include "psi4/libmoinfo/libmoinfo.h"
#include "mrcc.h"
#include "matrix.h"
#include "blas.h"
#include "psi4/libpsi4util/libpsi4util.h"

namespace psi {
namespace psimrcc {

// build_tau_intermediates("tau2[v][voo]{u}","3412");
// build_tau_intermediates("tau2[o][ovv]{u}","1234");

/**
 * @brief Computes the contractions
 * \f[ \tau_{ij}^{ab} = t_{ij}^{ab} + t_i^a t_j^b - t_i^b t_j^a \f]
 * \f[ \tilde{\tau}_{ij}^{ab} = t_{ij}^{ab} + \frac{1}{2} (t_i^a t_j^b - t_i^b t_j^a) \f]
 * \f[ \hat{\tau}_{ij}^{ab} = \frac{1}{2} t_{ij}^{ab} +  t_i^a t_j^b \f]
 * as described in J. Phys. Chem. vol. 94, pg. 4334 (1991).
 * See J. Phys. Chem. vol. 127, 024102 (2007) supplementary material for the spin-factored equations.
 */
void CCMRCC::build_tau_intermediates() {
    // t1t1[ov][ov]{u}, Ok
    wfn_->blas()->solve("t1t1_iame[ov][ov]{u} = #1432#   t1[o][v]{u} X t1[o][v]{u}");
    wfn_->blas()->solve("t1t1_IAME[OV][OV]{u} = #1432#   t1[O][V]{u} X t1[O][V]{u}");
    wfn_->blas()->solve("t1t1_iAMe[oV][Ov]{u} = #1432#   t1[o][v]{u} X t1[O][V]{u}");

    // tau[oo][vv]{u}, Ok
    wfn_->blas()->solve("tau[oo][vv]{u}  = t2[oo][vv]{u}");
    wfn_->blas()->solve("tau[oo][vv]{u} += #1324#   t1[o][v]{u} X t1[o][v]{u}");
    wfn_->blas()->solve("tau[oo][vv]{u} += #2314# - t1[o][v]{u} X t1[o][v]{u}");
    // tau[oO][vV]{u}, Ok
    wfn_->blas()->solve("tau[oO][vV]{u}  = t2[oO][vV]{u}");
    wfn_->blas()->solve("tau[oO][vV]{u} += #1324#   t1[o][v]{u} X t1[O][V]{u}");
    // tau[OO][VV]{u}, Ok
    wfn_->blas()->solve("tau[OO][VV]{u}  = t2[OO][VV]{u}");
    wfn_->blas()->solve("tau[OO][VV]{u} += #1324#   t1[O][V]{u} X t1[O][V]{u}");
    wfn_->blas()->solve("tau[OO][VV]{u} += #2314# - t1[O][V]{u} X t1[O][V]{u}");

    // tau[oO][Vv]{u}, Ok
    wfn_->blas()->solve("tau[oO][Vv]{u}  = #1243#   tau[oO][vV]{u}");

    // tau2[v][voo]{u}, Ok
    wfn_->blas()->solve("tau2[v][voo]{u}  = #3412# t2[oo][vv]{u}");
    wfn_->blas()->solve("tau2[v][voo]{u} += #3142# 1/2 t1[o][v]{u} X t1[o][v]{u}");
    wfn_->blas()->solve("tau2[v][voo]{u} += #4132# -1/2 t1[o][v]{u} X t1[o][v]{u}");

    // tau2[V][VOO]{u}, Ok
    wfn_->blas()->solve("tau2[V][VOO]{u}  = #3412# t2[OO][VV]{u}");
    wfn_->blas()->solve("tau2[V][VOO]{u} += #3142# 1/2 t1[O][V]{u} X t1[O][V]{u}");
    wfn_->blas()->solve("tau2[V][VOO]{u} += #4132# -1/2 t1[O][V]{u} X t1[O][V]{u}");

    // tau2[v][VoO]{u}, Ok
    wfn_->blas()->solve("tau2[v][VoO]{u}  = #3412# t2[oO][vV]{u}");
    wfn_->blas()->solve("tau2[v][VoO]{u} += #3142# 1/2 t1[o][v]{u} X t1[O][V]{u}");

    // tau2[V][vOo]{u}, Ok
    wfn_->blas()->solve("tau2[V][vOo]{u}  = #4321# t2[oO][vV]{u}");
    wfn_->blas()->solve("tau2[V][vOo]{u} += #4231# 1/2 t1[o][v]{u} X t1[O][V]{u}");

    // tau2[o][ovv]{u}, Ok
    wfn_->blas()->solve("tau2[o][ovv]{u}  = #1234# t2[oo][vv]{u}");
    wfn_->blas()->solve("tau2[o][ovv]{u} += #1324# 1/2 t1[o][v]{u} X t1[o][v]{u}");
    wfn_->blas()->solve("tau2[o][ovv]{u} += #2314# -1/2 t1[o][v]{u} X t1[o][v]{u}");

    // tau2[O][OVV]{u}, Ok
    wfn_->blas()->solve("tau2[O][OVV]{u}  = #1234# t2[OO][VV]{u}");
    wfn_->blas()->solve("tau2[O][OVV]{u} += #1324# 1/2 t1[O][V]{u} X t1[O][V]{u}");
    wfn_->blas()->solve("tau2[O][OVV]{u} += #2314# -1/2 t1[O][V]{u} X t1[O][V]{u}");

    // tau2[o][OvV]{u}, Ok
    wfn_->blas()->solve("tau2[o][OvV]{u}  = #1234# t2[oO][vV]{u}");
    wfn_->blas()->solve("tau2[o][OvV]{u} += #1324# 1/2 t1[o][v]{u} X t1[O][V]{u}");

    // tau2[O][oVv]{u}, Ok
    wfn_->blas()->solve("tau2[O][oVv]{u}  = #2143# t2[oO][vV]{u}");
    wfn_->blas()->solve("tau2[O][oVv]{u} += #2413# 1/2 t1[o][v]{u} X t1[O][V]{u}");

    /////////////////////////////////////////////
    // If you ask for tau3[pqrs] this will produce
    // tau3[psqr]
    /////////////////////////////////////////////

    // tau3[ov][ov]{u}
    wfn_->blas()->solve("tau3[ov][ov]{u}  = #1342# 1/2 t2[oo][vv]{u}");
    wfn_->blas()->solve("tau3[ov][ov]{u} += #1432# t1[o][v]{u} X t1[o][v]{u}");

    // tau3[OV][OV]{u}
    wfn_->blas()->solve("tau3[OV][OV]{u}  = #1342# 1/2 t2[OO][VV]{u}");
    wfn_->blas()->solve("tau3[OV][OV]{u} += #1432# t1[O][V]{u} X t1[O][V]{u}");

    // tau3[oV][vO]{u}
    wfn_->blas()->solve("tau3[oV][vO]{u}  = #1432# 1/2 t2[oO][vV]{u}");
    wfn_->blas()->solve("tau3[oV][vO]{u} += #1342# t1[o][v]{u} X t1[O][V]{u}");

    // tau3[Ov][Vo]{u}
    wfn_->blas()->solve("tau3[Ov][Vo]{u}  = #4123# 1/2 t2[oO][vV]{u}");
    wfn_->blas()->solve("tau3[Ov][Vo]{u} += #4213# t1[o][v]{u} X t1[O][V]{u}");

    // tau[oo][v>v]{u}, Ok
    // tau[OO][V>V]{u}, Ok
    wfn_->blas()->solve("tau[oo][v>v]{u}  = #1234# tau[oo][vv]{u}");
    wfn_->blas()->solve("tau[OO][V>V]{u}  = #1234# tau[OO][VV]{u}");

    wfn_->blas()->solve("tau[oO][v>=V]{u} = #1234# tau[oO][vV]{u}");
    wfn_->blas()->solve("tau[oO][V>=v]{u} = #1243# tau[oO][vV]{u}");
    wfn_->blas()->zero_right_four_diagonal("tau[oO][V>=v]{u}");
}

}  // namespace psimrcc
}  // namespace psi
