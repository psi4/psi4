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
#include "mrcc.h"
#include "matrix.h"
#include "blas.h"
#include "psi4/libpsi4util/libpsi4util.h"

namespace psi {
namespace psimrcc {

void CCMRCC::build_Z_intermediates() {
    // I am rewriting
    //   wfn_->blas()->append("Z_ijam[oov][o]{u} = #1234#  1/2 tau[oo][vv]{u} 2@2 <[vo]:[vv]>");
    // as
    wfn_->blas()->append("Z_ijam[oov][o]{u} = #1234#   tau[oo][vv]{u} 2@2 <[vo]|[vv]>");
    //

    wfn_->blas()->append("Z_iJaM[oOv][O]{u} = #1234#   tau[oO][vV]{u} 2@2 <[vo]|[vv]>");

    // I am rewriting
    //   wfn_->blas()->append("Z_iJAm[oOV][o]{u} = #1243# - tau[oO][vV]{u} 2@2 <[ov]|[vv]>");
    // as
    wfn_->blas()->append("Z_iJAm[oOV][o]{u} = #1234# - tau[oO][Vv]{u} 2@2 <[vo]|[vv]>");

    // I am rewriting
    // wfn_->blas()->append("Z_IJAM[OOV][O]{u} = #1234#  1/2 tau[OO][VV]{u} 2@2 <[vo]:[vv]>");
    // as
    wfn_->blas()->append("Z_IJAM[OOV][O]{u} = #1234#   tau[OO][VV]{u} 2@2 <[vo]|[vv]>");
    //
}

}  // namespace psimrcc
}  // namespace psi
