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

#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libmoinfo/libmoinfo.h"

#include "blas.h"
#include "mp2_ccsd.h"
#include "matrix.h"

namespace psi {
namespace psimrcc {

void MP2_CCSD::build_Z_intermediates() {
    wfn_->blas()->solve("Z_iJaM[aAa][O]{u} = #1234#   tau_oOvV[aA][vV]{u} 2@2 <[ao]|[vv]>");
    wfn_->blas()->solve("Z_iJAm[aAA][o]{u} = #1234# - tau_oOVv[aA][Vv]{u} 2@2 <[ao]|[vv]>");

    wfn_->blas()->solve("Z_iJaM[oAa][O]{u} = #1234#   tau_oOvV[oA][vV]{u} 2@2 <[ao]|[vv]>");
    wfn_->blas()->solve("Z_iJAm[oAA][o]{u} = #1234# - tau_oOVv[oA][Vv]{u} 2@2 <[ao]|[vv]>");

    wfn_->blas()->solve("Z_iJaM[aAv][O]{u} = #1234#   tau_oOvV[aA][vV]{u} 2@2 <[vo]|[vv]>");
}

}  // namespace psimrcc
}  // namespace psi
