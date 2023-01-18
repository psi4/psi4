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

void CCMRCC::build_W_intermediates() {
    build_W_mnij_intermediates();
    build_W_mNiJ_intermediates();
    build_W_MNIJ_intermediates();

    build_W_jbme_intermediates();
    build_W_JBme_intermediates();
    build_W_jBmE_intermediates();
    build_W_jbME_intermediates();
    build_W_JbMe_intermediates();
    build_W_JBME_intermediates();
}

void CCMRCC::build_W_mnij_intermediates() {
    wfn_->blas()->append("W_mnij[oo][oo]{u}  = <[oo]:[oo]>");
    wfn_->blas()->append("W_mnij[oo][oo]{u} += #1234# <[ooo]:[v]> 2@2 t1[o][v]{u}");
    wfn_->blas()->append("W_mnij[oo][oo]{u} += #1243# - <[ooo]:[v]> 2@2 t1[o][v]{u}");
    wfn_->blas()->append("W_mnij[oo][oo]{u} += 1/2 <[oo]:[vv]> 2@2 tau[oo][vv]{u}");
}

void CCMRCC::build_W_mNiJ_intermediates() {
    wfn_->blas()->append("W_mNiJ[oO][oO]{u}  = <[oo]|[oo]>");
    wfn_->blas()->append("W_mNiJ[oO][oO]{u} += #1234# <[ooo]|[v]> 2@2 t1[O][V]{u}");
    wfn_->blas()->append("W_mNiJ[oO][oO]{u} += #2143# <[ooo]|[v]> 2@2 t1[o][v]{u}");
    wfn_->blas()->append("W_mNiJ[oO][oO]{u} += <[oo]|[vv]> 2@2 tau[oO][vV]{u}");
}

void CCMRCC::build_W_MNIJ_intermediates() {
    wfn_->blas()->append("W_MNIJ[OO][OO]{u}  = <[oo]:[oo]>");
    wfn_->blas()->append("W_MNIJ[OO][OO]{u} += #1234# <[ooo]:[v]> 2@2 t1[O][V]{u}");
    wfn_->blas()->append("W_MNIJ[OO][OO]{u} += #1243# - <[ooo]:[v]> 2@2 t1[O][V]{u}");
    wfn_->blas()->append("W_MNIJ[OO][OO]{u} += 1/2 <[oo]:[vv]> 2@2 tau[OO][VV]{u}");
}

void CCMRCC::build_W_jbme_intermediates() {
    wfn_->blas()->append("W_jbme[ov][ov]{u}  = #3241# <[ov]:[vo]>");

    // This term uses an extra integral file
    // wfn_->blas()->append("W_jbme[ov][ov]{u} += #3241# <[ovv]:[v]> 2@2 t1[o][v]{u}");
    // I will rewrite it as two terms:
    wfn_->blas()->append("W_jbme[ov][ov]{u} += #3241#   <[v]|[ovv]> 1@2 t1[o][v]{u}");
    wfn_->blas()->append("W_jbme[ov][ov]{u} += #2431# - ([vvo]|[v]) 2@2 t1[o][v]{u}");
    //

    wfn_->blas()->append("W_jbme[ov][ov]{u} += #2314# - t1[o][v]{u} 1@1 <[o]:[oov]>");
    wfn_->blas()->append("W_jbme[ov][ov]{u} += - tau3[ov][ov]{u} 2@2 ([ov]:[ov])");
    wfn_->blas()->append("W_jbme[ov][ov]{u} += 1/2 t2[ov][OV]{u} 2@2 ([ov]|[ov])");
}

void CCMRCC::build_W_JBme_intermediates() {
    // Open-Shell
    wfn_->blas()->append("W_JBme[OV][ov]{o}  = #3241# <[ov]|[vo]>");
    // wfn_->blas()->append("W_JBme[OV][ov]{o} += #3241# <[ovv]|[v]> 2@2 t1[O][V]{o}");
    wfn_->blas()->append("W_JBme[OV][ov]{o} += #3241# <[v]|[ovv]> 1@2 t1[O][V]{o}");
    wfn_->blas()->append("W_JBme[OV][ov]{o} += #2314# - t1[O][V]{o} 1@1 <[o]|[oov]>");
    wfn_->blas()->append("W_JBme[OV][ov]{o} += - tau3[OV][OV]{o} 2@2 ([ov]|[ov])");
    wfn_->blas()->append("W_JBme[OV][ov]{o} += 1/2 t2[ov][OV]{o} 1@2 ([ov]:[ov])");
}

void CCMRCC::build_W_jBmE_intermediates() {
    wfn_->blas()->append("W_jBmE[oV][oV]{u}  = #3214# - <[ov]|[ov]>");
    wfn_->blas()->append("W_jBmE[oV][oV]{u} += #2431# - ([vvo]|[v]) 2@2 t1[o][v]{u}");
    wfn_->blas()->append("W_jBmE[oV][oV]{u} += #2341#   t1[O][V]{u} 1@1 <[o]|[ovo]>");
    wfn_->blas()->append("W_jBmE[oV][oV]{u} += tau3[oV][vO]{u} 2@2 <[ov]|[vo]>");
}

void CCMRCC::build_W_jbME_intermediates() {
    wfn_->blas()->append("W_jbME[ov][OV]{u}  = #3241# <[ov]|[vo]>");
    //  wfn_->blas()->append("W_jbME[ov][OV]{u} += #3241# <[ovv]|[v]> 2@2 t1[o][v]{u}");
    wfn_->blas()->append("W_jbME[ov][OV]{u} += #3241# <[v]|[ovv]> 1@2 t1[o][v]{u}");
    wfn_->blas()->append("W_jbME[ov][OV]{u} += #2314# - t1[o][v]{u} 1@1 <[o]|[oov]>");
    wfn_->blas()->append("W_jbME[ov][OV]{u} += - tau3[ov][ov]{u} 2@2 ([ov]|[ov])");
    wfn_->blas()->append("W_jbME[ov][OV]{u} += 1/2 t2[ov][OV]{u} 2@2 ([ov]:[ov])");
}

void CCMRCC::build_W_JbMe_intermediates() {
    // Open-shell
    wfn_->blas()->append("W_JbMe[Ov][Ov]{o}  = #3214# - <[ov]|[ov]>");
    wfn_->blas()->append("W_JbMe[Ov][Ov]{o} += #2431# - ([vvo]|[v]) 2@2 t1[O][V]{o}");
    wfn_->blas()->append("W_JbMe[Ov][Ov]{o} += #2341#   t1[o][v]{o} 1@1 <[o]|[ovo]>");
    wfn_->blas()->append("W_JbMe[Ov][Ov]{o} += tau3[Ov][Vo]{o} 2@2 <[ov]|[vo]>");
}

void CCMRCC::build_W_JBME_intermediates() {
    wfn_->blas()->append("W_JBME[OV][OV]{o}  = #3241# <[ov]:[vo]>");

    // This term uses an extra integral file:
    // wfn_->blas()->append("W_JBME[OV][OV]{o} += #3241# <[ovv]:[v]> 2@2 t1[O][V]{o}");
    // I will rewrite it as two terms:
    wfn_->blas()->append("W_JBME[OV][OV]{o} += #3241#   <[v]|[ovv]> 1@2 t1[O][V]{o}");
    wfn_->blas()->append("W_JBME[OV][OV]{o} += #2431# - ([vvo]|[v]) 2@2 t1[O][V]{o}");
    //

    wfn_->blas()->append("W_JBME[OV][OV]{o} += #2314# - t1[O][V]{o} 1@1 <[o]:[oov]>");
    wfn_->blas()->append("W_JBME[OV][OV]{o} += - tau3[OV][OV]{o} 2@2 ([ov]:[ov])");
    wfn_->blas()->append("W_JBME[OV][OV]{o} += 1/2 t2[ov][OV]{o} 1@2 ([ov]|[ov])");
}

}  // namespace psimrcc
}  // namespace psi
