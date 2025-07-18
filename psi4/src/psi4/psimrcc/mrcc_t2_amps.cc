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

/***************************************************************************
 *  PSIMRCC : Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/
#include "psi4/libmoinfo/libmoinfo.h"

#include "algebra_interface.h"
#include "blas.h"
#include "index.h"
#include "matrix.h"
#include "mrcc.h"

namespace psi {
namespace psimrcc {

// Many functions in this file have six-index local variables.
// This is a good place to look for large memory expenses.
// Because the variables are local, updating the wfn memory count is unnecessary.

void CCMRCC::build_t2_amplitudes() {
    build_t2_iJaB_amplitudes();
    build_t2_ijab_amplitudes();
    build_t2_IJAB_amplitudes();
}

void CCMRCC::build_t2_ijab_amplitudes() {
    if (wfn_->moinfo()->get_ref_size(UniqueOpenShellRefs) == 0) {
        wfn_->blas()->append("t2_eqns[oo][vv]{c}  = t2_eqns[oO][vV]{c}");
        wfn_->blas()->append("t2_eqns[oo][vv]{c} += #2134# - t2_eqns[oO][vV]{c}");
    } else {
        // Closed-shell
        wfn_->blas()->append("t2_eqns[oo][vv]{c}  = <[oo]:[vv]>");

        wfn_->blas()->append("t2_eqns[oo][vv]{c} += #3124# - t2[v][voo]{c} 1@2 F'_ae[v][v]{c}");
        wfn_->blas()->append("t2_eqns[oo][vv]{c} += #4123#   t2[v][voo]{c} 1@2 F'_ae[v][v]{c}");

        wfn_->blas()->append("t2_eqns[oo][vv]{c} += #1342#   t2[o][ovv]{c} 1@1 F'_mi[o][o]{c}");
        wfn_->blas()->append("t2_eqns[oo][vv]{c} += #2341# - t2[o][ovv]{c} 1@1 F'_mi[o][o]{c}");

        wfn_->blas()->append("t2_eqns[oo][vv]{c} += 1/2  W_mnij[oo][oo]{c} 1@1 tau[oo][vv]{c}");

        wfn_->blas()->append("t2_eqns[oo][v>v]{c} = tau[oo][v>v]{c} 2@2 <[v>v]:[v>v]>");

        wfn_->blas()->append("t2_eqns[oo][vv]{c} +>= #1234# t2_eqns[oo][v>v]{c}");

        wfn_->blas()->append("t2_eqns[oo][vv]{c} += #1234# - Z_ijam[oov][o]{c} 2@1 t1[o][v]{c}");
        wfn_->blas()->append("t2_eqns[oo][vv]{c} += #1243#   Z_ijam[oov][o]{c} 2@1 t1[o][v]{c}");

        wfn_->blas()->append("t2_eqns[oo][vv]{c} += #2413#   W_jbme[ov][ov]{c} 2@2 t2[ov][ov]{c}");
        wfn_->blas()->append("t2_eqns[oo][vv]{c} += #2314# - W_jbme[ov][ov]{c} 2@2 t2[ov][ov]{c}");

        wfn_->blas()->append("t2_eqns[oo][vv]{c} += #1423# - W_jbme[ov][ov]{c} 2@2 t2[ov][ov]{c}");
        wfn_->blas()->append("t2_eqns[oo][vv]{c} += #1324#   W_jbme[ov][ov]{c} 2@2 t2[ov][ov]{c}");

        //  wfn_->blas()->append("t2_eqns[oo][vv]{c} += #P-(34)P-(12)4213# - ([ov]:[vo]) 1@2 t1t1_iame[ov][ov]{c}");

        wfn_->blas()->append("t2_eqns[oo][vv]{c} += #4213# - ([ov]:[vo]) 1@2 t1t1_iame[ov][ov]{c}");
        wfn_->blas()->append("t2_eqns[oo][vv]{c} += #3214# + ([ov]:[vo]) 1@2 t1t1_iame[ov][ov]{c}");
        wfn_->blas()->append("t2_eqns[oo][vv]{c} += #4123# + ([ov]:[vo]) 1@2 t1t1_iame[ov][ov]{c}");
        wfn_->blas()->append("t2_eqns[oo][vv]{c} += #3124# - ([ov]:[vo]) 1@2 t1t1_iame[ov][ov]{c}");

        wfn_->blas()->append("t2_eqns[oo][vv]{c} += #2413#   W_jbME[ov][OV]{c} 2@2 t2[ov][OV]{c}");
        wfn_->blas()->append("t2_eqns[oo][vv]{c} += #2314# - W_jbME[ov][OV]{c} 2@2 t2[ov][OV]{c}");

        wfn_->blas()->append("t2_eqns[oo][vv]{c} += #1423# - W_jbME[ov][OV]{c} 2@2 t2[ov][OV]{c}");
        wfn_->blas()->append("t2_eqns[oo][vv]{c} += #1324#   W_jbME[ov][OV]{c} 2@2 t2[ov][OV]{c}");

        wfn_->blas()->append("t2_eqns[oo][vv]{c} += #1234#   t1[o][v]{c} 2@1 <[v]:[ovv]>");
        wfn_->blas()->append("t2_eqns[oo][vv]{c} += #2134# - t1[o][v]{c} 2@1 <[v]:[ovv]>");

        wfn_->blas()->append("t2_eqns[oo][vv]{c} += #3412# - t1[o][v]{c} 1@1 <[o]:[voo]>");
        wfn_->blas()->append("t2_eqns[oo][vv]{c} += #4312#   t1[o][v]{c} 1@1 <[o]:[voo]>");
    }

    // Open-shell
    wfn_->blas()->append("t2_eqns[oo][vv]{o}  = <[oo]:[vv]>");

    wfn_->blas()->append("t2_eqns[oo][vv]{o} += #3124# - t2[v][voo]{o} 1@2 F'_ae[v][v]{o}");
    wfn_->blas()->append("t2_eqns[oo][vv]{o} += #4123#   t2[v][voo]{o} 1@2 F'_ae[v][v]{o}");

    wfn_->blas()->append("t2_eqns[oo][vv]{o} += #1342#   t2[o][ovv]{o} 1@1 F'_mi[o][o]{o}");
    wfn_->blas()->append("t2_eqns[oo][vv]{o} += #2341# - t2[o][ovv]{o} 1@1 F'_mi[o][o]{o}");

    wfn_->blas()->append("t2_eqns[oo][vv]{o} += 1/2  W_mnij[oo][oo]{o} 1@1 tau[oo][vv]{o}");

    wfn_->blas()->append("t2_eqns[oo][v>v]{o} = tau[oo][v>v]{o} 2@2 <[v>v]:[v>v]>");

    wfn_->blas()->append("t2_eqns[oo][vv]{o} +>= #1234# t2_eqns[oo][v>v]{o}");

    wfn_->blas()->append("t2_eqns[oo][vv]{o} += #1234# - Z_ijam[oov][o]{o} 2@1 t1[o][v]{o}");
    wfn_->blas()->append("t2_eqns[oo][vv]{o} += #1243#   Z_ijam[oov][o]{o} 2@1 t1[o][v]{o}");

    wfn_->blas()->append("t2_eqns[oo][vv]{o} += #2413#   W_jbme[ov][ov]{o} 2@2 t2[ov][ov]{o}");
    wfn_->blas()->append("t2_eqns[oo][vv]{o} += #2314# - W_jbme[ov][ov]{o} 2@2 t2[ov][ov]{o}");

    wfn_->blas()->append("t2_eqns[oo][vv]{o} += #1423# - W_jbme[ov][ov]{o} 2@2 t2[ov][ov]{o}");
    wfn_->blas()->append("t2_eqns[oo][vv]{o} += #1324#   W_jbme[ov][ov]{o} 2@2 t2[ov][ov]{o}");

    wfn_->blas()->append("t2_eqns[oo][vv]{o} += #4213# - ([ov]:[vo]) 1@2 t1t1_iame[ov][ov]{o}");
    wfn_->blas()->append("t2_eqns[oo][vv]{o} += #3214# + ([ov]:[vo]) 1@2 t1t1_iame[ov][ov]{o}");
    wfn_->blas()->append("t2_eqns[oo][vv]{o} += #4123# + ([ov]:[vo]) 1@2 t1t1_iame[ov][ov]{o}");
    wfn_->blas()->append("t2_eqns[oo][vv]{o} += #3124# - ([ov]:[vo]) 1@2 t1t1_iame[ov][ov]{o}");

    wfn_->blas()->append("t2_eqns[oo][vv]{o} += #2413#   W_jbME[ov][OV]{o} 2@2 t2[ov][OV]{o}");
    wfn_->blas()->append("t2_eqns[oo][vv]{o} += #2314# - W_jbME[ov][OV]{o} 2@2 t2[ov][OV]{o}");

    wfn_->blas()->append("t2_eqns[oo][vv]{o} += #1423# - W_jbME[ov][OV]{o} 2@2 t2[ov][OV]{o}");
    wfn_->blas()->append("t2_eqns[oo][vv]{o} += #1324#   W_jbME[ov][OV]{o} 2@2 t2[ov][OV]{o}");

    wfn_->blas()->append("t2_eqns[oo][vv]{o} += #1234#   t1[o][v]{o} 2@1 <[v]:[ovv]>");
    wfn_->blas()->append("t2_eqns[oo][vv]{o} += #2134# - t1[o][v]{o} 2@1 <[v]:[ovv]>");

    wfn_->blas()->append("t2_eqns[oo][vv]{o} += #3412# - t1[o][v]{o} 1@1 <[o]:[voo]>");
    wfn_->blas()->append("t2_eqns[oo][vv]{o} += #4312#   t1[o][v]{o} 1@1 <[o]:[voo]>");
}

void CCMRCC::build_t2_iJaB_amplitudes() {
    // Closed-shell
    wfn_->blas()->append("t2_eqns[oO][vV]{c}  = <[oo]|[vv]>");

    wfn_->blas()->append("t2_eqns[oO][vV]{c} += #3214# t2[V][vOo]{c} 1@2 F'_ae[v][v]{c}");
    wfn_->blas()->append("t2_eqns[oO][vV]{c} += #4123# t2[v][VoO]{c} 1@2 F'_ae[v][v]{c}");

    wfn_->blas()->append("t2_eqns[oO][vV]{c} += #1432# - t2[O][oVv]{c} 1@1 F'_mi[o][o]{c}");
    wfn_->blas()->append("t2_eqns[oO][vV]{c} += #2341# - t2[o][OvV]{c} 1@1 F'_mi[o][o]{c}");

    wfn_->blas()->append("t2_eqns[oO][vV]{c} += W_mNiJ[oO][oO]{c} 1@1 tau[oO][vV]{c}");

    // wfn_->blas()->append("t2_eqns[oO][vV]{c} += tau[oO][vV]{c} 2@2 <[vv]|[vv]>");

    wfn_->blas()->append("t2_eqns[oO][vV]{c} += tau[oO][v>=V]{c} 2@2 <[vv]|[v>=v]>");
    wfn_->blas()->append("t2_eqns[oO][vV]{c} += #1243# tau[oO][V>=v]{c} 2@2 <[vv]|[v>=v]>");

    wfn_->blas()->append("t2_eqns[oO][vV]{c} += #1234#  - Z_iJaM[oOv][O]{c} 2@1 t1[O][V]{c}");
    wfn_->blas()->append("t2_eqns[oO][vV]{c} += #1243#    Z_iJAm[oOV][o]{c} 2@1 t1[o][v]{c}");

    wfn_->blas()->append("t2_eqns[oO][vV]{c} += #2413#   W_jbME[ov][OV]{c} 2@2 t2[ov][ov]{c}");
    wfn_->blas()->append("t2_eqns[oO][vV]{c} += #2413#   W_jbme[ov][ov]{c} 2@2 t2[ov][OV]{c}");

    wfn_->blas()->append("t2_eqns[oO][vV]{c} += #2314#   W_jBmE[oV][oV]{c} 2@2 t2[oV][Ov]{c}");
    wfn_->blas()->append("t2_eqns[oO][vV]{c} += #1423#   W_jBmE[oV][oV]{c} 2@1 t2[oV][Ov]{c}");

    wfn_->blas()->append("t2_eqns[oO][vV]{c} += #1324#   W_jbME[ov][OV]{c} 2@2 t2[OV][OV]{c}");

    wfn_->blas()->append("t2_eqns[oO][vV]{c} += #1324#   W_jbme[ov][ov]{c} 2@1 t2[ov][OV]{c}");

    wfn_->blas()->append("t2_eqns[oO][vV]{c} += #4213# - ([ov]|[vo]) 1@2 t1t1_iame[ov][ov]{c}");

    wfn_->blas()->append("t2_eqns[oO][vV]{c} += #2314# - <[ov]|[ov]> 1@2 t1t1_iAMe[oV][Ov]{c}");

    wfn_->blas()->append("t2_eqns[oO][vV]{c} += #1423# - <[ov]|[ov]> 1@1 t1t1_iAMe[oV][Ov]{c}");

    wfn_->blas()->append("t2_eqns[oO][vV]{c} += #3124# - ([ov]|[vo]) 1@2 t1t1_IAME[OV][OV]{c}");

    wfn_->blas()->append("t2_eqns[oO][vV]{c} += #1234#   t1[o][v]{c} 2@1 <[v]|[ovv]>");
    wfn_->blas()->append("t2_eqns[oO][vV]{c} += #2143#   t1[O][V]{c} 2@1 <[v]|[ovv]>");

    wfn_->blas()->append("t2_eqns[oO][vV]{c} += #3412# - t1[o][v]{c} 1@1 <[o]|[voo]>");
    wfn_->blas()->append("t2_eqns[oO][vV]{c} += #4321# - t1[O][V]{c} 1@1 <[o]|[voo]>");

    // Open-shell
    wfn_->blas()->append("t2_eqns[oO][vV]{o}  = <[oo]|[vv]>");

    wfn_->blas()->append("t2_eqns[oO][vV]{o} += #3214# t2[V][vOo]{o} 1@2 F'_AE[V][V]{o}");
    wfn_->blas()->append("t2_eqns[oO][vV]{o} += #4123# t2[v][VoO]{o} 1@2 F'_ae[v][v]{o}");

    wfn_->blas()->append("t2_eqns[oO][vV]{o} += #1432# - t2[O][oVv]{o} 1@1 F'_MI[O][O]{o}");
    wfn_->blas()->append("t2_eqns[oO][vV]{o} += #2341# - t2[o][OvV]{o} 1@1 F'_mi[o][o]{o}");

    wfn_->blas()->append("t2_eqns[oO][vV]{o} += W_mNiJ[oO][oO]{o} 1@1 tau[oO][vV]{o}");

    // wfn_->blas()->append("t2_eqns[oO][vV]{o} += tau[oO][vV]{o} 2@2 <[vv]|[vv]>");

    wfn_->blas()->append("t2_eqns[oO][vV]{o} += tau[oO][v>=V]{o} 2@2 <[vv]|[v>=v]>");
    wfn_->blas()->append("t2_eqns[oO][vV]{o} += #1243# tau[oO][V>=v]{o} 2@2 <[vv]|[v>=v]>");

    wfn_->blas()->append("t2_eqns[oO][vV]{o} += #1234#  - Z_iJaM[oOv][O]{o} 2@1 t1[O][V]{o}");
    wfn_->blas()->append("t2_eqns[oO][vV]{o} += #1243#    Z_iJAm[oOV][o]{o} 2@1 t1[o][v]{o}");

    wfn_->blas()->append("t2_eqns[oO][vV]{o} += #2413#   W_JBme[OV][ov]{o} 2@2 t2[ov][ov]{o}");
    wfn_->blas()->append("t2_eqns[oO][vV]{o} += #2413#   W_JBME[OV][OV]{o} 2@2 t2[ov][OV]{o}");

    wfn_->blas()->append("t2_eqns[oO][vV]{o} += #2314#   W_JbMe[Ov][Ov]{o} 2@2 t2[oV][Ov]{o}");
    wfn_->blas()->append("t2_eqns[oO][vV]{o} += #1423#   W_jBmE[oV][oV]{o} 2@1 t2[oV][Ov]{o}");

    wfn_->blas()->append("t2_eqns[oO][vV]{o} += #1324#   W_jbME[ov][OV]{o} 2@2 t2[OV][OV]{o}");

    wfn_->blas()->append("t2_eqns[oO][vV]{o} += #1324#   W_jbme[ov][ov]{o} 2@1 t2[ov][OV]{o}");

    wfn_->blas()->append("t2_eqns[oO][vV]{o} += #4213# - ([ov]|[vo]) 1@2 t1t1_iame[ov][ov]{o}");

    wfn_->blas()->append("t2_eqns[oO][vV]{o} += #2314# - <[ov]|[ov]> 1@2 t1t1_iAMe[oV][Ov]{o}");

    wfn_->blas()->append("t2_eqns[oO][vV]{o} += #1423# - <[ov]|[ov]> 1@1 t1t1_iAMe[oV][Ov]{o}");

    wfn_->blas()->append("t2_eqns[oO][vV]{o} += #3124# - ([ov]|[vo]) 1@2 t1t1_IAME[OV][OV]{o}");

    wfn_->blas()->append("t2_eqns[oO][vV]{o} += #1234#   t1[o][v]{o} 2@1 <[v]|[ovv]>");
    wfn_->blas()->append("t2_eqns[oO][vV]{o} += #2143#   t1[O][V]{o} 2@1 <[v]|[ovv]>");

    wfn_->blas()->append("t2_eqns[oO][vV]{o} += #3412# - t1[o][v]{o} 1@1 <[o]|[voo]>");
    wfn_->blas()->append("t2_eqns[oO][vV]{o} += #4321# - t1[O][V]{o} 1@1 <[o]|[voo]>");
}

void CCMRCC::build_t2_IJAB_amplitudes() {
    // Closed-shell
    wfn_->blas()->append("t2_eqns[OO][VV]{c}  = t2_eqns[oo][vv]{c}");

    // Open-shell
    wfn_->blas()->append("t2_eqns[OO][VV]{o}  = <[oo]:[vv]>");

    wfn_->blas()->append("t2_eqns[OO][VV]{o} += #3124# - t2[V][VOO]{o} 1@2 F'_AE[V][V]{o}");
    wfn_->blas()->append("t2_eqns[OO][VV]{o} += #4123#   t2[V][VOO]{o} 1@2 F'_AE[V][V]{o}");

    wfn_->blas()->append("t2_eqns[OO][VV]{o} += #1342#   t2[O][OVV]{o} 1@1 F'_MI[O][O]{o}");
    wfn_->blas()->append("t2_eqns[OO][VV]{o} += #2341# - t2[O][OVV]{o} 1@1 F'_MI[O][O]{o}");

    wfn_->blas()->append("t2_eqns[OO][VV]{o} += 1/2  W_MNIJ[OO][OO]{o} 1@1 tau[OO][VV]{o}");

    wfn_->blas()->append("t2_eqns[OO][V>V]{o} = tau[OO][V>V]{o} 2@2 <[v>v]:[v>v]>");

    wfn_->blas()->append("t2_eqns[OO][VV]{o} +>= #1234#  t2_eqns[OO][V>V]{o}");

    wfn_->blas()->append("t2_eqns[OO][VV]{o} += #1234# - Z_IJAM[OOV][O]{o} 2@1 t1[O][V]{o}");
    wfn_->blas()->append("t2_eqns[OO][VV]{o} += #1243#   Z_IJAM[OOV][O]{o} 2@1 t1[O][V]{o}");

    wfn_->blas()->append("t2_eqns[OO][VV]{o} += #2413#   W_JBME[OV][OV]{o} 2@2 t2[OV][OV]{o}");
    wfn_->blas()->append("t2_eqns[OO][VV]{o} += #2314# - W_JBME[OV][OV]{o} 2@2 t2[OV][OV]{o}");

    wfn_->blas()->append("t2_eqns[OO][VV]{o} += #1423# - W_JBME[OV][OV]{o} 2@2 t2[OV][OV]{o}");
    wfn_->blas()->append("t2_eqns[OO][VV]{o} += #1324#   W_JBME[OV][OV]{o} 2@2 t2[OV][OV]{o}");

    wfn_->blas()->append("t2_eqns[OO][VV]{o} += #4213# - ([ov]:[vo]) 1@2 t1t1_IAME[OV][OV]{o}");
    wfn_->blas()->append("t2_eqns[OO][VV]{o} += #3214# + ([ov]:[vo]) 1@2 t1t1_IAME[OV][OV]{o}");
    wfn_->blas()->append("t2_eqns[OO][VV]{o} += #4123# + ([ov]:[vo]) 1@2 t1t1_IAME[OV][OV]{o}");
    wfn_->blas()->append("t2_eqns[OO][VV]{o} += #3124# - ([ov]:[vo]) 1@2 t1t1_IAME[OV][OV]{o}");

    wfn_->blas()->append("t2_eqns[OO][VV]{o} += #2413#   W_JBme[OV][ov]{o} 2@1 t2[ov][OV]{o}");
    wfn_->blas()->append("t2_eqns[OO][VV]{o} += #2314# - W_JBme[OV][ov]{o} 2@1 t2[ov][OV]{o}");

    wfn_->blas()->append("t2_eqns[OO][VV]{o} += #1423# - W_JBme[OV][ov]{o} 2@1 t2[ov][OV]{o}");
    wfn_->blas()->append("t2_eqns[OO][VV]{o} += #1324#   W_JBme[OV][ov]{o} 2@1 t2[ov][OV]{o}");

    wfn_->blas()->append("t2_eqns[OO][VV]{o} += #1234#   t1[O][V]{o} 2@1 <[v]:[ovv]>");
    wfn_->blas()->append("t2_eqns[OO][VV]{o} += #2134# - t1[O][V]{o} 2@1 <[v]:[ovv]>");

    wfn_->blas()->append("t2_eqns[OO][VV]{o} += #3412# - t1[O][V]{o} 1@1 <[o]:[voo]>");
    wfn_->blas()->append("t2_eqns[OO][VV]{o} += #4312#   t1[O][V]{o} 1@1 <[o]:[voo]>");
}

void CCMRCC::build_t2_amplitudes_triples() {
    build_t2_ijab_amplitudes_triples_diagram1();
    build_t2_iJaB_amplitudes_triples_diagram1();
    build_t2_IJAB_amplitudes_triples_diagram1();

    build_t2_ijab_amplitudes_triples_diagram2();
    build_t2_iJaB_amplitudes_triples_diagram2();
    build_t2_IJAB_amplitudes_triples_diagram2();

    build_t2_ijab_amplitudes_triples_diagram3();
    build_t2_iJaB_amplitudes_triples_diagram3();
    build_t2_IJAB_amplitudes_triples_diagram3();
}

/**
 * @brief Computes the contraction
 * \f[ -\frac{1}{2} P(ij)\sum_{mne} t_{imn}^{abe} W_{mnje} \f]
 * as described in Fig. 1 of Phys. Chem. Chem. Phys. vol. 2, pg. 2047 (2000).
 * This algorithm hence performs
 * \f[ \frac{1}{2}\sum_{mne} t_{imn}^{abe} W_{mnje} \rightarrow  \{ -\bar{H}_{ij}^{ab},\bar{H}_{ji}^{ab} \} \f]
 * \f[ \sum_{mNE} t_{imN}^{abE} W_{mNjE} \rightarrow  \{ -\bar{H}_{ij}^{ab},\bar{H}_{ji}^{ab} \} \f]
 */
void CCMRCC::build_t2_ijab_amplitudes_triples_diagram1() {
    // Loop over references
    for (int ref = 0; ref < wfn_->moinfo()->get_nunique(); ref++) {
        int unique_ref = wfn_->moinfo()->get_ref_number(ref, UniqueRefs);

        // Grab the temporary matrices
        auto TijkabcMatTmp = wfn_->blas()->get_MatTmp("t3[ooo][vvv]", unique_ref, none);
        auto Tijkabc_matrix = TijkabcMatTmp->get_matrix();
        auto TijKabCMatTmp = wfn_->blas()->get_MatTmp("t3[ooO][vvV]", unique_ref, none);
        auto TijKabC_matrix = TijKabCMatTmp->get_matrix();

        auto WkijaMatTmp = wfn_->blas()->get_MatTmp("W_kija[o][oov]", unique_ref, none);
        auto WkiJAMatTmp = wfn_->blas()->get_MatTmp("W_kiJA[o][oOV]", unique_ref, none);
        auto Wkija_matrix = WkijaMatTmp->get_matrix();
        auto WkiJA_matrix = WkiJAMatTmp->get_matrix();

        auto HijabMatTmp = wfn_->blas()->get_MatTmp("t2_eqns[oo][vv]", unique_ref, none);

        // Grab the indexing for t3[iab][jkc]
        auto iab_indexing = wfn_->blas()->get_index("[ovv]");
        auto jkc_indexing = wfn_->blas()->get_index("[oov]");
        auto j_indexing = wfn_->blas()->get_index("[o]");

        auto& iab_tuples = iab_indexing->get_tuples();
        auto& jkc_tuples = jkc_indexing->get_tuples();

        // PART A: Sort T[ijk][abc]->T[iab][jkc]
        std::vector<double**> T_iabjkc(wfn_->nirrep(), nullptr);
        std::vector<double**> H_iabj(wfn_->nirrep(), nullptr);

        for (int h = 0; h < wfn_->moinfo()->get_nirreps(); h++) {
            // Allocate a block of T_iabjkc
            T_iabjkc[h] = block_matrix(iab_indexing->get_pairpi(h), jkc_indexing->get_pairpi(h));
            H_iabj[h] = block_matrix(iab_indexing->get_pairpi(h), j_indexing->get_pairpi(h));

            size_t iab_offset = iab_indexing->get_first(h);
            size_t jkc_offset = jkc_indexing->get_first(h);

            // AAA Contribution

            // Sort this block
            for (int iab = 0; iab < iab_indexing->get_pairpi(h); iab++) {
                int i = iab_tuples[iab_offset + iab][0];
                int a = iab_tuples[iab_offset + iab][1];
                int b = iab_tuples[iab_offset + iab][2];
                for (int jkc = 0; jkc < jkc_indexing->get_pairpi(h); jkc++) {
                    int j = jkc_tuples[jkc_offset + jkc][0];
                    int k = jkc_tuples[jkc_offset + jkc][1];
                    int c = jkc_tuples[jkc_offset + jkc][2];
                    T_iabjkc[h][iab][jkc] = TijkabcMatTmp->get_six_address_element(i, j, k, a, b, c);
                }
            }
            int m = iab_indexing->get_pairpi(h);
            int n = j_indexing->get_pairpi(h);
            int k = jkc_indexing->get_pairpi(h);
            if (m * n * k) {
                double alpha = 1.0;
                double beta = 0.0;
                C_DGEMM_22(m, n, k, alpha, &(T_iabjkc[h][0][0]), k, &(Wkija_matrix[h][0][0]), k, beta,
                           &(H_iabj[h][0][0]), n);
            }

            // Sort the result
            size_t j_offset = j_indexing->get_first(h);
            for (int iab = 0; iab < iab_indexing->get_pairpi(h); iab++) {
                int i = iab_tuples[iab_offset + iab][0];
                int a = iab_tuples[iab_offset + iab][1];
                int b = iab_tuples[iab_offset + iab][2];
                for (int j = 0; j < j_indexing->get_pairpi(h); j++) {
                    int abs_j = j + j_offset;
                    HijabMatTmp->add_four_address_element(i, abs_j, a, b, -0.5 * H_iabj[h][iab][j]);
                    HijabMatTmp->add_four_address_element(abs_j, i, a, b, 0.5 * H_iabj[h][iab][j]);
                }
            }

            // AAB Contribution

            // Sort this block
            for (int iab = 0; iab < iab_indexing->get_pairpi(h); iab++) {
                int i = iab_tuples[iab_offset + iab][0];
                int a = iab_tuples[iab_offset + iab][1];
                int b = iab_tuples[iab_offset + iab][2];
                for (int jkc = 0; jkc < jkc_indexing->get_pairpi(h); jkc++) {
                    int j = jkc_tuples[jkc_offset + jkc][0];
                    int k = jkc_tuples[jkc_offset + jkc][1];
                    int c = jkc_tuples[jkc_offset + jkc][2];
                    T_iabjkc[h][iab][jkc] = TijKabCMatTmp->get_six_address_element(i, j, k, a, b, c);
                }
            }
            m = iab_indexing->get_pairpi(h);
            n = j_indexing->get_pairpi(h);
            k = jkc_indexing->get_pairpi(h);
            if (m * n * k) {
                double alpha = 1.0;
                double beta = 0.0;
                C_DGEMM_22(m, n, k, alpha, &(T_iabjkc[h][0][0]), k, &(WkiJA_matrix[h][0][0]), k, beta,
                           &(H_iabj[h][0][0]), n);
            }

            // Sort the result
            j_offset = j_indexing->get_first(h);
            for (int iab = 0; iab < iab_indexing->get_pairpi(h); iab++) {
                int i = iab_tuples[iab_offset + iab][0];
                int a = iab_tuples[iab_offset + iab][1];
                int b = iab_tuples[iab_offset + iab][2];
                for (int j = 0; j < j_indexing->get_pairpi(h); j++) {
                    int abs_j = j + j_offset;
                    HijabMatTmp->add_four_address_element(i, abs_j, a, b, -H_iabj[h][iab][j]);
                    HijabMatTmp->add_four_address_element(abs_j, i, a, b, H_iabj[h][iab][j]);
                }
            }

            // Deallocate the memory for the block
            free_block(H_iabj[h]);
            free_block(T_iabjkc[h]);
        }
    }
}

/**
 * @brief Computes the contraction
 * \f[ -\frac{1}{2} P(ij)\sum_{mne} t_{imn}^{aBe} W_{mnJe} \f]
 * as described in Fig. 1 of Phys. Chem. Chem. Phys. vol. 2, pg. 2047 (2000).
 * This algorithm hence performs
 * \f[ \sum_{Mne} t_{inM}^{aeB} W_{MnJe} \rightarrow  \{ -\bar{H}_{iJ}^{aB} \} \f]
 * \f[ \sum_{mNe} t_{mNJ}^{aBE} W_{mNiE} \rightarrow  \{ \bar{H}_{iJ}^{aB} \} \f]
 * \f[ \frac{1}{2} \sum_{mne} t_{mnJ}^{aeB} W_{mnie} \rightarrow  \{ -\bar{H}_{iJ}^{aB} \} \f]
 * \f[ \frac{1}{2} \sum_{NME} t_{iNM}^{aBE} W_{MNJE} \rightarrow  \{ \bar{H}_{iJ}^{aB} \} \f]
 */
void CCMRCC::build_t2_iJaB_amplitudes_triples_diagram1() {
    // Loop over references
    for (int ref = 0; ref < wfn_->moinfo()->get_nunique(); ref++) {
        int unique_ref = wfn_->moinfo()->get_ref_number(ref, UniqueRefs);

        // Grab the temporary matrices
        CCMatTmp TijKabCMatTmp = wfn_->blas()->get_MatTmp("t3[ooO][vvV]", unique_ref, none);
        auto TijKabC_matrix = TijKabCMatTmp->get_matrix();
        CCMatTmp TiJKaBCMatTmp = wfn_->blas()->get_MatTmp("t3[oOO][vVV]", unique_ref, none);
        auto TiJKaBC_matrix = TiJKaBCMatTmp->get_matrix();

        CCMatTmp WkijaMatTmp = wfn_->blas()->get_MatTmp("W_kija[o][oov]", unique_ref, none);
        auto Wkija_matrix = WkijaMatTmp->get_matrix();
        CCMatTmp WKIjaMatTmp = wfn_->blas()->get_MatTmp("W_KIja[O][Oov]", unique_ref, none);
        auto WKIja_matrix = WKIjaMatTmp->get_matrix();
        CCMatTmp WkiJAMatTmp = wfn_->blas()->get_MatTmp("W_kiJA[o][oOV]", unique_ref, none);
        auto WkiJA_matrix = WkiJAMatTmp->get_matrix();
        CCMatTmp WKIJAMatTmp = wfn_->blas()->get_MatTmp("W_KIJA[O][OOV]", unique_ref, none);
        auto WKIJA_matrix = WKIJAMatTmp->get_matrix();

        CCMatTmp HiJaBMatTmp = wfn_->blas()->get_MatTmp("t2_eqns[oO][vV]", unique_ref, none);

        // Grab the indexing for t3[iab][jkc]
        auto jab_indexing = wfn_->blas()->get_index("[ovv]");
        auto kac_indexing = wfn_->blas()->get_index("[ovv]");
        auto iab_indexing = wfn_->blas()->get_index("[ovv]");
        auto iac_indexing = wfn_->blas()->get_index("[ovv]");
        auto kab_indexing = wfn_->blas()->get_index("[ovv]");
        auto kjb_indexing = wfn_->blas()->get_index("[oov]");
        auto kjc_indexing = wfn_->blas()->get_index("[oov]");
        auto ijc_indexing = wfn_->blas()->get_index("[oov]");
        auto ijb_indexing = wfn_->blas()->get_index("[oov]");
        auto j_indexing = wfn_->blas()->get_index("[o]");
        auto i_indexing = wfn_->blas()->get_index("[o]");

        auto& iab_tuples = iab_indexing->get_tuples();
        auto& jab_tuples = jab_indexing->get_tuples();
        auto& iac_tuples = iac_indexing->get_tuples();
        auto& kab_tuples = kab_indexing->get_tuples();
        auto& kac_tuples = kac_indexing->get_tuples();

        auto& kjb_tuples = kjb_indexing->get_tuples();
        auto& kjc_tuples = kjc_indexing->get_tuples();
        auto& ijc_tuples = ijc_indexing->get_tuples();
        auto& ijb_tuples = ijb_indexing->get_tuples();

        std::vector<double**> T_iackjb(wfn_->nirrep(), nullptr);
        std::vector<double**> H_iabj(wfn_->nirrep());

        for (int h = 0; h < wfn_->moinfo()->get_nirreps(); h++) {
            // Allocate a block of T_iabjkc
            T_iackjb[h] = block_matrix(iac_indexing->get_pairpi(h), kjb_indexing->get_pairpi(h));
            H_iabj[h] = block_matrix(iab_indexing->get_pairpi(h), j_indexing->get_pairpi(h));

            size_t iab_offset = iab_indexing->get_first(h);
            size_t iac_offset = iac_indexing->get_first(h);
            size_t kjb_offset = kjb_indexing->get_first(h);

            // AAA Contribution

            // PART A: Sort T[ijk][abc]->T[iac][kjb]
            for (int iac = 0; iac < iac_indexing->get_pairpi(h); iac++) {
                int i = iac_tuples[iac_offset + iac][0];
                int a = iac_tuples[iac_offset + iac][1];
                int c = iac_tuples[iac_offset + iac][2];
                for (int kjb = 0; kjb < kjb_indexing->get_pairpi(h); kjb++) {
                    int k = kjb_tuples[kjb_offset + kjb][0];
                    int j = kjb_tuples[kjb_offset + kjb][1];
                    int b = kjb_tuples[kjb_offset + kjb][2];
                    T_iackjb[h][iac][kjb] = TijKabCMatTmp->get_six_address_element(i, j, k, a, b, c);
                }
            }
            int m = iac_indexing->get_pairpi(h);
            int n = j_indexing->get_pairpi(h);
            int k = kjb_indexing->get_pairpi(h);
            if (m * n * k) {
                double alpha = 1.0;
                double beta = 0.0;
                C_DGEMM_22(m, n, k, alpha, &(T_iackjb[h][0][0]), k, &(WKIja_matrix[h][0][0]), k, beta,
                           &(H_iabj[h][0][0]), n);
            }

            // Sort the result
            size_t j_offset = j_indexing->get_first(h);
            for (int iab = 0; iab < iab_indexing->get_pairpi(h); iab++) {
                int i = iab_tuples[iab_offset + iab][0];
                int a = iab_tuples[iab_offset + iab][1];
                int b = iab_tuples[iab_offset + iab][2];
                for (int j = 0; j < j_indexing->get_pairpi(h); j++) {
                    int abs_j = j + j_offset;
                    HiJaBMatTmp->add_four_address_element(i, abs_j, a, b, -H_iabj[h][iab][j]);
                }
            }

            size_t kab_offset = kab_indexing->get_first(h);
            size_t ijc_offset = ijc_indexing->get_first(h);

            // ABB Contribution

            // Sort this block
            for (int kab = 0; kab < kab_indexing->get_pairpi(h); kab++) {
                int k = kab_tuples[kab_offset + kab][0];
                int a = kab_tuples[kab_offset + kab][1];
                int b = kab_tuples[kab_offset + kab][2];
                for (int ijc = 0; ijc < ijc_indexing->get_pairpi(h); ijc++) {
                    int i = ijc_tuples[ijc_offset + ijc][0];
                    int j = ijc_tuples[ijc_offset + ijc][1];
                    int c = ijc_tuples[ijc_offset + ijc][2];
                    T_iackjb[h][kab][ijc] = TiJKaBCMatTmp->get_six_address_element(i, j, k, a, b, c);
                }
            }

            m = kab_indexing->get_pairpi(h);
            n = i_indexing->get_pairpi(h);
            k = ijc_indexing->get_pairpi(h);
            if (m * n * k) {
                double alpha = 1.0;
                double beta = 0.0;
                C_DGEMM_22(m, n, k, alpha, &(T_iackjb[h][0][0]), k, &(WkiJA_matrix[h][0][0]), k, beta,
                           &(H_iabj[h][0][0]), n);
            }

            // Sort the result
            size_t jab_offset = jab_indexing->get_first(h);
            size_t i_offset = i_indexing->get_first(h);
            for (int jab = 0; jab < jab_indexing->get_pairpi(h); jab++) {
                int j = jab_tuples[jab_offset + jab][0];
                int a = jab_tuples[jab_offset + jab][1];
                int b = jab_tuples[jab_offset + jab][2];
                for (int i = 0; i < i_indexing->get_pairpi(h); i++) {
                    int abs_i = i + i_offset;
                    HiJaBMatTmp->add_four_address_element(abs_i, j, a, b, H_iabj[h][jab][i]);
                }
            }

            // t[mnJ][aeB] W[i][mne] Contribution

            size_t kac_offset = kac_indexing->get_first(h);
            size_t ijb_offset = ijb_indexing->get_first(h);

            // Sort this block
            for (int kac = 0; kac < kac_indexing->get_pairpi(h); kac++) {
                int k = kac_tuples[kac_offset + kac][0];
                int a = kac_tuples[kac_offset + kac][1];
                int c = kac_tuples[kac_offset + kac][2];
                for (int ijb = 0; ijb < ijb_indexing->get_pairpi(h); ijb++) {
                    int i = ijb_tuples[ijb_offset + ijb][0];
                    int j = ijb_tuples[ijb_offset + ijb][1];
                    int b = ijb_tuples[ijb_offset + ijb][2];
                    T_iackjb[h][kac][ijb] = TijKabCMatTmp->get_six_address_element(i, j, k, a, b, c);
                }
            }

            m = kac_indexing->get_pairpi(h);
            n = i_indexing->get_pairpi(h);
            k = ijb_indexing->get_pairpi(h);
            if (m * n * k) {
                double alpha = 1.0;
                double beta = 0.0;
                C_DGEMM_22(m, n, k, alpha, &(T_iackjb[h][0][0]), k, &(Wkija_matrix[h][0][0]), k, beta,
                           &(H_iabj[h][0][0]), n);
            }

            // Sort the result
            for (int jab = 0; jab < jab_indexing->get_pairpi(h); jab++) {
                int j = jab_tuples[jab_offset + jab][0];
                int a = jab_tuples[jab_offset + jab][1];
                int b = jab_tuples[jab_offset + jab][2];
                for (int i = 0; i < i_indexing->get_pairpi(h); i++) {
                    int abs_i = i + i_offset;
                    HiJaBMatTmp->add_four_address_element(abs_i, j, a, b, -0.5 * H_iabj[h][jab][i]);
                }
            }

            // t[iaB][MNE] W[J][MNE] Contribution
            size_t kjc_offset = kjc_indexing->get_first(h);

            // Sort this block
            for (int iab = 0; iab < iab_indexing->get_pairpi(h); iab++) {
                int i = iab_tuples[iab_offset + iab][0];
                int a = iab_tuples[iab_offset + iab][1];
                int b = iab_tuples[iab_offset + iab][2];
                for (int kjc = 0; kjc < kjc_indexing->get_pairpi(h); kjc++) {
                    int k = kjc_tuples[kjc_offset + kjc][0];
                    int j = kjc_tuples[kjc_offset + kjc][1];
                    int c = kjc_tuples[kjc_offset + kjc][2];
                    T_iackjb[h][iab][kjc] = TiJKaBCMatTmp->get_six_address_element(i, j, k, a, b, c);
                }
            }

            m = iab_indexing->get_pairpi(h);
            n = j_indexing->get_pairpi(h);
            k = kjc_indexing->get_pairpi(h);
            if (m * n * k) {
                double alpha = 1.0;
                double beta = 0.0;
                C_DGEMM_22(m, n, k, alpha, &(T_iackjb[h][0][0]), k, &(WKIJA_matrix[h][0][0]), k, beta,
                           &(H_iabj[h][0][0]), n);
            }

            // Sort the result
            for (int iab = 0; iab < iab_indexing->get_pairpi(h); iab++) {
                int i = iab_tuples[iab_offset + iab][0];
                int a = iab_tuples[iab_offset + iab][1];
                int b = iab_tuples[iab_offset + iab][2];
                for (int j = 0; j < j_indexing->get_pairpi(h); j++) {
                    int abs_j = j + j_offset;
                    HiJaBMatTmp->add_four_address_element(i, abs_j, a, b, 0.5 * H_iabj[h][iab][j]);
                }
            }
            // Deallocate the memory for the block
            free_block(T_iackjb[h]);
            free_block(H_iabj[h]);
        }
    }
}

/**
 * @brief Computes the contraction
 * \f[ -\frac{1}{2} P(IJ)\sum_{mne} t_{Imn}^{ABe} W_{mnJe} \f]
 * as described in Fig. 1 of Phys. Chem. Chem. Phys. vol. 2, pg. 2047 (2000).
 * This algorithm hence performs
 * \f[ \sum_{Mne} t_{nMI}^{eAB} W_{MnJe} \rightarrow  \{\bar{H}_{IJ}^{AB},-\bar{H}_{JI}^{AB} \} \f]
 * \f[ \frac{1}{2}\sum_{MNE} t_{IMN}^{ABE} W_{MNJE} \rightarrow  \{ -\bar{H}_{IJ}^{AB},\bar{H}_{JI}^{AB} \} \f]
 */
void CCMRCC::build_t2_IJAB_amplitudes_triples_diagram1() {
    // Loop over references
    for (int ref = 0; ref < wfn_->moinfo()->get_nunique(); ref++) {
        int unique_ref = wfn_->moinfo()->get_ref_number(ref, UniqueRefs);

        // Grab the temporary matrices
        auto TiJKaBCMatTmp = wfn_->blas()->get_MatTmp("t3[oOO][vVV]", unique_ref, none);
        auto TiJKaBC_matrix = TiJKaBCMatTmp->get_matrix();
        auto TIJKABCMatTmp = wfn_->blas()->get_MatTmp("t3[OOO][VVV]", unique_ref, none);
        auto TIJKABC_matrix = TIJKABCMatTmp->get_matrix();

        auto WKIjaMatTmp = wfn_->blas()->get_MatTmp("W_KIja[O][Oov]", unique_ref, none);
        auto WKIJAMatTmp = wfn_->blas()->get_MatTmp("W_KIJA[O][OOV]", unique_ref, none);
        auto WKIja_matrix = WKIjaMatTmp->get_matrix();
        auto WKIJA_matrix = WKIJAMatTmp->get_matrix();

        CCMatTmp HIJABMatTmp = wfn_->blas()->get_MatTmp("t2_eqns[OO][VV]", unique_ref, none);

        // Grab the indexing for t3[iab][jkc]
        auto iab_indexing = wfn_->blas()->get_index("[ovv]");
        auto kbc_indexing = wfn_->blas()->get_index("[ovv]");
        auto jia_indexing = wfn_->blas()->get_index("[oov]");
        auto jkc_indexing = wfn_->blas()->get_index("[oov]");
        auto j_indexing = wfn_->blas()->get_index("[o]");

        auto& iab_tuples = iab_indexing->get_tuples();
        auto& kbc_tuples = kbc_indexing->get_tuples();
        auto& jia_tuples = jia_indexing->get_tuples();
        auto& jkc_tuples = jkc_indexing->get_tuples();

        std::vector<double**> T_iabjkc(wfn_->nirrep(), nullptr);
        std::vector<double**> H_iabj(wfn_->nirrep(), nullptr);

        for (int h = 0; h < wfn_->moinfo()->get_nirreps(); h++) {
            // Allocate a block of T_iabjkc
            T_iabjkc[h] = block_matrix(kbc_indexing->get_pairpi(h), jia_indexing->get_pairpi(h));
            H_iabj[h] = block_matrix(kbc_indexing->get_pairpi(h), j_indexing->get_pairpi(h));

            size_t iab_offset = iab_indexing->get_first(h);
            size_t kbc_offset = kbc_indexing->get_first(h);
            size_t jia_offset = jia_indexing->get_first(h);
            size_t jkc_offset = jkc_indexing->get_first(h);

            // AAA Contribution

            // PART A: Sort T[ijk][abc]->T[kbc][jia]
            // Sort this block
            for (int kbc = 0; kbc < kbc_indexing->get_pairpi(h); kbc++) {
                int k = kbc_tuples[kbc_offset + kbc][0];
                int b = kbc_tuples[kbc_offset + kbc][1];
                int c = kbc_tuples[kbc_offset + kbc][2];
                for (int jia = 0; jia < jia_indexing->get_pairpi(h); jia++) {
                    int j = jia_tuples[jia_offset + jia][0];
                    int i = jia_tuples[jia_offset + jia][1];
                    int a = jia_tuples[jia_offset + jia][2];
                    T_iabjkc[h][kbc][jia] = TiJKaBCMatTmp->get_six_address_element(i, j, k, a, b, c);
                }
            }
            int m = kbc_indexing->get_pairpi(h);
            int n = j_indexing->get_pairpi(h);
            int k = jia_indexing->get_pairpi(h);
            if (m * n * k) {
                double alpha = 1.0;
                double beta = 0.0;
                C_DGEMM_22(m, n, k, alpha, &(T_iabjkc[h][0][0]), k, &(WKIja_matrix[h][0][0]), k, beta,
                           &(H_iabj[h][0][0]), n);
            }

            // Sort the result
            size_t j_offset = j_indexing->get_first(h);
            for (int iab = 0; iab < iab_indexing->get_pairpi(h); iab++) {
                int i = iab_tuples[iab_offset + iab][0];
                int a = iab_tuples[iab_offset + iab][1];
                int b = iab_tuples[iab_offset + iab][2];
                for (int j = 0; j < j_indexing->get_pairpi(h); j++) {
                    int abs_j = j + j_offset;
                    HIJABMatTmp->add_four_address_element(i, abs_j, a, b, H_iabj[h][iab][j]);
                    HIJABMatTmp->add_four_address_element(abs_j, i, a, b, -H_iabj[h][iab][j]);
                }
            }

            // BBB Contribution

            // Sort this block
            for (int iab = 0; iab < iab_indexing->get_pairpi(h); iab++) {
                int i = iab_tuples[iab_offset + iab][0];
                int a = iab_tuples[iab_offset + iab][1];
                int b = iab_tuples[iab_offset + iab][2];
                for (int jkc = 0; jkc < jkc_indexing->get_pairpi(h); jkc++) {
                    int j = jkc_tuples[jkc_offset + jkc][0];
                    int k = jkc_tuples[jkc_offset + jkc][1];
                    int c = jkc_tuples[jkc_offset + jkc][2];
                    T_iabjkc[h][iab][jkc] = TIJKABCMatTmp->get_six_address_element(i, j, k, a, b, c);
                }
            }
            m = iab_indexing->get_pairpi(h);
            n = j_indexing->get_pairpi(h);
            k = jkc_indexing->get_pairpi(h);
            if (m * n * k) {
                double alpha = 1.0;
                double beta = 0.0;
                C_DGEMM_22(m, n, k, alpha, &(T_iabjkc[h][0][0]), k, &(WKIJA_matrix[h][0][0]), k, beta,
                           &(H_iabj[h][0][0]), n);
            }

            // Sort the result
            for (int iab = 0; iab < iab_indexing->get_pairpi(h); iab++) {
                int i = iab_tuples[iab_offset + iab][0];
                int a = iab_tuples[iab_offset + iab][1];
                int b = iab_tuples[iab_offset + iab][2];
                for (int j = 0; j < j_indexing->get_pairpi(h); j++) {
                    int abs_j = j + j_offset;
                    HIJABMatTmp->add_four_address_element(i, abs_j, a, b, -0.5 * H_iabj[h][iab][j]);
                    HIJABMatTmp->add_four_address_element(abs_j, i, a, b, 0.5 * H_iabj[h][iab][j]);
                }
            }

            // Deallocate the memory for the block
            free_block(H_iabj[h]);
            free_block(T_iabjkc[h]);
        }
    }
}

/**
 * @brief Computes the contraction
 * \f[ -\frac{1}{2} P(ij)\sum_{mne} t_{imn}^{abe} W_{mnje} \f]
 * as described in Fig. 1 of Phys. Chem. Chem. Phys. vol. 2, pg. 2047 (2000).
 * This algorithm hence performs
 * \f[ \frac{1}{2}\sum_{mne} t_{imn}^{abe} W_{mnje} \rightarrow  \{ -\bar{H}_{ij}^{ab},\bar{H}_{ji}^{ab} \} \f]
 * \f[ \sum_{mNE} t_{imN}^{abE} W_{mNjE} \rightarrow  \{ -\bar{H}_{ij}^{ab},\bar{H}_{ji}^{ab} \} \f]
 */
void CCMRCC::build_t2_ijab_amplitudes_triples_diagram2() {
    // Loop over references
    for (int ref = 0; ref < wfn_->moinfo()->get_nunique(); ref++) {
        int unique_ref = wfn_->moinfo()->get_ref_number(ref, UniqueRefs);

        // Grab the temporary matrices
        auto TijkabcMatTmp = wfn_->blas()->get_MatTmp("t3[ooo][vvv]", unique_ref, none);
        auto Tijkabc_matrix = TijkabcMatTmp->get_matrix();
        auto TijKabCMatTmp = wfn_->blas()->get_MatTmp("t3[ooO][vvV]", unique_ref, none);
        auto TijKabC_matrix = TijKabCMatTmp->get_matrix();

        auto WaibcMatTmp = wfn_->blas()->get_MatTmp("W_aibc[v][ovv]", unique_ref, none);
        auto Waibc_matrix = WaibcMatTmp->get_matrix();
        auto WaIbCMatTmp = wfn_->blas()->get_MatTmp("W_aIbC[v][OvV]", unique_ref, none);
        auto WaIbC_matrix = WaIbCMatTmp->get_matrix();

        auto HijabMatTmp = wfn_->blas()->get_MatTmp("t2_eqns[oo][vv]", unique_ref, none);

        // Grab the indexing for t3[iab][jkc]
        auto ovv_indexing = wfn_->blas()->get_index("[ovv]");
        auto oov_indexing = wfn_->blas()->get_index("[oov]");
        auto v_indexing = wfn_->blas()->get_index("[v]");

        auto& ovv_tuples = ovv_indexing->get_tuples();
        auto& oov_tuples = oov_indexing->get_tuples();

        std::vector<double**> T_oovovv(wfn_->nirrep(), nullptr);
        std::vector<double**> H_ijab(wfn_->nirrep(), nullptr);

        for (int h = 0; h < wfn_->nirrep(); h++) {
            // Allocate a block of T_iabjkc
            T_oovovv[h] = block_matrix(oov_indexing->get_pairpi(h), ovv_indexing->get_pairpi(h));
            H_ijab[h] = block_matrix(oov_indexing->get_pairpi(h), v_indexing->get_pairpi(h));

            size_t ovv_offset = ovv_indexing->get_first(h);
            size_t oov_offset = oov_indexing->get_first(h);

            // AAA Contribution

            // Sort this block
            for (int kbc = 0; kbc < ovv_indexing->get_pairpi(h); kbc++) {
                int k = ovv_tuples[ovv_offset + kbc][0];
                int b = ovv_tuples[ovv_offset + kbc][1];
                int c = ovv_tuples[ovv_offset + kbc][2];
                for (int ija = 0; ija < oov_indexing->get_pairpi(h); ija++) {
                    int i = oov_tuples[oov_offset + ija][0];
                    int j = oov_tuples[oov_offset + ija][1];
                    int a = oov_tuples[oov_offset + ija][2];
                    T_oovovv[h][ija][kbc] = TijkabcMatTmp->get_six_address_element(i, j, k, a, b, c);
                }
            }
            int m = oov_indexing->get_pairpi(h);
            int n = v_indexing->get_pairpi(h);
            int k = ovv_indexing->get_pairpi(h);
            if (m * n * k) {
                double alpha = 1.0;
                double beta = 0.0;
                C_DGEMM_22(m, n, k, alpha, &(T_oovovv[h][0][0]), k, &(Waibc_matrix[h][0][0]), k, beta,
                           &(H_ijab[h][0][0]), n);
            }

            // Sort the result
            size_t v_offset = v_indexing->get_first(h);
            for (int ija = 0; ija < oov_indexing->get_pairpi(h); ija++) {
                int i = oov_tuples[oov_offset + ija][0];
                int j = oov_tuples[oov_offset + ija][1];
                int a = oov_tuples[oov_offset + ija][2];
                for (int b = 0; b < v_indexing->get_pairpi(h); b++) {
                    int abs_b = b + v_offset;
                    HijabMatTmp->add_four_address_element(i, j, a, abs_b, 0.5 * H_ijab[h][ija][b]);
                    HijabMatTmp->add_four_address_element(i, j, abs_b, a, -0.5 * H_ijab[h][ija][b]);
                }
            }

            // Sort this block
            for (int kbc = 0; kbc < ovv_indexing->get_pairpi(h); kbc++) {
                int k = ovv_tuples[ovv_offset + kbc][0];
                int b = ovv_tuples[ovv_offset + kbc][1];
                int c = ovv_tuples[ovv_offset + kbc][2];
                for (int ija = 0; ija < oov_indexing->get_pairpi(h); ija++) {
                    int i = oov_tuples[oov_offset + ija][0];
                    int j = oov_tuples[oov_offset + ija][1];
                    int a = oov_tuples[oov_offset + ija][2];
                    T_oovovv[h][ija][kbc] = TijKabCMatTmp->get_six_address_element(i, j, k, a, b, c);
                }
            }

            m = oov_indexing->get_pairpi(h);
            n = v_indexing->get_pairpi(h);
            k = ovv_indexing->get_pairpi(h);
            if (m * n * k) {
                double alpha = 1.0;
                double beta = 0.0;
                C_DGEMM_22(m, n, k, alpha, &(T_oovovv[h][0][0]), k, &(WaIbC_matrix[h][0][0]), k, beta,
                           &(H_ijab[h][0][0]), n);
            }

            // Sort the result
            v_offset = v_indexing->get_first(h);
            for (int ija = 0; ija < oov_indexing->get_pairpi(h); ija++) {
                int i = oov_tuples[oov_offset + ija][0];
                int j = oov_tuples[oov_offset + ija][1];
                int a = oov_tuples[oov_offset + ija][2];
                for (int b = 0; b < v_indexing->get_pairpi(h); b++) {
                    int abs_b = b + v_offset;
                    HijabMatTmp->add_four_address_element(i, j, a, abs_b, H_ijab[h][ija][b]);
                    HijabMatTmp->add_four_address_element(i, j, abs_b, a, -H_ijab[h][ija][b]);
                }
            }
            // Deallocate the memory for the block
            free_block(H_ijab[h]);
            free_block(T_oovovv[h]);
        }
    }
}

/**
 * @brief Computes the contraction
 * \f[ -\frac{1}{2} P(ij)\sum_{mne} t_{imn}^{abe} W_{mnje} \f]
 * as described in Fig. 1 of Phys. Chem. Chem. Phys. vol. 2, pg. 2047 (2000).
 * This algorithm hence performs
 * \f[ \frac{1}{2}\sum_{mne} t_{imn}^{abe} W_{mnje} \rightarrow  \{ -\bar{H}_{ij}^{ab},\bar{H}_{ji}^{ab} \} \f]
 * \f[ \sum_{mNE} t_{imN}^{abE} W_{mNjE} \rightarrow  \{ -\bar{H}_{ij}^{ab},\bar{H}_{ji}^{ab} \} \f]
 */
void CCMRCC::build_t2_iJaB_amplitudes_triples_diagram2() {
    // Loop over references
    for (int ref = 0; ref < wfn_->moinfo()->get_nunique(); ref++) {
        int unique_ref = wfn_->moinfo()->get_ref_number(ref, UniqueRefs);

        // Grab the temporary matrices
        auto TijKabCMatTmp = wfn_->blas()->get_MatTmp("t3[ooO][vvV]", unique_ref, none);
        auto TijKabC_matrix = TijKabCMatTmp->get_matrix();
        auto TiJKaBCMatTmp = wfn_->blas()->get_MatTmp("t3[oOO][vVV]", unique_ref, none);
        auto TiJKaBC_matrix = TiJKaBCMatTmp->get_matrix();

        auto WaibcMatTmp = wfn_->blas()->get_MatTmp("W_aibc[v][ovv]", unique_ref, none);
        auto Waibc_matrix = WaibcMatTmp->get_matrix();
        auto WaIbCMatTmp = wfn_->blas()->get_MatTmp("W_aIbC[v][OvV]", unique_ref, none);
        auto WaIbC_matrix = WaIbCMatTmp->get_matrix();
        auto WAiBcMatTmp = wfn_->blas()->get_MatTmp("W_AiBc[V][oVv]", unique_ref, none);
        auto WAiBc_matrix = WAiBcMatTmp->get_matrix();
        auto WAIBCMatTmp = wfn_->blas()->get_MatTmp("W_AIBC[V][OVV]", unique_ref, none);
        auto WAIBC_matrix = WAIBCMatTmp->get_matrix();

        auto HiJaBMatTmp = wfn_->blas()->get_MatTmp("t2_eqns[oO][vV]", unique_ref, none);

        // Grab the indexing for t3[iab][jkc]
        auto ovv_indexing = wfn_->blas()->get_index("[ovv]");
        auto oov_indexing = wfn_->blas()->get_index("[oov]");
        auto v_indexing = wfn_->blas()->get_index("[v]");

        auto& ovv_tuples = ovv_indexing->get_tuples();
        auto& oov_tuples = oov_indexing->get_tuples();

        std::vector<double**> T_oovovv(wfn_->nirrep(), nullptr);
        std::vector<double**> H_ijab(wfn_->nirrep(), nullptr);

        for (int h = 0; h < wfn_->moinfo()->get_nirreps(); h++) {
            // Allocate a block of T_iabjkc
            T_oovovv[h] = block_matrix(oov_indexing->get_pairpi(h), ovv_indexing->get_pairpi(h));
            H_ijab[h] = block_matrix(oov_indexing->get_pairpi(h), v_indexing->get_pairpi(h));

            size_t ovv_offset = ovv_indexing->get_first(h);
            size_t oov_offset = oov_indexing->get_first(h);

            // AAA Contribution

            // Sort this block
            for (int jab = 0; jab < ovv_indexing->get_pairpi(h); jab++) {
                int j = ovv_tuples[ovv_offset + jab][0];
                int a = ovv_tuples[ovv_offset + jab][1];
                int b = ovv_tuples[ovv_offset + jab][2];
                for (int ikc = 0; ikc < oov_indexing->get_pairpi(h); ikc++) {
                    int i = oov_tuples[oov_offset + ikc][0];
                    int k = oov_tuples[oov_offset + ikc][1];
                    int c = oov_tuples[oov_offset + ikc][2];
                    T_oovovv[h][ikc][jab] = TijKabCMatTmp->get_six_address_element(i, j, k, a, b, c);
                }
            }
            int m = oov_indexing->get_pairpi(h);
            int n = v_indexing->get_pairpi(h);
            int k = ovv_indexing->get_pairpi(h);
            if (m * n * k) {
                double alpha = 1.0;
                double beta = 0.0;
                C_DGEMM_22(m, n, k, alpha, &(T_oovovv[h][0][0]), k, &(Waibc_matrix[h][0][0]), k, beta,
                           &(H_ijab[h][0][0]), n);
            }

            // Sort the result
            size_t v_offset = v_indexing->get_first(h);
            for (int ijb = 0; ijb < oov_indexing->get_pairpi(h); ijb++) {
                int i = oov_tuples[oov_offset + ijb][0];
                int j = oov_tuples[oov_offset + ijb][1];
                int b = oov_tuples[oov_offset + ijb][2];
                for (int a = 0; a < v_indexing->get_pairpi(h); a++) {
                    int abs_a = a + v_offset;
                    HiJaBMatTmp->add_four_address_element(i, j, abs_a, b, 0.5 * H_ijab[h][ijb][a]);
                }
            }

            // Sort this block
            for (int kac = 0; kac < ovv_indexing->get_pairpi(h); kac++) {
                int k = ovv_tuples[ovv_offset + kac][0];
                int a = ovv_tuples[ovv_offset + kac][1];
                int c = ovv_tuples[ovv_offset + kac][2];
                for (int ijb = 0; ijb < oov_indexing->get_pairpi(h); ijb++) {
                    int i = oov_tuples[oov_offset + ijb][0];
                    int j = oov_tuples[oov_offset + ijb][1];
                    int b = oov_tuples[oov_offset + ijb][2];
                    T_oovovv[h][ijb][kac] = TiJKaBCMatTmp->get_six_address_element(i, j, k, a, b, c);
                }
            }
            m = oov_indexing->get_pairpi(h);
            n = v_indexing->get_pairpi(h);
            k = ovv_indexing->get_pairpi(h);
            if (m * n * k) {
                double alpha = 1.0;
                double beta = 0.0;
                C_DGEMM_22(m, n, k, alpha, &(T_oovovv[h][0][0]), k, &(WaIbC_matrix[h][0][0]), k, beta,
                           &(H_ijab[h][0][0]), n);
            }

            // Sort the result
            v_offset = v_indexing->get_first(h);
            for (int ijb = 0; ijb < oov_indexing->get_pairpi(h); ijb++) {
                int i = oov_tuples[oov_offset + ijb][0];
                int j = oov_tuples[oov_offset + ijb][1];
                int b = oov_tuples[oov_offset + ijb][2];
                for (int a = 0; a < v_indexing->get_pairpi(h); a++) {
                    int abs_a = a + v_offset;
                    HiJaBMatTmp->add_four_address_element(i, j, abs_a, b, H_ijab[h][ijb][a]);
                }
            }

            // Sort this block
            for (int jcb = 0; jcb < ovv_indexing->get_pairpi(h); jcb++) {
                int j = ovv_tuples[ovv_offset + jcb][0];
                int c = ovv_tuples[ovv_offset + jcb][1];
                int b = ovv_tuples[ovv_offset + jcb][2];
                for (int ika = 0; ika < oov_indexing->get_pairpi(h); ika++) {
                    int i = oov_tuples[oov_offset + ika][0];
                    int k = oov_tuples[oov_offset + ika][1];
                    int a = oov_tuples[oov_offset + ika][2];
                    T_oovovv[h][ika][jcb] = TijKabCMatTmp->get_six_address_element(i, j, k, a, b, c);
                }
            }
            m = oov_indexing->get_pairpi(h);
            n = v_indexing->get_pairpi(h);
            k = ovv_indexing->get_pairpi(h);
            if (m * n * k) {
                double alpha = 1.0;
                double beta = 0.0;
                C_DGEMM_22(m, n, k, alpha, &(T_oovovv[h][0][0]), k, &(WAiBc_matrix[h][0][0]), k, beta,
                           &(H_ijab[h][0][0]), n);
            }

            // Sort the result
            v_offset = v_indexing->get_first(h);
            for (int ija = 0; ija < oov_indexing->get_pairpi(h); ija++) {
                int i = oov_tuples[oov_offset + ija][0];
                int j = oov_tuples[oov_offset + ija][1];
                int a = oov_tuples[oov_offset + ija][2];
                for (int b = 0; b < v_indexing->get_pairpi(h); b++) {
                    int abs_b = b + v_offset;
                    HiJaBMatTmp->add_four_address_element(i, j, a, abs_b, H_ijab[h][ija][b]);
                }
            }

            // Sort this block
            for (int kbc = 0; kbc < ovv_indexing->get_pairpi(h); kbc++) {
                int k = ovv_tuples[ovv_offset + kbc][0];
                int b = ovv_tuples[ovv_offset + kbc][1];
                int c = ovv_tuples[ovv_offset + kbc][2];
                for (int ija = 0; ija < oov_indexing->get_pairpi(h); ija++) {
                    int i = oov_tuples[oov_offset + ija][0];
                    int j = oov_tuples[oov_offset + ija][1];
                    int a = oov_tuples[oov_offset + ija][2];
                    T_oovovv[h][ija][kbc] = TiJKaBCMatTmp->get_six_address_element(i, j, k, a, b, c);
                }
            }
            m = oov_indexing->get_pairpi(h);
            n = v_indexing->get_pairpi(h);
            k = ovv_indexing->get_pairpi(h);
            if (m * n * k) {
                double alpha = 1.0;
                double beta = 0.0;
                C_DGEMM_22(m, n, k, alpha, &(T_oovovv[h][0][0]), k, &(WAIBC_matrix[h][0][0]), k, beta,
                           &(H_ijab[h][0][0]), n);
            }

            // Sort the result
            v_offset = v_indexing->get_first(h);
            for (int ija = 0; ija < oov_indexing->get_pairpi(h); ija++) {
                int i = oov_tuples[oov_offset + ija][0];
                int j = oov_tuples[oov_offset + ija][1];
                int a = oov_tuples[oov_offset + ija][2];
                for (int b = 0; b < v_indexing->get_pairpi(h); b++) {
                    int abs_b = b + v_offset;
                    HiJaBMatTmp->add_four_address_element(i, j, a, abs_b, 0.5 * H_ijab[h][ija][b]);
                }
            }

            // Deallocate the memory for the block
            free_block(H_ijab[h]);
            free_block(T_oovovv[h]);
        }
    }
}

/**
 * @brief Computes the contraction
 * \f[ -\frac{1}{2} P(ij)\sum_{mne} t_{imn}^{abe} W_{mnje} \f]
 * as described in Fig. 1 of Phys. Chem. Chem. Phys. vol. 2, pg. 2047 (2000).
 * This algorithm hence performs
 * \f[ \frac{1}{2}\sum_{mne} t_{imn}^{abe} W_{mnje} \rightarrow  \{ -\bar{H}_{ij}^{ab},\bar{H}_{ji}^{ab} \} \f]
 * \f[ \sum_{mNE} t_{imN}^{abE} W_{mNjE} \rightarrow  \{ -\bar{H}_{ij}^{ab},\bar{H}_{ji}^{ab} \} \f]
 */
void CCMRCC::build_t2_IJAB_amplitudes_triples_diagram2() {
    // Loop over references
    for (int ref = 0; ref < wfn_->moinfo()->get_nunique(); ref++) {
        int unique_ref = wfn_->moinfo()->get_ref_number(ref, UniqueRefs);

        // Grab the temporary matrices
        auto TiJKaBCMatTmp = wfn_->blas()->get_MatTmp("t3[oOO][vVV]", unique_ref, none);
        auto TiJKaBC_matrix = TiJKaBCMatTmp->get_matrix();
        auto TIJKABCMatTmp = wfn_->blas()->get_MatTmp("t3[OOO][VVV]", unique_ref, none);
        auto TIJKABC_matrix = TIJKABCMatTmp->get_matrix();

        auto WAiBcMatTmp = wfn_->blas()->get_MatTmp("W_AiBc[V][oVv]", unique_ref, none);
        auto WAiBc_matrix = WAiBcMatTmp->get_matrix();
        auto WAIBCMatTmp = wfn_->blas()->get_MatTmp("W_AIBC[V][OVV]", unique_ref, none);
        auto WAIBC_matrix = WAIBCMatTmp->get_matrix();

        auto HIJABMatTmp = wfn_->blas()->get_MatTmp("t2_eqns[OO][VV]", unique_ref, none);

        // Grab the indexing for t3[iab][jkc]
        auto ovv_indexing = wfn_->blas()->get_index("[ovv]");
        auto oov_indexing = wfn_->blas()->get_index("[oov]");
        auto v_indexing = wfn_->blas()->get_index("[v]");

        auto& ovv_tuples = ovv_indexing->get_tuples();
        auto& oov_tuples = oov_indexing->get_tuples();

        std::vector<double**> T_oovovv(wfn_->nirrep(), nullptr);
        std::vector<double**> H_ijab(wfn_->nirrep());

        for (int h = 0; h < wfn_->moinfo()->get_nirreps(); h++) {
            // Allocate a block of T_iabjkc
            T_oovovv[h] = block_matrix(oov_indexing->get_pairpi(h), ovv_indexing->get_pairpi(h));
            H_ijab[h] = block_matrix(oov_indexing->get_pairpi(h), v_indexing->get_pairpi(h));

            size_t ovv_offset = ovv_indexing->get_first(h);
            size_t oov_offset = oov_indexing->get_first(h);

            // AAA Contribution

            // Sort this block
            for (int kbc = 0; kbc < ovv_indexing->get_pairpi(h); kbc++) {
                int k = ovv_tuples[ovv_offset + kbc][0];
                int b = ovv_tuples[ovv_offset + kbc][1];
                int c = ovv_tuples[ovv_offset + kbc][2];
                for (int ija = 0; ija < oov_indexing->get_pairpi(h); ija++) {
                    int i = oov_tuples[oov_offset + ija][0];
                    int j = oov_tuples[oov_offset + ija][1];
                    int a = oov_tuples[oov_offset + ija][2];
                    T_oovovv[h][ija][kbc] = TIJKABCMatTmp->get_six_address_element(i, j, k, a, b, c);
                }
            }
            int m = oov_indexing->get_pairpi(h);
            int n = v_indexing->get_pairpi(h);
            int k = ovv_indexing->get_pairpi(h);
            if (m * n * k) {
                double alpha = 1.0;
                double beta = 0.0;
                C_DGEMM_22(m, n, k, alpha, &(T_oovovv[h][0][0]), k, &(WAIBC_matrix[h][0][0]), k, beta,
                           &(H_ijab[h][0][0]), n);
            }

            // Sort the result
            size_t v_offset = v_indexing->get_first(h);
            for (int ija = 0; ija < oov_indexing->get_pairpi(h); ija++) {
                int i = oov_tuples[oov_offset + ija][0];
                int j = oov_tuples[oov_offset + ija][1];
                int a = oov_tuples[oov_offset + ija][2];
                for (int b = 0; b < v_indexing->get_pairpi(h); b++) {
                    int abs_b = b + v_offset;
                    HIJABMatTmp->add_four_address_element(i, j, a, abs_b, 0.5 * H_ijab[h][ija][b]);
                    HIJABMatTmp->add_four_address_element(i, j, abs_b, a, -0.5 * H_ijab[h][ija][b]);
                }
            }

            // Sort this block
            for (int ica = 0; ica < ovv_indexing->get_pairpi(h); ica++) {
                int i = ovv_tuples[ovv_offset + ica][0];
                int c = ovv_tuples[ovv_offset + ica][1];
                int a = ovv_tuples[ovv_offset + ica][2];
                for (int jkb = 0; jkb < oov_indexing->get_pairpi(h); jkb++) {
                    int j = oov_tuples[oov_offset + jkb][0];
                    int k = oov_tuples[oov_offset + jkb][1];
                    int b = oov_tuples[oov_offset + jkb][2];
                    T_oovovv[h][jkb][ica] = TiJKaBCMatTmp->get_six_address_element(i, j, k, a, b, c);
                }
            }

            m = oov_indexing->get_pairpi(h);
            n = v_indexing->get_pairpi(h);
            k = ovv_indexing->get_pairpi(h);
            if (m * n * k) {
                double alpha = 1.0;
                double beta = 0.0;
                C_DGEMM_22(m, n, k, alpha, &(T_oovovv[h][0][0]), k, &(WAiBc_matrix[h][0][0]), k, beta,
                           &(H_ijab[h][0][0]), n);
            }

            // Sort the result
            v_offset = v_indexing->get_first(h);
            for (int ija = 0; ija < oov_indexing->get_pairpi(h); ija++) {
                int i = oov_tuples[oov_offset + ija][0];
                int j = oov_tuples[oov_offset + ija][1];
                int a = oov_tuples[oov_offset + ija][2];
                for (int b = 0; b < v_indexing->get_pairpi(h); b++) {
                    int abs_b = b + v_offset;
                    HIJABMatTmp->add_four_address_element(i, j, a, abs_b, H_ijab[h][ija][b]);
                    HIJABMatTmp->add_four_address_element(i, j, abs_b, a, -H_ijab[h][ija][b]);
                }
            }
            // Deallocate the memory for the block
            free_block(H_ijab[h]);
            free_block(T_oovovv[h]);
        }
    }
}

/**
 * @brief Computes the contraction
 * \f[ \sum_{me} t_{ijm}^{abe} F''_{me} \f]
 * as described in Fig. 1 of Phys. Chem. Chem. Phys. vol. 2, pg. 2047 (2000).
 * This algorithm hence performs
 * \f[ \sum_{me} t_{ijm}^{abe} F''_{me} + \sum_{ME} t_{ijM}^{abE} F''_{ME} \rightarrow  \{ \bar{H}_{ij}^{ab} \} \f]
 */
void CCMRCC::build_t2_ijab_amplitudes_triples_diagram3() {
    // Loop over references
    for (int ref = 0; ref < wfn_->moinfo()->get_nunique(); ref++) {
        int unique_ref = wfn_->moinfo()->get_ref_number(ref, UniqueRefs);

        // Grab the temporary matrices
        auto HijabMatTmp = wfn_->blas()->get_MatTmp("t2_eqns[oo][vv]", unique_ref, none);
        auto TijkabcMatTmp = wfn_->blas()->get_MatTmp("t3[ooo][vvv]", unique_ref, none);
        auto TijKabCMatTmp = wfn_->blas()->get_MatTmp("t3[ooO][vvV]", unique_ref, none);
        auto FmeMatTmp = wfn_->blas()->get_MatTmp("F2_me[o][v]", unique_ref, none);
        auto FMEMatTmp = wfn_->blas()->get_MatTmp("F2_ME[O][V]", unique_ref, none);

        // Grab the indexing for t3[ijk][abc]
        auto& ij_tuples = HijabMatTmp->get_left()->get_tuples();
        auto& ab_tuples = HijabMatTmp->get_right()->get_tuples();
        auto& m_tuples = FmeMatTmp->get_left()->get_tuples();
        auto& e_tuples = FmeMatTmp->get_right()->get_tuples();

        auto Tijkabc_matrix = TijkabcMatTmp->get_matrix();
        auto TijKabC_matrix = TijKabCMatTmp->get_matrix();
        auto Hijab_matrix = HijabMatTmp->get_matrix();
        auto Fme_matrix = FmeMatTmp->get_matrix();
        auto FME_matrix = FMEMatTmp->get_matrix();
        auto ijkIndex = wfn_->blas()->get_index("[ooo]");
        auto abcIndex = wfn_->blas()->get_index("[vvv]");

        for (int h = 0; h < wfn_->moinfo()->get_nirreps(); h++) {
            size_t ij_offset = HijabMatTmp->get_left()->get_first(h);
            size_t ab_offset = HijabMatTmp->get_right()->get_first(h);
            for (int ab = 0; ab < HijabMatTmp->get_right_pairpi(h); ab++) {
                int a = ab_tuples[ab_offset + ab][0];
                int b = ab_tuples[ab_offset + ab][1];
                for (int ij = 0; ij < HijabMatTmp->get_left_pairpi(h); ij++) {
                    int i = ij_tuples[ij_offset + ij][0];
                    int j = ij_tuples[ij_offset + ij][1];
                    for (int m_sym = 0; m_sym < wfn_->moinfo()->get_nirreps(); m_sym++) {
                        size_t m_offset = FmeMatTmp->get_left()->get_first(m_sym);
                        size_t e_offset = FmeMatTmp->get_right()->get_first(m_sym);
                        for (int e = 0; e < FmeMatTmp->get_right_pairpi(m_sym); e++) {
                            int e_abs = e + e_offset;
                            size_t abe = abcIndex->get_tuple_rel_index(a, b, e_abs);
                            int abe_sym = abcIndex->get_tuple_irrep(a, b, e_abs);
                            for (int m = 0; m < FmeMatTmp->get_left_pairpi(m_sym); m++) {
                                int m_abs = m + m_offset;
                                size_t ijm = ijkIndex->get_tuple_rel_index(i, j, m_abs);
                                Hijab_matrix[h][ij][ab] += Tijkabc_matrix[abe_sym][ijm][abe] * Fme_matrix[m_sym][m][e];
                                Hijab_matrix[h][ij][ab] += TijKabC_matrix[abe_sym][ijm][abe] * FME_matrix[m_sym][m][e];
                            }
                        }
                    }
                }
            }
        }
    }
}

/**
 * @brief Computes the contraction
 * \f[ \sum_{me} t_{ijm}^{abe} F''_{me} \f]
 * as described in Fig. 1 of Phys. Chem. Chem. Phys. vol. 2, pg. 2047 (2000).
 * This algorithm hence performs
 * \f[ \sum_{me} t_{imJ}^{aeB} F''_{me} + \sum_{ME} t_{iMJ}^{aEB} F''_{ME} \rightarrow  \{ \bar{H}_{iJ}^{aB} \} \f]
 */
void CCMRCC::build_t2_iJaB_amplitudes_triples_diagram3() {
    // Loop over references
    for (int ref = 0; ref < wfn_->moinfo()->get_nunique(); ref++) {
        int unique_ref = wfn_->moinfo()->get_ref_number(ref, UniqueRefs);

        // Grab the temporary matrices
        auto HiJaBMatTmp = wfn_->blas()->get_MatTmp("t2_eqns[oO][vV]", unique_ref, none);
        auto TijKabCMatTmp = wfn_->blas()->get_MatTmp("t3[ooO][vvV]", unique_ref, none);
        auto TiJKaBCMatTmp = wfn_->blas()->get_MatTmp("t3[oOO][vVV]", unique_ref, none);
        auto FmeMatTmp = wfn_->blas()->get_MatTmp("F2_me[o][v]", unique_ref, none);
        auto FMEMatTmp = wfn_->blas()->get_MatTmp("F2_ME[O][V]", unique_ref, none);

        // Grab the indexing for t3[ijk][abc]
        auto& ij_tuples = HiJaBMatTmp->get_left()->get_tuples();
        auto& ab_tuples = HiJaBMatTmp->get_right()->get_tuples();
        auto& m_tuples = FmeMatTmp->get_left()->get_tuples();
        auto& e_tuples = FmeMatTmp->get_right()->get_tuples();

        auto TijKabC_matrix = TijKabCMatTmp->get_matrix();
        auto TiJKaBC_matrix = TiJKaBCMatTmp->get_matrix();
        auto HiJaB_matrix = HiJaBMatTmp->get_matrix();
        auto Fme_matrix = FmeMatTmp->get_matrix();
        auto FME_matrix = FMEMatTmp->get_matrix();
        auto ijkIndex = wfn_->blas()->get_index("[ooo]");
        auto abcIndex = wfn_->blas()->get_index("[vvv]");

        for (int h = 0; h < wfn_->moinfo()->get_nirreps(); h++) {
            size_t ij_offset = HiJaBMatTmp->get_left()->get_first(h);
            size_t ab_offset = HiJaBMatTmp->get_right()->get_first(h);
            for (int ab = 0; ab < HiJaBMatTmp->get_right_pairpi(h); ab++) {
                int a = ab_tuples[ab_offset + ab][0];
                int b = ab_tuples[ab_offset + ab][1];
                for (int ij = 0; ij < HiJaBMatTmp->get_left_pairpi(h); ij++) {
                    int i = ij_tuples[ij_offset + ij][0];
                    int j = ij_tuples[ij_offset + ij][1];
                    for (int m_sym = 0; m_sym < wfn_->moinfo()->get_nirreps(); m_sym++) {
                        size_t m_offset = FmeMatTmp->get_left()->get_first(m_sym);
                        size_t e_offset = FmeMatTmp->get_right()->get_first(m_sym);
                        for (int e = 0; e < FmeMatTmp->get_right_pairpi(m_sym); e++) {
                            int e_abs = e + e_offset;
                            size_t aeb = abcIndex->get_tuple_rel_index(a, e_abs, b);
                            int aeb_sym = abcIndex->get_tuple_irrep(a, e_abs, b);
                            for (int m = 0; m < FmeMatTmp->get_left_pairpi(m_sym); m++) {
                                int m_abs = m + m_offset;
                                size_t imj = ijkIndex->get_tuple_rel_index(i, m_abs, j);
                                HiJaB_matrix[h][ij][ab] += TijKabC_matrix[aeb_sym][imj][aeb] * Fme_matrix[m_sym][m][e];
                                HiJaB_matrix[h][ij][ab] += TiJKaBC_matrix[aeb_sym][imj][aeb] * FME_matrix[m_sym][m][e];
                            }
                        }
                    }
                }
            }
        }
    }
    //   wfn_->blas()->print("t2_test[oO][vV]{u}");
    //   wfn_->blas()->solve("ERROR{u} = 1000000.0 t2_test[oO][vV]{u} . t2_test[oO][vV]{u}");
    //   wfn_->blas()->print("ERROR{u}");
}

/**
 * @brief Computes the contraction
 * \f[ \sum_{me} t_{ijm}^{abe} F''_{me} \f]
 * as described in Fig. 1 of Phys. Chem. Chem. Phys. vol. 2, pg. 2047 (2000).
 * This algorithm hence performs
 * \f[ \sum_{me} t_{mIJ}^{eAB} F''_{me} + \sum_{ME} t_{MIJ}^{EAB} F''_{ME} \rightarrow  \{ \bar{H}_{IJ}^{AB} \} \f]
 */
void CCMRCC::build_t2_IJAB_amplitudes_triples_diagram3() {
    // Loop over references
    for (int ref = 0; ref < wfn_->moinfo()->get_nunique(); ref++) {
        int unique_ref = wfn_->moinfo()->get_ref_number(ref, UniqueRefs);

        // Grab the temporary matrices
        auto HIJABMatTmp = wfn_->blas()->get_MatTmp("t2_eqns[OO][VV]", unique_ref, none);
        auto TiJKaBCMatTmp = wfn_->blas()->get_MatTmp("t3[oOO][vVV]", unique_ref, none);
        auto TIJKABCMatTmp = wfn_->blas()->get_MatTmp("t3[OOO][VVV]", unique_ref, none);
        auto FmeMatTmp = wfn_->blas()->get_MatTmp("F2_me[o][v]", unique_ref, none);
        auto FMEMatTmp = wfn_->blas()->get_MatTmp("F2_ME[O][V]", unique_ref, none);

        // Grab the indexing for t3[ijk][abc]
        auto& ij_tuples = HIJABMatTmp->get_left()->get_tuples();
        auto& ab_tuples = HIJABMatTmp->get_right()->get_tuples();
        auto& m_tuples = FmeMatTmp->get_left()->get_tuples();
        auto& e_tuples = FmeMatTmp->get_right()->get_tuples();

        auto TiJKaBC_matrix = TiJKaBCMatTmp->get_matrix();
        auto TIJKABC_matrix = TIJKABCMatTmp->get_matrix();
        auto HIJAB_matrix = HIJABMatTmp->get_matrix();
        auto Fme_matrix = FmeMatTmp->get_matrix();
        auto FME_matrix = FMEMatTmp->get_matrix();
        auto ijkIndex = wfn_->blas()->get_index("[ooo]");
        auto abcIndex = wfn_->blas()->get_index("[vvv]");

        for (int h = 0; h < wfn_->moinfo()->get_nirreps(); h++) {
            size_t ij_offset = HIJABMatTmp->get_left()->get_first(h);
            size_t ab_offset = HIJABMatTmp->get_right()->get_first(h);
            for (int ab = 0; ab < HIJABMatTmp->get_right_pairpi(h); ab++) {
                int a = ab_tuples[ab_offset + ab][0];
                int b = ab_tuples[ab_offset + ab][1];
                for (int ij = 0; ij < HIJABMatTmp->get_left_pairpi(h); ij++) {
                    int i = ij_tuples[ij_offset + ij][0];
                    int j = ij_tuples[ij_offset + ij][1];
                    for (int m_sym = 0; m_sym < wfn_->moinfo()->get_nirreps(); m_sym++) {
                        size_t m_offset = FmeMatTmp->get_left()->get_first(m_sym);
                        size_t e_offset = FmeMatTmp->get_right()->get_first(m_sym);
                        for (int e = 0; e < FmeMatTmp->get_right_pairpi(m_sym); e++) {
                            int e_abs = e + e_offset;
                            size_t eab = abcIndex->get_tuple_rel_index(e_abs, a, b);
                            int eab_sym = abcIndex->get_tuple_irrep(e_abs, a, b);
                            for (int m = 0; m < FmeMatTmp->get_left_pairpi(m_sym); m++) {
                                int m_abs = m + m_offset;
                                size_t mij = ijkIndex->get_tuple_rel_index(m_abs, i, j);
                                HIJAB_matrix[h][ij][ab] += TiJKaBC_matrix[eab_sym][mij][eab] * Fme_matrix[m_sym][m][e];
                                HIJAB_matrix[h][ij][ab] += TIJKABC_matrix[eab_sym][mij][eab] * FME_matrix[m_sym][m][e];
                            }
                        }
                    }
                }
            }
        }
    }
    //   wfn_->blas()->print("t2_test[oO][vV]{u}");
    //   wfn_->blas()->solve("ERROR{u} = 1000000.0 t2_test[oO][vV]{u} . t2_test[oO][vV]{u}");
    //   wfn_->blas()->print("ERROR{u}");
}

/**
 * @brief Computes the contraction
 * \f[ -\frac{1}{2} P(ij)\sum_{mne} t_{imn}^{abe} W_{mnje} \f]
 * as described in Fig. 1 of Phys. Chem. Chem. Phys. vol. 2, pg. 2047 (2000).
 * This algorithm hence performs
 * \f[ \frac{1}{2}\sum_{mne} t_{imn}^{abe} W_{mnje} \rightarrow  \{ -\bar{H}_{ij}^{ab},\bar{H}_{ji}^{ab} \} \f]
 * \f[ \sum_{mNE} t_{imN}^{abE} W_{mNjE} \rightarrow  \{ -\bar{H}_{ij}^{ab},\bar{H}_{ji}^{ab} \} \f]
 */
/*
void CCMRCC::build_t2_ijab_amplitudes_triples_diagram1()
{
  // Loop over references
  for(int ref=0;ref<wfn_->moinfo()->get_nunique();ref++){
    int unique_ref  = wfn_->moinfo()->get_ref_number("u",ref);

    // Grab the temporary matrices
    CCMatTmp  TijkabcMatTmp = wfn_->blas()->get_MatTmp("t3[ooo][vvv]",unique_ref,none);
    double*** Tijkabc_matrix = TijkabcMatTmp->get_matrix();
    CCMatTmp  TijKabCMatTmp = wfn_->blas()->get_MatTmp("t3[ooO][vvV]",unique_ref,none);
    double*** TijKabC_matrix = TijKabCMatTmp->get_matrix();

    CCMatTmp  WkijaMatTmp = wfn_->blas()->get_MatTmp("W_kija[o][oov]",unique_ref,none);
    CCMatTmp  WkiJAMatTmp = wfn_->blas()->get_MatTmp("W_kiJA[o][oOV]",unique_ref,none);
    double*** Wkija_matrix = WkijaMatTmp->get_matrix();
    double*** WkiJA_matrix = WkiJAMatTmp->get_matrix();

    CCMatTmp  HijabMatTmp = wfn_->blas()->get_MatTmp("t2_eqns[oo][vv]",unique_ref,none);

    // Grab the indexing for t3[iab][jkc]
    CCIndex* iab_indexing = wfn_->blas()->get_index("[ovv]");
    CCIndex* jkc_indexing = wfn_->blas()->get_index("[oov]");
    CCIndex*   j_indexing = wfn_->blas()->get_index("[o]");


    short** iab_tuples = iab_indexing->get_tuples();
    short** jkc_tuples = jkc_indexing->get_tuples();

    // PART A: Sort T[ijk][abc]->T[iab][jkc]
    double ***T_iabjkc;
    double ***H_iabj;

    init_matrix<double**>(T_iabjkc,wfn_->moinfo()->get_nirreps());
    init_matrix<double**>(H_iabj,wfn_->moinfo()->get_nirreps());

    for(int h =0; h < wfn_->moinfo()->get_nirreps();h++){
      // Allocate a block of T_iabjkc
      init_matrix<double>(T_iabjkc[h],iab_indexing->get_pairpi(h),jkc_indexing->get_pairpi(h));
      init_matrix<double>(H_iabj[h],iab_indexing->get_pairpi(h),j_indexing->get_pairpi(h));

      size_t iab_offset = iab_indexing->get_first(h);
      size_t jkc_offset = jkc_indexing->get_first(h);

      // AAA Contribution

      // Sort this block
      for(int iab = 0;iab<iab_indexing->get_pairpi(h);iab++){
        int i = iab_tuples[iab_offset + iab][0];
        int a = iab_tuples[iab_offset + iab][1];
        int b = iab_tuples[iab_offset + iab][2];
        for(int jkc = 0;jkc<jkc_indexing->get_pairpi(h);jkc++){
          int j = jkc_tuples[jkc_offset + jkc][0];
          int k = jkc_tuples[jkc_offset + jkc][1];
          int c = jkc_tuples[jkc_offset + jkc][2];
          T_iabjkc[h][iab][jkc] = TijkabcMatTmp->get_six_address_element(i,j,k,a,b,c);
        }
      }
      int m = iab_indexing->get_pairpi(h);
      int n = j_indexing->get_pairpi(h);
      int k = jkc_indexing->get_pairpi(h);
      if(m*n*k){
        double alpha = 1.0;
        double beta  = 0.0;
        C_DGEMM_22(m, n, k, alpha,&(T_iabjkc[h][0][0]), k,
                   &(Wkija_matrix[h][0][0]), k, beta, &(H_iabj[h][0][0]),n);
      }

      // Sort the result
      size_t j_offset = j_indexing->get_first(h);
      for(int iab = 0;iab<iab_indexing->get_pairpi(h);iab++){
        int i = iab_tuples[iab_offset + iab][0];
        int a = iab_tuples[iab_offset + iab][1];
        int b = iab_tuples[iab_offset + iab][2];
        for(int j = 0;j<j_indexing->get_pairpi(h);j++){
          int abs_j = j + j_offset;
          HijabMatTmp->add_four_address_element(i,abs_j,a,b,-0.5*H_iabj[h][iab][j]);
          HijabMatTmp->add_four_address_element(abs_j,i,a,b,0.5*H_iabj[h][iab][j]);
        }
      }

      // AAB Contribution

      // Sort this block
      for(int iab = 0;iab<iab_indexing->get_pairpi(h);iab++){
        int i = iab_tuples[iab_offset + iab][0];
        int a = iab_tuples[iab_offset + iab][1];
        int b = iab_tuples[iab_offset + iab][2];
        for(int jkc = 0;jkc<jkc_indexing->get_pairpi(h);jkc++){
          int j = jkc_tuples[jkc_offset + jkc][0];
          int k = jkc_tuples[jkc_offset + jkc][1];
          int c = jkc_tuples[jkc_offset + jkc][2];
          T_iabjkc[h][iab][jkc] = TijKabCMatTmp->get_six_address_element(i,j,k,a,b,c);
        }
      }
      m = iab_indexing->get_pairpi(h);
      n = j_indexing->get_pairpi(h);
      k = jkc_indexing->get_pairpi(h);
      if(m*n*k){
        double alpha = 1.0;
        double beta  = 0.0;
        C_DGEMM_22(m, n, k, alpha,&(T_iabjkc[h][0][0]), k,
                   &(WkiJA_matrix[h][0][0]), k, beta, &(H_iabj[h][0][0]),n);
      }

      // Sort the result
      j_offset = j_indexing->get_first(h);
      for(int iab = 0;iab<iab_indexing->get_pairpi(h);iab++){
        int i = iab_tuples[iab_offset + iab][0];
        int a = iab_tuples[iab_offset + iab][1];
        int b = iab_tuples[iab_offset + iab][2];
        for(int j = 0;j<j_indexing->get_pairpi(h);j++){
          int abs_j = j + j_offset;
          HijabMatTmp->add_four_address_element(i,abs_j,a,b,-H_iabj[h][iab][j]);
          HijabMatTmp->add_four_address_element(abs_j,i,a,b,H_iabj[h][iab][j]);
        }
      }

      // Deallocate the memory for the block
      free_matrix<double>(H_iabj[h],iab_indexing->get_pairpi(h),j_indexing->get_pairpi(h));
      free_matrix<double>(T_iabjkc[h],iab_indexing->get_pairpi(h),jkc_indexing->get_pairpi(h));
    }
    free_matrix<double**>(H_iabj,wfn_->moinfo()->get_nirreps());
    free_matrix<double**>(T_iabjkc,wfn_->moinfo()->get_nirreps());
  }
//   wfn_->blas()->print("t2_test[oo][vv]{u}");
//   wfn_->blas()->solve("ERROR{u} = 1000000.0 t2_test[oo][vv]{u} . t2_test[oo][vv]{u}");
//   wfn_->blas()->print("ERROR{u}");
}
*/

}  // namespace psimrcc
}  // namespace psi
