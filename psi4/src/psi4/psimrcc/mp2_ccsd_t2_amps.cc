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

#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/liboptions/liboptions.h"

#include "mp2_ccsd.h"
#include "blas.h"

namespace psi {
namespace psimrcc {

void MP2_CCSD::build_mp2_t2_iJaB_amplitudes() {
    wfn_->blas()->solve("t2_eqns[oO][vV]{u}   = <[oo]|[vv]>");

    wfn_->blas()->solve("t2_eqns[oO][vV]{u} += #3214# t2[V][vOo]{u} 1@2 offdiagonal_F[v][v]{u}");
    wfn_->blas()->solve("t2_eqns[oO][vV]{u} += #4123# t2[v][VoO]{u} 1@2 offdiagonal_F[v][v]{u}");

    wfn_->blas()->solve("t2_eqns[oO][vV]{u} += #1432# - t2[O][oVv]{u} 1@1 offdiagonal_F[o][o]{u}");
    wfn_->blas()->solve("t2_eqns[oO][vV]{u} += #2341# - t2[o][OvV]{u} 1@1 offdiagonal_F[o][o]{u}");

    wfn_->blas()->solve("t2_delta[oO][vV]{u} = t2_eqns[oO][vV]{u} / d2[oO][vV]{u} - t2[oO][vV]{u}");

    // alpha,beta
    wfn_->blas()->solve("t2[oO][vV]{u}  = t2_eqns[oO][vV]{u} / d2[oO][vV]{u}");

    // alpha,alpha
    wfn_->blas()->solve("t2_eqns[oo][vv]{u}  = t2_eqns[oO][vV]{u}");
    wfn_->blas()->solve("t2_eqns[oo][vv]{u} += #2134# - t2_eqns[oO][vV]{u}");
    wfn_->blas()->solve("t2[oo][vv]{u}       = t2_eqns[oo][vv]{u} / d2[oo][vv]{u}");

    // beta, beta
    wfn_->blas()->solve("t2[OO][VV]{u}  = t2[oo][vv]{u}");
}

void MP2_CCSD::build_t2_iJaB_amplitudes() {
    START_TIMER("Building the T2_iJaB Amplitudes");

    // AAAA case (CCSD)
    wfn_->blas()->solve("HiJaB[aA][aA]{u}  = <[aa]|[aa]>");
    wfn_->blas()->solve("HiJaB[aA][aA]{u} += #3214# t2_VvOo[V][aAa]{u} 1@2 F'_AE[A][V]{u}");
    wfn_->blas()->solve("HiJaB[aA][aA]{u} += #4123# t2_vVoO[v][AaA]{u} 1@2 F'_ae[a][v]{u}");
    wfn_->blas()->solve("HiJaB[aA][aA]{u} += #1432# - t2_OoVv[O][aAa]{u} 1@1 F'_MI[O][A]{u}");
    wfn_->blas()->solve("HiJaB[aA][aA]{u} += #2341# - t2_oOvV[o][AaA]{u} 1@1 F'_mi[o][a]{u}");

    wfn_->blas()->solve("HiJaB[aA][aA]{u} += W_mNiJ[oO][aA]{u} 1@1 tau_oOvV[oO][aA]{u}");

    wfn_->blas()->solve("HiJaB[aA][aA]{u} += tau_oOvV[aA][vV]{u} 2@2 <[aa]|[vv]>");

    wfn_->blas()->solve("HiJaB[aA][aA]{u} += #1234#  - Z_iJaM[aAa][O]{u} 2@1 t1_OV[O][A]{u}");
    wfn_->blas()->solve("HiJaB[aA][aA]{u} += #1243#    Z_iJAm[aAA][o]{u} 2@1 t1_ov[o][a]{u}");

    wfn_->blas()->solve("HiJaB[aA][aA]{u} += #2413#   W_jbME[aa][OV]{u} 2@2 t2_ovov[aa][ov]{u}");
    wfn_->blas()->solve("HiJaB[aA][aA]{u} += #2413#   W_jbme[aa][ov]{u} 2@2 t2_ovOV[aa][OV]{u}");
    wfn_->blas()->solve("HiJaB[aA][aA]{u} += #2314#   W_jBmE[aA][oV]{u} 2@2 t2_oVOv[aA][Ov]{u}");
    wfn_->blas()->solve("HiJaB[aA][aA]{u} += #1423#   W_jBmE[aA][oV]{u} 2@1 t2_oVOv[oV][Aa]{u}");
    wfn_->blas()->solve("HiJaB[aA][aA]{u} += #1324#   W_jbME[aa][OV]{u} 2@2 t2_OVOV[AA][OV]{u}");
    wfn_->blas()->solve("HiJaB[aA][aA]{u} += #1324#   W_jbme[aa][ov]{u} 2@1 t2_ovOV[ov][AA]{u}");

    wfn_->blas()->solve("HiJaB[aA][aA]{u} += #4213# - ([ov]|[aa]) 1@2 t1t1_iame[aa][ov]{u}");
    wfn_->blas()->solve("HiJaB[aA][aA]{u} += #2314# - <[ov]|[aa]> 1@2 t1t1_iAMe[aA][Ov]{u}");
    wfn_->blas()->solve("HiJaB[aA][aA]{u} += #1423# - <[ov]|[aa]> 1@1 t1t1_iAMe[oV][Aa]{u}");
    wfn_->blas()->solve("HiJaB[aA][aA]{u} += #3124# - ([ov]|[aa]) 1@2 t1t1_IAME[AA][OV]{u}");

    wfn_->blas()->solve("HiJaB[aA][aA]{u} += #1234#   t1_ov[a][v]{u} 2@1 <[v]|[aaa]>");
    wfn_->blas()->solve("HiJaB[aA][aA]{u} += #2143#   t1_OV[A][V]{u} 2@1 <[v]|[aaa]>");
    wfn_->blas()->solve("HiJaB[aA][aA]{u} += #3412# - t1_ov[o][a]{u} 1@1 <[o]|[aaa]>");
    wfn_->blas()->solve("HiJaB[aA][aA]{u} += #4321# - t1_OV[O][A]{u} 1@1 <[o]|[aaa]>");

    if (options_.get_str("MP2_CCSD_METHOD") == "II") {
        // RAAA case (CCSD)
        wfn_->blas()->solve("HiJaB[oA][aA]{u}  = <[oa]|[aa]>");
        wfn_->blas()->solve("HiJaB[oA][aA]{u} += #3214# t2_VvOo[V][aAo]{u} 1@2 F'_AE[A][V]{u}");
        wfn_->blas()->solve("HiJaB[oA][aA]{u} += #4123# t2_vVoO[v][AoA]{u} 1@2 F'_ae[a][v]{u}");
        wfn_->blas()->solve("HiJaB[oA][aA]{u} += #1432# - t2_OoVv[O][oAa]{u} 1@1 F'_MI[O][A]{u}");
        wfn_->blas()->solve("HiJaB[oA][aA]{u} += #2341# - t2_oOvV[o][AaA]{u} 1@1 F'_mi[o][o]{u}");

        wfn_->blas()->solve("HiJaB[oA][aA]{u} += W_mNiJ[oO][oA]{u} 1@1 tau_oOvV[oO][aA]{u}");

        wfn_->blas()->solve("HiJaB[oA][aA]{u} += tau_oOvV[oA][vV]{u} 2@2 <[aa]|[vv]>");

        wfn_->blas()->solve("HiJaB[oA][aA]{u} += #1234#  - Z_iJaM[oAa][O]{u} 2@1 t1_OV[O][A]{u}");
        wfn_->blas()->solve("HiJaB[oA][aA]{u} += #1243#    Z_iJAm[oAA][o]{u} 2@1 t1_ov[o][a]{u}");

        wfn_->blas()->solve("HiJaB[oA][aA]{u} += #2413#   W_jbME[aa][OV]{u} 2@2 t2_ovov[oa][ov]{u}");
        wfn_->blas()->solve("HiJaB[oA][aA]{u} += #2413#   W_jbme[aa][ov]{u} 2@2 t2_ovOV[oa][OV]{u}");
        wfn_->blas()->solve("HiJaB[oA][aA]{u} += #2314#   W_jBmE[aA][oV]{u} 2@2 t2_oVOv[oA][Ov]{u}");
        wfn_->blas()->solve("HiJaB[oA][aA]{u} += #1423#   W_jBmE[oA][oV]{u} 2@1 t2_oVOv[oV][Aa]{u}");
        wfn_->blas()->solve("HiJaB[oA][aA]{u} += #1324#   W_jbME[oa][OV]{u} 2@2 t2_OVOV[AA][OV]{u}");
        wfn_->blas()->solve("HiJaB[oA][aA]{u} += #1324#   W_jbme[oa][ov]{u} 2@1 t2_ovOV[ov][AA]{u}");

        wfn_->blas()->solve("HiJaB[oA][aA]{u} += #4213# - ([ov]|[aa]) 1@2 t1t1_iame[oa][ov]{u}");
        wfn_->blas()->solve("HiJaB[oA][aA]{u} += #2314# - <[ov]|[aa]> 1@2 t1t1_iAMe[oA][Ov]{u}");
        wfn_->blas()->solve("HiJaB[oA][aA]{u} += #1423# - <[ov]|[oa]> 1@1 t1t1_iAMe[oV][Aa]{u}");
        wfn_->blas()->solve("HiJaB[oA][aA]{u} += #3124# - ([ov]|[ao]) 1@2 t1t1_IAME[AA][OV]{u}");

        wfn_->blas()->solve("HiJaB[oA][aA]{u} += #1234#   t1[o][v]{u} 2@1 <[v]|[aaa]>");
        wfn_->blas()->solve("HiJaB[oA][aA]{u} += #2143#   t1_OV[A][V]{u} 2@1 <[v]|[oaa]>");
        wfn_->blas()->solve("HiJaB[oA][aA]{u} += #3412# - t1_ov[o][a]{u} 1@1 <[o]|[aoa]>");
        wfn_->blas()->solve("HiJaB[oA][aA]{u} += #4321# - t1_OV[O][A]{u} 1@1 <[o]|[aao]>");

        wfn_->blas()->solve("HiJaB[aO][aA]{u} = #2143# HiJaB[oA][aA]{u}");

        // AARA case (CCSD)
        wfn_->blas()->solve("HiJaB[aA][vA]{u}  = <[aa]|[va]>");
        wfn_->blas()->solve("HiJaB[aA][vA]{u} += #3214# t2_VvOo[V][vAa]{u} 1@2 F'_AE[A][V]{u}");
        wfn_->blas()->solve("HiJaB[aA][vA]{u} += #4123# t2_vVoO[v][AaA]{u} 1@2 F'_ae[v][v]{u}");
        wfn_->blas()->solve("HiJaB[aA][vA]{u} += #1432# - t2_OoVv[O][aAv]{u} 1@1 F'_MI[O][A]{u}");
        wfn_->blas()->solve("HiJaB[aA][vA]{u} += #2341# - t2_oOvV[o][AvA]{u} 1@1 F'_mi[o][a]{u}");

        wfn_->blas()->solve("HiJaB[aA][vA]{u} += W_mNiJ[oO][aA]{u} 1@1 tau_oOvV[oO][vA]{u}");

        wfn_->blas()->solve("HiJaB[aA][vA]{u} += tau_oOvV[aA][vV]{u} 2@2 <[va]|[vv]>");

        wfn_->blas()->solve("HiJaB[aA][vA]{u} += #1234#  - Z_iJaM[aAv][O]{u} 2@1 t1_OV[O][A]{u}");
        wfn_->blas()->solve("HiJaB[aA][vA]{u} += #1243#    Z_iJAm[aAA][o]{u} 2@1 t1[o][v]{u}");

        wfn_->blas()->solve("HiJaB[aA][vA]{u} += #2413#   W_jbME[aa][OV]{u} 2@2 t2_ovov[av][ov]{u}");
        wfn_->blas()->solve("HiJaB[aA][vA]{u} += #2413#   W_jbme[aa][ov]{u} 2@2 t2_ovOV[av][OV]{u}");
        wfn_->blas()->solve("HiJaB[aA][vA]{u} += #2314#   W_jBmE[aV][oV]{u} 2@2 t2_oVOv[aA][Ov]{u}");
        wfn_->blas()->solve("HiJaB[aA][vA]{u} += #1423#   W_jBmE[aA][oV]{u} 2@1 t2_oVOv[oV][Av]{u}");
        wfn_->blas()->solve("HiJaB[aA][vA]{u} += #1324#   W_jbME[av][OV]{u} 2@2 t2_OVOV[AA][OV]{u}");
        wfn_->blas()->solve("HiJaB[aA][vA]{u} += #1324#   W_jbme[av][ov]{u} 2@1 t2_ovOV[ov][AA]{u}");

        wfn_->blas()->solve("HiJaB[aA][vA]{u} += #4213# - ([ov]|[aa]) 1@2 t1t1_iame[av][ov]{u}");
        wfn_->blas()->solve("HiJaB[aA][vA]{u} += #2314# - <[ov]|[av]> 1@2 t1t1_iAMe[aA][Ov]{u}");
        wfn_->blas()->solve("HiJaB[aA][vA]{u} += #1423# - <[ov]|[aa]> 1@1 t1t1_iAMe[oV][Av]{u}");
        wfn_->blas()->solve("HiJaB[aA][vA]{u} += #3124# - ([ov]|[va]) 1@2 t1t1_IAME[AA][OV]{u}");

        wfn_->blas()->solve("HiJaB[aA][vA]{u} += #1234#   t1_ov[a][v]{u} 2@1 <[v]|[ava]>");
        wfn_->blas()->solve("HiJaB[aA][vA]{u} += #2143#   t1_OV[A][V]{u} 2@1 <[v]|[aav]>");
        wfn_->blas()->solve("HiJaB[aA][vA]{u} += #3412# - t1[o][v]{u} 1@1 <[o]|[aaa]>");
        wfn_->blas()->solve("HiJaB[aA][vA]{u} += #4321# - t1_OV[O][A]{u} 1@1 <[o]|[vaa]>");

        wfn_->blas()->solve("HiJaB[aA][aV]{u} = #2143# HiJaB[aA][vA]{u}");

        wfn_->blas()->expand_spaces("HiJaB[oA][aA]{u}", "t2_eqns[oO][vV]{u}");
        wfn_->blas()->expand_spaces("HiJaB[aO][aA]{u}", "t2_eqns[oO][vV]{u}");

        wfn_->blas()->expand_spaces("HiJaB[aA][vA]{u}", "t2_eqns[oO][vV]{u}");
        wfn_->blas()->expand_spaces("HiJaB[aA][aV]{u}", "t2_eqns[oO][vV]{u}");
    }

    wfn_->blas()->expand_spaces("HiJaB[aA][aA]{u}", "t2_eqns[oO][vV]{u}");

    wfn_->blas()->solve("t2_delta[oO][vV]{u} = t2_eqns[oO][vV]{u} / d2[oO][vV]{u} - t2[oO][vV]{u}");

    wfn_->blas()->solve("t2[oO][vV]{u}  = t2_eqns[oO][vV]{u} / d2[oO][vV]{u}");
    END_TIMER("Building the T2_iJaB Amplitudes");
}

void MP2_CCSD::build_t2_ijab_amplitudes() {
    START_TIMER("Building the T2_ijab Amplitudes");
    wfn_->blas()->solve("t2_eqns[oo][vv]{u}  = t2_eqns[oO][vV]{u}");
    wfn_->blas()->solve("t2_eqns[oo][vv]{u} += #2134# - t2_eqns[oO][vV]{u}");
    wfn_->blas()->solve("t2[oo][vv]{u}  = t2_eqns[oo][vv]{u} / d2[oo][vv]{u}");
    END_TIMER("Building the T2_ijab Amplitudes");
}

void MP2_CCSD::build_t2_IJAB_amplitudes() {
    START_TIMER("Building the T2_IJAB Amplitudes");
    wfn_->blas()->solve("t2[OO][VV]{u}  = t2[oo][vv]{u}");
    END_TIMER("Building the T2_IJAB Amplitudes");
}

}  // namespace psimrcc
}  // namespace psi
