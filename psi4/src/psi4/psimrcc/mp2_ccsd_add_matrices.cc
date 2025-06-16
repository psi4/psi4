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

#include "blas.h"
#include "mp2_ccsd.h"

namespace psi {
namespace psimrcc {

void MP2_CCSD::add_matrices() {
    // O^4
    wfn_->blas()->add_Matrix("<[oo]:[oo]>");
    wfn_->blas()->add_Matrix("<[oo]|[oo]>");

    // O^3V
    wfn_->blas()->add_Matrix("([oo]:[ov])");
    wfn_->blas()->add_Matrix("([oo]|[ov])");
    wfn_->blas()->add_Matrix("<[o]:[voo]>");
    wfn_->blas()->add_Matrix("<[o]|[voo]>");
    wfn_->blas()->add_Matrix("<[ooo]|[v]>");
    wfn_->blas()->add_Matrix("<[o]:[oov]>");
    wfn_->blas()->add_Matrix("<[o]|[oov]>");
    wfn_->blas()->add_Matrix("<[o]|[ovo]>");

    // O^2V^2
    wfn_->blas()->add_Matrix("<[oo]:[vv]>");
    wfn_->blas()->add_Matrix("<[oo]|[vv]>");
    wfn_->blas()->add_Matrix("<[v]:[voo]>");
    wfn_->blas()->add_Matrix("<[v]|[voo]>");
    wfn_->blas()->add_Matrix("<[o]:[ovv]>");
    wfn_->blas()->add_Matrix("<[o]|[ovv]>");
    wfn_->blas()->add_Matrix("([ov]:[ov])");
    wfn_->blas()->add_Matrix("([ov]|[ov])");
    wfn_->blas()->add_Matrix("([ov]|[vo])");
    wfn_->blas()->add_Matrix("<[ov]|[ov]>");
    wfn_->blas()->add_Matrix("<[ov]|[vo]>");  // used only once

    // OV^3
    wfn_->blas()->add_Matrix("([ov]:[vv])");
    wfn_->blas()->add_Matrix("([ov]|[vv])");
    wfn_->blas()->add_Matrix("<[v]:[ovv]>");
    wfn_->blas()->add_Matrix("<[v]|[ovv]>");
    wfn_->blas()->add_Matrix("<[vo]|[vv]>");
    wfn_->blas()->add_Matrix("([vvo]|[v])");
    wfn_->blas()->add_Matrix("<[ovv]:[v]>");

    // Fock Matrix
    wfn_->blas()->add_Matrix("fock[o][o]{u}");
    wfn_->blas()->add_Matrix("fock[o][v]{u}");
    wfn_->blas()->add_Matrix("fock[v][v]{u}");
    wfn_->blas()->add_Matrix("fock[O][O]{u}");
    wfn_->blas()->add_Matrix("fock[O][V]{u}");
    wfn_->blas()->add_Matrix("fock[V][V]{u}");

    wfn_->blas()->add_Matrix("fock[v][o]{u}");
    wfn_->blas()->add_Matrix("fock[V][O]{u}");

    wfn_->blas()->add_Matrix("fock[oo]{u}");
    wfn_->blas()->add_Matrix("fock[ov]{u}");
    wfn_->blas()->add_Matrix("fock[vv]{u}");
    wfn_->blas()->add_Matrix("fock[OO]{u}");
    wfn_->blas()->add_Matrix("fock[OV]{u}");
    wfn_->blas()->add_Matrix("fock[VV]{u}");

    wfn_->blas()->add_Matrix("fock[vo]{u}");
    wfn_->blas()->add_Matrix("fock[VO]{u}");

    wfn_->blas()->add_Matrix("fock[ff]{u}");
    wfn_->blas()->add_Matrix("fock[FF]{u}");

    // Denominators
    wfn_->blas()->add_Matrix("d1[o][v]{u}");
    wfn_->blas()->add_Matrix("d1[O][V]{u}");
    wfn_->blas()->add_Matrix("d2[oo][vv]{u}");
    wfn_->blas()->add_Matrix("d2[oO][vV]{u}");
    wfn_->blas()->add_Matrix("d2[OO][VV]{u}");

    // Shifted denominators
    wfn_->blas()->add_Matrix("d'1[o][v]{u}");
    wfn_->blas()->add_Matrix("d'1[O][V]{u}");
    wfn_->blas()->add_Matrix("d'2[oo][vv]{u}");
    wfn_->blas()->add_Matrix("d'2[oO][vV]{u}");
    wfn_->blas()->add_Matrix("d'2[OO][VV]{u}");

    // Amplitudes
    wfn_->blas()->add_Matrix("t1[ov]{u}");
    wfn_->blas()->add_Matrix("t1[OV]{u}");
    wfn_->blas()->add_Matrix("t1[o][v]{u}");
    wfn_->blas()->add_Matrix("t1[O][V]{u}");

    wfn_->blas()->add_Matrix("t1_delta[o][v]{u}");

    wfn_->blas()->add_Matrix("t2[oO][vV]{u}");
    wfn_->blas()->add_Matrix("t2[oo][vv]{u}");
    wfn_->blas()->add_Matrix("t2[OO][VV]{u}");

    wfn_->blas()->add_Matrix("t2[ov][ov]{u}");
    wfn_->blas()->add_Matrix("t2[ov][OV]{u}");

    wfn_->blas()->add_Matrix("t2[o][ovv]{u}");
    wfn_->blas()->add_Matrix("t2[o][OvV]{u}");
    wfn_->blas()->add_Matrix("t2[O][oVv]{u}");
    wfn_->blas()->add_Matrix("t2[O][OVV]{u}");

    wfn_->blas()->add_Matrix("t2[v][voo]{u}");
    wfn_->blas()->add_Matrix("t2[v][VoO]{u}");
    wfn_->blas()->add_Matrix("t2[V][vOo]{u}");
    wfn_->blas()->add_Matrix("t2[V][VOO]{u}");

    wfn_->blas()->add_Matrix("tau[oo][vv]{u}");
    wfn_->blas()->add_Matrix("tau[oO][vV]{u}");
    wfn_->blas()->add_Matrix("tau[OO][VV]{u}");

    wfn_->blas()->add_Matrix("t2_delta[oO][vV]{u}");

    wfn_->blas()->add_Matrix("tau2[v][voo]{u}");
    wfn_->blas()->add_Matrix("tau2[v][VoO]{u}");
    wfn_->blas()->add_Matrix("tau2[V][VOO]{u}");
    wfn_->blas()->add_Matrix("tau2[V][vOo]{u}");

    wfn_->blas()->add_Matrix("tau2[o][ovv]{u}");
    wfn_->blas()->add_Matrix("tau2[o][OvV]{u}");
    wfn_->blas()->add_Matrix("tau2[O][oVv]{u}");
    wfn_->blas()->add_Matrix("tau2[O][OVV]{u}");

    wfn_->blas()->add_Matrix("W_mNiJ[oO][oO]{u}");

    // Similarity transformed Hamiltonian
    wfn_->blas()->add_Matrix("t1_eqns[a][a]{u}");
    wfn_->blas()->add_Matrix("t1_eqns[o][v]{u}");
    wfn_->blas()->add_Matrix("t1_eqns[O][V]{u}");
    wfn_->blas()->add_Matrix("t2_eqns[oo][vv]{u}");
    wfn_->blas()->add_Matrix("t2_eqns[oO][vV]{u}");
    wfn_->blas()->add_Matrix("t2_eqns[OO][VV]{u}");

    // F intermediates
    wfn_->blas()->add_Matrix("F_ae[v][v]{u}");
    wfn_->blas()->add_Matrix("F_AE[V][V]{u}");
    wfn_->blas()->add_Matrix("F_mi[o][o]{u}");
    wfn_->blas()->add_Matrix("F_MI[O][O]{u}");
    wfn_->blas()->add_Matrix("F_me[o][v]{u}");
    wfn_->blas()->add_Matrix("F_ME[O][V]{u}");
    wfn_->blas()->add_Matrix("F_me[ov]{u}");
    wfn_->blas()->add_Matrix("F_ME[OV]{u}");
    wfn_->blas()->add_Matrix("F'_mi[o][o]{u}");

    // MRPT2 Intermediates
    wfn_->blas()->add_Matrix("ERef{u}");
    wfn_->blas()->add_Matrix("EPT2{u}");
    wfn_->blas()->add_Matrix("Eaa{u}");
    wfn_->blas()->add_Matrix("Ebb{u}");
    wfn_->blas()->add_Matrix("Eaaaa{u}");
    wfn_->blas()->add_Matrix("Eabab{u}");
    wfn_->blas()->add_Matrix("Ebbbb{u}");

    // A^4
    wfn_->blas()->add_Matrix("<[aa]|[aa]>");
    wfn_->blas()->add_Matrix("<[aa]:[aa]>");

    // A^3O
    wfn_->blas()->add_Matrix("<[o]|[aaa]>");
    wfn_->blas()->add_Matrix("<[o]:[aaa]>");
    wfn_->blas()->add_Matrix("<[oa]|[aa]>");

    // A^3V
    wfn_->blas()->add_Matrix("<[v]|[aaa]>");
    wfn_->blas()->add_Matrix("<[v]:[aaa]>");
    wfn_->blas()->add_Matrix("<[aa]|[va]>");

    // A^2O^2
    wfn_->blas()->add_Matrix("<[oo]|[aa]>");
    wfn_->blas()->add_Matrix("<[oo]:[aa]>");
    wfn_->blas()->add_Matrix("<[o]|[aoa]>");
    wfn_->blas()->add_Matrix("<[o]|[aao]>");

    // A^2OV
    wfn_->blas()->add_Matrix("<[aa]|[ov]>");
    wfn_->blas()->add_Matrix("<[aa]|[vo]>");
    wfn_->blas()->add_Matrix("<[ov]|[aa]>");
    wfn_->blas()->add_Matrix("([ov]|[aa])");
    wfn_->blas()->add_Matrix("([ov]:[aa])");
    wfn_->blas()->add_Matrix("<[oa]|[va]>");
    wfn_->blas()->add_Matrix("<[oa]:[va]>");
    wfn_->blas()->add_Matrix("<[oa]|[av]>");
    wfn_->blas()->add_Matrix("<[v]|[oav]>");
    wfn_->blas()->add_Matrix("<[v]|[oaa]>");
    wfn_->blas()->add_Matrix("<[o]|[vaa]>");
    wfn_->blas()->add_Matrix("<[ov]|[av]>");

    // A^2V^2
    wfn_->blas()->add_Matrix("<[aa]|[vv]>");
    wfn_->blas()->add_Matrix("<[aa]:[vv]>");
    wfn_->blas()->add_Matrix("<[v]|[ava]>");
    wfn_->blas()->add_Matrix("<[v]|[aav]>");
    wfn_->blas()->add_Matrix("<[ov]:[va]>");

    // AO^2V
    wfn_->blas()->add_Matrix("<[a]:[voo]>");
    wfn_->blas()->add_Matrix("<[a]|[voo]>");
    wfn_->blas()->add_Matrix("<[o]|[oav]>");
    wfn_->blas()->add_Matrix("<[o]:[oav]>");
    wfn_->blas()->add_Matrix("<[o]|[ova]>");
    wfn_->blas()->add_Matrix("<[oav]:[v]>");
    wfn_->blas()->add_Matrix("<[ov]|[oa]>");
    wfn_->blas()->add_Matrix("([ov]|[ao])");
    wfn_->blas()->add_Matrix("<[oa]:[vo]>");
    wfn_->blas()->add_Matrix("<[oa]|[vo]>");
    wfn_->blas()->add_Matrix("<[oa]|[ov]>");

    // AOV^2
    wfn_->blas()->add_Matrix("<[a]:[ovv]>");
    wfn_->blas()->add_Matrix("<[a]|[ovv]>");
    wfn_->blas()->add_Matrix("<[ao]|[vv]>");
    wfn_->blas()->add_Matrix("([avo]|[v])");
    wfn_->blas()->add_Matrix("<[va]|[vv]>");
    wfn_->blas()->add_Matrix("([ov]|[va])");
    wfn_->blas()->add_Matrix("<[ov]|[va]>");

    // Effective Hamiltonian in the active space
    wfn_->blas()->add_Matrix("Hia[a][a]{u}");
    wfn_->blas()->add_Matrix("HIA[A][A]{u}");
    wfn_->blas()->add_Matrix("Hijab[aa][aa]{u}");
    wfn_->blas()->add_Matrix("HiJaB[aA][aA]{u}");
    wfn_->blas()->add_Matrix("HIJAB[AA][AA]{u}");

    wfn_->blas()->add_Matrix("HiJaB[oA][aA]{u}");
    wfn_->blas()->add_Matrix("HiJaB[aO][aA]{u}");
    wfn_->blas()->add_Matrix("HiJaB[aA][vA]{u}");
    wfn_->blas()->add_Matrix("HiJaB[aA][aV]{u}");

    wfn_->blas()->add_Matrix("t1_ov[a][v]{u}");
    wfn_->blas()->add_Matrix("t1_OV[A][V]{u}");
    wfn_->blas()->add_Matrix("t1_ov[o][a]{u}");
    wfn_->blas()->add_Matrix("t1_OV[O][A]{u}");

    // Newly added
    wfn_->blas()->add_Matrix("t2_oovv[ao][av]{u}");
    wfn_->blas()->add_Matrix("t2_oOvV[aO][aV]{u}");
    wfn_->blas()->add_Matrix("t2_oOvV[oA][vA]{u}");
    wfn_->blas()->add_Matrix("t2_OOVV[AO][AV]{u}");
    wfn_->blas()->add_Matrix("t2_oOvV[oA][aV]{u}");
    wfn_->blas()->add_Matrix("t2_oOvV[aO][vA]{u}");

    wfn_->blas()->add_Matrix("t2_oovv[a][ovv]{u}");
    wfn_->blas()->add_Matrix("t2_oOvV[a][OvV]{u}");
    wfn_->blas()->add_Matrix("t2_OoVv[A][oVv]{u}");
    wfn_->blas()->add_Matrix("t2_OOVV[A][OVV]{u}");
    wfn_->blas()->add_Matrix("t2_oovv[oo][aa]{u}");
    wfn_->blas()->add_Matrix("t2_oOvV[oO][aA]{u}");
    wfn_->blas()->add_Matrix("t2_OOVV[OO][AA]{u}");
    wfn_->blas()->add_Matrix("t2_oovv[aa][vv]{u}");
    wfn_->blas()->add_Matrix("t2_oOvV[aA][vV]{u}");
    wfn_->blas()->add_Matrix("t2_OOVV[AA][VV]{u}");

    wfn_->blas()->add_Matrix("t2_ovov[aa][ov]{u}");
    wfn_->blas()->add_Matrix("t2_ovov[oa][ov]{u}");
    wfn_->blas()->add_Matrix("t2_ovov[av][ov]{u}");

    wfn_->blas()->add_Matrix("t2_ovOV[aa][OV]{u}");
    wfn_->blas()->add_Matrix("t2_ovOV[oa][OV]{u}");
    wfn_->blas()->add_Matrix("t2_ovOV[ov][AA]{u}");
    wfn_->blas()->add_Matrix("t2_ovOV[av][OV]{u}");

    wfn_->blas()->add_Matrix("t2_oVOv[aA][Ov]{u}");
    wfn_->blas()->add_Matrix("t2_oVOv[oV][Aa]{u}");
    wfn_->blas()->add_Matrix("t2_oVOv[oA][Ov]{u}");
    wfn_->blas()->add_Matrix("t2_oVOv[oV][Av]{u}");
    wfn_->blas()->add_Matrix("t2_OVOV[AA][OV]{u}");

    wfn_->blas()->add_Matrix("t2_oovv[ao][vv]{u}");

    wfn_->blas()->add_Matrix("t2_vvoo[a][voo]{u}");
    wfn_->blas()->add_Matrix("t2_vVoO[a][VoO]{u}");
    wfn_->blas()->add_Matrix("t2_VvOo[A][vOo]{u}");
    wfn_->blas()->add_Matrix("t2_VVOO[A][VOO]{u}");

    wfn_->blas()->add_Matrix("t2_vvoo[v][aaa]{u}");
    wfn_->blas()->add_Matrix("t2_vVoO[v][AaA]{u}");
    wfn_->blas()->add_Matrix("t2_VvOo[V][aAa]{u}");
    wfn_->blas()->add_Matrix("t2_VVOO[V][AAA]{u}");

    wfn_->blas()->add_Matrix("t2_oovv[o][aaa]{u}");
    wfn_->blas()->add_Matrix("t2_oOvV[o][AaA]{u}");
    wfn_->blas()->add_Matrix("t2_OoVv[O][aAa]{u}");
    wfn_->blas()->add_Matrix("t2_OOVV[O][AAA]{u}");

    // T2 A^2O^2
    wfn_->blas()->add_Matrix("t2_OoVv[O][oAa]{u}");

    // T2 A^2OV
    wfn_->blas()->add_Matrix("t2_vVoO[v][AoA]{u}");
    wfn_->blas()->add_Matrix("t2_VvOo[V][aAo]{u}");
    wfn_->blas()->add_Matrix("t2_VvOo[V][vAa]{u}");

    wfn_->blas()->add_Matrix("fock[o][a]{u}");
    wfn_->blas()->add_Matrix("fock[O][A]{u}");
    wfn_->blas()->add_Matrix("fock[a][a]{u}");
    wfn_->blas()->add_Matrix("fock[a][v]{u}");
    wfn_->blas()->add_Matrix("fock[A][V]{u}");
    wfn_->blas()->add_Matrix("fock[A][A]{u}");
    wfn_->blas()->add_Matrix("fock[aa]{u}");
    wfn_->blas()->add_Matrix("fock[AA]{u}");

    wfn_->blas()->add_Matrix("F_ae[a][v]{u}");
    wfn_->blas()->add_Matrix("F_AE[A][V]{u}");
    wfn_->blas()->add_Matrix("F'_ae[a][v]{u}");
    wfn_->blas()->add_Matrix("F'_AE[A][V]{u}");

    // Off-diagonal Fock matrices
    wfn_->blas()->add_Matrix("offdiagonal_F[v][v]{u}");
    wfn_->blas()->add_Matrix("offdiagonal_F[o][o]{u}");

    wfn_->blas()->add_Matrix("F_ae[a][v]{u}");
    wfn_->blas()->add_Matrix("F_AE[A][V]{u}");
    wfn_->blas()->add_Matrix("F'_ae[v][v]{u}");
    wfn_->blas()->add_Matrix("F'_ae[a][v]{u}");
    wfn_->blas()->add_Matrix("F'_AE[A][V]{u}");
    wfn_->blas()->add_Matrix("F'_mi[o][a]{u}");
    wfn_->blas()->add_Matrix("F'_MI[O][A]{u}");

    wfn_->blas()->add_Matrix("tau_oOvV[oO][aA]{u}");
    wfn_->blas()->add_Matrix("tau_oOvV[aA][vV]{u}");
    wfn_->blas()->add_Matrix("tau_oOVv[aA][Vv]{u}");
    wfn_->blas()->add_Matrix("tau_oOvV[oA][vV]{u}");
    wfn_->blas()->add_Matrix("tau_oOvV[oO][vA]{u}");
    wfn_->blas()->add_Matrix("tau_oOVv[oA][Vv]{u}");

    wfn_->blas()->add_Matrix("tau3_ovov[aa][ov]{u}");
    wfn_->blas()->add_Matrix("tau3_ovov[oa][ov]{u}");
    wfn_->blas()->add_Matrix("tau3_ovov[av][ov]{u}");
    wfn_->blas()->add_Matrix("tau3_oVvO[aA][vO]{u}");
    wfn_->blas()->add_Matrix("tau3_oVvO[oA][vO]{u}");
    wfn_->blas()->add_Matrix("tau3_oVvO[aV][vO]{u}");

    wfn_->blas()->add_Matrix("t2_oovv[ao][va]{u}");
    wfn_->blas()->add_Matrix("t2_oovv[oo][va]{u}");

    wfn_->blas()->add_Matrix("t2_oOvV[aO][vA]{u}");
    wfn_->blas()->add_Matrix("t2_oOvV[oO][vA]{u}");
    wfn_->blas()->add_Matrix("t2_oOvV[o][AvA]{u}");
    wfn_->blas()->add_Matrix("t2_oOvV[aO][vV]{u}");

    wfn_->blas()->add_Matrix("t2_OoVv[O][aAv]{u}");

    wfn_->blas()->add_Matrix("W_mNiJ[oO][aA]{u}");
    wfn_->blas()->add_Matrix("W_mNiJ[oO][oA]{u}");

    wfn_->blas()->add_Matrix("W_jbme[aa][ov]{u}");
    wfn_->blas()->add_Matrix("W_jbme[oa][ov]{u}");
    wfn_->blas()->add_Matrix("W_jbme[av][ov]{u}");

    wfn_->blas()->add_Matrix("W_jbME[aa][OV]{u}");
    wfn_->blas()->add_Matrix("W_jbME[oa][OV]{u}");
    wfn_->blas()->add_Matrix("W_jbME[av][OV]{u}");

    wfn_->blas()->add_Matrix("W_jBmE[aA][oV]{u}");
    wfn_->blas()->add_Matrix("W_jBmE[oA][oV]{u}");
    wfn_->blas()->add_Matrix("W_jBmE[aV][oV]{u}");

    wfn_->blas()->add_Matrix("Z_iJaM[aAa][O]{u}");
    wfn_->blas()->add_Matrix("Z_iJAm[aAA][o]{u}");
    wfn_->blas()->add_Matrix("Z_iJaM[oAa][O]{u}");
    wfn_->blas()->add_Matrix("Z_iJAm[oAA][o]{u}");
    wfn_->blas()->add_Matrix("Z_iJaM[aAv][O]{u}");

    wfn_->blas()->add_Matrix("t1t1_iame[aa][ov]{u}");
    wfn_->blas()->add_Matrix("t1t1_iame[oa][ov]{u}");
    wfn_->blas()->add_Matrix("t1t1_iame[av][ov]{u}");
    wfn_->blas()->add_Matrix("t1t1_IAME[AA][OV]{u}");
    wfn_->blas()->add_Matrix("t1t1_iAMe[aA][Ov]{u}");
    wfn_->blas()->add_Matrix("t1t1_iAMe[oV][Aa]{u}");
    wfn_->blas()->add_Matrix("t1t1_iAMe[oA][Ov]{u}");
    wfn_->blas()->add_Matrix("t1t1_iAMe[oV][Av]{u}");
}

}  // namespace psimrcc
}  // namespace psi
