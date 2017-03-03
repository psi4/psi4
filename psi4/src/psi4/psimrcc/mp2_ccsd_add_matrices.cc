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

#include "blas.h"
#include "mp2_ccsd.h"

namespace psi{ namespace psimrcc{

void MP2_CCSD::add_matrices()
{
  // O^4
  blas->add_Matrix("<[oo]:[oo]>");
  blas->add_Matrix("<[oo]|[oo]>");

  // O^3V
  blas->add_Matrix("([oo]:[ov])");
  blas->add_Matrix("([oo]|[ov])");
  blas->add_Matrix("<[o]:[voo]>");
  blas->add_Matrix("<[o]|[voo]>");
  blas->add_Matrix("<[ooo]|[v]>");
  blas->add_Matrix("<[o]:[oov]>");
  blas->add_Matrix("<[o]|[oov]>");
  blas->add_Matrix("<[o]|[ovo]>");

  // O^2V^2
  blas->add_Matrix("<[oo]:[vv]>");
  blas->add_Matrix("<[oo]|[vv]>");
  blas->add_Matrix("<[v]:[voo]>");
  blas->add_Matrix("<[v]|[voo]>");
  blas->add_Matrix("<[o]:[ovv]>");
  blas->add_Matrix("<[o]|[ovv]>");
  blas->add_Matrix("([ov]:[ov])");
  blas->add_Matrix("([ov]|[ov])");
  blas->add_Matrix("([ov]|[vo])");
  blas->add_Matrix("<[ov]|[ov]>");
  blas->add_Matrix("<[ov]|[vo]>"); // used only once

  // OV^3
  blas->add_Matrix("([ov]:[vv])");
  blas->add_Matrix("([ov]|[vv])");
  blas->add_Matrix("<[v]:[ovv]>");
  blas->add_Matrix("<[v]|[ovv]>");
  blas->add_Matrix("<[vo]|[vv]>");
  blas->add_Matrix("([vvo]|[v])");
  blas->add_Matrix("<[ovv]:[v]>");

  // Fock Matrix
  blas->add_Matrix("fock[o][o]{u}");
  blas->add_Matrix("fock[o][v]{u}");
  blas->add_Matrix("fock[v][v]{u}");
  blas->add_Matrix("fock[O][O]{u}");
  blas->add_Matrix("fock[O][V]{u}");
  blas->add_Matrix("fock[V][V]{u}");

  blas->add_Matrix("fock[v][o]{u}");
  blas->add_Matrix("fock[V][O]{u}");

  blas->add_Matrix("fock[oo]{u}");
  blas->add_Matrix("fock[ov]{u}");
  blas->add_Matrix("fock[vv]{u}");
  blas->add_Matrix("fock[OO]{u}");
  blas->add_Matrix("fock[OV]{u}");
  blas->add_Matrix("fock[VV]{u}");

  blas->add_Matrix("fock[vo]{u}");
  blas->add_Matrix("fock[VO]{u}");

  blas->add_Matrix("fock[ff]{u}");
  blas->add_Matrix("fock[FF]{u}");


  // Denominators
  blas->add_Matrix("d1[o][v]{u}");
  blas->add_Matrix("d1[O][V]{u}");
  blas->add_Matrix("d2[oo][vv]{u}");
  blas->add_Matrix("d2[oO][vV]{u}");
  blas->add_Matrix("d2[OO][VV]{u}");

  // Shifted denominators
  blas->add_Matrix("d'1[o][v]{u}");
  blas->add_Matrix("d'1[O][V]{u}");
  blas->add_Matrix("d'2[oo][vv]{u}");
  blas->add_Matrix("d'2[oO][vV]{u}");
  blas->add_Matrix("d'2[OO][VV]{u}");

  // Amplitudes
  blas->add_Matrix("t1[ov]{u}");
  blas->add_Matrix("t1[OV]{u}");
  blas->add_Matrix("t1[o][v]{u}");
  blas->add_Matrix("t1[O][V]{u}");

  blas->add_Matrix("t1_delta[o][v]{u}");

  blas->add_Matrix("t2[oO][vV]{u}");
  blas->add_Matrix("t2[oo][vv]{u}");
  blas->add_Matrix("t2[OO][VV]{u}");

  blas->add_Matrix("t2[ov][ov]{u}");
  blas->add_Matrix("t2[ov][OV]{u}");

  blas->add_Matrix("t2[o][ovv]{u}");
  blas->add_Matrix("t2[o][OvV]{u}");
  blas->add_Matrix("t2[O][oVv]{u}");
  blas->add_Matrix("t2[O][OVV]{u}");

  blas->add_Matrix("t2[v][voo]{u}");
  blas->add_Matrix("t2[v][VoO]{u}");
  blas->add_Matrix("t2[V][vOo]{u}");
  blas->add_Matrix("t2[V][VOO]{u}");

  blas->add_Matrix("tau[oo][vv]{u}");
  blas->add_Matrix("tau[oO][vV]{u}");
  blas->add_Matrix("tau[OO][VV]{u}");

  blas->add_Matrix("t2_delta[oO][vV]{u}");

  blas->add_Matrix("tau2[v][voo]{u}");
  blas->add_Matrix("tau2[v][VoO]{u}");
  blas->add_Matrix("tau2[V][VOO]{u}");
  blas->add_Matrix("tau2[V][vOo]{u}");

  blas->add_Matrix("tau2[o][ovv]{u}");
  blas->add_Matrix("tau2[o][OvV]{u}");
  blas->add_Matrix("tau2[O][oVv]{u}");
  blas->add_Matrix("tau2[O][OVV]{u}");

  blas->add_Matrix("W_mNiJ[oO][oO]{u}");

  // Similarity transformed Hamiltonian
  blas->add_Matrix("t1_eqns[a][a]{u}");
  blas->add_Matrix("t1_eqns[o][v]{u}");
  blas->add_Matrix("t1_eqns[O][V]{u}");
  blas->add_Matrix("t2_eqns[oo][vv]{u}");
  blas->add_Matrix("t2_eqns[oO][vV]{u}");
  blas->add_Matrix("t2_eqns[OO][VV]{u}");

  // F intermediates
  blas->add_Matrix("F_ae[v][v]{u}");
  blas->add_Matrix("F_AE[V][V]{u}");
  blas->add_Matrix("F_mi[o][o]{u}");
  blas->add_Matrix("F_MI[O][O]{u}");
  blas->add_Matrix("F_me[o][v]{u}");
  blas->add_Matrix("F_ME[O][V]{u}");
  blas->add_Matrix("F_me[ov]{u}");
  blas->add_Matrix("F_ME[OV]{u}");
  blas->add_Matrix("F'_mi[o][o]{u}");

  // MRPT2 Intermediates
  blas->add_Matrix("ERef{u}");
  blas->add_Matrix("EPT2{u}");
  blas->add_Matrix("Eaa{u}");
  blas->add_Matrix("Ebb{u}");
  blas->add_Matrix("Eaaaa{u}");
  blas->add_Matrix("Eabab{u}");
  blas->add_Matrix("Ebbbb{u}");

  // A^4
  blas->add_Matrix("<[aa]|[aa]>");
  blas->add_Matrix("<[aa]:[aa]>");

  // A^3O
  blas->add_Matrix("<[o]|[aaa]>");
  blas->add_Matrix("<[o]:[aaa]>");
  blas->add_Matrix("<[oa]|[aa]>");

  // A^3V
  blas->add_Matrix("<[v]|[aaa]>");
  blas->add_Matrix("<[v]:[aaa]>");
  blas->add_Matrix("<[aa]|[va]>");

  // A^2O^2
  blas->add_Matrix("<[oo]|[aa]>");
  blas->add_Matrix("<[oo]:[aa]>");
  blas->add_Matrix("<[o]|[aoa]>");
  blas->add_Matrix("<[o]|[aao]>");

  // A^2OV
  blas->add_Matrix("<[aa]|[ov]>");
  blas->add_Matrix("<[aa]|[vo]>");
  blas->add_Matrix("<[ov]|[aa]>");
  blas->add_Matrix("([ov]|[aa])");
  blas->add_Matrix("([ov]:[aa])");
  blas->add_Matrix("<[oa]|[va]>");
  blas->add_Matrix("<[oa]:[va]>");
  blas->add_Matrix("<[oa]|[av]>");
  blas->add_Matrix("<[v]|[oav]>");
  blas->add_Matrix("<[v]|[oaa]>");
  blas->add_Matrix("<[o]|[vaa]>");
  blas->add_Matrix("<[ov]|[av]>");


  // A^2V^2
  blas->add_Matrix("<[aa]|[vv]>");
  blas->add_Matrix("<[aa]:[vv]>");
  blas->add_Matrix("<[v]|[ava]>");
  blas->add_Matrix("<[v]|[aav]>");
  blas->add_Matrix("<[ov]:[va]>");

  // AO^2V
  blas->add_Matrix("<[a]:[voo]>");
  blas->add_Matrix("<[a]|[voo]>");
  blas->add_Matrix("<[o]|[oav]>");
  blas->add_Matrix("<[o]:[oav]>");
  blas->add_Matrix("<[o]|[ova]>");
  blas->add_Matrix("<[oav]:[v]>");
  blas->add_Matrix("<[ov]|[oa]>");
  blas->add_Matrix("([ov]|[ao])");
  blas->add_Matrix("<[oa]:[vo]>");
  blas->add_Matrix("<[oa]|[vo]>");
  blas->add_Matrix("<[oa]|[ov]>");

  // AOV^2
  blas->add_Matrix("<[a]:[ovv]>");
  blas->add_Matrix("<[a]|[ovv]>");
  blas->add_Matrix("<[ao]|[vv]>");
  blas->add_Matrix("([avo]|[v])");
  blas->add_Matrix("<[va]|[vv]>");
  blas->add_Matrix("([ov]|[va])");
  blas->add_Matrix("<[ov]|[va]>");

  // Effective Hamiltonian in the active space
  blas->add_Matrix("Hia[a][a]{u}");
  blas->add_Matrix("HIA[A][A]{u}");
  blas->add_Matrix("Hijab[aa][aa]{u}");
  blas->add_Matrix("HiJaB[aA][aA]{u}");
  blas->add_Matrix("HIJAB[AA][AA]{u}");

  blas->add_Matrix("HiJaB[oA][aA]{u}");
  blas->add_Matrix("HiJaB[aO][aA]{u}");
  blas->add_Matrix("HiJaB[aA][vA]{u}");
  blas->add_Matrix("HiJaB[aA][aV]{u}");

  blas->add_Matrix("t1_ov[a][v]{u}");
  blas->add_Matrix("t1_OV[A][V]{u}");
  blas->add_Matrix("t1_ov[o][a]{u}");
  blas->add_Matrix("t1_OV[O][A]{u}");

  // Newly added
  blas->add_Matrix("t2_oovv[ao][av]{u}");
  blas->add_Matrix("t2_oOvV[aO][aV]{u}");
  blas->add_Matrix("t2_oOvV[oA][vA]{u}");
  blas->add_Matrix("t2_OOVV[AO][AV]{u}");
  blas->add_Matrix("t2_oOvV[oA][aV]{u}");
  blas->add_Matrix("t2_oOvV[aO][vA]{u}");

  blas->add_Matrix("t2_oovv[a][ovv]{u}");
  blas->add_Matrix("t2_oOvV[a][OvV]{u}");
  blas->add_Matrix("t2_OoVv[A][oVv]{u}");
  blas->add_Matrix("t2_OOVV[A][OVV]{u}");
  blas->add_Matrix("t2_oovv[oo][aa]{u}");
  blas->add_Matrix("t2_oOvV[oO][aA]{u}");
  blas->add_Matrix("t2_OOVV[OO][AA]{u}");
  blas->add_Matrix("t2_oovv[aa][vv]{u}");
  blas->add_Matrix("t2_oOvV[aA][vV]{u}");
  blas->add_Matrix("t2_OOVV[AA][VV]{u}");

  blas->add_Matrix("t2_ovov[aa][ov]{u}");
  blas->add_Matrix("t2_ovov[oa][ov]{u}");
  blas->add_Matrix("t2_ovov[av][ov]{u}");

  blas->add_Matrix("t2_ovOV[aa][OV]{u}");
  blas->add_Matrix("t2_ovOV[oa][OV]{u}");
  blas->add_Matrix("t2_ovOV[ov][AA]{u}");
  blas->add_Matrix("t2_ovOV[av][OV]{u}");

  blas->add_Matrix("t2_oVOv[aA][Ov]{u}");
  blas->add_Matrix("t2_oVOv[oV][Aa]{u}");
  blas->add_Matrix("t2_oVOv[oA][Ov]{u}");
  blas->add_Matrix("t2_oVOv[oV][Av]{u}");
  blas->add_Matrix("t2_OVOV[AA][OV]{u}");

  blas->add_Matrix("t2_oovv[ao][vv]{u}");

  blas->add_Matrix("t2_vvoo[a][voo]{u}");
  blas->add_Matrix("t2_vVoO[a][VoO]{u}");
  blas->add_Matrix("t2_VvOo[A][vOo]{u}");
  blas->add_Matrix("t2_VVOO[A][VOO]{u}");

  blas->add_Matrix("t2_vvoo[v][aaa]{u}");
  blas->add_Matrix("t2_vVoO[v][AaA]{u}");
  blas->add_Matrix("t2_VvOo[V][aAa]{u}");
  blas->add_Matrix("t2_VVOO[V][AAA]{u}");

  blas->add_Matrix("t2_oovv[o][aaa]{u}");
  blas->add_Matrix("t2_oOvV[o][AaA]{u}");
  blas->add_Matrix("t2_OoVv[O][aAa]{u}");
  blas->add_Matrix("t2_OOVV[O][AAA]{u}");

  // T2 A^2O^2
  blas->add_Matrix("t2_OoVv[O][oAa]{u}");

  // T2 A^2OV
  blas->add_Matrix("t2_vVoO[v][AoA]{u}");
  blas->add_Matrix("t2_VvOo[V][aAo]{u}");
  blas->add_Matrix("t2_VvOo[V][vAa]{u}");

  blas->add_Matrix("fock[o][a]{u}");
  blas->add_Matrix("fock[O][A]{u}");
  blas->add_Matrix("fock[a][a]{u}");
  blas->add_Matrix("fock[a][v]{u}");
  blas->add_Matrix("fock[A][V]{u}");
  blas->add_Matrix("fock[A][A]{u}");
  blas->add_Matrix("fock[aa]{u}");
  blas->add_Matrix("fock[AA]{u}");

  blas->add_Matrix("F_ae[a][v]{u}");
  blas->add_Matrix("F_AE[A][V]{u}");
  blas->add_Matrix("F'_ae[a][v]{u}");
  blas->add_Matrix("F'_AE[A][V]{u}");

  // Off-diagonal Fock matrices
  blas->add_Matrix("offdiagonal_F[v][v]{u}");
  blas->add_Matrix("offdiagonal_F[o][o]{u}");

  blas->add_Matrix("F_ae[a][v]{u}");
  blas->add_Matrix("F_AE[A][V]{u}");
  blas->add_Matrix("F'_ae[v][v]{u}");
  blas->add_Matrix("F'_ae[a][v]{u}");
  blas->add_Matrix("F'_AE[A][V]{u}");
  blas->add_Matrix("F'_mi[o][a]{u}");
  blas->add_Matrix("F'_MI[O][A]{u}");

  blas->add_Matrix("tau_oOvV[oO][aA]{u}");
  blas->add_Matrix("tau_oOvV[aA][vV]{u}");
  blas->add_Matrix("tau_oOVv[aA][Vv]{u}");
  blas->add_Matrix("tau_oOvV[oA][vV]{u}");
  blas->add_Matrix("tau_oOvV[oO][vA]{u}");
  blas->add_Matrix("tau_oOVv[oA][Vv]{u}");

  blas->add_Matrix("tau3_ovov[aa][ov]{u}");
  blas->add_Matrix("tau3_ovov[oa][ov]{u}");
  blas->add_Matrix("tau3_ovov[av][ov]{u}");
  blas->add_Matrix("tau3_oVvO[aA][vO]{u}");
  blas->add_Matrix("tau3_oVvO[oA][vO]{u}");
  blas->add_Matrix("tau3_oVvO[aV][vO]{u}");

  blas->add_Matrix("t2_oovv[ao][va]{u}");
  blas->add_Matrix("t2_oovv[oo][va]{u}");

  blas->add_Matrix("t2_oOvV[aO][vA]{u}");
  blas->add_Matrix("t2_oOvV[oO][vA]{u}");
  blas->add_Matrix("t2_oOvV[o][AvA]{u}");
  blas->add_Matrix("t2_oOvV[aO][vV]{u}");

  blas->add_Matrix("t2_OoVv[O][aAv]{u}");


  blas->add_Matrix("W_mNiJ[oO][aA]{u}");
  blas->add_Matrix("W_mNiJ[oO][oA]{u}");

  blas->add_Matrix("W_jbme[aa][ov]{u}");
  blas->add_Matrix("W_jbme[oa][ov]{u}");
  blas->add_Matrix("W_jbme[av][ov]{u}");

  blas->add_Matrix("W_jbME[aa][OV]{u}");
  blas->add_Matrix("W_jbME[oa][OV]{u}");
  blas->add_Matrix("W_jbME[av][OV]{u}");

  blas->add_Matrix("W_jBmE[aA][oV]{u}");
  blas->add_Matrix("W_jBmE[oA][oV]{u}");
  blas->add_Matrix("W_jBmE[aV][oV]{u}");

  blas->add_Matrix("Z_iJaM[aAa][O]{u}");
  blas->add_Matrix("Z_iJAm[aAA][o]{u}");
  blas->add_Matrix("Z_iJaM[oAa][O]{u}");
  blas->add_Matrix("Z_iJAm[oAA][o]{u}");
  blas->add_Matrix("Z_iJaM[aAv][O]{u}");

  blas->add_Matrix("t1t1_iame[aa][ov]{u}");
  blas->add_Matrix("t1t1_iame[oa][ov]{u}");
  blas->add_Matrix("t1t1_iame[av][ov]{u}");
  blas->add_Matrix("t1t1_IAME[AA][OV]{u}");
  blas->add_Matrix("t1t1_iAMe[aA][Ov]{u}");
  blas->add_Matrix("t1t1_iAMe[oV][Aa]{u}");
  blas->add_Matrix("t1t1_iAMe[oA][Ov]{u}");
  blas->add_Matrix("t1t1_iAMe[oV][Av]{u}");

}

}} /* End Namespaces */
