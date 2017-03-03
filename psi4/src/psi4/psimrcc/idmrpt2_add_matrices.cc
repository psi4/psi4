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

/***************************************************************************
 *  PSIMRCC : Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/

#include "blas.h"
#include "idmrpt2.h"

namespace psi{ namespace psimrcc{

void IDMRPT2::add_matrices()
{
  // O^4
  blas->add_Matrix("<[oo]:[oo]>");
  blas->add_Matrix("<[oo]|[oo]>");

  // O^2V^2
  blas->add_Matrix("<[oo]:[vv]>");
  blas->add_Matrix("<[oo]|[vv]>");

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

  blas->add_Matrix("t2[oO][vV]{u}");
  blas->add_Matrix("t2[oo][vv]{u}");
  blas->add_Matrix("t2[OO][VV]{u}");

  blas->add_Matrix("t2[o][ovv]{u}");
  blas->add_Matrix("t2[o][OvV]{u}");
  blas->add_Matrix("t2[O][oVv]{u}");
  blas->add_Matrix("t2[O][OVV]{u}");

  blas->add_Matrix("t2[v][voo]{u}");
  blas->add_Matrix("t2[v][VoO]{u}");
  blas->add_Matrix("t2[V][vOo]{u}");
  blas->add_Matrix("t2[V][VOO]{u}");

  // Similarity transformed Hamiltonian
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

  // A^3V
  blas->add_Matrix("<[v]|[aaa]>");
  blas->add_Matrix("<[v]:[aaa]>");

  // A^2O^2
  blas->add_Matrix("<[oo]|[aa]>");
  blas->add_Matrix("<[oo]:[aa]>");

  // A^2OV
  blas->add_Matrix("<[aa]|[ov]>");
//   blas->add_Matrix("<[ov]|[aa]>");
  blas->add_Matrix("([ov]|[aa])");
  blas->add_Matrix("([ov]:[aa])");

  // A^2V^2
  blas->add_Matrix("<[aa]|[vv]>");
  blas->add_Matrix("<[aa]:[vv]>");

  // AO^2V
  blas->add_Matrix("<[a]:[voo]>");
  blas->add_Matrix("<[a]|[voo]>");

  // AOV^2
  blas->add_Matrix("<[a]:[ovv]>");
  blas->add_Matrix("<[a]|[ovv]>");


  // Effective Hamiltonian in the active space
  blas->add_Matrix("Hia[a][a]{u}");
  blas->add_Matrix("HIA[A][A]{u}");
  blas->add_Matrix("Hijab[aa][aa]{u}");
  blas->add_Matrix("HiJaB[aA][aA]{u}");
  blas->add_Matrix("HIJAB[AA][AA]{u}");

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
  blas->add_Matrix("t2_OoVv[A][oVv]{o}");
  blas->add_Matrix("t2_OOVV[A][OVV]{o}");
  blas->add_Matrix("t2_oovv[oo][aa]{u}");
  blas->add_Matrix("t2_oOvV[oO][aA]{u}");
  blas->add_Matrix("t2_OOVV[OO][AA]{u}");
  blas->add_Matrix("t2_oovv[aa][vv]{u}");
  blas->add_Matrix("t2_oOvV[aA][vV]{u}");
  blas->add_Matrix("t2_OOVV[AA][VV]{u}");

  blas->add_Matrix("t2_ovov[aa][ov]{u}");
  blas->add_Matrix("t2_ovOV[aa][OV]{u}");
  blas->add_Matrix("t2_oVOv[aA][Ov]{u}");
  blas->add_Matrix("t2_OVOV[AA][OV]{u}");
  blas->add_Matrix("t2_ovOV[ov][AA]{u}");
  blas->add_Matrix("t2_oVOv[oV][Aa]{u}");

  blas->add_Matrix("t2_vvoo[a][voo]{u}");
  blas->add_Matrix("t2_vVoO[a][VoO]{u}");
  blas->add_Matrix("t2_VvOo[A][vOo]{o}");
  blas->add_Matrix("t2_VVOO[A][VOO]{o}");

  blas->add_Matrix("t2_vvoo[v][aaa]{u}");
  blas->add_Matrix("t2_vVoO[v][AaA]{u}");
  blas->add_Matrix("t2_VvOo[V][aAa]{u}");
  blas->add_Matrix("t2_VVOO[V][AAA]{o}");

  blas->add_Matrix("t2_oovv[o][aaa]{u}");
  blas->add_Matrix("t2_oOvV[o][AaA]{u}");
  blas->add_Matrix("t2_OoVv[O][aAa]{u}");
  blas->add_Matrix("t2_OOVV[O][AAA]{u}");


  blas->add_Matrix("fock[o][a]{u}");
  blas->add_Matrix("fock[O][A]{o}");
  blas->add_Matrix("fock[a][a]{u}");
  blas->add_Matrix("fock[a][v]{u}");
  blas->add_Matrix("fock[A][V]{o}");
  blas->add_Matrix("fock[A][A]{o}");
  blas->add_Matrix("fock[aa]{u}");
  blas->add_Matrix("fock[AA]{u}");
  blas->add_Matrix("fock[ff]{u}");
  blas->add_Matrix("fock[FF]{u}");

  blas->add_Matrix("factor_mk{u}");
  blas->add_Matrix("neg_factor_mk{u}");
}

}} /* End Namespaces */
