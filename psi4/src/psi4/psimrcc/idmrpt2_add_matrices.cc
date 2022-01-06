/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

#include "blas.h"
#include "idmrpt2.h"

namespace psi {
namespace psimrcc {

void IDMRPT2::add_matrices() {
    // O^4
    wfn_->blas()->add_Matrix("<[oo]:[oo]>");
    wfn_->blas()->add_Matrix("<[oo]|[oo]>");

    // O^2V^2
    wfn_->blas()->add_Matrix("<[oo]:[vv]>");
    wfn_->blas()->add_Matrix("<[oo]|[vv]>");

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

    wfn_->blas()->add_Matrix("t2[oO][vV]{u}");
    wfn_->blas()->add_Matrix("t2[oo][vv]{u}");
    wfn_->blas()->add_Matrix("t2[OO][VV]{u}");

    wfn_->blas()->add_Matrix("t2[o][ovv]{u}");
    wfn_->blas()->add_Matrix("t2[o][OvV]{u}");
    wfn_->blas()->add_Matrix("t2[O][oVv]{u}");
    wfn_->blas()->add_Matrix("t2[O][OVV]{u}");

    wfn_->blas()->add_Matrix("t2[v][voo]{u}");
    wfn_->blas()->add_Matrix("t2[v][VoO]{u}");
    wfn_->blas()->add_Matrix("t2[V][vOo]{u}");
    wfn_->blas()->add_Matrix("t2[V][VOO]{u}");

    // Similarity transformed Hamiltonian
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

    // A^3V
    wfn_->blas()->add_Matrix("<[v]|[aaa]>");
    wfn_->blas()->add_Matrix("<[v]:[aaa]>");

    // A^2O^2
    wfn_->blas()->add_Matrix("<[oo]|[aa]>");
    wfn_->blas()->add_Matrix("<[oo]:[aa]>");

    // A^2OV
    wfn_->blas()->add_Matrix("<[aa]|[ov]>");
    //   wfn_->blas()->add_Matrix("<[ov]|[aa]>");
    wfn_->blas()->add_Matrix("([ov]|[aa])");
    wfn_->blas()->add_Matrix("([ov]:[aa])");

    // A^2V^2
    wfn_->blas()->add_Matrix("<[aa]|[vv]>");
    wfn_->blas()->add_Matrix("<[aa]:[vv]>");

    // AO^2V
    wfn_->blas()->add_Matrix("<[a]:[voo]>");
    wfn_->blas()->add_Matrix("<[a]|[voo]>");

    // AOV^2
    wfn_->blas()->add_Matrix("<[a]:[ovv]>");
    wfn_->blas()->add_Matrix("<[a]|[ovv]>");

    // Effective Hamiltonian in the active space
    wfn_->blas()->add_Matrix("Hia[a][a]{u}");
    wfn_->blas()->add_Matrix("HIA[A][A]{u}");
    wfn_->blas()->add_Matrix("Hijab[aa][aa]{u}");
    wfn_->blas()->add_Matrix("HiJaB[aA][aA]{u}");
    wfn_->blas()->add_Matrix("HIJAB[AA][AA]{u}");

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
    wfn_->blas()->add_Matrix("t2_OoVv[A][oVv]{o}");
    wfn_->blas()->add_Matrix("t2_OOVV[A][OVV]{o}");
    wfn_->blas()->add_Matrix("t2_oovv[oo][aa]{u}");
    wfn_->blas()->add_Matrix("t2_oOvV[oO][aA]{u}");
    wfn_->blas()->add_Matrix("t2_OOVV[OO][AA]{u}");
    wfn_->blas()->add_Matrix("t2_oovv[aa][vv]{u}");
    wfn_->blas()->add_Matrix("t2_oOvV[aA][vV]{u}");
    wfn_->blas()->add_Matrix("t2_OOVV[AA][VV]{u}");

    wfn_->blas()->add_Matrix("t2_ovov[aa][ov]{u}");
    wfn_->blas()->add_Matrix("t2_ovOV[aa][OV]{u}");
    wfn_->blas()->add_Matrix("t2_oVOv[aA][Ov]{u}");
    wfn_->blas()->add_Matrix("t2_OVOV[AA][OV]{u}");
    wfn_->blas()->add_Matrix("t2_ovOV[ov][AA]{u}");
    wfn_->blas()->add_Matrix("t2_oVOv[oV][Aa]{u}");

    wfn_->blas()->add_Matrix("t2_vvoo[a][voo]{u}");
    wfn_->blas()->add_Matrix("t2_vVoO[a][VoO]{u}");
    wfn_->blas()->add_Matrix("t2_VvOo[A][vOo]{o}");
    wfn_->blas()->add_Matrix("t2_VVOO[A][VOO]{o}");

    wfn_->blas()->add_Matrix("t2_vvoo[v][aaa]{u}");
    wfn_->blas()->add_Matrix("t2_vVoO[v][AaA]{u}");
    wfn_->blas()->add_Matrix("t2_VvOo[V][aAa]{u}");
    wfn_->blas()->add_Matrix("t2_VVOO[V][AAA]{o}");

    wfn_->blas()->add_Matrix("t2_oovv[o][aaa]{u}");
    wfn_->blas()->add_Matrix("t2_oOvV[o][AaA]{u}");
    wfn_->blas()->add_Matrix("t2_OoVv[O][aAa]{u}");
    wfn_->blas()->add_Matrix("t2_OOVV[O][AAA]{u}");

    wfn_->blas()->add_Matrix("fock[o][a]{u}");
    wfn_->blas()->add_Matrix("fock[O][A]{o}");
    wfn_->blas()->add_Matrix("fock[a][a]{u}");
    wfn_->blas()->add_Matrix("fock[a][v]{u}");
    wfn_->blas()->add_Matrix("fock[A][V]{o}");
    wfn_->blas()->add_Matrix("fock[A][A]{o}");
    wfn_->blas()->add_Matrix("fock[aa]{u}");
    wfn_->blas()->add_Matrix("fock[AA]{u}");
    wfn_->blas()->add_Matrix("fock[ff]{u}");
    wfn_->blas()->add_Matrix("fock[FF]{u}");

    wfn_->blas()->add_Matrix("factor_mk{u}");
    wfn_->blas()->add_Matrix("neg_factor_mk{u}");
}

}  // namespace psimrcc
}  // namespace psi
