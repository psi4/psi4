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
#include "psi4/libpsi4util/libpsi4util.h"

#include "mp2_ccsd.h"
#include "matrix.h"
#include "blas.h"

namespace psi {
namespace psimrcc {

void MP2_CCSD::build_F_intermediates() {
    build_F_ae_intermediates();
    build_F_AE_intermediates();

    build_F_mi_intermediates();
    build_F_MI_intermediates();

    build_F_me_intermediates();
    build_F_ME_intermediates();

    build_F_prime_ae_intermediates();
    build_F_prime_AE_intermediates();

    build_F_prime_mi_intermediates();
    build_F_prime_MI_intermediates();
}

void MP2_CCSD::build_offdiagonal_F() {
    wfn_->blas()->solve("offdiagonal_F[v][v]{u} = fock[v][v]{u}");
    wfn_->blas()->solve_zero_two_diagonal("offdiagonal_F[v][v]{u}");

    wfn_->blas()->solve("offdiagonal_F[o][o]{u} = fock[o][o]{u}");
    wfn_->blas()->solve_zero_two_diagonal("offdiagonal_F[o][o]{u}");
}

void MP2_CCSD::build_F_ae_intermediates() {
    // Add the VV Fock matrix with the diagonal terms zeroed
    wfn_->blas()->solve("F_ae[v][v]{u} = fock[v][v]{u}");
    wfn_->blas()->solve_zero_two_diagonal("F_ae[v][v]{u}");

    wfn_->blas()->solve("F_ae[v][v]{u} += -1/2 t1[o][v]{u} 1@1 fock[o][v]{u}");

    wfn_->blas()->solve("F_ae[v][v]{u} += #12# ([ov]:[vv]) 1@1 t1[ov]{u}");
    wfn_->blas()->solve("F_ae[v][v]{u} += #12# ([ov]|[vv]) 1@1 t1[OV]{u} ");

    wfn_->blas()->solve("F_ae[v][v]{u} += -1/2 tau2[v][voo]{u} 2@2 <[v]:[voo]>");
    wfn_->blas()->solve("F_ae[v][v]{u} += - tau2[v][VoO]{u} 2@2 <[v]|[voo]>");

    wfn_->blas()->reduce_spaces("F_ae[a][v]{u}", "F_ae[v][v]{u}");
}

void MP2_CCSD::build_F_AE_intermediates() {
    // Add the VV Fock matrix with the diagonal terms zeroed
    wfn_->blas()->solve("F_AE[V][V]{u} = fock[V][V]{u}");

    wfn_->blas()->solve_zero_two_diagonal("F_AE[V][V]{u}");

    wfn_->blas()->solve("F_AE[V][V]{u} += -1/2 t1[O][V]{u} 1@1 fock[O][V]{u}");

    wfn_->blas()->solve("F_AE[V][V]{u} += #12# ([ov]:[vv]) 1@1 t1[OV]{u}");
    wfn_->blas()->solve("F_AE[V][V]{u} += #12# ([ov]|[vv]) 1@1 t1[ov]{u} ");

    wfn_->blas()->solve("F_AE[V][V]{u} += -1/2 tau2[V][VOO]{u} 2@2 <[v]:[voo]>");
    wfn_->blas()->solve("F_AE[V][V]{u} += - tau2[V][vOo]{u} 2@2 <[v]|[voo]>");

    wfn_->blas()->reduce_spaces("F_AE[A][V]{u}", "F_AE[V][V]{u}");
}

void MP2_CCSD::build_F_mi_intermediates() {
    // Add the VV Fock matrix with the diagonal terms zeroed
    wfn_->blas()->solve("F_mi[o][o]{u} = fock[o][o]{u}");
    wfn_->blas()->solve_zero_two_diagonal("F_mi[o][o]{u}");

    wfn_->blas()->solve("F_mi[o][o]{u} += 1/2 fock[o][v]{u} 2@2 t1[o][v]{u}");

    wfn_->blas()->solve("F_mi[o][o]{u} += #12# ([oo]:[ov]) 2@1 t1[ov]{u}");
    wfn_->blas()->solve("F_mi[o][o]{u} += #12# ([oo]|[ov]) 2@1 t1[OV]{u} ");

    wfn_->blas()->solve("F_mi[o][o]{u} += 1/2  <[o]:[ovv]> 2@2 tau2[o][ovv]{u}");
    wfn_->blas()->solve("F_mi[o][o]{u} +=      <[o]|[ovv]> 2@2 tau2[o][OvV]{u} ");
}

void MP2_CCSD::build_F_MI_intermediates() {
    // Open-shell
    // Add the VV Fock matrix with the diagonal terms zeroed
    wfn_->blas()->solve("F_MI[O][O]{u} = fock[O][O]{u}");

    wfn_->blas()->solve_zero_two_diagonal("F_MI[O][O]{u}");

    wfn_->blas()->solve("F_MI[O][O]{u} += 1/2 fock[O][V]{u} 2@2 t1[O][V]{u}");

    wfn_->blas()->solve("F_MI[O][O]{u} += #12# ([oo]:[ov]) 2@1 t1[OV]{u}");
    wfn_->blas()->solve("F_MI[O][O]{u} += #12# ([oo]|[ov]) 2@1 t1[ov]{u} ");

    wfn_->blas()->solve("F_MI[O][O]{u} += 1/2  <[o]:[ovv]> 2@2 tau2[O][OVV]{u}");
    wfn_->blas()->solve("F_MI[O][O]{u} +=      <[o]|[ovv]> 2@2 tau2[O][oVv]{u} ");
}

void MP2_CCSD::build_F_prime_mi_intermediates() {
    wfn_->blas()->solve("F'_mi[o][o]{u} = F_mi[o][o]{u}");
    wfn_->blas()->solve("F'_mi[o][o]{u} += #12# 1/2 F_me[o][v]{u} 2@2 t1[o][v]{u}");
    wfn_->blas()->reduce_spaces("F'_mi[o][a]{u}", "F'_mi[o][o]{u}");
}

void MP2_CCSD::build_F_prime_MI_intermediates() {
    wfn_->blas()->reduce_spaces("F'_MI[O][A]{u}", "F_MI[O][O]{u}");
    wfn_->blas()->solve("F'_MI[O][A]{u} += #12# 1/2 F_ME[O][V]{u} 2@2 t1_OV[A][V]{u}");
}

void MP2_CCSD::build_F_me_intermediates() {
    // Add the VV Fock matrix with the diagonal terms zeroed
    wfn_->blas()->solve("F_me[o][v]{u} = fock[o][v]{u}");

    wfn_->blas()->solve("F_me[o][v]{u} += #12# ([ov]:[ov]) 2@1 t1[ov]{u}");
    wfn_->blas()->solve("F_me[o][v]{u} += #12# ([ov]|[ov]) 2@1 t1[OV]{u} ");

    wfn_->blas()->solve("F_me[ov]{u} = #12# F_me[o][v]{u}");
}

void MP2_CCSD::build_F_ME_intermediates() {
    // Add the VV Fock matrix with the diagonal terms zeroed
    wfn_->blas()->solve("F_ME[O][V]{u} = fock[O][V]{u}");

    wfn_->blas()->solve("F_ME[O][V]{u} += #12# ([ov]:[ov]) 2@1 t1[OV]{u}");
    wfn_->blas()->solve("F_ME[O][V]{u} += #12# ([ov]|[ov]) 2@1 t1[ov]{u} ");

    wfn_->blas()->solve("F_ME[OV]{u} = #12# F_ME[O][V]{u}");
}

void MP2_CCSD::build_F_prime_ae_intermediates() {
    // Closed-Shell + Open-Shell Spin-Adapted Form
    // Add the VV Fock matrix with the diagonal terms zeroed
    //   wfn_->blas()->solve("F'_ae[a][v]{u}  = F_ae[a][v]{u}");
    //   wfn_->blas()->solve("F'_ae[a][v]{u} += #12# -1/2 t1_ov[o][a]{u} 1@1 F_me[o][v]{u}");

    wfn_->blas()->solve("F'_ae[v][v]{u}  = F_ae[v][v]{u}");
    wfn_->blas()->solve("F'_ae[v][v]{u} += #12# -1/2 t1[o][v]{u} 1@1 F_me[o][v]{u}");

    wfn_->blas()->reduce_spaces("F'_ae[a][v]{u}", "F'_ae[v][v]{u}");
}

void MP2_CCSD::build_F_prime_AE_intermediates() {
    // Open-Shell
    // Add the VV Fock matrix with the diagonal terms zeroed
    wfn_->blas()->solve("F'_AE[A][V]{u}  = F_AE[A][V]{u}");
    wfn_->blas()->solve("F'_AE[A][V]{u} += #12# -1/2 t1_OV[O][A]{u} 1@1 F_ME[O][V]{u}");
}

}  // namespace psimrcc
}  // namespace psi
