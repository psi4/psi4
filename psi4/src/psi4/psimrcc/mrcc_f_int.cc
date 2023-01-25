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

// clang-format off
/**
 * @brief Computes the contractions
 * \f[ F_{ae} = (1-\delta_{ae})f_{ae} - \frac{1}{2} \sum_m f_{me} t_m^a + \sum_{mf} t_m^f <ma||fe> - \frac{1}{2} \sum_{mnf} \tilde{\tau}_{in}^{af} <mn||ef> \f]
 * \f[ F_{mi} = (1-\delta_{mi})f_{mi} + \frac{1}{2} \sum_e f_{me} t_i^e + \sum_{en} t_n^e <mn||ie> + \frac{1}{2} \sum_{nef} \tilde{\tau}_{in}^{ef} <mn||ef> \f]
 * \f[ F_{me} = f_{me} + \sum_{nf} t_n^f <mn||ef> \f]
 * \f[ F'_{ae} = F_{ae} - \frac{1}{2} \sum_m t_m^a F_{me} \f]
 * \f[ F'_{mi} = F_{mi} + \frac{1}{2} \sum_e t_i^e F_{me} \f]
 * as described in J. Phys. Chem. vol. 94, pg. 4334 (1991).
 * See J. Phys. Chem. vol. 127, 024102 (2007) supplementary material for the spin-factored equations.
 */
// clang-format on
void CCMRCC::build_F_intermediates() {
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

    if (triples_type == ccsd_t) {
        build_F2_me_intermediates();
        build_F2_ME_intermediates();
    }
}

void CCMRCC::build_F_ae_intermediates() {
    // Closed-shell
    // Add the VV Fock matrix with the diagonal terms zeroed
    wfn_->blas()->append("F_ae[v][v]{c} = fock[v][v]{c}");
    wfn_->blas()->append_zero_two_diagonal("F_ae[v][v]{c}");

    wfn_->blas()->append("F_ae[v][v]{c} += -1/2 t1[o][v]{c} 1@1 fock[o][v]{c}");

    wfn_->blas()->append("F_ae[v][v]{c} += #12# ([ov]:[vv]) 1@1 t1[ov]{c}");
    wfn_->blas()->append("F_ae[v][v]{c} += #12# ([ov]|[vv]) 1@1 t1[ov]{c} ");

    wfn_->blas()->append("F_ae[v][v]{c} += -1/2 tau2[v][voo]{c} 2@2 <[v]:[voo]>");
    wfn_->blas()->append("F_ae[v][v]{c} += - tau2[v][VoO]{c} 2@2 <[v]|[voo]>");

    // Open-shell
    // Add the VV Fock matrix with the diagonal terms zeroed
    wfn_->blas()->append("F_ae[v][v]{o} = fock[v][v]{o}");
    wfn_->blas()->append_zero_two_diagonal("F_ae[v][v]{o}");

    wfn_->blas()->append("F_ae[v][v]{o} += -1/2 t1[o][v]{o} 1@1 fock[o][v]{o}");
    wfn_->blas()->append("F_ae[v][v]{o} += #12# ([ov]:[vv]) 1@1 t1[ov]{o}");
    wfn_->blas()->append("F_ae[v][v]{o} += #12# ([ov]|[vv]) 1@1 t1[OV]{o} ");

    wfn_->blas()->append("F_ae[v][v]{o} += -1/2 tau2[v][voo]{o} 2@2 <[v]:[voo]>");
    wfn_->blas()->append("F_ae[v][v]{o} += - tau2[v][VoO]{o} 2@2 <[v]|[voo]>");
}

void CCMRCC::build_F_AE_intermediates() {
    // Open-shell
    // Add the VV Fock matrix with the diagonal terms zeroed
    wfn_->blas()->append("F_AE[V][V]{o} = fock[V][V]{o}");

    wfn_->blas()->append_zero_two_diagonal("F_AE[V][V]{o}");

    wfn_->blas()->append("F_AE[V][V]{o} += -1/2 t1[O][V]{o} 1@1 fock[O][V]{o}");

    wfn_->blas()->append("F_AE[V][V]{o} += #12# ([ov]:[vv]) 1@1 t1[OV]{o}");
    wfn_->blas()->append("F_AE[V][V]{o} += #12# ([ov]|[vv]) 1@1 t1[ov]{o} ");

    wfn_->blas()->append("F_AE[V][V]{o} += -1/2 tau2[V][VOO]{o} 2@2 <[v]:[voo]>");
    wfn_->blas()->append("F_AE[V][V]{o} += - tau2[V][vOo]{o} 2@2 <[v]|[voo]>");
}

void CCMRCC::build_F_mi_intermediates() {
    // Closed-shell
    // Add the VV Fock matrix with the diagonal terms zeroed
    wfn_->blas()->append("F_mi[o][o]{c} = fock[o][o]{c}");
    wfn_->blas()->append_zero_two_diagonal("F_mi[o][o]{c}");

    wfn_->blas()->append("F_mi[o][o]{c} += 1/2 fock[o][v]{c} 2@2 t1[o][v]{c}");

    wfn_->blas()->append("F_mi[o][o]{c} += #12# ([oo]:[ov]) 2@1 t1[ov]{c}");
    wfn_->blas()->append("F_mi[o][o]{c} += #12# ([oo]|[ov]) 2@1 t1[ov]{c} ");

    wfn_->blas()->append("F_mi[o][o]{c} += 1/2  <[o]:[ovv]> 2@2 tau2[o][ovv]{c}");
    wfn_->blas()->append("F_mi[o][o]{c} +=      <[o]|[ovv]> 2@2 tau2[o][OvV]{c} ");

    // Open-shell
    // Add the VV Fock matrix with the diagonal terms zeroed
    wfn_->blas()->append("F_mi[o][o]{o} = fock[o][o]{o}");
    wfn_->blas()->append_zero_two_diagonal("F_mi[o][o]{o}");

    wfn_->blas()->append("F_mi[o][o]{o} += 1/2 fock[o][v]{o} 2@2 t1[o][v]{o}");

    wfn_->blas()->append("F_mi[o][o]{o} += #12# ([oo]:[ov]) 2@1 t1[ov]{o}");
    wfn_->blas()->append("F_mi[o][o]{o} += #12# ([oo]|[ov]) 2@1 t1[OV]{o} ");

    wfn_->blas()->append("F_mi[o][o]{o} += 1/2  <[o]:[ovv]> 2@2 tau2[o][ovv]{o}");
    wfn_->blas()->append("F_mi[o][o]{o} +=      <[o]|[ovv]> 2@2 tau2[o][OvV]{o} ");
}

void CCMRCC::build_F_MI_intermediates() {
    // Open-shell
    // Add the VV Fock matrix with the diagonal terms zeroed
    wfn_->blas()->append("F_MI[O][O]{o} = fock[O][O]{o}");

    wfn_->blas()->append_zero_two_diagonal("F_MI[O][O]{o}");

    wfn_->blas()->append("F_MI[O][O]{o} += 1/2 fock[O][V]{o} 2@2 t1[O][V]{o}");

    wfn_->blas()->append("F_MI[O][O]{o} += #12# ([oo]:[ov]) 2@1 t1[OV]{o}");
    wfn_->blas()->append("F_MI[O][O]{o} += #12# ([oo]|[ov]) 2@1 t1[ov]{o} ");

    wfn_->blas()->append("F_MI[O][O]{o} += 1/2  <[o]:[ovv]> 2@2 tau2[O][OVV]{o}");
    wfn_->blas()->append("F_MI[O][O]{o} +=      <[o]|[ovv]> 2@2 tau2[O][oVv]{o} ");
}

void CCMRCC::build_F_me_intermediates() {
    // Closed-Shell
    // Add the VV Fock matrix with the diagonal terms zeroed
    wfn_->blas()->append("F_me[o][v]{c} = fock[o][v]{c}");

    wfn_->blas()->append("F_me[o][v]{c} += #12# ([ov]:[ov]) 2@1 t1[ov]{c}");
    wfn_->blas()->append("F_me[o][v]{c} += #12# ([ov]|[ov]) 2@1 t1[ov]{c} ");

    wfn_->blas()->append("F_me[ov]{c} = #12# F_me[o][v]{c}");

    // Open-Shell
    // Add the VV Fock matrix with the diagonal terms zeroed
    wfn_->blas()->append("F_me[o][v]{o} = fock[o][v]{o}");

    wfn_->blas()->append("F_me[o][v]{o} += #12# ([ov]:[ov]) 2@1 t1[ov]{o}");
    wfn_->blas()->append("F_me[o][v]{o} += #12# ([ov]|[ov]) 2@1 t1[OV]{o} ");

    wfn_->blas()->append("F_me[ov]{o} = #12# F_me[o][v]{o}");
}

void CCMRCC::build_F_ME_intermediates() {
    if (triples_type >= ccsd_t) {
        wfn_->blas()->append("F_ME[O][V]{c} = fock[O][V]{c}");

        wfn_->blas()->append("F_ME[O][V]{c} += #12# ([ov]:[ov]) 2@1 t1[OV]{c}");
        wfn_->blas()->append("F_ME[O][V]{c} += #12# ([ov]|[ov]) 2@1 t1[OV]{c} ");
    }

    // Open-Shell
    // Add the VV Fock matrix with the diagonal terms zeroed
    wfn_->blas()->append("F_ME[O][V]{o} = fock[O][V]{o}");

    wfn_->blas()->append("F_ME[O][V]{o} += #12# ([ov]:[ov]) 2@1 t1[OV]{o}");
    wfn_->blas()->append("F_ME[O][V]{o} += #12# ([ov]|[ov]) 2@1 t1[ov]{o} ");

    wfn_->blas()->append("F_ME[OV]{o} = #12# F_ME[O][V]{o}");
}

void CCMRCC::build_F_prime_ae_intermediates() {
    // Closed-Shell + Open-Shell Spin-Adapted Form
    // Add the VV Fock matrix with the diagonal terms zeroed
    wfn_->blas()->append("F'_ae[v][v]{u}  = F_ae[v][v]{u}");
    wfn_->blas()->append("F'_ae[v][v]{u} += #12# -1/2 t1[o][v]{u} 1@1 F_me[o][v]{u}");
}

void CCMRCC::build_F_prime_AE_intermediates() {
    // Open-Shell
    // Add the VV Fock matrix with the diagonal terms zeroed
    wfn_->blas()->append("F'_AE[V][V]{o}  = F_AE[V][V]{o}");
    wfn_->blas()->append("F'_AE[V][V]{o} += #12# -1/2 t1[O][V]{o} 1@1 F_ME[O][V]{o}");
}

void CCMRCC::build_F_prime_mi_intermediates() {
    // Closed-Shell + Open-Shell Spin-Adapted Form
    // Add the VV Fock matrix with the diagonal terms zeroed
    wfn_->blas()->append("F'_mi[o][o]{u}  = F_mi[o][o]{u}");
    wfn_->blas()->append("F'_mi[o][o]{u} += #12# 1/2 F_me[o][v]{u} 2@2 t1[o][v]{u}");
}

void CCMRCC::build_F_prime_MI_intermediates() {
    // Open-Shell
    // Add the VV Fock matrix with the diagonal terms zeroed
    wfn_->blas()->append("F'_MI[O][O]{o}  = F_MI[O][O]{o}");
    wfn_->blas()->append("F'_MI[O][O]{o} += #12# 1/2 F_ME[O][V]{o} 2@2 t1[O][V]{o}");
}

void CCMRCC::build_F2_me_intermediates() {
    //  // Closed-Shell
    //  // Add the VV Fock matrix with the diagonal terms zeroed
    //  if(triples_type==ccsdt_1a){
    //    wfn_->blas()->solve("F2_me[o][v]{c} = fock[o][v]{c}");
    //  }else{
    //    wfn_->blas()->solve("F2_me[o][v]{c} = F_me[o][v]{c}");
    //  }
    //
    //  // Open-Shell
    //  // Add the VV Fock matrix with the diagonal terms zeroed
    //  if(triples_type==ccsdt_1a){
    //    wfn_->blas()->solve("F2_me[o][v]{o} = fock[o][v]{o}");
    //  }else{
    //    wfn_->blas()->solve("F2_me[o][v]{o} = F_me[o][v]{o}");
    //  }
    //  );
}

void CCMRCC::build_F2_ME_intermediates() {
    //  // Closed-Shell
    //  if(triples_type==ccsdt_1a){
    //    wfn_->blas()->solve("F2_ME[O][V]{c} = fock[O][V]{c}");
    //  }else{
    //    wfn_->blas()->solve("F2_ME[O][V]{c} = F_ME[O][V]{c}");
    //  }
    //
    //  // Open-Shell
    //  if(triples_type==ccsdt_1a){
    //    wfn_->blas()->solve("F2_ME[O][V]{o} = fock[O][V]{o}");
    //  }else{
    //    wfn_->blas()->solve("F2_ME[O][V]{o} = F_ME[O][V]{o}");
    //  }
}

}  // namespace psimrcc
}  // namespace psi
