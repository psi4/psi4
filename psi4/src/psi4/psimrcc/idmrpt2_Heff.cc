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
 *  PSIMRCC
 *  Copyright (C) 2007 by Francesco Evangelista and Andrew Simmonett
 *  frank@ccc.uga.edu   andysim@ccc.uga.edu
 *  A multireference coupled cluster code
 ***************************************************************************/
#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/libpsi4util/libpsi4util.h"

#include "blas.h"
#include "idmrpt2.h"
#include "matrix.h"

namespace psi {
namespace psimrcc {

void IDMRPT2::build_Heff_mrpt2_offdiagonal() {
    build_Heff_uv();
    build_Heff_UV();
    build_Heff_uVxY();
    build_Heff_uvxy();
    build_Heff_UVXY();

    intvec occ_to_act = wfn_->moinfo()->get_occ_to_actv();
    intvec vir_to_act = wfn_->moinfo()->get_vir_to_actv();

    for (int i = 0; i < wfn_->moinfo()->get_ref_size(AllRefs); ++i) {
        int i_unique = wfn_->moinfo()->get_ref_number(i);
        // Find the off_diagonal elements for reference i
        // Loop over reference j (in a safe way)
        for (int j = 0; j < wfn_->moinfo()->get_ref_size(AllRefs); j++) {
            if (i != j) {
                std::vector<std::pair<int, int> > alpha_internal_excitation =
                    wfn_->moinfo()->get_alpha_internal_excitation(i, j);
                std::vector<std::pair<int, int> > beta_internal_excitation =
                    wfn_->moinfo()->get_beta_internal_excitation(i, j);
                double sign_internal_excitation = wfn_->moinfo()->get_sign_internal_excitation(i, j);

                double element = 0.0;
                if (i == i_unique) {
                    // Set alpha-alpha single excitations
                    if ((alpha_internal_excitation.size() == 1) && (beta_internal_excitation.size() == 0))
                        element = sign_internal_excitation *
                                  wfn_->blas()
                                      ->get_MatTmp("Hia[a][a]", i_unique, none)
                                      ->get_two_address_element(occ_to_act[alpha_internal_excitation[0].first],
                                                                vir_to_act[alpha_internal_excitation[0].second]);

                    // Set beta-beta single excitations
                    if ((alpha_internal_excitation.size() == 0) && (beta_internal_excitation.size() == 1))
                        element = sign_internal_excitation *
                                  wfn_->blas()
                                      ->get_MatTmp("HIA[A][A]", i_unique, none)
                                      ->get_two_address_element(occ_to_act[beta_internal_excitation[0].first],
                                                                vir_to_act[beta_internal_excitation[0].second]);

                    // Set (alpha,alpha)->(alpha,alpha) double excitations
                    if ((alpha_internal_excitation.size() == 2) && (beta_internal_excitation.size() == 0))
                        element = sign_internal_excitation *
                                  wfn_->blas()
                                      ->get_MatTmp("Hijab[aa][aa]", i_unique, none)
                                      ->get_four_address_element(occ_to_act[alpha_internal_excitation[0].first],
                                                                 occ_to_act[alpha_internal_excitation[1].first],
                                                                 vir_to_act[alpha_internal_excitation[0].second],
                                                                 vir_to_act[alpha_internal_excitation[1].second]);

                    // Set (alpha,beta)->(alpha,beta) double excitations
                    if ((alpha_internal_excitation.size() == 1) && (beta_internal_excitation.size() == 1))
                        element = sign_internal_excitation *
                                  wfn_->blas()
                                      ->get_MatTmp("HiJaB[aA][aA]", i_unique, none)
                                      ->get_four_address_element(occ_to_act[alpha_internal_excitation[0].first],
                                                                 occ_to_act[beta_internal_excitation[0].first],
                                                                 vir_to_act[alpha_internal_excitation[0].second],
                                                                 vir_to_act[beta_internal_excitation[0].second]);

                    // Set (beta,beta)->(beta,beta) double excitations
                    if ((alpha_internal_excitation.size() == 0) && (beta_internal_excitation.size() == 2))
                        element = sign_internal_excitation *
                                  wfn_->blas()
                                      ->get_MatTmp("HIJAB[AA][AA]", i_unique, none)
                                      ->get_four_address_element(occ_to_act[beta_internal_excitation[0].first],
                                                                 occ_to_act[beta_internal_excitation[1].first],
                                                                 vir_to_act[beta_internal_excitation[0].second],
                                                                 vir_to_act[beta_internal_excitation[1].second]);
                } else {
                    // Set alpha-alpha single excitations
                    if ((alpha_internal_excitation.size() == 1) && (beta_internal_excitation.size() == 0))
                        element = sign_internal_excitation *
                                  wfn_->blas()
                                      ->get_MatTmp("HIA[A][A]", i_unique, none)
                                      ->get_two_address_element(occ_to_act[alpha_internal_excitation[0].first],
                                                                vir_to_act[alpha_internal_excitation[0].second]);

                    // Set beta-beta single excitations
                    if ((alpha_internal_excitation.size() == 0) && (beta_internal_excitation.size() == 1))
                        element = sign_internal_excitation *
                                  wfn_->blas()
                                      ->get_MatTmp("Hia[a][a]", i_unique, none)
                                      ->get_two_address_element(occ_to_act[beta_internal_excitation[0].first],
                                                                vir_to_act[beta_internal_excitation[0].second]);

                    // Set (alpha,alpha)->(alpha,alpha) double excitations
                    if ((alpha_internal_excitation.size() == 2) && (beta_internal_excitation.size() == 0))
                        element = sign_internal_excitation *
                                  wfn_->blas()
                                      ->get_MatTmp("HIJAB[AA][AA]", i_unique, none)
                                      ->get_four_address_element(occ_to_act[alpha_internal_excitation[0].first],
                                                                 occ_to_act[alpha_internal_excitation[1].first],
                                                                 vir_to_act[alpha_internal_excitation[0].second],
                                                                 vir_to_act[alpha_internal_excitation[1].second]);

                    // Set (alpha,beta)->(alpha,beta) double excitations
                    if ((alpha_internal_excitation.size() == 1) && (beta_internal_excitation.size() == 1))
                        element = sign_internal_excitation *
                                  wfn_->blas()
                                      ->get_MatTmp("HiJaB[aA][aA]", i_unique, none)
                                      ->get_four_address_element(occ_to_act[beta_internal_excitation[0].first],
                                                                 occ_to_act[alpha_internal_excitation[0].first],
                                                                 vir_to_act[beta_internal_excitation[0].second],
                                                                 vir_to_act[alpha_internal_excitation[0].second]);

                    // Set (beta,beta)->(beta,beta) double excitations
                    if ((alpha_internal_excitation.size() == 0) && (beta_internal_excitation.size() == 2))
                        element = sign_internal_excitation *
                                  wfn_->blas()
                                      ->get_MatTmp("Hijab[aa][aa]", i_unique, none)
                                      ->get_four_address_element(occ_to_act[beta_internal_excitation[0].first],
                                                                 occ_to_act[beta_internal_excitation[1].first],
                                                                 vir_to_act[beta_internal_excitation[0].second],
                                                                 vir_to_act[beta_internal_excitation[1].second]);
                }
                Heff_mrpt2[j][i] = element;
            }
        }
    }
}

}  // namespace psimrcc
}  // namespace psi
