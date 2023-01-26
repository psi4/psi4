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

/**
 *  @file updater.cc
 *  @ingroup (PSIMRCC)
 *  @brief Contains methods for updating the CC equations
 */

#include <vector>
//#include <string>
//
#include <cstdio>
#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/liboptions/liboptions.h"

#include "blas.h"
#include "matrix.h"
#include "matrixtmp.h"
#include "updater.h"

namespace psi {

namespace psimrcc {

Updater::Updater(std::shared_ptr<PSIMRCCWfn> wfn, Options &options) : options_(options), wfn_(wfn) {}

Updater::~Updater() {}

void Updater::zero_internal_amps() {
    if (options_.get_bool("ZERO_INTERNAL_AMPS")) {
        // Zero internal amplitudes for unique reference i
        for (int i = 0; i < wfn_->moinfo()->get_nunique(); i++) {
            int unique_i = wfn_->moinfo()->get_ref_number(i, UniqueRefs);
            // Loop over reference j
            for (int j = 0; j < wfn_->moinfo()->get_ref_size(AllRefs); j++) {
                std::vector<std::pair<int, int> > alpha_internal_excitation =
                    wfn_->moinfo()->get_alpha_internal_excitation(unique_i, j);
                std::vector<std::pair<int, int> > beta_internal_excitation =
                    wfn_->moinfo()->get_beta_internal_excitation(unique_i, j);

                // Zero alpha-alpha single excitations
                if ((alpha_internal_excitation.size() == 1) && (beta_internal_excitation.size() == 0)) {
                    wfn_->blas()
                        ->get_MatTmp("t1[o][v]", unique_i, none)
                        ->set_two_address_element(alpha_internal_excitation[0].first,
                                                  alpha_internal_excitation[0].second, 0.0);
                }

                // Zero beta-beta single excitations
                if ((alpha_internal_excitation.size() == 0) && (beta_internal_excitation.size() == 1))
                    wfn_->blas()
                        ->get_MatTmp("t1[O][V]", unique_i, none)
                        ->set_two_address_element(beta_internal_excitation[0].first, beta_internal_excitation[0].second,
                                                  0.0);

                // Zero (alpha,alpha)->(alpha,alpha) double excitations (all permutations)
                if ((alpha_internal_excitation.size() == 2) && (beta_internal_excitation.size() == 0)) {
                    wfn_->blas()
                        ->get_MatTmp("t2[oo][vv]", unique_i, none)
                        ->set_four_address_element(
                            alpha_internal_excitation[0].first, alpha_internal_excitation[1].first,
                            alpha_internal_excitation[0].second, alpha_internal_excitation[1].second, 0.0);
                    wfn_->blas()
                        ->get_MatTmp("t2[oo][vv]", unique_i, none)
                        ->set_four_address_element(
                            alpha_internal_excitation[0].first, alpha_internal_excitation[1].first,
                            alpha_internal_excitation[1].second, alpha_internal_excitation[0].second, 0.0);
                    wfn_->blas()
                        ->get_MatTmp("t2[oo][vv]", unique_i, none)
                        ->set_four_address_element(
                            alpha_internal_excitation[1].first, alpha_internal_excitation[0].first,
                            alpha_internal_excitation[0].second, alpha_internal_excitation[1].second, 0.0);
                    wfn_->blas()
                        ->get_MatTmp("t2[oo][vv]", unique_i, none)
                        ->set_four_address_element(
                            alpha_internal_excitation[1].first, alpha_internal_excitation[0].first,
                            alpha_internal_excitation[1].second, alpha_internal_excitation[0].second, 0.0);
                }

                // Zero (alpha,beta)->(alpha,beta) double excitations
                if ((alpha_internal_excitation.size() == 1) && (beta_internal_excitation.size() == 1)) {
                    wfn_->blas()
                        ->get_MatTmp("t2[oO][vV]", unique_i, none)
                        ->set_four_address_element(
                            alpha_internal_excitation[0].first, beta_internal_excitation[0].first,
                            alpha_internal_excitation[0].second, beta_internal_excitation[0].second, 0.0);
                }

                // Zero (beta,beta)->(beta,beta) double excitations (all permutations)
                if ((alpha_internal_excitation.size() == 0) && (beta_internal_excitation.size() == 2)) {
                    wfn_->blas()
                        ->get_MatTmp("t2[OO][VV]", unique_i, none)
                        ->set_four_address_element(beta_internal_excitation[0].first, beta_internal_excitation[1].first,
                                                   beta_internal_excitation[0].second,
                                                   beta_internal_excitation[1].second, 0.0);
                    wfn_->blas()
                        ->get_MatTmp("t2[OO][VV]", unique_i, none)
                        ->set_four_address_element(beta_internal_excitation[0].first, beta_internal_excitation[1].first,
                                                   beta_internal_excitation[1].second,
                                                   beta_internal_excitation[0].second, 0.0);
                    wfn_->blas()
                        ->get_MatTmp("t2[OO][VV]", unique_i, none)
                        ->set_four_address_element(beta_internal_excitation[1].first, beta_internal_excitation[0].first,
                                                   beta_internal_excitation[0].second,
                                                   beta_internal_excitation[1].second, 0.0);
                    wfn_->blas()
                        ->get_MatTmp("t2[OO][VV]", unique_i, none)
                        ->set_four_address_element(beta_internal_excitation[1].first, beta_internal_excitation[0].first,
                                                   beta_internal_excitation[1].second,
                                                   beta_internal_excitation[0].second, 0.0);
                }
            }
        }
    } else {
        outfile->Printf(
            "\n  Warning: the internal amplitudes are not zeroed.\n  This is not proper Mk-MRCC. Size-extensivity "
            "might be lost\n");
    }
}

void Updater::zero_t1_internal_amps() {
    if (options_.get_bool("ZERO_INTERNAL_AMPS")) {
        // Zero internal amplitudes for unique reference i
        for (int i = 0; i < wfn_->moinfo()->get_nunique(); i++) {
            int unique_i = wfn_->moinfo()->get_ref_number(i, UniqueRefs);
            // Loop over reference j
            for (int j = 0; j < wfn_->moinfo()->get_ref_size(AllRefs); j++) {
                std::vector<std::pair<int, int> > alpha_internal_excitation =
                    wfn_->moinfo()->get_alpha_internal_excitation(unique_i, j);
                std::vector<std::pair<int, int> > beta_internal_excitation =
                    wfn_->moinfo()->get_beta_internal_excitation(unique_i, j);

                // Zero alpha-alpha single excitations
                if ((alpha_internal_excitation.size() == 1) && (beta_internal_excitation.size() == 0))
                    wfn_->blas()
                        ->get_MatTmp("t1[o][v]", unique_i, none)
                        ->set_two_address_element(alpha_internal_excitation[0].first,
                                                  alpha_internal_excitation[0].second, 0.0);

                // Zero beta-beta single excitations
                if ((alpha_internal_excitation.size() == 0) && (beta_internal_excitation.size() == 1))
                    wfn_->blas()
                        ->get_MatTmp("t1[O][V]", unique_i, none)
                        ->set_two_address_element(beta_internal_excitation[0].first, beta_internal_excitation[0].second,
                                                  0.0);
            }
        }
    } else {
        outfile->Printf(
            "\n  Warning: the internal amplitudes are not zeroed.\n  This is not proper Mk-MRCC. Size-extensivity "
            "might be lost\n");
    }
}

void Updater::zero_internal_delta_amps() {
    if (options_.get_bool("ZERO_INTERNAL_AMPS")) {
        // Zero internal amplitudes for unique reference i
        for (int i = 0; i < wfn_->moinfo()->get_nunique(); i++) {
            int unique_i = wfn_->moinfo()->get_ref_number(i, UniqueRefs);
            // Loop over reference j
            for (int j = 0; j < wfn_->moinfo()->get_ref_size(AllRefs); j++) {
                std::vector<std::pair<int, int> > alpha_internal_excitation =
                    wfn_->moinfo()->get_alpha_internal_excitation(unique_i, j);
                std::vector<std::pair<int, int> > beta_internal_excitation =
                    wfn_->moinfo()->get_beta_internal_excitation(unique_i, j);

                // Zero alpha-alpha single excitations
                if ((alpha_internal_excitation.size() == 1) && (beta_internal_excitation.size() == 0))
                    wfn_->blas()
                        ->get_MatTmp("t1_delta[o][v]", unique_i, none)
                        ->set_two_address_element(alpha_internal_excitation[0].first,
                                                  alpha_internal_excitation[0].second, 0.0);

                // Zero beta-beta single excitations
                if ((alpha_internal_excitation.size() == 0) && (beta_internal_excitation.size() == 1))
                    wfn_->blas()
                        ->get_MatTmp("t1_delta[O][V]", unique_i, none)
                        ->set_two_address_element(beta_internal_excitation[0].first, beta_internal_excitation[0].second,
                                                  0.0);

                // Zero (alpha,alpha)->(alpha,alpha) double excitations (all permutations)
                if ((alpha_internal_excitation.size() == 2) && (beta_internal_excitation.size() == 0)) {
                    wfn_->blas()
                        ->get_MatTmp("t2_delta[oo][vv]", unique_i, none)
                        ->set_four_address_element(
                            alpha_internal_excitation[0].first, alpha_internal_excitation[1].first,
                            alpha_internal_excitation[0].second, alpha_internal_excitation[1].second, 0.0);
                    wfn_->blas()
                        ->get_MatTmp("t2_delta[oo][vv]", unique_i, none)
                        ->set_four_address_element(
                            alpha_internal_excitation[0].first, alpha_internal_excitation[1].first,
                            alpha_internal_excitation[1].second, alpha_internal_excitation[0].second, 0.0);
                    wfn_->blas()
                        ->get_MatTmp("t2_delta[oo][vv]", unique_i, none)
                        ->set_four_address_element(
                            alpha_internal_excitation[1].first, alpha_internal_excitation[0].first,
                            alpha_internal_excitation[0].second, alpha_internal_excitation[1].second, 0.0);
                    wfn_->blas()
                        ->get_MatTmp("t2_delta[oo][vv]", unique_i, none)
                        ->set_four_address_element(
                            alpha_internal_excitation[1].first, alpha_internal_excitation[0].first,
                            alpha_internal_excitation[1].second, alpha_internal_excitation[0].second, 0.0);
                }

                // Zero (alpha,beta)->(alpha,beta) double excitations
                if ((alpha_internal_excitation.size() == 1) && (beta_internal_excitation.size() == 1)) {
                    wfn_->blas()
                        ->get_MatTmp("t2_delta[oO][vV]", unique_i, none)
                        ->set_four_address_element(
                            alpha_internal_excitation[0].first, beta_internal_excitation[0].first,
                            alpha_internal_excitation[0].second, beta_internal_excitation[0].second, 0.0);
                }

                // Zero (beta,beta)->(beta,beta) double excitations (all permutations)
                if ((alpha_internal_excitation.size() == 0) && (beta_internal_excitation.size() == 2)) {
                    wfn_->blas()
                        ->get_MatTmp("t2_delta[OO][VV]", unique_i, none)
                        ->set_four_address_element(beta_internal_excitation[0].first, beta_internal_excitation[1].first,
                                                   beta_internal_excitation[0].second,
                                                   beta_internal_excitation[1].second, 0.0);
                    wfn_->blas()
                        ->get_MatTmp("t2_delta[OO][VV]", unique_i, none)
                        ->set_four_address_element(beta_internal_excitation[0].first, beta_internal_excitation[1].first,
                                                   beta_internal_excitation[1].second,
                                                   beta_internal_excitation[0].second, 0.0);
                    wfn_->blas()
                        ->get_MatTmp("t2_delta[OO][VV]", unique_i, none)
                        ->set_four_address_element(beta_internal_excitation[1].first, beta_internal_excitation[0].first,
                                                   beta_internal_excitation[0].second,
                                                   beta_internal_excitation[1].second, 0.0);
                    wfn_->blas()
                        ->get_MatTmp("t2_delta[OO][VV]", unique_i, none)
                        ->set_four_address_element(beta_internal_excitation[1].first, beta_internal_excitation[0].first,
                                                   beta_internal_excitation[1].second,
                                                   beta_internal_excitation[0].second, 0.0);
                }
            }
        }
    } else {
        outfile->Printf(
            "\n  Warning: the internal amplitudes are not zeroed.\n  This is not proper Mk-MRCC. Size-extensivity "
            "might be lost\n");
    }
}

}  // namespace psimrcc
}  // namespace psi
