/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#include "mrccsd_t.h"

namespace psi {
namespace psimrcc {

void MRCCSD_T::compute_ooo_contribution_to_Heff_restricted(int i, int j, int k, int mu, BlockMatrix* T3) {
    // Find the off_diagonal elements for reference mu
    // Loop over reference nu (in a safe way)
    for (int nu = 0; nu < nrefs; nu++) {
        if (nu != mu) {
            std::vector<std::pair<int, int> > alpha_internal_excitation =
                wfn_->moinfo()->get_alpha_internal_excitation(mu, nu);
            std::vector<std::pair<int, int> > beta_internal_excitation =
                wfn_->moinfo()->get_beta_internal_excitation(mu, nu);
            double sign_internal_excitation = wfn_->moinfo()->get_sign_internal_excitation(mu, nu);

            // Set (alpha)->(alpha) single excitations
            if ((alpha_internal_excitation.size() == 1) && (beta_internal_excitation.size() == 0)) {
                d_h_eff[nu][mu] += sign_internal_excitation * compute_A_ooo_contribution_to_Heff_restricted(
                                                                  alpha_internal_excitation[0].first,
                                                                  alpha_internal_excitation[0].second, i, j, k, mu, T3);
            }
        }
    }
}

void MRCCSD_T::compute_ooO_contribution_to_Heff_restricted(int i, int j, int k, int mu, BlockMatrix* T3) {
    // Find the off_diagonal elements for reference mu
    // Loop over reference nu (in a safe way)
    for (int nu = 0; nu < nrefs; nu++) {
        if (nu != mu) {
            std::vector<std::pair<int, int> > alpha_internal_excitation =
                wfn_->moinfo()->get_alpha_internal_excitation(mu, nu);
            std::vector<std::pair<int, int> > beta_internal_excitation =
                wfn_->moinfo()->get_beta_internal_excitation(mu, nu);
            double sign_internal_excitation = wfn_->moinfo()->get_sign_internal_excitation(mu, nu);

            // Set (alpha)->(alpha) single excitations
            if ((alpha_internal_excitation.size() == 1) && (beta_internal_excitation.size() == 0)) {
                d_h_eff[nu][mu] += sign_internal_excitation * compute_A_ooO_contribution_to_Heff_restricted(
                                                                  alpha_internal_excitation[0].first,
                                                                  alpha_internal_excitation[0].second, i, j, k, mu, T3);
            }
            // Set (beta)->(beta) single excitations
            if ((alpha_internal_excitation.size() == 0) && (beta_internal_excitation.size() == 1)) {
                d_h_eff[nu][mu] += sign_internal_excitation * compute_B_ooO_contribution_to_Heff_restricted(
                                                                  beta_internal_excitation[0].first,
                                                                  beta_internal_excitation[0].second, i, j, k, mu, T3);
            }
            // Set (alpha,beta)->(alpha,beta) double excitations
            if ((alpha_internal_excitation.size() == 1) && (beta_internal_excitation.size() == 1)) {
                d_h_eff[nu][mu] += sign_internal_excitation * compute_AB_ooO_contribution_to_Heff_restricted(
                                                                  alpha_internal_excitation[0].first,
                                                                  beta_internal_excitation[0].first,
                                                                  alpha_internal_excitation[0].second,
                                                                  beta_internal_excitation[0].second, i, j, k, mu, T3);
            }
        }
    }
}

void MRCCSD_T::compute_oOO_contribution_to_Heff_restricted(int i, int j, int k, int mu, BlockMatrix* T3) {
    // Find the off_diagonal elements for reference mu
    // Loop over reference nu (in a safe way)
    for (int nu = 0; nu < nrefs; nu++) {
        if (nu != mu) {
            std::vector<std::pair<int, int> > alpha_internal_excitation =
                wfn_->moinfo()->get_alpha_internal_excitation(mu, nu);
            std::vector<std::pair<int, int> > beta_internal_excitation =
                wfn_->moinfo()->get_beta_internal_excitation(mu, nu);
            double sign_internal_excitation = wfn_->moinfo()->get_sign_internal_excitation(mu, nu);

            // Set (alpha)->(alpha) single excitations
            if ((alpha_internal_excitation.size() == 1) && (beta_internal_excitation.size() == 0)) {
                d_h_eff[nu][mu] += sign_internal_excitation * compute_A_oOO_contribution_to_Heff_restricted(
                                                                  alpha_internal_excitation[0].first,
                                                                  alpha_internal_excitation[0].second, i, j, k, mu, T3);
            }
            // Set (beta)->(beta) single excitations
            if ((alpha_internal_excitation.size() == 0) && (beta_internal_excitation.size() == 1)) {
                d_h_eff[nu][mu] += sign_internal_excitation * compute_B_oOO_contribution_to_Heff_restricted(
                                                                  beta_internal_excitation[0].first,
                                                                  beta_internal_excitation[0].second, i, j, k, mu, T3);
            }
            // Set (alpha,beta)->(alpha,beta) double excitations
            if ((alpha_internal_excitation.size() == 1) && (beta_internal_excitation.size() == 1)) {
                d_h_eff[nu][mu] += sign_internal_excitation * compute_AB_oOO_contribution_to_Heff_restricted(
                                                                  alpha_internal_excitation[0].first,
                                                                  beta_internal_excitation[0].first,
                                                                  alpha_internal_excitation[0].second,
                                                                  beta_internal_excitation[0].second, i, j, k, mu, T3);
            }
        }
    }
}

void MRCCSD_T::compute_OOO_contribution_to_Heff_restricted(int i, int j, int k, int mu, BlockMatrix* T3) {
    // Find the off_diagonal elements for reference mu
    // Loop over reference nu (in a safe way)
    for (int nu = 0; nu < nrefs; nu++) {
        if (nu != mu) {
            std::vector<std::pair<int, int> > alpha_internal_excitation =
                wfn_->moinfo()->get_alpha_internal_excitation(mu, nu);
            std::vector<std::pair<int, int> > beta_internal_excitation =
                wfn_->moinfo()->get_beta_internal_excitation(mu, nu);
            double sign_internal_excitation = wfn_->moinfo()->get_sign_internal_excitation(mu, nu);

            // Set (beta)->(beta) single excitations
            if ((alpha_internal_excitation.size() == 0) && (beta_internal_excitation.size() == 1)) {
                d_h_eff[nu][mu] += sign_internal_excitation * compute_B_OOO_contribution_to_Heff_restricted(
                                                                  beta_internal_excitation[0].first,
                                                                  beta_internal_excitation[0].second, i, j, k, mu, T3);
            }
        }
    }
}

}  // namespace psimrcc
}  // namespace psi
