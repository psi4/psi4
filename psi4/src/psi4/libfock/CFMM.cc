/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

#include "jk.h"
#include "SplitJK.h"

#include <algorithm>

using namespace psi;

namespace psi {

DirectCFMM::DirectCFMM(std::shared_ptr<BasisSet> primary, Options& options) : SplitJK(primary, options) {
    cfmmtree_ = std::make_shared<CFMMTree>(primary_, nullptr, options_);
}

void DirectCFMM::build_G_component(std::vector<std::shared_ptr<Matrix>>& D, std::vector<std::shared_ptr<Matrix>>& J,
    std::vector<std::shared_ptr<TwoBodyAOInt> >& eri_computers) {

    timer_on("CFMM: J");

    cfmmtree_->build_J(eri_computers, D, J);
    num_computed_shells_ = cfmmtree_->num_computed_shells();

    timer_off("CFMM: J");
}

void DirectCFMM::print_header() const {
    if (print_) {
        outfile->Printf("  ==> Continuous Fast Multipole Method (CFMM) <==\n\n");
        outfile->Printf("    Primary Basis: %11s\n", primary_->name().c_str());
        outfile->Printf("    Max Multipole Order: %11d\n", cfmmtree_->lmax());
        outfile->Printf("    Max Tree Depth: %11d\n", cfmmtree_->nlevels());
        outfile->Printf("\n");
    }
}

size_t DirectCFMM::num_computed_shells() {
    return num_computed_shells_;
}

DFCFMM::DFCFMM(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, 
               Options& options) : DirectDFJ(primary, auxiliary, options) {

    df_cfmm_tree_ = std::make_shared<CFMMTree>(primary_, auxiliary_, options_);

    const int nshell_aux = auxiliary_->nshell();
    J_metric_shell_diag_.resize(nshell_aux, 0.0);
    for (int shell = 0; shell < nshell_aux; shell++) {
        const int bf_start = auxiliary_->shell(shell).function_index();
        const int bf_end = bf_start + auxiliary_->shell(shell).nfunction();
        for (int bf = bf_start; bf < bf_end; bf++) {
            J_metric_shell_diag_[shell] =
                std::max(J_metric_shell_diag_[shell], J_metric_->get(bf, bf));
        }
    }
}

void DFCFMM::build_G_component(std::vector<std::shared_ptr<Matrix>>& D, std::vector<std::shared_ptr<Matrix>>& J,
    std::vector<std::shared_ptr<TwoBodyAOInt> >& eri_computers) {

    timer_on("DFCFMM: J");

    int naux = auxiliary_->nbf();

    if (gamma_.size() != D.size()) {
        gamma_.clear();
        gamma_.resize(D.size());
        for (int i = 0; i < D.size(); i++) {
            gamma_[i] = std::make_shared<Matrix>(naux, 1);
        }
    }

    // Build gammaP = (P|uv)Duv
    df_cfmm_tree_->df_set_contraction(ContractionType::DF_AUX_PRI);
    df_cfmm_tree_->build_J(eri_computers, D, gamma_, J_metric_shell_diag_);
    const size_t computed_triplets_first_contraction = df_cfmm_tree_->num_computed_shells();

    // Solve for gammaQ => (P|Q)*gammaQ = gammaP
    for (int i = 0; i < D.size(); i++) {
        SharedMatrix Jmet_copy = J_metric_->clone();
        std::vector<int> ipiv(naux);

        C_DGESV(naux, 1, Jmet_copy->pointer()[0], naux, ipiv.data(), gamma_[i]->pointer()[0], naux);
    }

    // Build Juv = (uv|Q) * gammaQ
    df_cfmm_tree_->df_set_contraction(ContractionType::DF_PRI_AUX);
    df_cfmm_tree_->build_J(eri_computers, gamma_, J, J_metric_shell_diag_);
    const size_t computed_triplets_second_contraction = df_cfmm_tree_->num_computed_shells();

    num_computed_shells_ = computed_triplets_first_contraction + computed_triplets_second_contraction;

    timer_off("DFCFMM: J");
}

void DFCFMM::print_header() const {
    if (print_) {
        outfile->Printf("  ==> CFMM-Accelerated Direct Density Fitted J <==\n\n");
        outfile->Printf("    Primary Basis: %11s\n", primary_->name().c_str());
        outfile->Printf("    Auxiliary Basis: %11s\n", auxiliary_->name().c_str());
        outfile->Printf("    Max Multipole Order: %11d\n", df_cfmm_tree_->lmax());
        outfile->Printf("    Max Tree Depth: %11d\n", df_cfmm_tree_->nlevels());
        outfile->Printf("\n");
    }
}

size_t DFCFMM::num_computed_shells() {
    return num_computed_shells_;
}

}  // namespace psi
