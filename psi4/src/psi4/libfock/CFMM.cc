#include "jk.h"
#include "SplitJK.h"

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
}

void DFCFMM::build_G_component(std::vector<std::shared_ptr<Matrix>>& D, std::vector<std::shared_ptr<Matrix>>& J,
    std::vector<std::shared_ptr<TwoBodyAOInt> >& eri_computers) {

    timer_on("DFCFMM: J");

    int naux = auxiliary_->nbf();
    int nshell_aux = auxiliary_->nshell();

    if (gamma.size() == 0) {
        gamma.resize(D.size());
        for (int i = 0; i < D.size(); i++) {
            gamma[i] = std::make_shared<Matrix>(naux, 1);
        }
    }

    // diagonal shell maxima of J_metric_ for screening
    std::vector<double> J_metric_shell_diag(nshell_aux, 0.0);
    for (size_t s = 0; s < nshell_aux; s++) {
        int bf_start = auxiliary_->shell(s).function_index();
        int bf_end = bf_start + auxiliary_->shell(s).nfunction();
        for (size_t bf = bf_start; bf < bf_end; bf++) {
            J_metric_shell_diag[s] = std::max(J_metric_shell_diag[s], J_metric_->get(bf, bf));
        }
    }

    // Build gammaP = (P|uv)Duv
    df_cfmm_tree_->df_set_contraction(ContractionType::DF_AUX_PRI);
    df_cfmm_tree_->build_J(eri_computers, D, gamma, J_metric_shell_diag);
    const size_t computed_triplets_first_contraction = df_cfmm_tree_->num_computed_shells();

    // Solve for gammaQ => (P|Q)*gammaQ = gammaP
    for (int i = 0; i < D.size(); i++) {
        SharedMatrix Jmet_copy = J_metric_->clone();
        std::vector<int> ipiv(naux);

        C_DGESV(naux, 1, Jmet_copy->pointer()[0], naux, ipiv.data(), gamma[i]->pointer()[0], naux);
    }

    // Build Juv = (uv|Q) * gammaQ
    df_cfmm_tree_->df_set_contraction(ContractionType::DF_PRI_AUX);
    df_cfmm_tree_->build_J(eri_computers, gamma, J, J_metric_shell_diag);
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

}
