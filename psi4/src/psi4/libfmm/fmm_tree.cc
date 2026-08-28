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

#include "psi4/libfmm/multipoles_helper.h"
#include "psi4/libfmm/fmm_tree.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector3.h"
#include "psi4/libmints/gshell.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/onebody.h"
#include "psi4/libmints/multipoles.h"
#include "psi4/libmints/overlap.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <memory>
#include <unordered_map>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi {

// The tree traversal and the separation of direct near-field interactions from
// multipole far-field interactions follow the continuous fast multipole method
// (CFMM): White et al., Chem. Phys. Lett. 230, 8 (1994),
// doi:10.1016/0009-2614(94)01128-1, and White et al., Chem. Phys. Lett. 253,
// 268 (1996), doi:10.1016/0009-2614(96)00175-3.

ShellPair::ShellPair(std::shared_ptr<BasisSet>& bs1, std::shared_ptr<BasisSet>& bs2, std::pair<int, int> pair_index,
                     std::shared_ptr<HarmonicCoefficients>& mpole_coefs, double cfmm_extent_tol) {
    bs1_ = bs1;
    bs2_ = bs2;
    pair_index_ = pair_index;

    const GaussianShell& P_shell = bs1_->shell(pair_index.first);
    const GaussianShell& Q_shell = bs2_->shell(pair_index.second);

    Vector3 P_center = P_shell.center();
    Vector3 Q_center = Q_shell.center();

    center_ = Vector3(0.0, 0.0, 0.0);
    exp_ = INFINITY;

    int nprim_p = P_shell.nprimitive();
    int nprim_q = Q_shell.nprimitive();
    for (int pp = 0; pp < nprim_p; pp++) {
        double P_exponent = P_shell.exp(pp);
        for (int qp = 0; qp < nprim_q; qp++) {
            double Q_exponent = Q_shell.exp(qp);

            const double product_exponent = std::abs(P_exponent + Q_exponent);
            Vector3 product_center =
                (P_exponent * P_center + Q_exponent * Q_center) / product_exponent;

            center_ += product_center;
            exp_ = std::min(exp_, product_exponent);
        }
    }
    center_ /= (nprim_p * nprim_q);
    // The Gaussian-product center is obtained from the Gaussian product
    // theorem. The radial extent is the radius at which the diffuse envelope
    // exp(-exp_ * r^2 / 2) reaches cfmm_extent_tol; this extent controls the
    // CFMM near-/far-field partition described by White et al. (1994).
    extent_ = std::sqrt(-2.0 * std::log(cfmm_extent_tol) / exp_);

    mpole_coefs_ = mpole_coefs;

}

void ShellPair::calculate_mpoles(Vector3 box_center, std::shared_ptr<OneBodyAOInt>& s_ints,
                                 std::shared_ptr<OneBodyAOInt>& mpole_ints, int lmax) {

    int P = pair_index_.first;
    int Q = pair_index_.second;

    // Calculate the overlap integrals (Order 0 multipole integrals)
    s_ints->compute_shell(P, Q);
    const double* sbuffer = s_ints->buffers()[0];

    // Calculate the multipole integrals
    mpole_ints->compute_shell(P, Q);
    const double* mbuffer = mpole_ints->buffers()[0];

    const GaussianShell& P_shell = bs1_->shell(P);
    const GaussianShell& Q_shell = bs2_->shell(Q);

    int p_start = P_shell.start();
    int num_p = P_shell.nfunction();

    int q_start = Q_shell.start();
    int num_q = Q_shell.nfunction();

    for (int p = p_start; p < p_start + num_p; p++) {
        int dp = p - p_start;
        for (int q = q_start; q < q_start + num_q; q++) {
            int dq = q - q_start;

            auto pair_multipoles = std::make_shared<RealSolidHarmonics>(lmax, box_center, SolidHarmonicsType::Regular);

            pair_multipoles->add(0, 0, sbuffer[dp * num_q + dq]);

            // Convert Cartesian AO multipole integrals to real regular solid
            // harmonics. The moment definition is Eq. (49) of Reine et al.,
            // WIREs Comput. Mol. Sci. 2, 290 (2012), doi:10.1002/wcms.78;
            // normalization and phase conventions follow Helgaker, Jorgensen,
            // and Olsen, Molecular Electronic-Structure Theory (Wiley, 2000),
            // doi:10.1002/9781119019572, pp. 412-414.
            int running_index = 0;
            for (int l = 1; l <= lmax; l++) {
                int l_ncart = ncart(l);
                for (int m = -l; m <= l; m++) {
                    int mu = m_addr(m);
                    const auto& mpole_terms = mpole_coefs_->get_terms(l, mu);

                    int cartesian_component = 0;
                    for (int ii = 0; ii <= l; ii++) {
                        int a = l - ii;
                        for (int jj = 0; jj <= ii; jj++) {
                            int b = ii - jj;
                            int c = jj;
                            int cartesian_index = a * l_ncart * l_ncart + b * l_ncart + c;

                            const auto term = mpole_terms.find(cartesian_index);
                            if (term != mpole_terms.end()) {
                                double coefficient = term->second;
                                int integral_index = cartesian_component + running_index;
                                pair_multipoles->add(
                                    l, mu,
                                    pow(-1.0, (double)l + 1) * coefficient *
                                        mbuffer[integral_index * num_p * num_q + dp * num_q + dq]);
                            }
                            cartesian_component += 1;
                        } // end jj
                    } // end ii
                } // end m loop
                running_index += l_ncart;
            } // end l
            mpoles_.push_back(pair_multipoles);
        } // end q
    } // end p

}

CFMMBox::CFMMBox(std::shared_ptr<CFMMBox> parent, std::vector<std::shared_ptr<ShellPair>> primary_shell_pairs,
                std::vector<std::shared_ptr<ShellPair>> auxiliary_shell_pairs, Vector3 origin, double length,
                int level, int lmax, int ws) {

    parent_ = parent;
    primary_shell_pairs_ = primary_shell_pairs;
    auxiliary_shell_pairs_ = auxiliary_shell_pairs;
    origin_ = origin;
    center_ = origin_ + 0.5 * Vector3(length, length, length);
    length_ = length;
    level_ = level;
    lmax_ = lmax;
    ws_ = ws;
    ws_max_ = ws;

}

void CFMMBox::make_children() {

    int nchild = (level_ > 0) ? 16 : 8;
    std::vector<std::vector<std::shared_ptr<ShellPair>>> child_shell_pair_primary_buffer(nchild);
    std::vector<std::vector<std::shared_ptr<ShellPair>>> child_shell_pair_auxiliary_buffer(nchild);

    // Max WS at the child's level
    int child_level_max_ws = std::max(2, (int) std::pow(2, level_+1));
    int diffuse_child_max_ws = child_level_max_ws;

    // Fill order (ws,z,y,x) (0)000 (0)001 (0)010 (0)011 (0)100 (0)101 (0)110 (0)111
    // (1)000 (1)001 (1)010 (1)011 (1)100 (1)101 (1)110 (1)111
    for (std::shared_ptr<ShellPair> shell_pair : primary_shell_pairs_) {
        Vector3 sp_center = shell_pair->get_center();
        double x = sp_center[0];
        double y = sp_center[1];
        double z = sp_center[2];
        double extent = shell_pair->get_extent();
        int ws = std::max(2, 2 * (int)std::ceil(extent / length_));

        int xbit = (x < center_[0]) ? 0 : 1;
        int ybit = (y < center_[1]) ? 0 : 1;
        int zbit = (z < center_[2]) ? 0 : 1;
        int rbit = (level_ == 0 || ws < 2 * ws_) ? 0 : 1;

        int boxind = 8 * rbit + 4 * zbit + 2 * ybit + 1 * xbit;
        child_shell_pair_primary_buffer[boxind].push_back(shell_pair);

        int child_ws = std::max(2, (int)std::ceil(2.0 * extent / length_));
        if (child_ws > diffuse_child_max_ws) diffuse_child_max_ws = child_ws;
    }

    // auxiliary_shell_pairs_ is empty for a direct primary-basis contraction.
    for (std::shared_ptr<ShellPair> shell_pair : auxiliary_shell_pairs_) {
        Vector3 sp_center = shell_pair->get_center();
        double x = sp_center[0];
        double y = sp_center[1];
        double z = sp_center[2];
        double extent = shell_pair->get_extent();
        int ws = std::max(2, 2 * (int)std::ceil(extent / length_));

        int xbit = (x < center_[0]) ? 0 : 1;
        int ybit = (y < center_[1]) ? 0 : 1;
        int zbit = (z < center_[2]) ? 0 : 1;
        int rbit = (level_ == 0 || ws < 2 * ws_) ? 0 : 1;

        int boxind = 8 * rbit + 4 * zbit + 2 * ybit + 1 * xbit;
        child_shell_pair_auxiliary_buffer[boxind].push_back(shell_pair);

        int child_ws = std::max(2, (int)std::ceil(2.0 * extent / length_));
        if (child_ws > diffuse_child_max_ws) diffuse_child_max_ws = child_ws;
    }

    // Make the children
    for (int boxind = 0; boxind < nchild; boxind++) {
        int xbit = boxind % 2;
        int ybit = (boxind / 2) % 2;
        int zbit = (boxind / 4) % 2;
        int rbit = (boxind / 8) % 2;
        Vector3 new_origin = origin_ + Vector3(xbit * 0.5 * length_, ybit * 0.5 * length_, zbit * 0.5 * length_);
        int child_ws = 2 * ws_ - 2 + 2 * rbit;
        children_.push_back(std::make_shared<CFMMBox>(this->get(), child_shell_pair_primary_buffer[boxind], child_shell_pair_auxiliary_buffer[boxind],
                                                      new_origin, 0.5 * length_, level_ + 1, lmax_, child_ws));
        if (child_ws == child_level_max_ws) children_[boxind]->set_ws_max(diffuse_child_max_ws);
    }
}

void CFMMBox::set_regions() {

    // Creates a temporary parent shared pointer
    std::shared_ptr<CFMMBox> parent = parent_.lock();

    // Parent is not a nullpointer
    if (parent) {
        constexpr auto axis_sign = [](double coordinate) constexpr noexcept {
            return static_cast<double>((coordinate >= 1.0e-8) - (coordinate <= -1.0e-8));
        };

        // Near field or local far fields are from children of parents
        // and children of parent's near field
        for (std::shared_ptr<CFMMBox> parent_nf : parent->near_field_) {
            for (std::shared_ptr<CFMMBox> child : parent_nf->children_) {
                if (child->nshell_pair() == 0) continue;
                // WS Max formulation takes the most diffuse branch into account
                int ref_ws = (ws_max_ + child->ws_max_) / 2;

                Vector3 Rab = child->center_ - center_;
                double dx = axis_sign(Rab[0]);
                double dy = axis_sign(Rab[1]);
                double dz = axis_sign(Rab[2]);

                Rab = Rab - length_ * Vector3(dx, dy, dz);
                double rab = std::sqrt(Rab.dot(Rab));

                if (rab <= length_ * ref_ws) {
                    near_field_.push_back(child);
                } else {
                    local_far_field_.push_back(child);
                }
            }
        }
    } else {
        near_field_.push_back(this->get());
    }
}

void CFMMBox::compute_far_field_contribution(std::shared_ptr<CFMMBox> lff_box) {
    for (int N = 0; N < Vff_.size(); N++) {
        std::shared_ptr<RealSolidHarmonics> far_field = lff_box->mpoles_[N]->far_field_vector(center_);
        Vff_[N]->add(far_field);
    }
}

void CFMMBox::add_parent_far_field_contribution() {
    // Temporary parent shared pointer object
    std::shared_ptr<CFMMBox> parent = parent_.lock();

    if (parent) {
        // Downward pass: translate the parent's irregular local expansion to
        // the child center, then accumulate it in the child's local expansion.
        for (int N = 0; N < Vff_.size(); N++) {
            auto parent_contribution = parent->Vff_[N]->translate(center_);
            Vff_[N]->add(parent_contribution);
        }
    }
}

void CFMMBox::compute_multipoles(const std::vector<SharedMatrix>& D, ContractionType contraction_type) {

    if (mpoles_.size() == 0) {
        mpoles_.resize(D.size());
        Vff_.resize(D.size());
    }

    // Create multipoles and far field vectors for each density matrix
    for (int N = 0; N < D.size(); N++) {
        mpoles_[N] = std::make_shared<RealSolidHarmonics>(lmax_, center_, SolidHarmonicsType::Regular);
        Vff_[N] = std::make_shared<RealSolidHarmonics>(lmax_, center_, SolidHarmonicsType::Irregular);
    }

    bool is_primary = (contraction_type == ContractionType::DF_AUX_PRI || contraction_type == ContractionType::DIRECT);
    std::vector<std::shared_ptr<ShellPair>>& ref_shell_pairs = (is_primary) ? primary_shell_pairs_ : auxiliary_shell_pairs_;
    if (ref_shell_pairs.empty()) return;

    std::shared_ptr<BasisSet> bs1 = ref_shell_pairs[0]->bs1();
    std::shared_ptr<BasisSet> bs2 = ref_shell_pairs[0]->bs2();
    int nbf = (is_primary) ? bs1->nbf() : 1;

    // Leaf-level upward pass: contract shell-pair moments with the density.
    for (const auto& sp : ref_shell_pairs) {

        std::vector<std::shared_ptr<RealSolidHarmonics>>& sp_mpoles = sp->get_mpoles();

        std::pair<int, int> PQ = sp->get_shell_pair_index();
        int P = PQ.first;
        int Q = PQ.second;

        double prefactor = (P == Q || !is_primary) ? 1.0 : 2.0;

        const GaussianShell& Pshell = bs1->shell(P);
        const GaussianShell& Qshell = bs2->shell(Q);

        int p_start = Pshell.start();
        int num_p = Pshell.nfunction();

        int q_start = Qshell.start();
        int num_q = Qshell.nfunction();

        for (int N = 0; N < D.size(); N++) {
            double* Dp = D[N]->pointer()[0];
            for (int p = p_start; p < p_start + num_p; p++) {
                int dp = p - p_start;
                for (int q = q_start; q < q_start + num_q; q++) {
                    int dq = q - q_start;

                    std::shared_ptr<RealSolidHarmonics> basis_multipole = sp_mpoles[dp * num_q + dq]->copy();

                    basis_multipole->scale(prefactor * Dp[p * nbf + q]);
                    mpoles_[N]->add(basis_multipole);
                } // end q
            } // end p
        } // end N
    }
}

void CFMMBox::compute_mpoles_from_children() {

    int nmat = 0;
    // Upward pass: translate each child's regular expansion to the parent
    // center and accumulate it. This is the multipole-to-multipole step of the
    // CFMM hierarchy described by White et al. (1994, 1996).
    for (std::shared_ptr<CFMMBox> child : children_) {
        nmat = std::max(nmat, (int)child->mpoles_.size());
    }

    if (mpoles_.size() == 0) {
        mpoles_.resize(nmat);
        Vff_.resize(nmat);
    }

    // Create multipoles and far field vectors for each density matrix
    for (int N = 0; N < nmat; N++) {
        mpoles_[N] = std::make_shared<RealSolidHarmonics>(lmax_, center_, SolidHarmonicsType::Regular);
        Vff_[N] = std::make_shared<RealSolidHarmonics>(lmax_, center_, SolidHarmonicsType::Irregular);
    }

    for (std::shared_ptr<CFMMBox> child : children_) {
        if (child->nshell_pair() == 0) continue;
        for (int N = 0; N < nmat; N++) {
            std::shared_ptr<RealSolidHarmonics> child_multipoles = child->mpoles_[N]->translate(center_);
            mpoles_[N]->add(child_multipoles);
        }
    }
}

CFMMTree::CFMMTree(std::shared_ptr<BasisSet> primary, std::shared_ptr<BasisSet> auxiliary, Options& options)
                    : primary_(primary), auxiliary_(auxiliary), options_(options) {

    if (primary && auxiliary) {
        contraction_type_ = ContractionType::DF_AUX_PRI;
    } else if (primary) {
        contraction_type_ = ContractionType::DIRECT;
    } else if (auxiliary) {
        contraction_type_ = ContractionType::METRIC;
    } else {
        throw PSIEXCEPTION("CFMMTree requires a primary and/or auxiliary basis set.");
    }

    molecule_ = primary_ ? primary_->molecule() : auxiliary_->molecule();
    nlevels_ = options_.get_int("CFMM_GRAIN");
    if (nlevels_ <= 2) {
        throw PSIEXCEPTION("CFMM_GRAIN must be at least 3.");
    } else if (nlevels_ >= 6) {
        throw PSIEXCEPTION("CFMM_GRAIN must not exceed 5.");
    }
    lmax_ = options_.get_int("CFMM_ORDER");
    if (lmax_ < 0) {
        throw PSIEXCEPTION("CFMM_ORDER must be nonnegative.");
    }

    nthread_ = 1;
#ifdef _OPENMP
    nthread_ = Process::environment.get_n_threads();
#endif

    print_ = options_.get_int("PRINT");
    bench_ = options_.get_int("BENCH");

    density_screening_ = (options_.get_str("SCREENING") == "DENSITY");
    ints_tolerance_ =
        (options_.get_str("SCREENING") == "NONE") ? 0.0 : options_.get_double("INTS_TOLERANCE");

    int num_boxes = (nlevels_ == 1) ? 1 : (0.5 * std::pow(16, nlevels_) + 7) / 15;
    tree_.resize(num_boxes);

    mpole_coefs_ = std::make_shared<HarmonicCoefficients>(lmax_, SolidHarmonicsType::Regular);
    double cfmm_extent_tol = options.get_double("CFMM_EXTENT_TOLERANCE");
    if (cfmm_extent_tol <= 0.0 || cfmm_extent_tol >= 1.0) {
        throw PSIEXCEPTION("CFMM_EXTENT_TOLERANCE must be greater than zero and less than one.");
    }

    std::shared_ptr<BasisSet> zero = BasisSet::zero_ao_basis_set();

    std::shared_ptr<IntegralFactory> factory;
    if (contraction_type_ == ContractionType::DIRECT) {
        factory = std::make_shared<IntegralFactory>(primary_);
    } else if (contraction_type_ == ContractionType::DF_AUX_PRI) {
        factory = std::make_shared<IntegralFactory>(auxiliary_, zero, primary_, primary_);
    } else if (contraction_type_ == ContractionType::METRIC) {
        factory = std::make_shared<IntegralFactory>(auxiliary_, zero, auxiliary_, zero);
    }

    std::shared_ptr<TwoBodyAOInt> shellpair_int = std::shared_ptr<TwoBodyAOInt>(factory->eri());

    const auto& ints_bra_shell_pairs = shellpair_int->shell_pairs_bra();
    size_t bra_nshell_pairs = ints_bra_shell_pairs.size();

    if (contraction_type_ == ContractionType::DIRECT) primary_shell_pairs_.resize(bra_nshell_pairs);
    else auxiliary_shell_pairs_.resize(bra_nshell_pairs);

#pragma omp parallel for num_threads(nthread_)
    for (size_t pair_index = 0; pair_index < bra_nshell_pairs; pair_index++) {
        const auto& pair = ints_bra_shell_pairs[pair_index];
        if (contraction_type_ == ContractionType::DIRECT) {
            primary_shell_pairs_[pair_index] = std::make_shared<ShellPair>(primary_, primary_, pair, mpole_coefs_, cfmm_extent_tol);
        } else {
            auxiliary_shell_pairs_[pair_index] = std::make_shared<ShellPair>(auxiliary_, zero, pair, mpole_coefs_, cfmm_extent_tol);
        }
    }

    if (contraction_type_ == ContractionType::DF_AUX_PRI) {
        const auto& ints_ket_shell_pairs = shellpair_int->shell_pairs_ket();
        size_t ket_nshell_pairs = ints_ket_shell_pairs.size();
        primary_shell_pairs_.resize(ket_nshell_pairs);

#pragma omp parallel for num_threads(nthread_)
        for (size_t pair_index = 0; pair_index < ket_nshell_pairs; pair_index++) {
            const auto& pair = ints_ket_shell_pairs[pair_index];
            primary_shell_pairs_[pair_index] = std::make_shared<ShellPair>(primary_, primary_, pair, mpole_coefs_, cfmm_extent_tol);
        }
    }

    timer_on("CFMMTree: Setup");

    make_root_node();
    make_children();
    sort_leaf_boxes();
    setup_regions();
    setup_local_far_field_task_pairs();
    setup_shellpair_info();
    if (contraction_type_ == ContractionType::DIRECT || contraction_type_ == ContractionType::DF_AUX_PRI)
        calculate_shellpair_multipoles(true);
    if (contraction_type_ == ContractionType::METRIC || contraction_type_ == ContractionType::DF_AUX_PRI)
        calculate_shellpair_multipoles(false);

    timer_off("CFMMTree: Setup");

    if (print_ >= 2) print_out();
}

void CFMMTree::df_set_contraction(ContractionType contraction_type) {
    if (contraction_type_ != ContractionType::DF_PRI_AUX && contraction_type_ != ContractionType::DF_AUX_PRI) {
        throw PSIEXCEPTION("Cannot reset the contraction type of a non-three-index DF CFMM tree.");
    }
    if (contraction_type == ContractionType::DIRECT || contraction_type == ContractionType::METRIC) {
        throw PSIEXCEPTION("Cannot reset a DF CFMM tree to a non-DF contraction type.");
    }
    contraction_type_ = contraction_type;
}

void CFMMTree::sort_leaf_boxes() {

    // Starting and ending leaf node box indices
    int start = (nlevels_ == 1) ? 0 : (0.5 * std::pow(16, nlevels_-1) + 7) / 15;
    int end = (nlevels_ == 1) ? 1 : (0.5 * std::pow(16, nlevels_) + 7) / 15;

    for (int bi = start; bi < end; bi++) {
        std::shared_ptr<CFMMBox> box = tree_[bi];
        if (box->nshell_pair() > 0) sorted_leaf_boxes_.push_back(box);
    }

    auto box_compare = [](const std::shared_ptr<CFMMBox> &a,
                                const std::shared_ptr<CFMMBox> &b) { return a->nshell_pair() > b->nshell_pair(); };

    std::sort(sorted_leaf_boxes_.begin(), sorted_leaf_boxes_.end(), box_compare);

}

void CFMMTree::make_root_node() {
    double min_dim = molecule_->x(0);
    double max_dim = molecule_->x(0);

    for (int atom = 0; atom < molecule_->natom(); atom++) {
        double x = molecule_->x(atom);
        double y = molecule_->y(atom);
        double z = molecule_->z(atom);
        min_dim = std::min(x, min_dim);
        min_dim = std::min(y, min_dim);
        min_dim = std::min(z, min_dim);
        max_dim = std::max(x, max_dim);
        max_dim = std::max(y, max_dim);
        max_dim = std::max(z, max_dim);
    }

    max_dim += 0.1; // Add a small buffer to the box

    Vector3 origin = Vector3(min_dim, min_dim, min_dim);
    double length = (max_dim - min_dim);

    tree_[0] = std::make_shared<CFMMBox>(nullptr, primary_shell_pairs_, auxiliary_shell_pairs_, origin, length, 0, lmax_, 2);
}

void CFMMTree::make_children() {

#pragma omp parallel num_threads(nthread_)
    {
        for (int level = 0; level <= nlevels_ - 2; level += 1) {
            int start, end;
            if (level == 0) {
                start = 0;
                end = 1;
            } else {
                start = (0.5 * std::pow(16, level) + 7) / 15;
                end = (0.5 * std::pow(16, level+1) + 7) / 15;
            }

#pragma omp for
            for (int bi = start; bi < end; bi++) {
                tree_[bi]->make_children();
                auto children = tree_[bi]->get_children();

                for (int ci = 0; ci < children.size(); ci++) {
                    int ti = (level == 0) ? ci + 1 : bi * 16 - 7 + ci;
                    tree_[ti] = children[ci];
                }
            }
        }
    }

}

void CFMMTree::setup_regions() {

#pragma omp parallel num_threads(nthread_)
    {
        for (int level = 0; level <= nlevels_ - 1; level += 1) {
            int start, end;
            if (level == 0) {
                start = 0;
                end = 1;
            } else {
                start = (0.5 * std::pow(16, level) + 7) / 15;
                end = (0.5 * std::pow(16, level+1) + 7) / 15;
            }

#pragma omp for
            for (int bi = start; bi < end; bi++) {
                if (tree_[bi]->nshell_pair() == 0) continue;
                tree_[bi]->set_regions();
            }
        }
    }

}

void CFMMTree::setup_shellpair_info() {

    int primary_task_index = 0;
    int auxiliary_task_index = 0;
    for (int i = 0; i < sorted_leaf_boxes_.size(); i++) {
        std::shared_ptr<CFMMBox> curr = sorted_leaf_boxes_[i];
        auto& primary_shellpairs = curr->get_primary_shell_pairs();
        auto& auxiliary_shellpairs = curr->get_auxiliary_shell_pairs();
        auto& nf_boxes = curr->near_field_boxes();

        for (auto& primary_sp : primary_shellpairs) {
            auto PQ = primary_sp->get_shell_pair_index();
            int P = PQ.first;
            int Q = PQ.second;

            primary_shellpair_tasks_.emplace_back(P, Q);
            primary_shellpair_list_.push_back(primary_sp);
            primary_shellpair_to_box_.push_back(curr);
            primary_shellpair_to_nf_boxes_.push_back({});

            for (int nfi = 0; nfi < nf_boxes.size(); nfi++) {
                std::shared_ptr<CFMMBox> neighbor = nf_boxes[nfi];
                if (neighbor->nshell_pair() == 0) continue;
                primary_shellpair_to_nf_boxes_[primary_task_index].push_back(neighbor);
            }
            primary_task_index += 1;
        }

        for (auto& auxiliary_sp : auxiliary_shellpairs) {
            auto RS = auxiliary_sp->get_shell_pair_index();
            int R = RS.first;
            int S = RS.second;

            auxiliary_shellpair_tasks_.emplace_back(R, S);
            auxiliary_shellpair_list_.push_back(auxiliary_sp);
            auxiliary_shellpair_to_box_.push_back(curr);
            auxiliary_shellpair_to_nf_boxes_.push_back({});

            for (int nfi = 0; nfi < nf_boxes.size(); nfi++) {
                std::shared_ptr<CFMMBox> neighbor = nf_boxes[nfi];
                if (neighbor->nshell_pair() == 0) continue;
                auxiliary_shellpair_to_nf_boxes_[auxiliary_task_index].push_back(neighbor);
            }
            auxiliary_task_index += 1;
        }
    }

}

bool CFMMTree::shell_significant(int P, int Q, int R, int S, std::vector<std::shared_ptr<TwoBodyAOInt>>& ints,
                             const std::vector<SharedMatrix>& D) {
    // Density screening multiplies the Schwarz shell-quartet bound by the
    // largest relevant density block. This is the direct-SCF screening strategy
    // of Haser and Ahlrichs, J. Comput. Chem. 10, 104 (1989),
    // doi:10.1002/jcc.540100111, expressed through Psi4's integral engine.
    if (density_screening_) {
        double D_PQ = 0.0;
        double D_RS = 0.0;

        double prefactor = (D.size() == 1) ? 4.0 : 2.0;

        for (int i = 0; i < D.size(); i++) {
            D_PQ += ints[0]->shell_pair_max_density(i, P, Q);
            D_RS += ints[0]->shell_pair_max_density(i, R, S);
        }

        double screen_val = prefactor * std::max(D_PQ, D_RS) * std::sqrt(ints[0]->shell_ceiling2(P, Q, R, S));

        if (screen_val >= ints_tolerance_) return true;
        else return false;

    } else {
        return ints[0]->shell_significant(P, Q, R, S);
    }
}

void CFMMTree::setup_local_far_field_task_pairs() {

    // First access is level, second access is the list of local far field per box for that box
    lff_task_pairs_per_level_.resize(nlevels_);

    // Build the task pairs
    for (int level = 0; level < nlevels_; level++) {
        int start, end;
        if (level == 0) {
            start = 0;
            end = 1;
        } else {
            start = (0.5 * std::pow(16, level) + 7) / 15;
            end = (0.5 * std::pow(16, level+1) + 7) / 15;
        }

        for (int bi = start; bi < end; bi++) {
            std::shared_ptr<CFMMBox> box = tree_[bi];
            if (box->nshell_pair() == 0) continue;
            for (auto& lff : box->local_far_field_boxes()) {
                if (lff->nshell_pair() == 0) continue;
                lff_task_pairs_per_level_[level].emplace_back(box, lff);
            }
        }
    }
}

void CFMMTree::calculate_shellpair_multipoles(bool is_primary) {

    timer_on("CFMMTree: Shell-Pair Multipole Ints");

    // Build multipole integrals for the selected primary or auxiliary basis.
    std::vector<std::shared_ptr<OneBodyAOInt>> sints;
    std::vector<std::shared_ptr<OneBodyAOInt>> mpints;

    std::shared_ptr<IntegralFactory> int_factory;

    if (is_primary) {
        int_factory = std::make_shared<IntegralFactory>(primary_);
    } else {
        auto zero = BasisSet::zero_ao_basis_set();
        int_factory = std::make_shared<IntegralFactory>(auxiliary_, zero, auxiliary_, zero);
    }

    for (int thread = 0; thread < nthread_; thread++) {
        sints.push_back(std::shared_ptr<OneBodyAOInt>(int_factory->ao_overlap()));
        mpints.push_back(std::shared_ptr<OneBodyAOInt>(int_factory->ao_multipoles(lmax_)));
    }

    std::vector<std::pair<int, int>>& shellpair_tasks = (is_primary) ? primary_shellpair_tasks_ : auxiliary_shellpair_tasks_;
    std::vector<std::shared_ptr<ShellPair>>& shellpair_list = (is_primary) ? primary_shellpair_list_ : auxiliary_shellpair_list_;
    std::vector<std::shared_ptr<CFMMBox>>& shellpair_to_box = (is_primary) ? primary_shellpair_to_box_ : auxiliary_shellpair_to_box_;

#pragma omp parallel for num_threads(nthread_) schedule(guided)
    for (int task = 0; task < shellpair_tasks.size(); task++) {
        std::shared_ptr<ShellPair> shellpair = shellpair_list[task];
        std::shared_ptr<CFMMBox> box = shellpair_to_box[task];

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        mpints[thread]->set_origin(box->center());
        shellpair->calculate_mpoles(box->center(), sints[thread], mpints[thread], lmax_);
    }

    timer_off("CFMMTree: Shell-Pair Multipole Ints");

}

void CFMMTree::calculate_multipoles(const std::vector<SharedMatrix>& D) {
    timer_on("CFMMTree: Box Multipoles");

    // Build density-contracted regular expansions in the leaf boxes.
#pragma omp parallel num_threads(nthread_)
    {
#pragma omp for
        for (int bi = 0; bi < sorted_leaf_boxes_.size(); bi++) {
            sorted_leaf_boxes_[bi]->compute_multipoles(D, contraction_type_);
        }

        // Upward pass: aggregate leaf moments into successively coarser boxes.
        for (int level = nlevels_ - 2; level >= 0; level -= 1) {
            int start, end;
            if (level == 0) {
                start = 0;
                end = 1;
            } else {
                start = (0.5 * std::pow(16, level) + 7) / 15;
                end = (0.5 * std::pow(16, level+1) + 7) / 15;
            }

#pragma omp for
            for (int bi = start; bi < end; bi++) {
                if (tree_[bi]->nshell_pair() == 0) continue;
                tree_[bi]->compute_mpoles_from_children();
            }
        }
    }

    timer_off("CFMMTree: Box Multipoles");
}

void CFMMTree::compute_far_field() {

    timer_on("CFMMTree: Far Field Vector");

#pragma omp parallel num_threads(nthread_)
    {
        // Interaction pass: convert the regular multipoles of each local
        // far-field source box into an irregular local expansion at the target.
        // Downward pass: propagate the parent local expansion to its children.
        for (int level = 0; level < nlevels_; level++) {
            const auto& all_box_pairs = lff_task_pairs_per_level_[level];
#pragma omp for
            for (int box_pair = 0; box_pair < all_box_pairs.size(); box_pair++) {
                std::shared_ptr<CFMMBox> box1 = all_box_pairs[box_pair].first;
                std::shared_ptr<CFMMBox> box2 = all_box_pairs[box_pair].second;
                box1->compute_far_field_contribution(box2);
            }

            int start, end;
            if (level == 0) {
                start = 0;
                end = 1;
            } else {
                start = (0.5 * std::pow(16, level) + 7) / 15;
                end = (0.5 * std::pow(16, level+1) + 7) / 15;
            }

#pragma omp for
            for (int bi = start; bi < end; bi++) {
                if (tree_[bi]->nshell_pair() == 0) continue;
                tree_[bi]->add_parent_far_field_contribution();
            }
        }
    }

    timer_off("CFMMTree: Far Field Vector");

}

void CFMMTree::build_nf_J(std::vector<std::shared_ptr<TwoBodyAOInt>>& ints,
                          const std::vector<SharedMatrix>& D, std::vector<SharedMatrix>& J,
                          const std::vector<double>& metric_shell_diagonal_max) {
    if (contraction_type_ == ContractionType::DIRECT) build_nf_direct_J(ints, D, J);
    else if (contraction_type_ == ContractionType::DF_AUX_PRI)
        build_nf_gamma_P(ints, D, J, metric_shell_diagonal_max);
    else if (contraction_type_ == ContractionType::DF_PRI_AUX)
        build_nf_df_J(ints, D, J, metric_shell_diagonal_max);
    else if (contraction_type_ == ContractionType::METRIC)
        build_nf_metric(ints, D, J);
}

void CFMMTree::build_nf_direct_J(std::vector<std::shared_ptr<TwoBodyAOInt>>& ints,
                          const std::vector<SharedMatrix>& D, std::vector<SharedMatrix>& J) {

    timer_on("CFMMTree: Near Field Direct J");

    // Explicitly evaluate the near-field part of
    //     J_pq = sum_rs (pq|rs) D_rs.
    // The complementary, well-separated part is supplied by the CFMM
    // multipole expansion; see White et al. (1994, 1996).

    int nshell = primary_->nshell();
    // Maximum space (r_nbf * s_nbf) to allocate per task
    size_t max_alloc = 0;

    // A map of the function (num_r * num_s) offsets per shell-pair in a box pair
    std::unordered_map<int, int> offsets;

    int start = (nlevels_ == 1) ? 0 : (0.5 * std::pow(16, nlevels_-1) + 7) / 15;
    int end = (nlevels_ == 1) ? 1 : (0.5 * std::pow(16, nlevels_) + 7) / 15;

    for (int bi = start; bi < end; bi++) {
        auto& RSshells = tree_[bi]->get_primary_shell_pairs();
        int RSoff = 0;
        for (int RSind = 0; RSind < RSshells.size(); RSind++) {
            std::pair<int, int> RS = RSshells[RSind]->get_shell_pair_index();
            int R = RS.first;
            int S = RS.second;
            offsets[R * nshell + S] = RSoff;
            int Rfunc = primary_->shell(R).nfunction();
            int Sfunc = primary_->shell(S).nfunction();
            RSoff += Rfunc * Sfunc;
        }
        max_alloc = std::max((size_t) RSoff, max_alloc);
    }

    // Thread-local accumulation buffers exploit eightfold ERI permutational symmetry.
    std::vector<std::vector<std::vector<double>>> JT;
    for (int thread = 0; thread < nthread_; thread++) {
        std::vector<std::vector<double>> J2;
        for (size_t N = 0; N <D.size(); N++) {
            std::vector<double> temp(2 * max_alloc);
            J2.push_back(temp);
        }
        JT.push_back(J2);
    }

    // Number of computed shell quartets for optional benchmark reporting.
    size_t computed_shells = 0L;

#pragma omp parallel for num_threads(nthread_) schedule(guided) reduction(+ : computed_shells)
    for (int task = 0; task < primary_shellpair_tasks_.size(); task++) {
        int P = primary_shellpair_tasks_[task].first;
        int Q = primary_shellpair_tasks_[task].second;

        const GaussianShell& Pshell = primary_->shell(P);
        const GaussianShell& Qshell = primary_->shell(Q);

        int p_start = Pshell.start();
        int num_p = Pshell.nfunction();

        int q_start = Qshell.start();
        int num_q = Qshell.nfunction();

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        for (const auto& nf_box : primary_shellpair_to_nf_boxes_[task]) {
            auto& RSshells = nf_box->get_primary_shell_pairs();

            bool touched = false;
            for (const auto& RSsh : RSshells) {
                auto RS = RSsh->get_shell_pair_index();
                int R = RS.first;
                int S = RS.second;

                if (R * nshell + S > P * nshell + Q) continue;
                if (!shell_significant(P, Q, R, S, ints, D)) continue;

                if (ints[thread]->compute_shell(P, Q, R, S) == 0) continue;
                computed_shells++;

                const GaussianShell& Rshell = primary_->shell(R);
                const GaussianShell& Sshell = primary_->shell(S);

                int r_start = Rshell.start();
                int num_r = Rshell.nfunction();

                int s_start = Sshell.start();
                int num_s = Sshell.nfunction();

                double prefactor = 1.0;
                if (P != Q) prefactor *= 2;
                if (R != S) prefactor *= 2;
                if (P == R && Q == S) prefactor *= 0.5;

                int RSoff = offsets.at(R * nshell + S);

                const double* pqrs = ints[thread]->buffer();

                for (int N = 0; N < D.size(); N++) {
                    double** Dp = D[N]->pointer();
                    double* JTp = JT[thread][N].data();
                    const double* pqrs2 = pqrs;

                    if (!touched) {
                        std::memset(static_cast<void*>(&JTp[0L * max_alloc]), 0, max_alloc * sizeof(double));
                        std::memset(static_cast<void*>(&JTp[1L * max_alloc]), 0, max_alloc * sizeof(double));
                    }

                    // Contraction into box shell pairs to improve parallel performance
                    double* J1p = &JTp[0L * max_alloc];
                    double* J2p = &JTp[1L * max_alloc];

                    for (int p = p_start; p < p_start + num_p; p++) {
                        int dp = p - p_start;
                        for (int q = q_start; q < q_start + num_q; q++) {
                            int dq = q - q_start;
                            for (int r = r_start; r < r_start + num_r; r++) {
                                int dr = r - r_start;
                                for (int s = s_start; s < s_start + num_s; s++) {
                                    int ds = s - s_start;

                                    int pq = dp * num_q + dq;
                                    int rs = RSoff + dr * num_s + ds;

                                    J1p[pq] += prefactor * (*pqrs2) * Dp[r][s];
                                    J2p[rs] += prefactor * (*pqrs2) * Dp[p][q];
                                    pqrs2++;
                                } // end s
                            } // end r
                        } // end q
                    } // end p
                } // end N
                touched = true;
            } // end RSshells
            if (!touched) continue;

            // Accumulate the thread-local shell blocks into J.
            for (int N = 0; N < D.size(); N++) {
                double** Jp = J[N]->pointer();
                double* JTp = JT[thread][N].data();

                double* J1p = &JTp[0L * max_alloc];
                double* J2p = &JTp[1L * max_alloc];

                for (int p = p_start; p < p_start + num_p; p++) {
                    int dp = p - p_start;
                    for (int q = q_start; q < q_start + num_q; q++) {
                        int dq = q - q_start;

                        int pq = dp * num_q + dq;
#pragma omp atomic
                        Jp[p][q] += J1p[pq];
                    }
                }

                for (const auto& RSsh : RSshells) {
                    std::pair<int, int> RS = RSsh->get_shell_pair_index();
                    int R = RS.first;
                    int S = RS.second;

                    int RSoff = offsets.at(R * nshell + S);

                    const GaussianShell& Rshell = primary_->shell(R);
                    const GaussianShell& Sshell = primary_->shell(S);

                    int r_start = Rshell.start();
                    int num_r = Rshell.nfunction();

                    int s_start = Sshell.start();
                    int num_s = Sshell.nfunction();

                    for (int r = r_start; r < r_start + num_r; r++) {
                        int dr = r - r_start;
                        for (int s = s_start; s < s_start + num_s; s++) {
                            int ds = s - s_start;
                            int rs = RSoff + dr * num_s + ds;
#pragma omp atomic
                            Jp[r][s] += J2p[rs];
                        }
                    }
                }
            }
        } // end nf_box
    } // end primary shell-pair tasks

    num_computed_shells_ = computed_shells;

    if (bench_) {
        auto mode = std::ostream::app;
        auto printer = PsiOutStream("bench.dat", mode);
        size_t ntri = nshell * (nshell + 1L) / 2L;
        size_t possible_shells = ntri * (ntri + 1L) / 2L;
        double computed_fraction = ((double) computed_shells) / possible_shells;
        printer.Printf("CFMM Near Field: Computed %20zu Shell Quartets out of %20zu, (%11.3E ratio)\n",
                        computed_shells, possible_shells, computed_fraction);
    }

    timer_off("CFMMTree: Near Field Direct J");
}

void CFMMTree::build_nf_gamma_P(std::vector<std::shared_ptr<TwoBodyAOInt>>& ints,
                                const std::vector<SharedMatrix>& densities,
                                std::vector<SharedMatrix>& gamma_p,
                                const std::vector<double>& metric_shell_diagonal_max) {
    timer_on("DF CFMM: Near Field Gamma P");

    // First density-fitting contraction:
    //     gamma_P = sum_uv (P|uv) D_uv.
    // Together with build_nf_df_J, this implements the two three-index
    // contractions used by DF-CFMM; see Sodt, Subotnik, and Head-Gordon,
    // J. Chem. Phys. 125, 194109 (2006), doi:10.1063/1.2370949, and Lazarski,
    // Burow, and Sierka, J. Chem. Theory Comput. 11, 3029 (2015),
    // doi:10.1021/acs.jctc.5b00252.
    // gamma_P is a scratch intermediate, so it must be reset for every build.
    for (size_t ind = 0; ind < densities.size(); ind++) {
        gamma_p[ind]->zero();
    }

    int pri_nshell = primary_->nshell();

    const bool screen_integrals = ints_tolerance_ > 0.0;

    // The density-weighted three-center bound requires one maximum diagonal
    // metric value for every auxiliary shell.
    if (screen_integrals && metric_shell_diagonal_max.size() != static_cast<size_t>(auxiliary_->nshell())) {
        throw PSIEXCEPTION(
            "CFMMTree::build_nf_gamma_P requires one metric-shell bound per auxiliary shell when screening is "
            "enabled.");
    }

    // Maximum density magnitude in each primary shell-pair block UV.
    Matrix density_shell_pair_max(pri_nshell, pri_nshell);
    auto density_shell_pair_maxp = density_shell_pair_max.pointer();

    if (screen_integrals) {
        for (size_t U = 0; U < pri_nshell; U++) {
            int u_start = primary_->shell(U).start();
            int num_u = primary_->shell(U).nfunction();

            for (size_t V = 0; V < pri_nshell; V++) {
                int v_start = primary_->shell(V).start();
                int num_v = primary_->shell(V).nfunction();

                for (size_t i = 0; i < densities.size(); i++) {
                    auto Dp = densities[i]->pointer();
                    for (size_t u = u_start; u < u_start + num_u; u++) {
                        for (size_t v = v_start; v < v_start + num_v; v++) {
                            density_shell_pair_maxp[U][V] =
                                std::max(density_shell_pair_maxp[U][V], std::abs(Dp[u][v]));
                        }
                    }
                }
            }
        }
    }

    size_t computed_shells = 0L;

#pragma omp parallel for num_threads(nthread_) schedule(guided) reduction(+ : computed_shells)
    for (int task = 0; task < auxiliary_shellpair_tasks_.size(); task++) {

        int P = auxiliary_shellpair_tasks_[task].first;
        const GaussianShell& Pshell = auxiliary_->shell(P);

        int p_start = Pshell.start();
        int num_p = Pshell.nfunction();

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        for (const auto& nf_box : auxiliary_shellpair_to_nf_boxes_[task]) {
            auto& UVshells = nf_box->get_primary_shell_pairs();

            for (const auto& UVsh : UVshells) {
                auto UV = UVsh->get_shell_pair_index();
                int U = UV.first;
                int V = UV.second;

                if (screen_integrals) {
                    double screen_val = density_shell_pair_maxp[U][V] * density_shell_pair_maxp[U][V] *
                                        metric_shell_diagonal_max[P] * ints[thread]->shell_pair_value(U, V);
                    if (screen_val < ints_tolerance_ * ints_tolerance_) continue;
                }
                computed_shells++;

                int u_start = primary_->shell(U).start();
                int num_u = primary_->shell(U).nfunction();

                int v_start = primary_->shell(V).start();
                int num_v = primary_->shell(V).nfunction();

                ints[thread]->compute_shell(P, 0, U, V);
                const double* buffer = ints[thread]->buffer();

                double prefactor = 2.0;
                if (U == V) prefactor *= 0.5;

                for (int i = 0; i < densities.size(); i++) {
                    double** Dp = densities[i]->pointer();
                    double* three_index_integrals = const_cast<double*>(buffer);
                    double* gamma_data = gamma_p[i]->pointer()[0];

                    std::vector<double> density_block(num_u * num_v, 0.0);
                    double* density_block_data = density_block.data();

                    for (int u = u_start; u < u_start + num_u; u++) {
                        for (int v = v_start; v < v_start + num_v; v++) {
                            *density_block_data = Dp[u][v];
                            density_block_data++;
                        }
                    }
                    C_DGEMV('N', num_p, num_u * num_v, prefactor, three_index_integrals, num_u * num_v,
                            density_block.data(), 1, 1.0, &(gamma_data[p_start]), 1);

                } // end i
            } // UV shells
        } // NF Boxes
    } // task

    num_computed_shells_ = computed_shells;

    timer_off("DF CFMM: Near Field Gamma P");
}

void CFMMTree::build_nf_df_J(std::vector<std::shared_ptr<TwoBodyAOInt>>& ints,
                             const std::vector<SharedMatrix>& auxiliary_coefficients,
                             std::vector<SharedMatrix>& J,
                             const std::vector<double>& metric_shell_diagonal_max) {
    timer_on("DF CFMM: Near Field J");

    // Second density-fitting contraction:
    //     J_uv = sum_Q (uv|Q) gamma_Q.
    // The near-field integral contraction here is completed by the far-field
    // multipole contribution in build_ff_J.

    int pri_nshell = primary_->nshell();
    int aux_nshell = auxiliary_->nshell();
    int nmat = auxiliary_coefficients.size();

    int max_nbf_per_shell = 0;
    for (int P = 0; P < pri_nshell; P++) {
        max_nbf_per_shell = std::max(max_nbf_per_shell, primary_->shell(P).nfunction());
    }

    // Thread-local J buffers used by the BLAS contraction.
    std::vector<std::vector<SharedMatrix>> JT;

    for (int thread = 0; thread < nthread_; thread++) {
        std::vector<SharedMatrix> J2;
        for (size_t ind = 0; ind < nmat; ind++) {
            J2.push_back(std::make_shared<Matrix>(max_nbf_per_shell, max_nbf_per_shell));
            J2[ind]->zero();
        }
        JT.push_back(J2);
    }

    const bool screen_integrals = ints_tolerance_ > 0.0;
    if (screen_integrals && metric_shell_diagonal_max.size() != static_cast<size_t>(aux_nshell)) {
        throw PSIEXCEPTION(
            "CFMMTree::build_nf_df_J requires one metric-shell bound per auxiliary shell when screening is enabled.");
    }

    // Maximum auxiliary-coefficient magnitude in each auxiliary shell.
    std::vector<double> auxiliary_shell_max(aux_nshell, 0.0);
    if (screen_integrals) {
        for (size_t i = 0; i < auxiliary_coefficients.size(); i++) {
            for (int P = 0; P < aux_nshell; P++) {
                int p_start = auxiliary_->shell(P).start();
                int num_p = auxiliary_->shell(P).nfunction();
                for (int p = p_start; p < p_start + num_p; p++) {
                    auxiliary_shell_max[P] =
                        std::max(auxiliary_shell_max[P], std::abs(auxiliary_coefficients[i]->get(p, 0)));
                }
            }
        }
    }

    size_t computed_shells = 0L;

#pragma omp parallel for num_threads(nthread_) schedule(guided) reduction(+ : computed_shells)
    for (int task = 0; task < primary_shellpair_tasks_.size(); task++) {

        int U = primary_shellpair_tasks_[task].first;
        int V = primary_shellpair_tasks_[task].second;

        const GaussianShell& Ushell = primary_->shell(U);
        const GaussianShell& Vshell = primary_->shell(V);

        int u_start = Ushell.start();
        int num_u = Ushell.nfunction();

        int v_start = Vshell.start();
        int num_v = Vshell.nfunction();

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        double prefactor = 2.0;
        if (U == V) prefactor *= 0.5;

        for (const auto& nf_box : primary_shellpair_to_nf_boxes_[task]) {
            auto& Qshells = nf_box->get_auxiliary_shell_pairs();

            for (const auto& Qsh : Qshells) {

                int Q = Qsh->get_shell_pair_index().first;

                if (screen_integrals) {
                    double screen_val = auxiliary_shell_max[Q] * auxiliary_shell_max[Q] *
                                        metric_shell_diagonal_max[Q] * ints[thread]->shell_pair_value(U, V);
                    if (screen_val < ints_tolerance_ * ints_tolerance_) continue;
                }
                computed_shells++;

                int q_start = auxiliary_->shell(Q).start();
                int num_q = auxiliary_->shell(Q).nfunction();

                ints[thread]->compute_shell(Q, 0, U, V);

                const double* buffer = ints[thread]->buffer();

                for (int i = 0; i < auxiliary_coefficients.size(); i++) {
                    double* j_buffer = JT[thread][i]->pointer()[0];
                    double* auxiliary_data = auxiliary_coefficients[i]->pointer()[0];
                    double* three_index_integrals = const_cast<double*>(buffer);

                    C_DGEMV('T', num_q, num_u * num_v, prefactor, three_index_integrals, num_u * num_v,
                            &(auxiliary_data[q_start]), 1, 1.0, j_buffer, 1);

                } // end i
            } // end Qsh
        } // end nf box

        // Accumulate the thread-local buffer into the output matrix.

        for (int i = 0; i < auxiliary_coefficients.size(); i++) {
            double* JTp = JT[thread][i]->pointer()[0];
            double** Jp = J[i]->pointer();
            for (int u = u_start; u < u_start + num_u; u++) {
                int du = u - u_start;
                for (int v = v_start; v < v_start + num_v; v++) {
                    int dv = v - v_start;

                    Jp[u][v] += JTp[du * num_v + dv];
                }
            }
            JT[thread][i]->zero();
        }
    } // end task

    num_computed_shells_ = computed_shells;

    timer_off("DF CFMM: Near Field J");
}

void CFMMTree::build_nf_metric(std::vector<std::shared_ptr<TwoBodyAOInt>>& ints,
                               const std::vector<SharedMatrix>& D, std::vector<SharedMatrix>& J) {
    // This path is deliberately retained as the integration point for a
    // future CFMM-accelerated auxiliary Coulomb-metric contraction.
    throw PSIEXCEPTION("The CFMM auxiliary-metric contraction is not implemented.");
}

void CFMMTree::build_ff_J(std::vector<SharedMatrix>& J) {

    timer_on("CFMMTree: Far Field J");

    bool is_primary = (contraction_type_ == ContractionType::DF_PRI_AUX || contraction_type_ == ContractionType::DIRECT);
    int nbf = (is_primary) ? primary_->nbf() : 1;

    std::shared_ptr<BasisSet>& ref_basis = (is_primary) ? primary_ : auxiliary_;
    const auto& shellpair_tasks = (is_primary) ? primary_shellpair_tasks_ : auxiliary_shellpair_tasks_;
    const auto& shellpair_list = (is_primary) ? primary_shellpair_list_ : auxiliary_shellpair_list_;
    const auto& shellpair_to_box = (is_primary) ? primary_shellpair_to_box_ : auxiliary_shellpair_to_box_;

#pragma omp parallel for num_threads(nthread_) schedule(guided)
    for (int task = 0; task < shellpair_tasks.size(); task++) {
        const auto& shellpair = shellpair_list[task];
        const auto& box = shellpair_to_box[task];
        const auto& Vff = box->far_field_vector();
        const auto& shellpair_mpoles = shellpair->get_mpoles();

        int P = shellpair_tasks[task].first;
        int Q = shellpair_tasks[task].second;

        double prefactor = (!is_primary || P == Q) ? 1.0 : 2.0;

        const GaussianShell& Pshell = shellpair->bs1()->shell(P);
        const GaussianShell& Qshell = shellpair->bs2()->shell(Q);

        int p_start = Pshell.start();
        int num_p = Pshell.nfunction();

        int q_start = Qshell.start();
        int num_q = Qshell.nfunction();

        for (int p = p_start; p < p_start + num_p; p++) {
            int dp = p - p_start;
            for (int q = q_start; q < q_start + num_q; q++) {
                int dq = q - q_start;
                for (int N = 0; N < J.size(); N++) {
                    double* Jp = J[N]->pointer()[0];
                    // Far field multipole contributions
#pragma omp atomic
                    Jp[p * nbf + q] += prefactor * Vff[N]->dot(shellpair_mpoles[dp * num_q + dq]);
                } // end N
            } // end q
        } // end p
    }

    timer_off("CFMMTree: Far Field J");
}

void CFMMTree::build_J(std::vector<std::shared_ptr<TwoBodyAOInt>>& ints,
                        const std::vector<SharedMatrix>& D, std::vector<SharedMatrix>& J,
                        const std::vector<double>& metric_shell_diagonal_max) {

    timer_on("CFMMTree: J");

    num_computed_shells_ = 0L;

    // J is additive here. CompositeJK controls whether a full build begins
    // from zero or an incremental build retains the accumulated matrix. The
    // gamma_P scratch intermediate is reset locally in build_nf_gamma_P.

    // Update the densities
    if (density_screening_ && contraction_type_ == ContractionType::DIRECT) {
        for (int thread = 0; thread < nthread_; thread++) {
            ints[thread]->update_density(D);
        }
    }

    // Compute multipoles and far field
    calculate_multipoles(D);
    compute_far_field();

    // Compute near field J and far field J
    build_nf_J(ints, D, J, metric_shell_diagonal_max);
    build_ff_J(J);

    // Hermitize the square primary-basis result after accumulation.
    if (contraction_type_ == ContractionType::DIRECT || (contraction_type_ == ContractionType::DF_PRI_AUX)) {
        for (int ind = 0; ind < D.size(); ind++) {
            J[ind]->hermitivitize();
        }
    }

    timer_off("CFMMTree: J");
}

void CFMMTree::print_out() {
    for (int bi = 0; bi < tree_.size(); bi++) {
        std::shared_ptr<CFMMBox> box = tree_[bi];
        const auto& shell_pairs = box->get_primary_shell_pairs();
        int nshells = shell_pairs.size();
        int level = box->get_level();
        int ws = box->get_ws();
        if (nshells > 0) {
            outfile->Printf("  BOX INDEX: %d, LEVEL: %d, WS: %d, NSHELLS: %d\n", bi, level, ws, nshells);
        }
    }
}

}  // namespace psi
