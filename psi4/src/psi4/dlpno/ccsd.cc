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

#include "dlpno.h"
#include "sparse.h"

#include "psi4/lib3index/3index.h"
#include "psi4/libdiis/diismanager.h"
#include "psi4/libfock/cubature.h"
#include "psi4/libfock/points.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/local.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/orthog.h"
#include "psi4/libmints/twobody.h"
#include "psi4/libmints/vector.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libqt/qt.h"

#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi {
namespace dlpno {

DLPNOCCSD::DLPNOCCSD(SharedWavefunction ref_wfn, Options& options) : DLPNOBase(ref_wfn, options) {}
DLPNOCCSD::~DLPNOCCSD() {}

inline SharedMatrix DLPNOCCSD::S_PNO(const int ij, const int mn) {
    int i, j, m, n;
    std::tie(i, j) = ij_to_i_j_[ij];
    std::tie(m, n) = ij_to_i_j_[mn];
    
    
    const int m_ij = lmopair_to_lmos_dense_[ij][m], n_ij = lmopair_to_lmos_dense_[ij][n];
    if (m_ij == -1 || n_ij == -1) {
        outfile->Printf("Invalid PNO Pairs (%d, %d) and (%d, %d)\n", i, j, m, n);
        throw PSIEXCEPTION("Invalid PNO pairs!");
    }

    const int nlmo_ij = lmopair_to_lmos_[ij].size();

    int mn_ij; 
    if (m_ij < n_ij) {
        mn_ij = n_ij * nlmo_ij + m_ij;
    } else {
        mn_ij = m_ij * nlmo_ij + n_ij;
    }

    if (i < j) {
        const int ji = ij_to_ji_[ij];
        return S_pno_ij_mn_[ji][mn_ij];
    } else {
        return S_pno_ij_mn_[ij][mn_ij];
    }
}

void DLPNOCCSD::compute_pno_overlaps() {
    const int n_lmo_pairs = ij_to_i_j_.size();
    S_pno_ij_mn_.resize(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        int i, j;
        std::tie(i, j) = ij_to_i_j_[ij];

        const int npno_ij = n_pno_[ij];
        const int nlmo_ij = lmopair_to_lmos_[ij].size();

        if (npno_ij == 0 || i < j) continue;

        S_pno_ij_mn_[ij].resize(nlmo_ij * nlmo_ij);

        for (int mn_ij = 0; mn_ij < nlmo_ij * nlmo_ij; mn_ij++) {
            const int m_ij = mn_ij / nlmo_ij, n_ij = mn_ij % nlmo_ij;
            const int m = lmopair_to_lmos_[ij][m_ij], n = lmopair_to_lmos_[ij][n_ij];
            const int mn = i_j_to_ij_[m][n];

            if (mn == -1 || n_pno_[mn] == 0 || m_ij < n_ij) continue;

            S_pno_ij_mn_[ij][mn_ij] = submatrix_rows_and_cols(*S_pao_, lmopair_to_paos_[ij], lmopair_to_paos_[mn]);
            S_pno_ij_mn_[ij][mn_ij] = linalg::triplet(X_pno_[ij], S_pno_ij_mn_[ij][mn_ij], X_pno_[mn], true, false, false);
        }
    }
}

void DLPNOCCSD::estimate_memory() {
    outfile->Printf("  ==> DLPNO-CCSD Memory Estimate <== \n\n");

    int n_lmo_pairs = ij_to_i_j_.size();

    size_t pno_overlap_memory = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : pno_overlap_memory)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        int i, j;
        std::tie(i, j) = ij_to_i_j_[ij];

        const int npno_ij = n_pno_[ij];
        const int nlmo_ij = lmopair_to_lmos_[ij].size();
        if (npno_ij == 0 || i < j) continue;

        for (int mn_ij = 0; mn_ij < nlmo_ij * nlmo_ij; mn_ij++) {
            const int m_ij = mn_ij / nlmo_ij, n_ij = mn_ij % nlmo_ij;
            const int m = lmopair_to_lmos_[ij][m_ij], n = lmopair_to_lmos_[ij][n_ij];
            const int mn = i_j_to_ij_[m][n];
            if (mn == -1 || n_pno_[mn] == 0 || m_ij < n_ij) continue;

            pno_overlap_memory += n_pno_[ij] * n_pno_[mn];
        }
    }

    size_t oooo = 0;
    size_t ooov = 0;
    size_t oovv = 0;
    size_t ovvv = 0;
    size_t qvv = 0;

#pragma omp parallel for schedule(dynamic) reduction(+ : oooo, ooov, oovv, ovvv, qvv)
    for (int ij = 0; ij < n_lmo_pairs; ij++) {
        int i, j;
        std::tie(i, j) = ij_to_i_j_[ij];

        const int naux_ij = lmopair_to_ribfs_[ij].size();
        const int nlmo_ij = lmopair_to_lmos_[ij].size();
        const int npno_ij = n_pno_[ij];

        oooo += nlmo_ij * nlmo_ij;
        ooov += 2 * nlmo_ij * npno_ij;
        oovv += 4 * npno_ij * npno_ij;
        ovvv += 3 * npno_ij * npno_ij * npno_ij;
        if (i >= j) qvv += naux_ij * npno_ij * npno_ij;
    }

    const size_t total_df_memory = qij_memory_ + qia_memory_ + qab_memory_;
    const size_t total_pno_int_memory = oooo + ooov + oovv + ovvv + qvv;
    const size_t total_memory = total_df_memory + pno_overlap_memory + total_pno_int_memory;

    // 2^30 bytes per GiB
    outfile->Printf("    (q | i j) integrals    : %.3f [GiB]\n", qij_memory_ * pow(2.0, -30) * sizeof(double));
    outfile->Printf("    (q | i a) integrals    : %.3f [GiB]\n", qia_memory_ * pow(2.0, -30) * sizeof(double));
    outfile->Printf("    (q | a b) integrals    : %.3f [GiB]\n", qab_memory_ * pow(2.0, -30) * sizeof(double));
    outfile->Printf("    (m i | n j) integrals  : %.3f [GiB]\n", oooo * pow(2.0, -30) * sizeof(double));
    outfile->Printf("    1-virtual PNO integrals: %.3f [GiB]\n", ooov * pow(2.0, -30) * sizeof(double));
    outfile->Printf("    2-virtual PNO integrals: %.3f [GiB]\n", oovv * pow(2.0, -30) * sizeof(double));
    outfile->Printf("    3-virtual PNO integrals: %.3f [GiB]\n", ovvv * pow(2.0, -30) * sizeof(double));
    outfile->Printf("    4-virtual PNO integrals: %.3f [GiB]\n", qvv * pow(2.0, -30) * sizeof(double));
    outfile->Printf("    PNO/PNO overlaps       : %.3f [GiB]\n\n", pno_overlap_memory * pow(2.0, -30) * sizeof(double));
    outfile->Printf("    Total Memory Given     : %.3f [GiB]\n", memory_ * pow(2.0, -30));
    outfile->Printf("    Total Memory Required  : %.3f [GiB]\n\n", total_memory * pow(2.0, -30) * sizeof(double));
    
    if (total_memory * sizeof(double) < memory_) {
        outfile->Printf("    Storing virtual/virtual integrals...\n\n");
        virtual_storage_ = CORE;
    } else {
        outfile->Printf("    Computing virtual/virtual integrals as needed...\n\n");
        virtual_storage_ = DIRECT;
    }
}

std::vector<double> DLPNOCCSD::compute_pair_energies() {
    int n_lmo_pairs = ij_to_i_j_.size();
    std::vector<double> e_ijs(n_lmo_pairs);

    // Step 1: compute SC-LMP2 pair energies
#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        int i, j;
        std::tie(i, j) = ij_to_i_j_[ij];
        int ji = ij_to_ji_[ij];

        if (i > j) continue;

        //                                                   //
        // ==> Assemble (ia|jb) for pair ij in PAO basis <== //
        //                                                   //

        // number of PAOs in the pair domain (before removing linear dependencies)
        int npao_ij = lmopair_to_paos_[ij].size();  // X_pao_ij->rowspi(0);

        // number of auxiliary basis in the domain
        int naux_ij = lmopair_to_ribfs_[ij].size();

        auto i_qa = std::make_shared<Matrix>("Three-index Integrals", naux_ij, npao_ij);
        auto j_qa = std::make_shared<Matrix>("Three-index Integrals", naux_ij, npao_ij);

        for (int q_ij = 0; q_ij < naux_ij; q_ij++) {
            int q = lmopair_to_ribfs_[ij][q_ij];
            int centerq = ribasis_->function_to_center(q);
            for (int a_ij = 0; a_ij < npao_ij; a_ij++) {
                int a = lmopair_to_paos_[ij][a_ij];
                i_qa->set(q_ij, a_ij, qia_[q]->get(riatom_to_lmos_ext_dense_[centerq][i], riatom_to_paos_ext_dense_[centerq][a]));
                j_qa->set(q_ij, a_ij, qia_[q]->get(riatom_to_lmos_ext_dense_[centerq][j], riatom_to_paos_ext_dense_[centerq][a]));
            }
        }

        auto A_solve = submatrix_rows_and_cols(*full_metric_, lmopair_to_ribfs_[ij], lmopair_to_ribfs_[ij]);
        C_DGESV_wrapper(A_solve, i_qa);

        auto K_pao_ij = linalg::doublet(i_qa, j_qa, true, false);

        //                                      //
        // ==> Canonicalize PAOs of pair ij <== //
        //                                      //

        auto S_pao_ij = submatrix_rows_and_cols(*S_pao_, lmopair_to_paos_[ij], lmopair_to_paos_[ij]);
        auto F_pao_ij = submatrix_rows_and_cols(*F_pao_, lmopair_to_paos_[ij], lmopair_to_paos_[ij]);

        SharedMatrix X_pao_ij;  // canonical transformation of this domain's PAOs to
        SharedVector e_pao_ij;  // energies of the canonical PAOs
        std::tie(X_pao_ij, e_pao_ij) = orthocanonicalizer(S_pao_ij, F_pao_ij);

        // S_pao_ij = linalg::triplet(X_pao_ij, S_pao_ij, X_pao_ij, true, false, false);
        F_pao_ij = linalg::triplet(X_pao_ij, F_pao_ij, X_pao_ij, true, false, false);
        K_pao_ij = linalg::triplet(X_pao_ij, K_pao_ij, X_pao_ij, true, false, false);

        // number of PAOs in the domain after removing linear dependencies
        int npao_can_ij = X_pao_ij->colspi(0);
        auto T_pao_ij = K_pao_ij->clone();
        for (int a = 0; a < npao_can_ij; ++a) {
            for (int b = 0; b < npao_can_ij; ++b) {
                T_pao_ij->set(a, b, T_pao_ij->get(a, b) /
                                        (-e_pao_ij->get(b) + -e_pao_ij->get(a) + F_lmo_->get(i, i) + F_lmo_->get(j, j)));
            }
        }

        size_t nvir_ij = K_pao_ij->rowspi(0);

        auto Tt_pao_ij = T_pao_ij->clone();
        Tt_pao_ij->scale(2.0);
        Tt_pao_ij->subtract(T_pao_ij->transpose());

        // mp2 energy of this LMO pair before transformation to PNOs
        double e_ij = K_pao_ij->vector_dot(Tt_pao_ij);

        e_ijs[ij] = e_ij;
        if (i != j) e_ijs[ji] = e_ij;
    }

    return e_ijs;
}

double DLPNOCCSD::filter_pairs(const std::vector<double>& e_ijs, const std::vector<std::vector<int>>& strong_pairs,
                                double tolerance) {
    int natom = molecule_->natom();
    int nbf = basisset_->nbf();
    int naocc = i_j_to_ij_.size();
    int n_lmo_pairs = ij_to_i_j_.size();
    int naux = ribasis_->nbf();
    int npao = C_pao_->colspi(0);  // same as nbf

    // Step 2. Split up strong and weak pairs based on e_ijs
    int strong_pair_count = 0, weak_pair_count = 0;
    double delta_e = 0.0;

    ij_to_i_j_strong_.clear();
    ij_to_i_j_weak_.clear();

    i_j_to_ij_strong_.clear();
    i_j_to_ij_weak_.clear();

    i_j_to_ij_strong_.resize(naocc);
    i_j_to_ij_weak_.resize(naocc);

    for (size_t i = 0; i < naocc; i++) {
        i_j_to_ij_strong_[i].resize(naocc, -1);
        i_j_to_ij_weak_[i].resize(naocc, -1);

        for (size_t j = 0; j < naocc; j++) {
            int ij = i_j_to_ij_[i][j];
            if (ij == -1) continue;

            if (std::fabs(e_ijs[ij]) >= tolerance && strong_pairs[i][j] != -1) { // Strong Pair
                i_j_to_ij_strong_[i][j] = strong_pair_count;
                ij_to_i_j_strong_.push_back(std::make_pair(i,j));
                ++strong_pair_count;
            } else { // Weak Pair
                i_j_to_ij_weak_[i][j] = weak_pair_count;
                ij_to_i_j_weak_.push_back(std::make_pair(i,j));
                delta_e += e_ijs[ij];
                ++weak_pair_count;
            }
        }
    }

    ij_to_ji_strong_.clear();
    ij_to_ji_weak_.clear();

    for (size_t ij = 0; ij < ij_to_i_j_strong_.size(); ++ij) {
        size_t i, j;
        std::tie(i, j) = ij_to_i_j_strong_[ij];
        ij_to_ji_strong_.push_back(i_j_to_ij_strong_[j][i]);
    }
    
    for (size_t ij = 0; ij < ij_to_i_j_weak_.size(); ++ij) {
        size_t i, j;
        std::tie(i, j) = ij_to_i_j_weak_[ij];
        ij_to_ji_weak_.push_back(i_j_to_ij_weak_[j][i]);
    }

    return delta_e;
}

void DLPNOCCSD::reset_sparsity() {
    int natom = molecule_->natom();
    int nbf = basisset_->nbf();
    int nshell = basisset_->nshell();
    int naux = ribasis_->nbf();
    int naocc = nalpha_ - nfrzc();
    int npao = C_pao_->colspi(0);  // same as nbf

    auto bf_to_atom = std::vector<int>(nbf);
    auto ribf_to_atom = std::vector<int>(naux);

    for (size_t i = 0; i < nbf; ++i) {
        bf_to_atom[i] = basisset_->function_to_center(i);
    }

    for (size_t i = 0; i < naux; ++i) {
        ribf_to_atom[i] = ribasis_->function_to_center(i);
    }

    // map from LMO to local DF domain (aux basis functions)
    // locality determined via mulliken charges

    lmo_to_ribfs_.clear();
    lmo_to_riatoms_.clear();

    lmo_to_ribfs_.resize(naocc);
    lmo_to_riatoms_.resize(naocc);

    for (size_t i = 0; i < naocc; ++i) {
        // atomic mulliken populations for this orbital
        std::vector<double> mkn_pop(natom, 0.0);

        auto P_i = reference_wavefunction_->S()->clone();

        for (size_t u = 0; u < nbf; u++) {
            P_i->scale_row(0, u, C_lmo_->get(u, i));
            P_i->scale_column(0, u, C_lmo_->get(u, i));
        }

        for (size_t u = 0; u < nbf; u++) {
            int centerU = basisset_->function_to_center(u);
            double p_uu = P_i->get(u, u);

            for (size_t v = 0; v < nbf; v++) {
                int centerV = basisset_->function_to_center(v);
                double p_vv = P_i->get(v, v);

                // off-diag pops (p_uv) split between u and v prop to diag pops
                double p_uv = P_i->get(u, v);
                mkn_pop[centerU] += p_uv * ((p_uu) / (p_uu + p_vv));
                mkn_pop[centerV] += p_uv * ((p_vv) / (p_uu + p_vv));
            }
        }

        // if non-zero mulliken pop on atom, include atom in the LMO's fitting domain
        for (size_t a = 0; a < natom; a++) {
            if (fabs(mkn_pop[a]) > T_CUT_MKN_) {
                lmo_to_riatoms_[i].push_back(a);

                // each atom's aux orbitals are all-or-nothing for each LMO
                for (int u : atom_to_ribf_[a]) {
                    lmo_to_ribfs_[i].push_back(u);
                }
            }
        }
    }

    // map from LMO to local virtual domain (PAOs)
    // locality determined via differential overlap integrals

    lmo_to_paos_.resize(naocc);
    lmo_to_paoatoms_.resize(naocc);
    for (size_t i = 0; i < naocc; ++i) {
        // PAO domains determined by differential overlap integral
        std::vector<int> lmo_to_paos_temp;
        for (size_t u = 0; u < nbf; ++u) {
            if (fabs(DOI_iu_->get(i, u)) > T_CUT_DO_) {
                lmo_to_paos_temp.push_back(u);
            }
        }

        // if any PAO on an atom is in the list, we take all of the PAOs on that atom
        lmo_to_paos_[i] = contract_lists(lmo_to_paos_temp, atom_to_bf_);

        // contains the same information as previous map
        lmo_to_paoatoms_[i] = block_list(lmo_to_paos_[i], bf_to_atom);
    }
}

void DLPNOCCSD::recompute_pair_domains() {
    int natom = molecule_->natom();
    int nbf = basisset_->nbf();
    int naocc = i_j_to_ij_.size();
    int n_lmo_pairs = ij_to_i_j_.size();
    int naux = ribasis_->nbf();
    int npao = C_pao_->colspi(0);  // same as nbf

    // Step 3. Recompute Sparse Maps based on new info
    ij_to_i_j_ = ij_to_i_j_strong_;
    i_j_to_ij_ = i_j_to_ij_strong_;
    ij_to_ji_ = ij_to_ji_strong_;

    n_lmo_pairs = ij_to_i_j_.size();

    // map from (LMO, LMO) pair to local auxiliary and virtual domains
    // LMO pair domains are the union of LMO domains

    lmopair_to_paos_.resize(n_lmo_pairs);
    lmopair_to_paoatoms_.resize(n_lmo_pairs);
    lmopair_to_ribfs_.resize(n_lmo_pairs);
    lmopair_to_riatoms_.resize(n_lmo_pairs);

#pragma omp parallel for
    for (size_t ij = 0; ij < n_lmo_pairs; ++ij) {
        size_t i, j;
        std::tie(i, j) = ij_to_i_j_[ij];

        lmopair_to_paos_[ij] = merge_lists(lmo_to_paos_[i], lmo_to_paos_[j]);
        lmopair_to_paoatoms_[ij] = merge_lists(lmo_to_paoatoms_[i], lmo_to_paoatoms_[j]);

        lmopair_to_ribfs_[ij] = merge_lists(lmo_to_ribfs_[i], lmo_to_ribfs_[j]);
        lmopair_to_riatoms_[ij] = merge_lists(lmo_to_riatoms_[i], lmo_to_riatoms_[j]);
    }

    // Create a list of lmos that "interact" with a lmo_pair by differential overlap
    lmopair_to_lmos_.clear();
    lmopair_to_lmos_.resize(n_lmo_pairs);
    lmopair_to_lmos_dense_.resize(n_lmo_pairs);

#pragma omp parallel for
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        int i, j;
        std::tie(i, j) = ij_to_i_j_[ij];
        lmopair_to_lmos_[ij].clear();
        lmopair_to_lmos_dense_[ij] = std::vector<int>(naocc, -1);

        int m_ij = 0;
        for (int m = 0; m < naocc; ++m) {
            int im = i_j_to_ij_[i][m];
            int jm = i_j_to_ij_[j][m];

            if (im != -1 && jm != -1) {
                lmopair_to_lmos_[ij].push_back(m);
                lmopair_to_lmos_dense_[ij][m] = m_ij;
                m_ij++;
            }
        }
    }

    print_aux_pair_domains();
    print_pao_pair_domains();

    // determine maps to extended LMO domains, which are the union of an LMO's domain with domains
    //   of all interacting LMOs
    lmo_to_riatoms_ext_ = extend_maps(lmo_to_riatoms_, ij_to_i_j_);

    riatom_to_lmos_ext_ = invert_map(lmo_to_riatoms_ext_, natom);
    riatom_to_paos_ext_ = chain_maps(riatom_to_lmos_ext_, lmo_to_paos_);

    // We'll use these maps to screen the the local MO transform (first index):
    //   (mn|Q) * C_mi -> (in|Q)
    riatom_to_atoms1_ = chain_maps(riatom_to_lmos_ext_, lmo_to_atoms_);
    riatom_to_shells1_ = chain_maps(riatom_to_atoms1_, atom_to_shell_);
    riatom_to_bfs1_ = chain_maps(riatom_to_atoms1_, atom_to_bf_);

    // We'll use these maps to screen the projected AO transform (second index):
    //   (mn|Q) * C_nu -> (mu|Q)
    riatom_to_atoms2_ = chain_maps(riatom_to_lmos_ext_, chain_maps(lmo_to_paos_, pao_to_atoms_));
    riatom_to_shells2_ = chain_maps(riatom_to_atoms2_, atom_to_shell_);
    riatom_to_bfs2_ = chain_maps(riatom_to_atoms2_, atom_to_bf_);

    // Need dense versions of previous maps for quick lookup

    // riatom_to_lmos_ext_dense_[riatom][lmo] is the index of lmo in riatom_to_lmos_ext_[riatom]
    //   (if present), else -1
    riatom_to_lmos_ext_dense_.resize(natom);
    // riatom_to_paos_ext_dense_[riatom][pao] is the index of pao in riatom_to_paos_ext_[riatom]
    //   (if present), else -1
    riatom_to_paos_ext_dense_.resize(natom);

    // riatom_to_atoms1_dense_(1,2)[riatom][a] is true if the orbitals basis functions of atom A
    //   are needed for the (LMO,PAO) transform
    riatom_to_atoms1_dense_.resize(natom);
    riatom_to_atoms2_dense_.resize(natom);

    for (int a_ri = 0; a_ri < natom; a_ri++) {
        riatom_to_lmos_ext_dense_[a_ri] = std::vector<int>(naocc, -1);
        riatom_to_paos_ext_dense_[a_ri] = std::vector<int>(npao, -1);
        riatom_to_atoms1_dense_[a_ri] = std::vector<bool>(natom, false);
        riatom_to_atoms2_dense_[a_ri] = std::vector<bool>(natom, false);

        for (int i_ind = 0; i_ind < riatom_to_lmos_ext_[a_ri].size(); i_ind++) {
            int i = riatom_to_lmos_ext_[a_ri][i_ind];
            riatom_to_lmos_ext_dense_[a_ri][i] = i_ind;
        }
        for (int u_ind = 0; u_ind < riatom_to_paos_ext_[a_ri].size(); u_ind++) {
            int u = riatom_to_paos_ext_[a_ri][u_ind];
            riatom_to_paos_ext_dense_[a_ri][u] = u_ind;
        }
        for (int a_bf : riatom_to_atoms1_[a_ri]) {
            riatom_to_atoms1_dense_[a_ri][a_bf] = true;
        }
        for (int a_bf : riatom_to_atoms2_[a_ri]) {
            riatom_to_atoms2_dense_[a_ri][a_bf] = true;
        }
    }

    lmopair_lmo_to_riatom_lmo_.clear();
    lmopair_pao_to_riatom_pao_.clear();

    lmopair_lmo_to_riatom_lmo_.resize(n_lmo_pairs);
    lmopair_pao_to_riatom_pao_.resize(n_lmo_pairs);
#pragma omp parallel for
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        int naux_ij = lmopair_to_ribfs_[ij].size();
        lmopair_lmo_to_riatom_lmo_[ij].resize(naux_ij);
        lmopair_pao_to_riatom_pao_[ij].resize(naux_ij);

        for (int q_ij = 0; q_ij < naux_ij; q_ij++) {
            int q = lmopair_to_ribfs_[ij][q_ij];
            int centerq = ribasis_->function_to_center(q);

            int nlmo_ij = lmopair_to_lmos_[ij].size();
            int npao_ij = lmopair_to_paos_[ij].size();
            lmopair_lmo_to_riatom_lmo_[ij][q_ij].resize(nlmo_ij);
            lmopair_pao_to_riatom_pao_[ij][q_ij].resize(npao_ij);

            for (int m_ij = 0; m_ij < nlmo_ij; m_ij++) {
                int m = lmopair_to_lmos_[ij][m_ij];
                int m_sparse = riatom_to_lmos_ext_dense_[centerq][m];
                lmopair_lmo_to_riatom_lmo_[ij][q_ij][m_ij] = m_sparse;
            }

            for (int a_ij = 0; a_ij < npao_ij; a_ij++) {
                int a = lmopair_to_paos_[ij][a_ij];
                int a_sparse = riatom_to_paos_ext_dense_[centerq][a];
                lmopair_pao_to_riatom_pao_[ij][q_ij][a_ij] = a_sparse;
            }
        }
    }
}

void DLPNOCCSD::ccsd_pair_prescreening() {
    int natom = molecule_->natom();
    int nbf = basisset_->nbf();
    int naocc = i_j_to_ij_.size();
    int n_lmo_pairs = ij_to_i_j_.size();
    int naux = ribasis_->nbf();
    int npao = C_pao_->colspi(0);  // same as nbf

    outfile->Printf("\n  ==> Determining Strong and Weak Pairs (Crude Prescreening Step) <==\n\n");
    const std::vector<double>& e_ijs_crude = compute_pair_energies();
    de_lmp2_ = filter_pairs(e_ijs_crude, i_j_to_ij_, T_CUT_PAIRS_MP2_);

    int n_strong_pairs = ij_to_i_j_strong_.size();
    int n_weak_pairs = ij_to_i_j_weak_.size();

    outfile->Printf("    Number of Eliminated Pairs   = %d\n", n_weak_pairs);
    outfile->Printf("    Number of Remaining Pairs    = %d\n", n_strong_pairs);
    outfile->Printf("    Strong Pairs / Total Pairs   = (%.2f %%)\n", (100.0 * n_strong_pairs) / (naocc * naocc));
    outfile->Printf("    SC-LMP2 Weak Pair Correction = %.12f\n\n", de_lmp2_);

    // Reset Sparsity After 
    T_CUT_MKN_ *= 0.01;
    T_CUT_DO_ *= 0.5;

    reset_sparsity();
    recompute_pair_domains();

    // Recompute 3-center integrals (with refined pair domains)
    compute_qia();

    outfile->Printf("\n  ==> Determining Strong and Weak Pairs (Refined Prescreening Step) <==\n\n");
    const std::vector<double>& e_ijs = compute_pair_energies();
    const auto i_j_to_ij_strong_copy = i_j_to_ij_strong_;
    de_lmp2_ += filter_pairs(e_ijs, i_j_to_ij_strong_copy, T_CUT_PAIRS_);

    n_strong_pairs = ij_to_i_j_strong_.size();
    n_weak_pairs = ij_to_i_j_weak_.size();

    outfile->Printf("    Number of Eliminated Pairs   = %d\n", n_weak_pairs);
    outfile->Printf("    Number of Remaining Pairs    = %d\n", n_strong_pairs);
    outfile->Printf("    Strong Pairs / Total Pairs   = (%.2f %%)\n", (100.0 * n_strong_pairs) / (naocc * naocc));
    outfile->Printf("    SC-LMP2 Weak Pair Correction = %.12f\n\n", de_lmp2_);

    recompute_pair_domains();
}

void DLPNOCCSD::compute_cc_integrals() {
    outfile->Printf("    Computing CC integrals...\n\n");

    int n_lmo_pairs = ij_to_i_j_.size();
    // 0 virtual
    K_mnij_.resize(n_lmo_pairs);
    // 1 virtual
    K_bar_.resize(n_lmo_pairs);
    L_bar_.resize(n_lmo_pairs);
    // 2 virtual
    J_ijab_.resize(n_lmo_pairs);
    L_iajb_.resize(n_lmo_pairs);
    M_iajb_.resize(n_lmo_pairs);
    // 3 virtual
    K_tilde_chem_.resize(n_lmo_pairs);
    K_tilde_phys_.resize(n_lmo_pairs);
    L_tilde_.resize(n_lmo_pairs);
    // 4 virtual
    if (virtual_storage_ == CORE) Qab_ij_.resize(n_lmo_pairs);

    size_t qvv_memory = 0;
    size_t qvv_svd_memory = 0;

#pragma omp parallel for schedule(dynamic, 1) reduction(+ : qvv_memory) reduction(+ : qvv_svd_memory)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        int i, j;
        std::tie(i, j) = ij_to_i_j_[ij];
        const int ji = ij_to_ji_[ij];

        // number of PNOs in the pair domain
        const int npno_ij = n_pno_[ij];
        if (npno_ij == 0) continue;

        // number of LMOs in the pair domain
        const int nlmo_ij = lmopair_to_lmos_[ij].size();
        // number of PAOs in the pair domain (before removing linear dependencies)
        const int npao_ij = lmopair_to_paos_[ij].size();
        // number of auxiliary functions in the pair domain
        const int naux_ij = lmopair_to_ribfs_[ij].size();

        auto q_pair = std::make_shared<Matrix>(naux_ij, 1);

        auto q_io = std::make_shared<Matrix>(naux_ij, nlmo_ij);
        auto q_jo = std::make_shared<Matrix>(naux_ij, nlmo_ij);

        auto q_jv = std::make_shared<Matrix>(naux_ij, npno_ij);
        auto q_vv = std::make_shared<Matrix>(naux_ij, npno_ij * npno_ij);

        for (int q_ij = 0; q_ij < naux_ij; q_ij++) {
            const int q = lmopair_to_ribfs_[ij][q_ij];
            const int centerq = ribasis_->function_to_center(q);

            const int i_sparse = riatom_to_lmos_ext_dense_[centerq][i];
            const int j_sparse = riatom_to_lmos_ext_dense_[centerq][j];
            const std::vector<int> i_slice(1, i_sparse);
            const std::vector<int> j_slice(1, j_sparse);

            q_pair->set(q_ij, 0, (*qij_[q])(i_sparse, j_sparse));
            
            auto q_io_tmp = submatrix_rows_and_cols(*qij_[q], i_slice, 
                                lmopair_lmo_to_riatom_lmo_[ij][q_ij]);
            C_DCOPY(nlmo_ij, &(*q_io_tmp)(0,0), 1, &(*q_io)(q_ij, 0), 1);

            auto q_jo_tmp = submatrix_rows_and_cols(*qij_[q], j_slice, 
                                lmopair_lmo_to_riatom_lmo_[ij][q_ij]);
            C_DCOPY(nlmo_ij, &(*q_jo_tmp)(0,0), 1, &(*q_jo)(q_ij, 0), 1);

            auto q_jv_tmp = submatrix_rows_and_cols(*qia_[q], j_slice,
                                lmopair_pao_to_riatom_pao_[ij][q_ij]);
            q_jv_tmp = linalg::doublet(q_jv_tmp, X_pno_[ij], false, false);
            C_DCOPY(npno_ij, &(*q_jv_tmp)(0,0), 1, &(*q_jv)(q_ij, 0), 1);

            SharedMatrix q_vv_tmp;
            if (T_CUT_EIG_ > 0.0) {
                SharedMatrix P;
                SharedVector D;
                std::tie(P, D) = qab_svd_[q];
                auto qab_temp1 = linalg::doublet(X_pno_[ij], submatrix_rows(*P, lmopair_pao_to_riatom_pao_[ij][q_ij]), true, false);
                auto qab_temp2 = qab_temp1->clone();
                // auto Dtemp = std::make_shared<Matrix>(D->dim(), D->dim());
                // Dtemp->set_diagonal(D);
                for (int i = 0; i < D->dim(); ++i) qab_temp1->scale_column(0, i, D->get(i));

                q_vv_tmp = linalg::doublet(qab_temp1, qab_temp2, false, true);
            } else {
                q_vv_tmp = submatrix_rows_and_cols(*qab_[q], lmopair_pao_to_riatom_pao_[ij][q_ij],
                                lmopair_pao_to_riatom_pao_[ij][q_ij]);
                q_vv_tmp = linalg::triplet(X_pno_[ij], q_vv_tmp, X_pno_[ij], true, false, false);
            }
            C_DCOPY(npno_ij * npno_ij, &(*q_vv_tmp)(0,0), 1, &(*q_vv)(q_ij, 0), 1);
        }

        auto A_solve = submatrix_rows_and_cols(*full_metric_, lmopair_to_ribfs_[ij], lmopair_to_ribfs_[ij]);
        A_solve->power(0.5, 1.0e-14);

        C_DGESV_wrapper(A_solve->clone(), q_pair);
        C_DGESV_wrapper(A_solve->clone(), q_io);
        C_DGESV_wrapper(A_solve->clone(), q_jo);
        C_DGESV_wrapper(A_solve->clone(), q_jv);
        C_DGESV_wrapper(A_solve, q_vv);

        K_mnij_[ij] = linalg::doublet(q_io, q_jo, true, false);
        K_bar_[ij] = linalg::doublet(q_io, q_jv, true, false);
        J_ijab_[ij] = linalg::doublet(q_pair, q_vv, true, false);
        J_ijab_[ij]->reshape(npno_ij, npno_ij);

        K_tilde_chem_[ji] = linalg::doublet(q_jv, q_vv, true, false);
        K_tilde_phys_[ji] = std::make_shared<Matrix>(npno_ij, npno_ij * npno_ij);
        auto K_tilde_temp = std::make_shared<Matrix>(npno_ij, npno_ij * npno_ij);

        for (int a_ij = 0; a_ij < npno_ij; a_ij++) {
            for (int e_ij = 0; e_ij < npno_ij; e_ij++) {
                for (int f_ij = 0; f_ij < npno_ij; f_ij++) {
                    (*K_tilde_phys_[ji])(a_ij, e_ij * npno_ij + f_ij) = 
                                    (*K_tilde_chem_[ji])(e_ij, a_ij * npno_ij + f_ij);
                    (*K_tilde_temp)(a_ij, e_ij * npno_ij + f_ij) = 
                                    (*K_tilde_chem_[ji])(f_ij, a_ij * npno_ij + e_ij);
                }
            }
        }
        L_tilde_[ji] = K_tilde_chem_[ji]->clone();
        L_tilde_[ji]->scale(2.0);
        L_tilde_[ji]->subtract(K_tilde_temp);

        if (virtual_storage_ == CORE && i <= j) {
            // SVD Decomposition of DF-ERIs
            // DOI: 10.1063/1.4905005
            auto [U, S, V] = q_vv->svd_temps();
            q_vv->svd(U, S, V);

            int nsvd_ij = 0;
            std::vector<int> slice_indices;
            while (nsvd_ij < S->dim() && S->get(nsvd_ij) >= T_CUT_SVD_) {
                U->scale_column(0, nsvd_ij, S->get(nsvd_ij));
                slice_indices.push_back(nsvd_ij);
                nsvd_ij += 1;
            }

            // U(Q_ij, r_ij) S(r_ij) V(r_ij, a_ij * e_ij)
            // U(Q_ij, s_ij) S(s_ij) V(s_ij, b_ij * f_ij)
            U = submatrix_cols(*U, slice_indices);
            auto B_rs = linalg::doublet(U, U, true, false);
            B_rs->power(0.5, 1.0e-14);
            
            V = linalg::doublet(B_rs, submatrix_rows(*V, slice_indices));

            Qab_ij_[ij].resize(nsvd_ij);
            for (int q_ij = 0; q_ij < nsvd_ij; q_ij++) {
                Qab_ij_[ij][q_ij] = std::make_shared<Matrix>(npno_ij, npno_ij);
                C_DCOPY(npno_ij * npno_ij, &(*V)(q_ij, 0), 1, &(*Qab_ij_[ij][q_ij])(0,0), 1);
            }

            qvv_memory += naux_ij * npno_ij * npno_ij;
            qvv_svd_memory += nsvd_ij * npno_ij * npno_ij;
        }

        // L_iajb
        L_iajb_[ij] = K_iajb_[ij]->clone();
        L_iajb_[ij]->scale(2.0);
        L_iajb_[ij]->subtract(K_iajb_[ij]->transpose());

        // Lt_iajb
        M_iajb_[ij] = K_iajb_[ij]->clone();
        M_iajb_[ij]->scale(2.0);
        M_iajb_[ij]->subtract(J_ijab_[ij]);
    }

    // Antisymmetrize K_mbij integrals
#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        int i, j;
        std::tie(i, j) = ij_to_i_j_[ij];
        const int ji = ij_to_ji_[ij];

        // number of PNOs in the pair domain
        const int npno_ij = n_pno_[ij];
        if (npno_ij == 0) continue;

        L_bar_[ij] = K_bar_[ij]->clone();
        L_bar_[ij]->scale(2.0);
        L_bar_[ij]->subtract(K_bar_[ji]);
    }

    outfile->Printf("    SVD Memory/DF Memory : %8.4f \n\n", static_cast<double>(qvv_svd_memory) / qvv_memory);

}

SharedMatrix DLPNOCCSD::compute_Fmi(const std::vector<SharedMatrix>& tau_tilde) {
    timer_on("Compute Fmi");

    int n_lmo_pairs = ij_to_i_j_.size();
    int naocc = nalpha_ - nfrzc();

    // Equation 40, Term 1
    auto Fmi = F_lmo_->clone();

#pragma omp parallel for schedule(dynamic, 1)
    for (int mi = 0; mi < n_lmo_pairs; ++mi) {
        int m, i;
        std::tie(m, i) = ij_to_i_j_[mi];

        if (n_pno_[mi] == 0) continue;

        for (int n_mi = 0; n_mi < lmopair_to_lmos_[mi].size(); n_mi++) {
            int n = lmopair_to_lmos_[mi][n_mi];
            int mn = i_j_to_ij_[m][n], nn = i_j_to_ij_[n][n], in = i_j_to_ij_[i][n];
            int nm = ij_to_ji_[mn];

            if (n_pno_[mn] == 0 || n_pno_[in] == 0) continue;

            int i_mn = lmopair_to_lmos_dense_[mn][i];

            // Equation 40, Term 3
            std::vector<int> i_mn_slice(1, i_mn);
            auto l_mn_temp = submatrix_rows(*L_bar_[mn], i_mn_slice);
            auto S_nn_mn = S_PNO(nn, mn);
            l_mn_temp = linalg::doublet(S_nn_mn, l_mn_temp, false, true);
            (*Fmi)(m,i) += l_mn_temp->vector_dot(T_ia_[n]);

            // Equation 40, Term 4
            auto S_in_mn = S_PNO(in, mn);
            l_mn_temp = linalg::triplet(S_in_mn, L_iajb_[mn], S_in_mn, false, false, true);
            (*Fmi)(m,i) += l_mn_temp->vector_dot(tau_tilde[in]);
        }
    }

    timer_off("Compute Fmi");

    return Fmi;
}

std::vector<SharedMatrix> DLPNOCCSD::compute_Fbe(const std::vector<SharedMatrix>& tau_tilde) {
    timer_on("Compute Fbe");

    int n_lmo_pairs = ij_to_i_j_.size();
    int naocc = nalpha_ - nfrzc();

    std::vector<SharedMatrix> Fbe(n_lmo_pairs);

#pragma omp parallel for
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        int i, j;
        std::tie(i, j) = ij_to_i_j_[ij];

        int npno_ij = n_pno_[ij];
        if (npno_ij == 0) continue;

        // Equation 39, Term 1
        Fbe[ij] = std::make_shared<Matrix>("Fbe", npno_ij, npno_ij);
        Fbe[ij]->zero();
        Fbe[ij]->set_diagonal(e_pno_[ij]);

        for (int m_ij = 0; m_ij < lmopair_to_lmos_[ij].size(); m_ij++) {
            int m = lmopair_to_lmos_[ij][m_ij];
            int im = i_j_to_ij_[i][m], mm = i_j_to_ij_[m][m], mj = i_j_to_ij_[m][j];
            int npno_mm = n_pno_[mm], npno_mj = n_pno_[mj];

            if (n_pno_[mj] != 0) {
                auto S_ij_mj = S_PNO(ij, mj);
                auto S_mm_mj = S_PNO(mm, mj);
                auto T_m_temp = linalg::doublet(S_mm_mj, T_ia_[m], true, false);

                auto Fbe_mj = linalg::doublet(T_m_temp, L_tilde_[mj], true, false);
                Fbe_mj->reshape(n_pno_[mj], n_pno_[mj]);

                Fbe[ij]->add(linalg::triplet(S_ij_mj, Fbe_mj, S_ij_mj, false, false, true));
            }

            for (int n_ij = 0; n_ij < lmopair_to_lmos_[ij].size(); n_ij++) {
                int n = lmopair_to_lmos_[ij][n_ij];
                int mn = i_j_to_ij_[m][n];

                if (mn != -1 && n_pno_[mn] != 0) {
                    auto S_ij_mn = S_PNO(ij, mn);
                    auto tau_L_temp = linalg::triplet(tau_tilde[mn], L_iajb_[mn], S_ij_mn, false, true, true);
                    Fbe[ij]->subtract(linalg::doublet(S_ij_mn, tau_L_temp, false, false));
                }
            }
        }

    }

    timer_off("Compute Fbe");

    return Fbe;
}

std::vector<SharedMatrix> DLPNOCCSD::compute_Fme() {
    timer_on("Compute Fme");

    int n_lmo_pairs = ij_to_i_j_.size();
    int naocc = nalpha_ - nfrzc();

    std::vector<SharedMatrix> Fme(n_lmo_pairs);

#pragma omp parallel for
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        int i, j;
        std::tie(i, j) = ij_to_i_j_[ij];

        int nlmo_ij = lmopair_to_lmos_[ij].size();
        int npno_ij = n_pno_[ij];

        if (npno_ij == 0) continue;

        Fme[ij] = std::make_shared<Matrix>("Fme", nlmo_ij, npno_ij);
        Fme[ij]->zero();

        for (int m_ij = 0; m_ij < lmopair_to_lmos_[ij].size(); m_ij++) {
            int m = lmopair_to_lmos_[ij][m_ij];
            int mm = i_j_to_ij_[m][m], mj = i_j_to_ij_[m][j];

            for (int n_ij = 0; n_ij < lmopair_to_lmos_[ij].size(); n_ij++) {
                int n = lmopair_to_lmos_[ij][n_ij];
                int mn = i_j_to_ij_[m][n], nn = i_j_to_ij_[n][n];

                if (mn != -1 && n_pno_[mn] != 0) {
                    auto S_mn_nn = S_PNO(mn, nn);
                    auto T_n_temp = linalg::doublet(S_mn_nn, T_ia_[n], false, false);

                    auto S_mn_ij = S_PNO(mn, ij);
                    auto F_me_temp = linalg::triplet(S_mn_ij, L_iajb_[mn], T_n_temp, true, false, false);
                    C_DAXPY(npno_ij, 1.0, &(*F_me_temp)(0,0), 1, &(*Fme[ij])(m_ij, 0), 1);
                }
            }
        }
    }

    timer_off("Compute Fme");

    return Fme;
}

std::vector<SharedMatrix> DLPNOCCSD::compute_Wmnij(const std::vector<SharedMatrix>& tau) {
    timer_on("Compute Wmnij");
    int n_lmo_pairs = ij_to_i_j_.size();
    int naocc = nalpha_ - nfrzc();

    std::vector<SharedMatrix> Wmnij(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int mn = 0; mn < n_lmo_pairs; ++mn) {
        int m, n;
        std::tie(m, n) = ij_to_i_j_[mn];
        int nm = ij_to_ji_[mn];

        int npno_mn = n_pno_[mn];
        if (npno_mn == 0) continue;

        Wmnij[mn] = K_mnij_[mn]->clone();

        int nlmo_mn = lmopair_to_lmos_[mn].size();
        auto T_i_mn = std::make_shared<Matrix>(nlmo_mn, npno_mn);
        T_i_mn->zero();

        for (int i_mn = 0; i_mn < lmopair_to_lmos_[mn].size(); i_mn++) {
            int i = lmopair_to_lmos_[mn][i_mn];
            int ii = i_j_to_ij_[i][i], im = i_j_to_ij_[i][m], in = i_j_to_ij_[i][n];

            auto S_ii_mn = S_PNO(ii, mn);
            auto T_temp = linalg::doublet(S_ii_mn, T_ia_[i], true, false);
            C_DCOPY(npno_mn, &(*T_temp)(0,0), 1, &(*T_i_mn)(i_mn, 0), 1);
        }

        Wmnij[mn]->add(linalg::doublet(K_bar_[mn], T_i_mn, false, true));
        Wmnij[mn]->add(linalg::doublet(T_i_mn, K_bar_[nm], false, true));

        for (int ij_mn = 0; ij_mn < nlmo_mn * nlmo_mn; ++ij_mn) {
            int i_mn = ij_mn / nlmo_mn, j_mn = ij_mn % nlmo_mn;
            int i = lmopair_to_lmos_[mn][i_mn], j = lmopair_to_lmos_[mn][j_mn];
            int ij = i_j_to_ij_[i][j];
            
            if (ij != -1 && n_pno_[ij] != 0) {
                auto S_mn_ij = S_PNO(mn, ij);
                auto K_temp = linalg::triplet(S_mn_ij, K_iajb_[mn], S_mn_ij, true, false, false);
                (*Wmnij[mn])(i_mn, j_mn) += K_temp->vector_dot(tau[ij]);
            }
        }
    }

    timer_off("Compute Wmnij");

    return Wmnij;
}

std::vector<SharedMatrix> DLPNOCCSD::compute_Wmbej(const std::vector<SharedMatrix>& tau_bar) {
    timer_on("Compute Wmbej");
    int n_lmo_pairs = ij_to_i_j_.size();
    int naocc = nalpha_ - nfrzc();

    std::vector<SharedMatrix> Wmbej(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int mj = 0; mj < n_lmo_pairs; ++mj) {
        int m, j;
        std::tie(m, j) = ij_to_i_j_[mj];
        int mm = i_j_to_ij_[m][m], jj = i_j_to_ij_[j][j], jm = ij_to_ji_[mj];
        int npno_mm = n_pno_[mm], npno_mj = n_pno_[mj];

        if (npno_mj == 0) continue;

        Wmbej[mj] = K_iajb_[jm]->clone();
        
        auto S_mm_mj = S_PNO(mm, mj);
        auto S_mj_jj = S_PNO(mj, jj);
        auto tia_temp = linalg::doublet(S_mj_jj, T_ia_[j], false, false);

        auto K_temp1 = std::make_shared<Matrix>(npno_mj, npno_mj);
        K_temp1 = K_tilde_chem_[mj]->clone();
        K_temp1->reshape(npno_mj * npno_mj, npno_mj);
        K_temp1 = linalg::doublet(K_temp1, tia_temp, false, false);
        K_temp1->reshape(npno_mj, npno_mj);
        auto eye = std::make_shared<Matrix>(npno_mj, npno_mj);
        eye->identity();
        K_temp1 = linalg::doublet(eye, K_temp1, false, true);
        
        Wmbej[mj]->add(K_temp1);

        for (int n_mj = 0; n_mj < lmopair_to_lmos_[mj].size(); n_mj++) {
            int n = lmopair_to_lmos_[mj][n_mj];
            int mn = i_j_to_ij_[m][n], nn = i_j_to_ij_[n][n], jn = i_j_to_ij_[j][n];
            int nm = ij_to_ji_[mn], nj = ij_to_ji_[jn];

            auto S_nn_mj = S_PNO(nn, mj);
            auto t_n_temp = linalg::doublet(S_nn_mj, T_ia_[n], true, false);
            C_DGER(npno_mj, npno_mj, -1.0, &(*t_n_temp)(0, 0), 1, &(*K_bar_[jm])(n_mj, 0), 1, &(*Wmbej[mj])(0, 0), npno_mj);

            if (n_pno_[mn] == 0 || n_pno_[nj] == 0) continue;

            auto S_mn_jn = S_PNO(mn, jn);
            auto S_mj_mn = S_PNO(mj, mn);
            auto S_jn_mj = S_PNO(jn, mj);
            auto tau_temp = linalg::triplet(S_mn_jn, tau_bar[jn], S_jn_mj, false, false, false);
            auto K_mn_temp = linalg::doublet(S_mj_mn, K_iajb_[mn], false, false);
            Wmbej[mj]->subtract(linalg::doublet(tau_temp, K_mn_temp, true, true));

            auto T_nj_temp = linalg::triplet(S_mn_jn, T_iajb_[nj], S_jn_mj, false, false, false);
            auto L_mn_temp = linalg::doublet(S_mj_mn, L_iajb_[mn], false, false);
            auto TL_temp = linalg::doublet(T_nj_temp, L_mn_temp, true, true);
            TL_temp->scale(0.5);
            Wmbej[mj]->add(TL_temp);
        }
    }

    timer_off("Compute Wmbej");

    return Wmbej;
}

std::vector<SharedMatrix> DLPNOCCSD::compute_Wmbje(const std::vector<SharedMatrix>& tau_bar) {

    timer_on("Compute Wmbje");

    int n_lmo_pairs = ij_to_i_j_.size();
    int naocc = nalpha_ - nfrzc();

    std::vector<SharedMatrix> Wmbje(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int mj = 0; mj < n_lmo_pairs; ++mj) {
        int m, j;
        std::tie(m, j) = ij_to_i_j_[mj];
        int mm = i_j_to_ij_[m][m], jj = i_j_to_ij_[j][j], jm = ij_to_ji_[mj];
        int npno_mm = n_pno_[mm], npno_mj = n_pno_[mj];

        if (npno_mj == 0) continue;

        Wmbje[mj] = J_ijab_[mj]->clone();
        Wmbje[mj]->scale(-1.0);
        
        auto S_mm_mj = S_PNO(mm, mj);
        auto S_mj_jj = S_PNO(mj, jj);
        auto tia_temp = linalg::doublet(S_mj_jj, T_ia_[j], false, false);

        auto K_temp1 = K_tilde_chem_[mj]->clone();
        K_temp1 = linalg::doublet(tia_temp, K_temp1, true, false);
        K_temp1->reshape(npno_mj, npno_mj);
        auto eye = std::make_shared<Matrix>(npno_mj, npno_mj);
        eye->identity();
        K_temp1 = linalg::doublet(eye, K_temp1, false, false);

        Wmbje[mj]->subtract(K_temp1);

        for (int n_mj = 0; n_mj < lmopair_to_lmos_[mj].size(); n_mj++) {
            int n = lmopair_to_lmos_[mj][n_mj];
            int mn = i_j_to_ij_[m][n], nn = i_j_to_ij_[n][n], jn = i_j_to_ij_[j][n];
            int nm = ij_to_ji_[mn], nj = ij_to_ji_[jn];

            auto S_nn_mj = S_PNO(nn, mj);
            auto t_n_temp = linalg::doublet(S_nn_mj, T_ia_[n], true, false);

            if (n_pno_[nj] == 0 || n_pno_[mn] == 0) continue;

            auto S_mn_mj = S_PNO(mn, mj);
            auto S_mn_jn = S_PNO(mn, jn);
            auto S_jn_jm = S_PNO(jn, jm);

            int j_mn = lmopair_to_lmos_dense_[mn][j];
            std::vector<int> j_mn_slice(1, j_mn);
            auto K_mn_temp = linalg::doublet(submatrix_rows(*K_bar_[mn], j_mn_slice), S_mn_mj)->transpose();
            C_DGER(npno_mj, npno_mj, 1.0, &(*t_n_temp)(0,0), 1, &(*K_mn_temp)(0, 0), 1, &(*Wmbje[mj])(0, 0), npno_mj);

            auto tau_temp = linalg::triplet(S_mn_jn, tau_bar[jn], S_jn_jm, false, false, false);
            K_mn_temp = linalg::doublet(K_iajb_[mn], S_mn_mj, false, false);
            Wmbje[mj]->add(linalg::doublet(tau_temp, K_mn_temp, true, false));
        }
    }

    timer_off("Compute Wmbje");

    return Wmbje;
}

void DLPNOCCSD::lccsd_iterations() {

    int n_lmo_pairs = ij_to_i_j_.size();
    int naocc = nalpha_ - nfrzc();

    outfile->Printf("\n  ==> Local CCSD <==\n\n");
    outfile->Printf("    E_CONVERGENCE = %.2e\n", options_.get_double("E_CONVERGENCE"));
    outfile->Printf("    R_CONVERGENCE = %.2e\n\n", options_.get_double("R_CONVERGENCE"));
    outfile->Printf("                     Corr. Energy    Delta E     Max R1     Max R2\n");

    // => Initialize Singles Amplitudes <= //

    T_ia_.resize(naocc);
#pragma omp parallel for
    for (int i = 0; i < naocc; ++i) {
        int ii = i_j_to_ij_[i][i];
        T_ia_[i] = std::make_shared<Matrix>(n_pno_[ii], 1);
        T_ia_[i]->zero();
    }

    // => Initialize Dressed Doubles Amplitudes <= //

    std::vector<SharedMatrix> tau(n_lmo_pairs);
    std::vector<SharedMatrix> tau_tilde(n_lmo_pairs);
    std::vector<SharedMatrix> tau_bar(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        if (n_pno_[ij] == 0) continue;

        tau[ij] = T_iajb_[ij]->clone();
        tau_tilde[ij] = T_iajb_[ij]->clone();
        tau_bar[ij] = T_iajb_[ij]->clone();
        tau_bar[ij]->scale(0.5);
    }

    // => Initialize Residuals <= //

    std::vector<SharedMatrix> R_ia(naocc);
    std::vector<SharedMatrix> Rn_iajb(n_lmo_pairs);
    std::vector<SharedMatrix> R_iajb(n_lmo_pairs);

    int iteration = 0, max_iteration = options_.get_int("DLPNO_MAXITER");
    double e_curr = 0.0, e_prev = 0.0, r1_curr = 0.0, r2_curr = 0.0;
    bool e_converged = false, r_converged = false;

    DIISManager diis1(options_.get_int("DIIS_MAX_VECS"), "LCCSD DIIS (T1)", DIISManager::RemovalPolicy::LargestError, DIISManager::StoragePolicy::InCore);
    DIISManager diis2(options_.get_int("DIIS_MAX_VECS"), "LCCSD DIIS (T2)", DIISManager::RemovalPolicy::LargestError, DIISManager::StoragePolicy::InCore);

    while (!(e_converged && r_converged)) {
        // RMS of residual per single LMO, for assesing convergence
        std::vector<double> R_ia_rms(naocc, 0.0);
        // RMS of residual per LMO pair, for assessing convergence
        std::vector<double> R_iajb_rms(n_lmo_pairs, 0.0);

        // Build one-particle intermediates
        auto Fmi = compute_Fmi(tau_tilde);
        auto Fbe = compute_Fbe(tau_tilde);
        auto Fme = compute_Fme();

        // Build two-particle W intermediates
        auto Wmnij = compute_Wmnij(tau);
        auto Wmbej = compute_Wmbej(tau_bar);
        auto Wmbje = compute_Wmbje(tau_bar);

        // Calculate singles residuals from current amplitudes
#pragma omp parallel for
        for (int i = 0; i < naocc; i++) {
            int ii = i_j_to_ij_[i][i];
            int npno_ii = n_pno_[ii];

            // Madriaga Eq. 34, Term 2
            R_ia[i] = linalg::doublet(Fbe[ii], T_ia_[i], false, false);
            /*
            for (int a_ii = 0; a_ii < npno_ii; a_ii++) {
                (*R_ia[i])(a_ii, 0) += e_pno_[ii]->get(a_ii) * (*T_ia_[i])(a_ii, 0);
            }
            */

            for (int m = 0; m < naocc; m++) {
                int im = i_j_to_ij_[i][m], mm = i_j_to_ij_[m][m], mi = i_j_to_ij_[m][i];
                int npno_mm = n_pno_[mm], npno_mi = n_pno_[mi];

                if (im == -1 || n_pno_[im] == 0) continue;

                // Madriaga Eq. 34, Term 5
                auto S_mm_im = S_PNO(mm, im);
                auto temp_t1 = linalg::doublet(S_mm_im, T_ia_[m], true, false);
                auto S_im_ii = S_PNO(im, ii);
                R_ia[i]->add(linalg::triplet(S_im_ii, M_iajb_[im], temp_t1, true, false, false));

                // Madriaga Eq. 34, Term 3
                auto S_mm_ii = S_PNO(mm, ii);
                auto T_m_temp = linalg::doublet(S_mm_ii, T_ia_[m], true, false);
                T_m_temp->scale(Fmi->get(m,i));
                R_ia[i]->subtract(T_m_temp);

                // Madriaga Eq. 34, Term 4
                int m_im = lmopair_to_lmos_dense_[im][m];
                std::vector<int> m_im_slice(1, m_im);
                auto Fe_im = submatrix_rows(*Fme[im], m_im_slice);
                R_ia[i]->add(linalg::triplet(S_im_ii, Tt_iajb_[im], Fe_im, true, false, true));

                // Madriaga Eq. 34, Term 6
                auto T_mi_mm = Tt_iajb_[mi]->clone();
                T_mi_mm->reshape(npno_mi * npno_mi, 1);
                // (npno_ii, npno_mm) (npno_mm, npno_mm * npno_mm) (npno_mm * npno_mm, 1)
                R_ia[i]->add(linalg::triplet(S_im_ii, K_tilde_phys_[mi], T_mi_mm, true, false, false));

                // Madriaga Eq. 34, Term 7
                for (int n_im = 0; n_im < lmopair_to_lmos_[im].size(); n_im++) {
                    int n = lmopair_to_lmos_[im][n_im];
                    int mn = i_j_to_ij_[m][n], nm = i_j_to_ij_[n][m];
                    if (n_pno_[mn] == 0) continue;

                    int i_mn = lmopair_to_lmos_dense_[mn][i];
                    std::vector<int> i_mn_slice(1, i_mn);
                    auto L_temp = submatrix_rows(*L_bar_[mn], i_mn_slice);

                    auto S_ii_mn = S_PNO(ii, mn);
                    R_ia[i]->subtract(linalg::triplet(S_ii_mn, T_iajb_[mn], L_temp, false, false, true));
                }
            }
            R_ia_rms[i] = R_ia[i]->rms();
        }

        // Calculate residuals from current amplitudes
#pragma omp parallel for schedule(dynamic, 1)
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            int i, j;
            std::tie(i, j) = ij_to_i_j_[ij];
            int ii = i_j_to_ij_[i][i], jj = i_j_to_ij_[j][j], ji = ij_to_ji_[ij];
            int npno_ij = n_pno_[ij], npno_ii = n_pno_[ii], npno_jj = n_pno_[jj];

            if (npno_ij == 0) continue;

            // Buffers for R2 (Save memory)
            SharedMatrix r2_temp;

            // Madriaga Eq. 35, Term 1
            Rn_iajb[ij] = K_iajb_[ij]->clone();
            Rn_iajb[ij]->scale(0.5);

            // Madriaga Eq. 35, Term 2a
            Rn_iajb[ij]->add(linalg::doublet(T_iajb_[ij], Fbe[ij], false, true));
            /*
            for (int a_ij = 0; a_ij < npno_ij; a_ij++) {
                for (int b_ij = 0; b_ij < npno_ij; b_ij++) {
                    (*Rn_iajb[ij])(a_ij, b_ij) += e_pno_[ij]->get(b_ij) * (*T_iajb_[ij])(a_ij, b_ij);
                }
            }
            */

            // Madriaga Eq. 35, Term 5
            if (virtual_storage_ == CORE) { // High Memory Algorithm
                int nsvd_ij = (i > j) ? Qab_ij_[ji].size() : Qab_ij_[ij].size();
                for (int q_ij = 0; q_ij < nsvd_ij; q_ij++) {
                    if (i > j) r2_temp = linalg::triplet(Qab_ij_[ji][q_ij], tau[ij], Qab_ij_[ji][q_ij]);
                    else r2_temp = linalg::triplet(Qab_ij_[ij][q_ij], tau[ij], Qab_ij_[ij][q_ij]);
                    r2_temp->scale(0.5);
                    Rn_iajb[ij]->add(r2_temp);
                }
            } else { // Low Memory Algorithm
                int naux_ij = lmopair_to_ribfs_[ij].size();
                auto q_ab = std::make_shared<Matrix>(naux_ij, npno_ij * npno_ij);
                for (int q_ij = 0; q_ij < naux_ij; q_ij++) {
                    int q = lmopair_to_ribfs_[ij][q_ij];
                    int centerq = ribasis_->function_to_center(q);

                    auto ab_temp = submatrix_rows_and_cols(*qab_[q], lmopair_pao_to_riatom_pao_[ij][q_ij], 
                                                        lmopair_pao_to_riatom_pao_[ij][q_ij]);
                    ab_temp = linalg::triplet(X_pno_[ij], ab_temp, X_pno_[ij], true, false, false);
                    C_DCOPY(npno_ij * npno_ij, &(*ab_temp)(0,0), 1, &(*q_ab)(q_ij, 0), 1);
                }
                auto A_solve = submatrix_rows_and_cols(*full_metric_, lmopair_to_ribfs_[ij], lmopair_to_ribfs_[ij]);
                A_solve->power(0.5, 1.0e-14);
                C_DGESV_wrapper(A_solve, q_ab);

                for (int q_ij = 0; q_ij < lmopair_to_ribfs_[ij].size(); q_ij++) {
                    std::vector<int> q_ij_slice(1, q_ij);
                    r2_temp = submatrix_rows(*q_ab, q_ij_slice);
                    r2_temp->reshape(npno_ij, npno_ij);
                    r2_temp = linalg::triplet(r2_temp, tau[ij], r2_temp);
                    r2_temp->scale(0.5);
                    Rn_iajb[ij]->add(r2_temp);
                }
                q_ab = nullptr;
            }

            // Madriaga Eq. 35, Term 12
            auto S_ij_ii = S_PNO(ij, ii);
            r2_temp = linalg::doublet(S_ij_ii, T_ia_[i], false, false);
            r2_temp = linalg::doublet(r2_temp, K_tilde_phys_[ji], true, false);
            r2_temp->reshape(npno_ij, npno_ij);
            auto eye = std::make_shared<Matrix>(npno_ij, npno_ij);
            eye->identity();
            Rn_iajb[ij]->add(linalg::doublet(r2_temp, eye, true, false));

            for (int m_ij = 0; m_ij < lmopair_to_lmos_[ij].size(); m_ij++) {
                int m = lmopair_to_lmos_[ij][m_ij];
                int im = i_j_to_ij_[i][m], mm = i_j_to_ij_[m][m], mj = i_j_to_ij_[m][j];
                int mi = ij_to_ji_[im];
                int npno_mj = n_pno_[mj];

                // Shared Intermediates
                auto S_ij_im = S_PNO(ij, im);
                auto S_ij_mj = S_PNO(ij, mj);
                auto S_ij_mm = S_PNO(ij, mm);
                auto temp_t1 = linalg::doublet(S_ij_mm, T_ia_[m], false, false);

                // Madriaga Eq. 35, Term 2b
                std::vector<int> m_ij_slice(1, m_ij);
                r2_temp = submatrix_rows(*Fme[ij], m_ij_slice);
                r2_temp = linalg::doublet(T_iajb_[ij], r2_temp, false, true);
                C_DGER(npno_ij, npno_ij, -0.5, &(*r2_temp)(0,0), 1, &(*temp_t1)(0,0), 1, &(*Rn_iajb[ij])(0,0), npno_ij);
                
                if (n_pno_[mj] != 0) {
                    // Madriaga Eq. 35, Term 6 (Zmbij term)
                    r2_temp = linalg::triplet(S_ij_mj, tau[ij], S_ij_mj, true, false, false);
                    r2_temp->reshape(npno_mj * npno_mj, 1);
                    r2_temp = linalg::doublet(K_tilde_phys_[mj], r2_temp, false, false);
                    r2_temp = linalg::doublet(S_ij_mj, r2_temp, false, false);
                    C_DGER(npno_ij, npno_ij, -1.0, &(*temp_t1)(0,0), 1, &(*r2_temp)(0,0), 1, &(*Rn_iajb[ij])(0,0), npno_ij);

                    // Madriaga Eq. 35, Term 10
                    auto S_ii_mj = S_PNO(ii, mj);
                    auto T_i_temp = linalg::doublet(S_ii_mj, T_ia_[i], true, false);
                    r2_temp = linalg::triplet(T_i_temp, K_iajb_[mj], S_ij_mj, true, false, true);
                    C_DGER(npno_ij, npno_ij, -1.0, &(*temp_t1)(0,0), 1, &(*r2_temp)(0,0), 1, &(*Rn_iajb[ij])(0,0), npno_ij);

                    // Madriaga Eq. 35, Term 11
                    r2_temp = linalg::triplet(S_ij_mj, J_ijab_[mj], T_i_temp, false, false, false);
                    C_DGER(npno_ij, npno_ij, -1.0, &(*r2_temp)(0,0), 1, &(*temp_t1)(0,0), 1, &(*Rn_iajb[ij])(0,0), npno_ij);
                }

                if (n_pno_[im] != 0) {
                    // Madriaga Eq. 35, Term 3
                    auto S_im_ij = S_PNO(im, ij);
                    r2_temp = linalg::triplet(S_im_ij, T_iajb_[im], S_im_ij, true, false, false);
                    int m_jj = lmopair_to_lmos_dense_[jj][m];
                    std::vector<int> m_jj_slice(1, m_jj);
                    double tf_dot = T_ia_[j]->vector_dot(submatrix_rows(*Fme[jj], m_jj_slice)->transpose());
                    r2_temp->scale((*Fmi)(m,j) + 0.5 * tf_dot);
                    Rn_iajb[ij]->subtract(r2_temp);
                }

                // Madriaga Eq. 35, Term 13
                r2_temp = submatrix_rows(*K_bar_[ij], m_ij_slice)->transpose();
                C_DGER(npno_ij, npno_ij, -1.0, &(*temp_t1)(0,0), 1, &(*r2_temp)(0,0), 1, &(*Rn_iajb[ij])(0,0), npno_ij);

                if (n_pno_[mi] != 0 && n_pno_[mj] != 0) {
                    auto S_im_ij = S_PNO(im, ij);
                    auto S_mj_mi = S_PNO(mj, mi);
                    auto Wmbej_mj = linalg::triplet(S_ij_mj, Wmbej[mj], S_mj_mi, false, false, false);
                    auto Wmbje_mj = linalg::triplet(S_ij_mj, Wmbje[mj], S_mj_mi, false, false, false);
                    auto Wmbje_mi = linalg::triplet(S_ij_im, Wmbje[mi], S_mj_mi, false, false, true);

                    // Madriaga Eq. 35, Term 7
                    r2_temp = T_iajb_[im]->clone();
                    r2_temp->subtract(T_iajb_[im]->transpose());
                    Rn_iajb[ij]->add(linalg::triplet(S_im_ij, r2_temp, Wmbej_mj, true, false, true));

                    // Madriaga Eq. 35, Term 8
                    r2_temp = Wmbej_mj->clone();
                    r2_temp->add(Wmbje_mj);
                    Rn_iajb[ij]->add(linalg::triplet(S_im_ij, T_iajb_[im], r2_temp, true, false, true));

                    // Madriaga Eq. 35, Term 9
                    Rn_iajb[ij]->add(linalg::triplet(S_ij_mj, T_iajb_[mj], Wmbje_mi, false, false, true));
                }

                // Madriaga Eq. 35, Term 4
                for (int n_ij = 0; n_ij < lmopair_to_lmos_[ij].size(); n_ij++) {
                    int n = lmopair_to_lmos_[ij][n_ij];
                    int mn = i_j_to_ij_[m][n];

                    if (mn == -1 || n_pno_[mn] == 0) continue;

                    auto S_ij_mn = S_PNO(ij, mn);

                    int i_mn = lmopair_to_lmos_dense_[mn][i], j_mn = lmopair_to_lmos_dense_[mn][j];
                    r2_temp = linalg::triplet(S_ij_mn, tau[mn], S_ij_mn, false, false, true);
                    r2_temp->scale(0.5 * Wmnij[mn]->get(i_mn, j_mn));
                    Rn_iajb[ij]->add(r2_temp);
                }
            }
        }

        // Compute Doubles Residual from Non-Symmetrized Doubles Residual
#pragma omp parallel for schedule(dynamic, 1)
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            int ji = ij_to_ji_[ij];
            R_iajb[ij] = std::make_shared<Matrix>("R_iajb", n_pno_[ij], n_pno_[ij]);

            if (n_pno_[ij] == 0) continue;

            R_iajb[ij] = Rn_iajb[ij]->clone();
            R_iajb[ij]->add(Rn_iajb[ji]->transpose());

            R_iajb_rms[ij] = R_iajb[ij]->rms();
        }

        // Update Singles Amplitude
#pragma omp parallel for
        for (int i = 0; i < naocc; ++i) {
            int ii = i_j_to_ij_[i][i];
            for (int a_ij = 0; a_ij < n_pno_[ii]; ++a_ij) {
                (*T_ia_[i])(a_ij, 0) -= (*R_ia[i])(a_ij, 0) / (e_pno_[ii]->get(a_ij) - F_lmo_->get(i,i));
            }
        }

        // Update Doubles Amplitude
#pragma omp parallel for schedule(dynamic, 1)
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            int i, j;
            std::tie(i, j) = ij_to_i_j_[ij];
            int ii = i_j_to_ij_[i][i], jj = i_j_to_ij_[j][j];

            if (n_pno_[ij] == 0) continue;

            for (int a_ij = 0; a_ij < n_pno_[ij]; ++a_ij) {
                for (int b_ij = 0; b_ij < n_pno_[ij]; ++b_ij) {
                    (*T_iajb_[ij])(a_ij, b_ij) -= (*R_iajb[ij])(a_ij, b_ij) / 
                                (e_pno_[ij]->get(a_ij) + e_pno_[ij]->get(b_ij) - F_lmo_->get(i,i) - F_lmo_->get(j,j));
                }
            }
        }

        // DIIS Extrapolation
        auto T_ia_flat = flatten_mats(T_ia_);
        auto R_ia_flat = flatten_mats(R_ia);

        auto T_iajb_flat = flatten_mats(T_iajb_);
        auto R_iajb_flat = flatten_mats(R_iajb);

        if (iteration == 0) {
            diis1.set_error_vector_size(R_ia_flat.get());
            diis1.set_vector_size(T_ia_flat.get());

            diis2.set_error_vector_size(R_iajb_flat.get());
            diis2.set_vector_size(T_iajb_flat.get());
        }

        diis1.add_entry(R_ia_flat.get(), T_ia_flat.get());
        diis1.extrapolate(T_ia_flat.get());

        diis2.add_entry(R_iajb_flat.get(), T_iajb_flat.get());
        diis2.extrapolate(T_iajb_flat.get());

        copy_flat_mats(T_ia_flat, T_ia_);
        copy_flat_mats(T_iajb_flat, T_iajb_);

        // Update Special Doubles Amplitudes
#pragma omp parallel for schedule(dynamic, 1)
        for (int ij = 0; ij < n_lmo_pairs; ij++) {
            int i, j;
            std::tie(i, j) = ij_to_i_j_[ij];
            int ii = i_j_to_ij_[i][i], jj = i_j_to_ij_[j][j];

            if (n_pno_[ij] == 0) continue;

            Tt_iajb_[ij] = T_iajb_[ij]->clone();
            Tt_iajb_[ij]->scale(2.0);
            Tt_iajb_[ij]->subtract(T_iajb_[ij]->transpose());

            auto S_ii_ij = S_PNO(ii, ij);
            auto S_jj_ij = S_PNO(jj, ij);
            auto tia_temp = linalg::doublet(S_ii_ij, T_ia_[i], true, false);
            auto tjb_temp = linalg::doublet(S_jj_ij, T_ia_[j], true, false);

            for (int a_ij = 0; a_ij < n_pno_[ij]; ++a_ij) {
                for (int b_ij = 0; b_ij < n_pno_[ij]; ++b_ij) {
                    double t1_cont = tia_temp->get(a_ij, 0) * tjb_temp->get(b_ij, 0);
                    double t2_cont = T_iajb_[ij]->get(a_ij, b_ij);

                    tau[ij]->set(a_ij, b_ij, t2_cont + t1_cont);
                    tau_tilde[ij]->set(a_ij, b_ij, t2_cont + 0.5 * t1_cont);
                    tau_bar[ij]->set(a_ij, b_ij, 0.5 * t2_cont + t1_cont);
                }
            }
        }

        // evaluate convergence using current amplitudes and residuals
        e_prev = e_curr;
        // Compute LCCSD energy
        e_curr = 0.0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : e_curr)
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            if (n_pno_[ij] == 0) continue;

            e_curr += tau[ij]->vector_dot(L_iajb_[ij]);
        }
        double r_curr1 = *max_element(R_ia_rms.begin(), R_ia_rms.end());
        double r_curr2 = *max_element(R_iajb_rms.begin(), R_iajb_rms.end());

        r_converged = (fabs(r_curr1) < options_.get_double("R_CONVERGENCE"));
        r_converged &= (fabs(r_curr2) < options_.get_double("R_CONVERGENCE"));
        e_converged = (fabs(e_curr - e_prev) < options_.get_double("E_CONVERGENCE"));

        outfile->Printf("  @LCCSD iter %3d: %16.12f %10.3e %10.3e %10.3e\n", iteration, e_curr, e_curr - e_prev, r_curr1, r_curr2);

        iteration++;

        if (iteration > max_iteration) {
            throw PSIEXCEPTION("Maximum DLPNO iterations exceeded.");
        }
    }

    e_lccsd_ = e_curr;
}

double DLPNOCCSD::compute_energy() {

    timer_on("DLPNO-CCSD");

    print_header();

    timer_on("Setup Orbitals");
    setup_orbitals();
    timer_off("Setup Orbitals");

    timer_on("Overlap Ints");
    compute_overlap_ints();
    timer_off("Overlap Ints");

    timer_on("Dipole Ints");
    compute_dipole_ints();
    timer_off("Dipole Ints");

    // Adjust parameters for "crude" prescreening
    T_CUT_MKN_ *= 100;
    T_CUT_DO_ *= 2;

    timer_on("Sparsity");
    prep_sparsity();
    timer_off("Sparsity");

    timer_on("Determine Strong/Weak Pairs");
    compute_qia();
    compute_metric();
    ccsd_pair_prescreening();
    timer_off("Determine Strong/Weak Pairs");

    timer_on("DF Ints");
    print_integral_sparsity();
    compute_qij();
    compute_qia();
    compute_qab();
    timer_off("DF Ints");

    timer_on("PNO Transform");
    pno_transform();
    timer_off("PNO Transform");

    timer_on("PNO Overlaps");
    compute_pno_overlaps();
    timer_off("PNO Overlaps");

    timer_on("CC Integrals");
    estimate_memory();
    compute_cc_integrals();
    qij_.clear();
    qia_.clear();
    if (virtual_storage_ != DIRECT) {
        qab_.clear();
        qab_svd_.clear();
    }
    timer_off("CC Integrals");

    timer_on("LCCSD");
    lccsd_iterations();
    timer_off("LCCSD");

    print_results();

    timer_off("DLPNO-CCSD");

    double e_scf = reference_wavefunction_->energy();
    double e_ccsd_corr = e_lccsd_ + de_lmp2_ + de_dipole_ + de_pno_total_;
    double e_ccsd_total = e_scf + e_ccsd_corr;

    set_scalar_variable("CCSD CORRELATION ENERGY", e_ccsd_corr);
    set_scalar_variable("CURRENT CORRELATION ENERGY", e_ccsd_corr);
    set_scalar_variable("CCSD TOTAL ENERGY", e_ccsd_total);
    set_scalar_variable("CURRENT ENERGY", e_ccsd_total);

    return e_ccsd_total;
}

void DLPNOCCSD::print_integral_sparsity() {
    // statistics for number of (MN|K) shell triplets we need to compute

    int nbf = basisset_->nbf();
    int nshell = basisset_->nshell();
    int naux = ribasis_->nbf();
    int naocc = nalpha_ - nfrzc();

    size_t triplets = 0;          // computed (MN|K) triplets with no screening
    size_t triplets_lmo = 0;      // computed (MN|K) triplets with only LMO screening
    size_t triplets_pao = 0;      // computed (MN|K) triplets with only PAO screening
    size_t triplets_lmo_lmo = 0;  // computed (MN|K) triplets with LMO and LMO screening
    size_t triplets_lmo_pao = 0;  // computed (MN|K) triplets with LMO and PAO screening
    size_t triplets_pao_pao = 0;  // computed (MN|K) triplets with PAO and PAO screening

    for (size_t atom = 0; atom < riatom_to_shells1_.size(); atom++) {
        size_t nshellri_atom = atom_to_rishell_[atom].size();
        triplets += nshell * nshell * nshellri_atom;
        triplets_lmo += riatom_to_shells1_[atom].size() * nshell * nshellri_atom;
        triplets_pao += nshell * riatom_to_shells2_[atom].size() * nshellri_atom;
        triplets_lmo_lmo += riatom_to_shells1_[atom].size() * riatom_to_shells1_[atom].size() * nshellri_atom;
        triplets_lmo_pao += riatom_to_shells1_[atom].size() * riatom_to_shells2_[atom].size() * nshellri_atom;
        triplets_pao_pao += riatom_to_shells2_[atom].size() * riatom_to_shells2_[atom].size() * nshellri_atom;
    }
    size_t screened_total = 3 * triplets - triplets_lmo_lmo - triplets_lmo_pao - triplets_pao_pao;
    size_t screened_lmo = triplets - triplets_lmo;
    size_t screened_pao = triplets - triplets_pao;

    // statistics for the number of (iu|Q) integrals we're left with after the transformation

    size_t total_integrals = (size_t)naocc * nbf * naux + naocc * naocc * naux + nbf * nbf * naux;
    size_t actual_integrals = 0;

    qij_memory_ = 0;
    qia_memory_ = 0;
    qab_memory_ = 0;

    for (size_t atom = 0; atom < riatom_to_shells1_.size(); atom++) {
        qij_memory_ +=
            riatom_to_lmos_ext_[atom].size() * riatom_to_lmos_ext_[atom].size() * atom_to_ribf_[atom].size();
        qia_memory_ +=
            riatom_to_lmos_ext_[atom].size() * riatom_to_paos_ext_[atom].size() * atom_to_ribf_[atom].size();
        qab_memory_ +=
            riatom_to_paos_ext_[atom].size() * riatom_to_paos_ext_[atom].size() * atom_to_ribf_[atom].size();
    }

    actual_integrals = qij_memory_ + qia_memory_ + qab_memory_;

    // number of doubles * (2^3 bytes / double) * (1 GiB / 2^30 bytes)
    double total_memory = total_integrals * pow(2.0, -27);
    double actual_memory = actual_integrals * pow(2.0, -27);
    double screened_memory = total_memory - actual_memory;

    outfile->Printf("\n");
    outfile->Printf("    Coefficient sparsity in AO -> LMO transform: %6.2f %% \n", screened_lmo * 100.0 / triplets);
    outfile->Printf("    Coefficient sparsity in AO -> PAO transform: %6.2f %% \n", screened_pao * 100.0 / triplets);
    outfile->Printf("    Coefficient sparsity in combined transforms: %6.2f %% \n", screened_total * 100.0 / (3.0 * triplets));
    outfile->Printf("\n");
    outfile->Printf("    Storing transformed LMO/LMO, LMO/PAO, and PAO/PAO integrals in sparse format.\n");
    outfile->Printf("    Required memory: %.3f GiB (%.2f %% reduction from dense format) \n", actual_memory,
                    screened_memory * 100.0 / total_memory);
}

void DLPNOCCSD::print_header() {
    outfile->Printf("   --------------------------------------------\n");
    outfile->Printf("                    DLPNO-CCSD                 \n");
    outfile->Printf("                   by Andy Jiang               \n");
    outfile->Printf("   --------------------------------------------\n\n");
    outfile->Printf("  DLPNO convergence set to %s.\n\n", options_.get_str("PNO_CONVERGENCE").c_str());
    outfile->Printf("  Detailed DLPNO thresholds and cutoffs:\n");
    outfile->Printf("    T_CUT_DO     = %6.3e \n", T_CUT_DO_);
    outfile->Printf("    T_CUT_PNO    = %6.3e \n", T_CUT_PNO_);
    outfile->Printf("    T_CUT_PAIRS  = %6.3e \n", T_CUT_PAIRS_);
    outfile->Printf("    T_CUT_MKN    = %6.3e \n", T_CUT_MKN_);
    outfile->Printf("    T_CUT_EIG    = %6.3e \n", T_CUT_EIG_);
    outfile->Printf("    T_CUT_SVD    = %6.3e \n", T_CUT_SVD_);
    outfile->Printf("    DIAG_SCALE   = %6.3e \n", T_CUT_PNO_DIAG_SCALE_);
    outfile->Printf("    T_CUT_DO_ij  = %6.3e \n", options_.get_double("T_CUT_DO_ij"));
    outfile->Printf("    T_CUT_PRE    = %6.3e \n", options_.get_double("T_CUT_PRE"));
    outfile->Printf("    T_CUT_DO_PRE = %6.3e \n", options_.get_double("T_CUT_DO_PRE"));
    outfile->Printf("    T_CUT_CLMO   = %6.3e \n", options_.get_double("T_CUT_CLMO"));
    outfile->Printf("    T_CUT_CPAO   = %6.3e \n", options_.get_double("T_CUT_CPAO"));
    outfile->Printf("    S_CUT        = %6.3e \n", options_.get_double("S_CUT"));
    outfile->Printf("    F_CUT        = %6.3e \n", options_.get_double("F_CUT"));
    outfile->Printf("\n");
}

void DLPNOCCSD::print_results() {
    outfile->Printf("  \n");
    outfile->Printf("  Total DLPNO-CCSD Correlation Energy: %16.12f \n", e_lccsd_ + de_lmp2_ + de_pno_total_ + de_dipole_);
    outfile->Printf("    CCSD Correlation Energy:           %16.12f \n", e_lccsd_);
    outfile->Printf("    Weak Pair MP2 Correction:          %16.12f \n", de_lmp2_);
    outfile->Printf("    LMO Truncation Correction:         %16.12f \n", de_dipole_);
    outfile->Printf("    PNO Truncation Correction:         %16.12f \n", de_pno_total_);
    outfile->Printf("\n\n  @Total DLPNO-CCSD Energy: %16.12f \n", variables_["SCF TOTAL ENERGY"] + e_lccsd_ + de_lmp2_ + de_pno_total_ + de_dipole_);
}

}  // namespace dlpno
}  // namespace psi