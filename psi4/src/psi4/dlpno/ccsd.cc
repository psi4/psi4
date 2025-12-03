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

#include <ctime>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi {
namespace dlpno {

DLPNOCCSD::DLPNOCCSD(SharedWavefunction ref_wfn, Options& options) : DLPNO(ref_wfn, options) {}
DLPNOCCSD::~DLPNOCCSD() {}

inline SharedMatrix DLPNOCCSD::S_PNO(const int ij, const int mn) {
    int i, j, m, n;
    std::tie(i, j) = ij_to_i_j_[ij];
    std::tie(m, n) = ij_to_i_j_[mn];

    int ji = ij_to_ji_[ij];
    int nm = ij_to_ji_[mn];

    if (i == m) { // S(ij, mn) -> S(ij, in) -> S(ji, ni)
        return S_pno_ij_kj_[ji][n];
    } else if (i == n) { // S(ij, mn) -> S(ij, mi) -> S(ji, mi)
        return S_pno_ij_kj_[ji][m];
    } else if (j == m) { // S(ij, mn) -> S(ij, jn) -> S(ij, nj)
        return S_pno_ij_kj_[ij][n];
    } else if (j == n) { // S(ij, mn) -> S(ij, mj)
        return S_pno_ij_kj_[ij][m];
    } else if (m == n) {
        const int ij_idx = (i > j) ? ji : ij;
        return S_pno_ij_nn_[ij_idx][m];
    } else {
        int i, j, m, n;
        std::tie(i, j) = ij_to_i_j_[ij];
        std::tie(m, n) = ij_to_i_j_[mn];

        const int m_ij = lmopair_to_lmos_dense_[ij][m], n_ij = lmopair_to_lmos_dense_[ij][n];
        if (m_ij == -1 || n_ij == -1 || low_memory_overlap_) {
            auto S_ij_mn = submatrix_rows_and_cols(*S_pao_, lmopair_to_paos_[ij], lmopair_to_paos_[mn]);
            return linalg::triplet(X_pno_[ij], S_ij_mn, X_pno_[mn], true, false, false);
        }
        
        const int nlmo_ij = lmopair_to_lmos_[ij].size();

        int mn_ij; 
        if (m_ij > n_ij) {
            mn_ij = n_ij * nlmo_ij + m_ij;
        } else {
            mn_ij = m_ij * nlmo_ij + n_ij;
        }
        
        if (i > j) {
            const int ji = ij_to_ji_[ij];
            return S_pno_ij_mn_[ji][mn_ij];
        } else {
            return S_pno_ij_mn_[ij][mn_ij];
        }
    }
}

inline std::vector<SharedMatrix> DLPNOCCSD::QIA_PNO(const int ij) {
    auto &[i, j] = ij_to_i_j_[ij];
    int pair_idx = (i > j) ? ij_to_ji_[ij] : ij;

    if (write_qia_pno_) {
        int naux_ij = lmopair_to_ribfs_[ij].size();
        int nlmo_ij = lmopair_to_lmos_[ij].size();
        int npno_ij = n_pno_[ij];

        std::stringstream toc_entry;
        toc_entry << "QIA (PNO) " << pair_idx;

        auto q_ov = std::make_shared<Matrix>(toc_entry.str(), naux_ij, nlmo_ij * npno_ij);
#pragma omp critical
        q_ov->load(psio_, PSIF_DLPNO_QIA_PNO, psi::Matrix::SubBlocks);

        std::vector<SharedMatrix> Qma_pno(naux_ij);
        for (int q_ij = 0; q_ij < naux_ij; q_ij++) {
            Qma_pno[q_ij] = std::make_shared<Matrix>(nlmo_ij, npno_ij);
            ::memcpy(&(*Qma_pno[q_ij])(0, 0), &(*q_ov)(q_ij, 0), nlmo_ij * npno_ij * sizeof(double));
        }
        return Qma_pno;
    } else {
        return Qma_ij_[pair_idx];
    }
}

inline std::vector<SharedMatrix> DLPNOCCSD::QAB_PNO(const int ij) {
    auto &[i, j] = ij_to_i_j_[ij];
    int pair_idx = (i > j) ? ij_to_ji_[ij] : ij;

    if (write_qab_pno_) {
        int naux_ij = lmopair_to_ribfs_[ij].size();
        int npno_ij = n_pno_[ij];

        std::stringstream toc_entry;
        toc_entry << "QAB (PNO) " << pair_idx;

        auto q_vv = std::make_shared<Matrix>(toc_entry.str(), naux_ij, npno_ij * npno_ij);
#pragma omp critical
        q_vv->load(psio_, PSIF_DLPNO_QAB_PNO, psi::Matrix::ThreeIndexLowerTriangle);

        std::vector<SharedMatrix> Qab_pno(naux_ij);
        for (int q_ij = 0; q_ij < naux_ij; q_ij++) {
            Qab_pno[q_ij] = std::make_shared<Matrix>(npno_ij, npno_ij);
            ::memcpy(&(*Qab_pno[q_ij])(0, 0), &(*q_vv)(q_ij, 0), npno_ij * npno_ij * sizeof(double));
        }
        return Qab_pno;
    } else {
        return Qab_ij_[pair_idx];
    }
}

void DLPNOCCSD::compute_pno_overlaps() {

    const int naocc = i_j_to_ij_.size();
    const int n_lmo_pairs = ij_to_i_j_.size();
    
    S_pno_ij_kj_.resize(n_lmo_pairs);
    S_pno_ij_nn_.resize(n_lmo_pairs);
    if (!low_memory_overlap_) S_pno_ij_mn_.resize(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        int i, j;
        std::tie(i, j) = ij_to_i_j_[ij];

        S_pno_ij_kj_[ij].resize(naocc);

        const int npno_ij = n_pno_[ij];
        const int nlmo_ij = lmopair_to_lmos_[ij].size();

        for (int k = 0; k < naocc; ++k) {
            int kj = i_j_to_ij_[k][j];

            if (kj == -1) continue;

            S_pno_ij_kj_[ij][k] = submatrix_rows_and_cols(*S_pao_, lmopair_to_paos_[ij], lmopair_to_paos_[kj]);
            S_pno_ij_kj_[ij][k] = linalg::triplet(X_pno_[ij], S_pno_ij_kj_[ij][k], X_pno_[kj], true, false, false);
        }

        if (i > j) continue;

        S_pno_ij_nn_[ij].resize(naocc);

        for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
            int k = lmopair_to_lmos_[ij][k_ij];
            int kk = i_j_to_ij_[k][k];

            S_pno_ij_nn_[ij][k] = submatrix_rows_and_cols(*S_pao_, lmopair_to_paos_[ij], lmopair_to_paos_[kk]);
            S_pno_ij_nn_[ij][k] = linalg::triplet(X_pno_[ij], S_pno_ij_nn_[ij][k], X_pno_[kk], true, false, false);
        }

        if (!low_memory_overlap_) {
            S_pno_ij_mn_[ij].resize(nlmo_ij * nlmo_ij);

            for (int mn_ij = 0; mn_ij < nlmo_ij * nlmo_ij; ++mn_ij) {
                const int m_ij = mn_ij / nlmo_ij, n_ij = mn_ij % nlmo_ij;
                const int m = lmopair_to_lmos_[ij][m_ij], n = lmopair_to_lmos_[ij][n_ij];
                const int mn = i_j_to_ij_[m][n];

                // If these are true, the quantity has already been formed in a previous intermediate
                if (i == m || i == n || j == m || j == n || m == n) continue;
                if (mn == -1 || m_ij > n_ij) continue;

                S_pno_ij_mn_[ij][mn_ij] = submatrix_rows_and_cols(*S_pao_, lmopair_to_paos_[ij], lmopair_to_paos_[mn]);
                S_pno_ij_mn_[ij][mn_ij] = linalg::triplet(X_pno_[ij], S_pno_ij_mn_[ij][mn_ij], X_pno_[mn], true, false, false);
            }
        }
    }
}

void DLPNOCCSD::estimate_memory() {
    outfile->Printf("  ==> DLPNO-CCSD Memory Requirements <== \n\n");

    int naocc = i_j_to_ij_.size();
    int n_lmo_pairs = ij_to_i_j_.size();

    // Estimate memory cost of storing overlap integrals
    size_t low_overlap_memory = 0, pno_overlap_memory = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : low_overlap_memory, pno_overlap_memory)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];

        const int npno_ij = n_pno_[ij];
        const int nlmo_ij = lmopair_to_lmos_[ij].size();

        // These account for the memory costs of cheaper overlap integrals S(ij, kj) and S(ij, kk)
        for (int k = 0; k < naocc; ++k) {
            int kj = i_j_to_ij_[k][j];
            int ik = i_j_to_ij_[i][k];

            if (kj != -1) {
                pno_overlap_memory += n_pno_[ij] * n_pno_[kj];
                low_overlap_memory += n_pno_[ij] * n_pno_[kj];
            }

            if (i <= j && lmopair_to_lmos_dense_[ij][k] != -1) {
                int kk = i_j_to_ij_[k][k];
                pno_overlap_memory += n_pno_[ij] * n_pno_[kk];
                low_overlap_memory += n_pno_[ij] * n_pno_[kk];
            }
        }

        // These account for the memory costs of the expensive S(ij, kl) PNO overlaps
        if (i <= j) {
            for (int mn_ij = 0; mn_ij < nlmo_ij * nlmo_ij; mn_ij++) {
                const int m_ij = mn_ij / nlmo_ij, n_ij = mn_ij % nlmo_ij;
                const int m = lmopair_to_lmos_[ij][m_ij], n = lmopair_to_lmos_[ij][n_ij];
                const int mn = i_j_to_ij_[m][n];
                // If these are true, the quantity has already been formed in a previous intermediate
                if (i == m || i == n || j == m || j == n || m == n) continue;
                if (mn == -1 || m_ij > n_ij) continue;

                pno_overlap_memory += n_pno_[ij] * n_pno_[mn];
            } // end mn_ij
        } // end if
    } // end ij

    size_t oo, ov, vv, vv_non_proj, vvv, qo, qv, qov, qvv;

    // oo => n_lmo_pairs * (nlmo_{ij}, nlmo_{ij})-like quantities: \beta_{ij}^{kl} (1 case over strong pairs, restricted indexing)

    // ov => n_lmo_pairs * (nlmo_{ij}, npno_{ij})-like quantities: K_mibj, J_ijmb, L_mibj (3 cases over all pairs)
    // F_kc (1 case over strong pairs, restricted indexing)

    // vv => n_lmo_pairs * (npno_{ij}, npno_{ij})-like quantities: K_iajb, T_iajb, Tt_iajb, L_iajb, gamma_{ki}^{ac},
    // delta_{ik}^{ac} (6 cases over all pairs)
    // F_ab_t1, F_ab (2 cases over strong pairs only, restricted indexing)

    // vv_non_proj => (i k | a_{ij} c_{kj}) and (i a_{ij} | k c_{kj}) integrals (2 cases over strong pairs)

    // vvv => (i a_{ij} | b_{ij} c_{ij}) (1 case over all pairs)

    // qo => B^{Q}_{k_{ij} i} (2 cases over strong pairs)
    // qv => B^{Q}_{a_{ij} i} (2 cases over strong pairs)
    // qov => B^{Q_{ij}}_{k_{ij} a_{ij}} (1 case over strong pairs, restricted indexing)
    // qvv => B^{Q_{ij}}_{a_{ij} b_{ij}} (1 case over strong pairs, restricted indexing)

#pragma omp parallel for schedule(dynamic) reduction(+ : oo, ov, vv, vv_non_proj, vvv, qo, qv, qov, qvv)
    for (int ij = 0; ij < n_lmo_pairs; ij++) {
        auto &[i, j] = ij_to_i_j_[ij];

        const int naux_ij = lmopair_to_ribfs_[ij].size();
        const int nlmo_ij = lmopair_to_lmos_[ij].size();
        const int npno_ij = n_pno_[ij];
        const bool is_strong_pair = i_j_to_ij_strong_[i][j] != -1;
        
        ov += 3 * nlmo_ij * npno_ij; // 3 cases over all pairs
        vv += 6 * npno_ij * npno_ij; // 6 cases over all pairs
        vvv += npno_ij * npno_ij * npno_ij; // 1 case over all pairs

        // Tensors of the other types are only computed over strong pairs
        if (is_strong_pair) {
            // Memory cost of the special non-projected J_{ik}^{ac} and K_{ik}^{ac} integrals
            for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
                int k = lmopair_to_lmos_[ij][k_ij];
                int ik = i_j_to_ij_[i][k], jk = i_j_to_ij_[j][k];
                vv_non_proj += 2 * n_pno_[ij] * n_pno_[jk]; // 2 cases
            } // end k_ij

            qo += 2 * naux_ij * nlmo_ij; // 2 cases over all strong pairs
            qv += 2 * naux_ij * npno_ij; // 2 cases over all strong pairs

            // Only over unique pairs
            if (i <= j) {
                oo += nlmo_ij * nlmo_ij; // 1 case over strong pairs, restricted indexing
                ov += nlmo_ij * npno_ij; // 1 case over strong pairs, restricted indexing
                vv += 2 * npno_ij * npno_ij; // 2 cases over strong pairs, restricted indexing
                qov += naux_ij * nlmo_ij * npno_ij; // 1 case over strong pairs, restricted indexing
                qvv += naux_ij * npno_ij * npno_ij; // 1 case over strong pairs, restricted indexing
            } // end if

        } // end if
    } // end ij

    // Estimate the amount of buffer space required
    size_t thread_buffer_a = 0, thread_buffer_b = 0;

    for (int ij = 0; ij < n_lmo_pairs; ij++) {
        const int naux_ij = lmopair_to_ribfs_[ij].size();
        const int nlmo_ij = lmopair_to_lmos_[ij].size();
        const int npao_ij = lmopair_to_paos_[ij].size();
        const int npno_ij = n_pno_[ij];

        // Determine size of extended_pao_domain
        std::vector<int> extended_pao_domain;
        
        extended_pao_domain = lmopair_to_paos_[ij];
        for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
            int k = lmopair_to_lmos_[ij][k_ij];
            extended_pao_domain = merge_lists(extended_pao_domain, lmo_to_paos_[k]);
        }
        const int npao_ext_ij = extended_pao_domain.size();

        size_t buff_a_est = naux_ij * npao_ext_ij * (nlmo_ij + npno_ij); // Buffer size used in computed projected integrals
        buff_a_est += naux_ij * npno_ij * (nlmo_ij + npno_ij); // Buffer size for storing computing qov and qvv
        thread_buffer_a = std::max(thread_buffer_a, buff_a_est);

        size_t buff_b_est = 2 * naux_ij * npno_ij * (nlmo_ij + npno_ij); // Buffer size for qov and qvv (worst case 2 copies of each)
        thread_buffer_b = std::max(thread_buffer_b, buff_b_est);
    }

    low_memory_overlap_ = options_.get_bool("LOW_MEMORY_OVERLAP");
    if (low_memory_overlap_) pno_overlap_memory = low_overlap_memory;

    write_qia_pno_ = options_.get_bool("WRITE_QIA_PNO");
    if (write_qia_pno_) qov = 0;

    write_qab_pno_ = options_.get_bool("WRITE_QAB_PNO");
    if (write_qab_pno_) qvv = 0;

// Thread and OMP Parallel info
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif

    const size_t total_df_memory = qij_memory_ + qia_memory_ + qab_memory_;
    const size_t total_pno_int_memory = oo + ov + vv + vv_non_proj + vvv + qo + qv + qov + qvv;

    size_t memory_integrals = total_df_memory + total_pno_int_memory + thread_buffer_a * nthreads;
    size_t memory_ccsd = total_df_memory + total_pno_int_memory + pno_overlap_memory + thread_buffer_b * nthreads;

    // 1 GB = 1000^3 = 10^9 Bytes
    const double DOUBLES_TO_GB = pow(10.0, -9) * sizeof(double);
    const double WORDS_TO_GB = pow(10.0, -9);

    outfile->Printf("    *** Common Quantities ***\n");
    outfile->Printf("    (q | i j) [AUX, LMO]          : %8.3f [GB]\n", qij_memory_ * DOUBLES_TO_GB);
    outfile->Printf("    (q | i a) [AUX, LMO, PAO]     : %8.3f [GB]\n", qia_memory_ * DOUBLES_TO_GB);
    outfile->Printf("    (q | a b) [AUX, PAO]          : %8.3f [GB]\n", qab_memory_ * DOUBLES_TO_GB);
    outfile->Printf("    (k_{ij}, l_{ij})-like         : %8.3f [GB]\n", oo * DOUBLES_TO_GB);
    outfile->Printf("    (k_{ij}, c_{ij})-like         : %8.3f [GB]\n", ov * DOUBLES_TO_GB);
    outfile->Printf("    (a_{ij}, b_{ij})-like         : %8.3f [GB]\n", vv * DOUBLES_TO_GB);
    outfile->Printf("    (a_{ij}, c_{kj})-like         : %8.3f [GB]\n", vv_non_proj * DOUBLES_TO_GB);
    outfile->Printf("    (i a_{ij} | b_{ij}, c_{ij})   : %8.3f [GB]\n", vvv * DOUBLES_TO_GB);
    outfile->Printf("    (Q_{ij} | k_{ij} i)           : %8.3f [GB]\n", qo * DOUBLES_TO_GB);
    outfile->Printf("    (Q_{ij} | a_{ij} i)           : %8.3f [GB]\n", qv * DOUBLES_TO_GB);
    outfile->Printf("    (Q_{ij} | m_{ij} a_{ij})      : %8.3f [GB]\n", qov * DOUBLES_TO_GB);
    outfile->Printf("    (Q_{ij} | a_{ij} b_{ij})      : %8.3f [GB]\n\n", qvv * DOUBLES_TO_GB);

    outfile->Printf("    *** Maximum ERI buffer space of %.3f [GB] per pair, over %d threads...\n\n", thread_buffer_a * DOUBLES_TO_GB, nthreads);

    outfile->Printf("    Buffer Space for ERIs         : %8.3f [GB]\n", thread_buffer_a * nthreads * DOUBLES_TO_GB);
    outfile->Printf("    Total Memory Required (ERIs)  : %8.3f [GB]\n\n", memory_integrals * DOUBLES_TO_GB);

    outfile->Printf("    *** Quantities for LCCSD Iterations ***\n");
    outfile->Printf("    PNO overlaps                  : %8.3f [GB]\n\n", pno_overlap_memory * DOUBLES_TO_GB);

    outfile->Printf("    *** Maximum CCSD buffer space of %.3f [GB] per pair, over %d threads...\n\n", thread_buffer_b * DOUBLES_TO_GB, nthreads);

    outfile->Printf("    Buffer Space for Iterations   : %8.3f [GB]\n", thread_buffer_b * nthreads * DOUBLES_TO_GB);
    outfile->Printf("    Total Memory Required (LCCSD) : %8.3f [GB]\n\n", memory_ccsd * DOUBLES_TO_GB);
    
    outfile->Printf("    Total Memory Given            : %8.3f [GB]\n\n", memory_ * WORDS_TO_GB);

    // Memory checks!!!
    bool memory_changed = false;

    if (std::max(memory_ccsd, memory_integrals) * sizeof(double) > 0.9 * memory_) {
        outfile->Printf("  Total Required Memory is more than 90%% of Available Memory!\n");
        outfile->Printf("    Attempting to switch to semi-direct low memory PNO overlap algorithm...\n");

        memory_ccsd += (low_overlap_memory - pno_overlap_memory);
        low_memory_overlap_ = true;
        memory_changed = true;
        pno_overlap_memory = low_overlap_memory;
        outfile->Printf("    Required Memory Reduced to %.3f [GB]\n\n", std::max(memory_ccsd, memory_integrals) * DOUBLES_TO_GB);
    }

    if (std::max(memory_ccsd, memory_integrals) * sizeof(double) > 0.9 * memory_) {
        outfile->Printf("  Total Required Memory is (still) more than 90%% of Available Memory!\n");
        outfile->Printf("    Attempting to switch to disk IO for (Q_{ij}|m_{ij} a_{ij}) integrals...\n");

        memory_ccsd -= qov;
        memory_integrals -= qov;
        write_qia_pno_ = true;
        memory_changed = true;
        qov = 0L;
        outfile->Printf("    Required Memory Reduced to %.3f [GB]\n\n", std::max(memory_ccsd, memory_integrals) * DOUBLES_TO_GB);
    }

    if (std::max(memory_ccsd, memory_integrals) * sizeof(double) > 0.9 * memory_) {
        outfile->Printf("  Total Required Memory is (still) more than 90%% of Available Memory!\n");
        outfile->Printf("    Attempting to switch to disk IO for (Q_{ij}|a_{ij} b_{ij}) integrals...\n");

        memory_ccsd -= qvv;
        memory_integrals -= qvv;
        write_qab_pno_ = true;
        memory_changed = true;
        qvv = 0L;
        outfile->Printf("    Required Memory Reduced to %.3f [GB]\n\n", std::max(memory_ccsd, memory_integrals) * DOUBLES_TO_GB);
    }

    if (std::max(memory_ccsd, memory_integrals) * sizeof(double) > 0.9 * memory_) {
        outfile->Printf("  Total Required Memory is (still) more than 90%% of Available Memory!\n");
        outfile->Printf("    We exhausted all of our options!!! This computation cannot continue...\n");

        throw PSIEXCEPTION("   Too little memory given for DLPNO-CCSD Algorithm!");
    }

    if (memory_changed) {
        outfile->Printf("  ==> (Updated) DLPNO-CCSD Memory Requirements <== \n\n");

        outfile->Printf("    *** Common Quantities ***\n");
        outfile->Printf("    (q | i j) [AUX, LMO]          : %8.3f [GB]\n", qij_memory_ * DOUBLES_TO_GB);
        outfile->Printf("    (q | i a) [AUX, LMO, PAO]     : %8.3f [GB]\n", qia_memory_ * DOUBLES_TO_GB);
        outfile->Printf("    (q | a b) [AUX, PAO]          : %8.3f [GB]\n", qab_memory_ * DOUBLES_TO_GB);
        outfile->Printf("    (k_{ij}, l_{ij})-like         : %8.3f [GB]\n", oo * DOUBLES_TO_GB);
        outfile->Printf("    (k_{ij}, c_{ij})-like         : %8.3f [GB]\n", ov * DOUBLES_TO_GB);
        outfile->Printf("    (a_{ij}, b_{ij})-like         : %8.3f [GB]\n", vv * DOUBLES_TO_GB);
        outfile->Printf("    (a_{ij}, c_{kj})-like         : %8.3f [GB]\n", vv_non_proj * DOUBLES_TO_GB);
        outfile->Printf("    (i a_{ij} | b_{ij}, c_{ij})   : %8.3f [GB]\n", vvv * DOUBLES_TO_GB);
        outfile->Printf("    (Q_{ij} | k_{ij} i)           : %8.3f [GB]\n", qo * DOUBLES_TO_GB);
        outfile->Printf("    (Q_{ij} | a_{ij} i)           : %8.3f [GB]\n", qv * DOUBLES_TO_GB);
        outfile->Printf("    (Q_{ij} | m_{ij} a_{ij})      : %8.3f [GB]\n", qov * DOUBLES_TO_GB);
        outfile->Printf("    (Q_{ij} | a_{ij} b_{ij})      : %8.3f [GB]\n\n", qvv * DOUBLES_TO_GB);

        outfile->Printf("    *** Maximum ERI buffer space of %.3f [GB] per pair, over %d threads...\n\n", thread_buffer_a * DOUBLES_TO_GB, nthreads);

        outfile->Printf("    Buffer Space for ERIs         : %8.3f [GB]\n", thread_buffer_a * nthreads * DOUBLES_TO_GB);
        outfile->Printf("    Total Memory Required (ERIs)  : %8.3f [GB]\n\n", memory_integrals * DOUBLES_TO_GB);

        outfile->Printf("    *** Quantities for LCCSD Iterations ***\n");
        outfile->Printf("    PNO overlaps                  : %8.3f [GB]\n\n", pno_overlap_memory * DOUBLES_TO_GB);

        outfile->Printf("    *** Maximum CCSD buffer space of %.3f [GB] per pair, over %d threads...\n\n", thread_buffer_b * DOUBLES_TO_GB, nthreads);

        outfile->Printf("    Buffer Space for Iterations   : %8.3f [GB]\n", thread_buffer_b * nthreads * DOUBLES_TO_GB);
        outfile->Printf("    Total Memory Required (LCCSD) : %8.3f [GB]\n\n", memory_ccsd * DOUBLES_TO_GB);
        
        outfile->Printf("    Total Memory Given            : %8.3f [GB]\n\n", memory_ * WORDS_TO_GB);
    }

    if (low_memory_overlap_) {
        outfile->Printf("    Using low memory PNO overlap algorithm...  \n\n");
    } else {
        outfile->Printf("    Using high memory PNO overlap algorithm... \n\n");
    }

    if (write_qia_pno_) {
        outfile->Printf("    Writing (Q_{ij}|m_{ij} a_{ij}) integrals to disk...\n\n");
    } else {
        outfile->Printf("    Storing (Q_{ij}|m_{ij} a_{ij}) integrals in RAM... \n\n");
    }

    if (write_qab_pno_) {
        outfile->Printf("    Writing (Q_{ij}|a_{ij} b_{ij}) integrals to disk...\n\n");
    } else {
        outfile->Printf("    Storing (Q_{ij}|a_{ij} b_{ij}) integrals in RAM... \n\n");
    }
}

template<bool crude> std::vector<double> DLPNOCCSD::compute_pair_energies() {
    /*
        If crude, runs semicanonical (non-iterative) MP2
        If non-crude, computes PNOs (through PNO transform), runs full iterative LMP2
    */

    int nbf = basisset_->nbf();
    int naocc = i_j_to_ij_.size();
    int n_lmo_pairs = ij_to_i_j_.size();

    const int MIN_PNOS = options_.get_int("MIN_PNOS");

    std::vector<double> e_ijs(n_lmo_pairs);

    if constexpr (crude) {
        outfile->Printf("\n  ==> Semi-Canonical MP2 Pair Prescreening <==\n\n");
    }

    if constexpr (!crude) {
        outfile->Printf("\n  ==> Forming Pair Natural Orbitals (for LMP2) <==\n");

        K_iajb_.resize(n_lmo_pairs);
        T_iajb_.resize(n_lmo_pairs);
        Tt_iajb_.resize(n_lmo_pairs);
        X_pno_.resize(n_lmo_pairs);
        e_pno_.resize(n_lmo_pairs);

        n_pno_.resize(n_lmo_pairs);
        occ_pno_.resize(n_lmo_pairs);
        trace_pno_.resize(n_lmo_pairs);
        e_ratio_pno_.resize(n_lmo_pairs);
        de_pno_.resize(n_lmo_pairs);
    }

    // Step 1: compute SC-LMP2 pair energies
#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];
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

        S_pao_ij = linalg::triplet(X_pao_ij, S_pao_ij, X_pao_ij, true, false, false);
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

        // Compute semi-canonical MP2 pair energy using PAOs
        size_t nvir_ij = K_pao_ij->rowspi(0);

        auto Tt_pao_ij = T_pao_ij->clone();
        Tt_pao_ij->scale(2.0);
        Tt_pao_ij->subtract(T_pao_ij->transpose());

        // MP2 energy of this LMO pair before transformation to PNOs
        double e_ij_initial = K_pao_ij->vector_dot(Tt_pao_ij);

        e_ijs[ij] = e_ij_initial;
        if (i < j) {
            e_ijs[ji] = e_ij_initial;
        }

        // Compute PNOs in the non-crude prescreening step :)
        if constexpr (!crude) {
            //                                           //
            // ==> Canonical PAOs  to Canonical PNOs <== //
            //                                           //

            // PNOs defined in (DOI: 10.1063/1.3086717), EQ 17 through EQ 24

            // Construct pair density from amplitudes
            auto D_ij = linalg::doublet(Tt_pao_ij, T_pao_ij, false, true);
            D_ij->add(linalg::doublet(Tt_pao_ij, T_pao_ij, true, false));
            if (i == j) D_ij->scale(0.5);

            // Diagonalization of pair density gives PNOs (in basis of the LMO's virtual domain) and PNO occ numbers
            auto X_pno_ij = std::make_shared<Matrix>("eigenvectors", nvir_ij, nvir_ij);
            Vector pno_occ("eigenvalues", nvir_ij);
            D_ij->diagonalize(*X_pno_ij, pno_occ, descending);

            double t_cut_scale = 1.0;
            if (i == j) t_cut_scale *= T_CUT_PNO_DIAG_SCALE_;
            if (i < ncore_ || j < ncore_) t_cut_scale *= T_CUT_PNO_CORE_SCALE_;

            double occ_total = 0.0;
            for (size_t a = 0; a < nvir_ij; ++a) {
                occ_total += pno_occ.get(a);
            }

            double e_pno = 0.0;
            double occ_pno = 0.0;

            int nvir_ij_final = 0;
            std::vector<int> a_curr;
            auto K_pno_init = linalg::triplet(X_pno_ij, K_pao_ij, X_pno_ij, true, false, false);
            auto Tt_pno_init = linalg::triplet(X_pno_ij, Tt_pao_ij, X_pno_ij, true, false, false);

            for (size_t a = 0; a < nvir_ij; ++a) {
                if (fabs(pno_occ.get(a)) >= t_cut_scale * T_CUT_PNO_MP2_ || occ_pno / occ_total < T_CUT_TRACE_MP2_ ||
                        std::fabs(e_pno) < T_CUT_ENERGY_MP2_ * std::fabs(e_ij_initial) || a < MIN_PNOS) {
                    // Energy criteria
                    e_pno = submatrix_rows_and_cols(*K_pno_init, a_curr, a_curr)->vector_dot(submatrix_rows_and_cols(*Tt_pno_init, a_curr, a_curr));

                    // Trace criteria
                    occ_pno += pno_occ.get(a);

                    // Add PNO to list
                    a_curr.push_back(a);
                    nvir_ij_final++;
                }
            }

            // Make sure there is at least one PNO per pair :)
            nvir_ij_final = std::max(1, nvir_ij_final);

            Dimension zero(1);
            Dimension dim_final(1);
            dim_final.fill(nvir_ij_final);

            // This transformation gives orbitals that are orthonormal but not canonical
            X_pno_ij = X_pno_ij->get_block({zero, X_pno_ij->rowspi()}, {zero, dim_final});
            pno_occ = pno_occ.get_block({zero, dim_final});

            SharedMatrix pno_canon;
            SharedVector e_pno_ij;
            std::tie(pno_canon, e_pno_ij) = canonicalizer(X_pno_ij, F_pao_ij);

            // This transformation gives orbitals that are orthonormal and canonical
            X_pno_ij = linalg::doublet(X_pno_ij, pno_canon, false, false);

            auto K_pno_ij = linalg::triplet(X_pno_ij, K_pao_ij, X_pno_ij, true, false, false);
            auto T_pno_ij = linalg::triplet(X_pno_ij, T_pao_ij, X_pno_ij, true, false, false);
            auto Tt_pno_ij = linalg::triplet(X_pno_ij, Tt_pao_ij, X_pno_ij, true, false, false);

            // mp2 energy of this LMO pair after transformation to PNOs and truncation
            double e_ij_trunc = K_pno_ij->vector_dot(Tt_pno_ij);

            // truncation error
            double de_pno_ij = e_ij_initial - e_ij_trunc;

            X_pno_ij = linalg::doublet(X_pao_ij, X_pno_ij, false, false);

            // Set values for relavant PNO-related quantities
            K_iajb_[ij] = K_pno_ij;
            T_iajb_[ij] = T_pno_ij;
            Tt_iajb_[ij] = Tt_pno_ij;
            X_pno_[ij] = X_pno_ij;
            e_pno_[ij] = e_pno_ij;
            n_pno_[ij] = X_pno_ij->colspi(0);
            occ_pno_[ij] = pno_occ.get(n_pno_[ij] - 1);
            trace_pno_[ij] = occ_pno / occ_total;
            e_ratio_pno_[ij] = e_ij_trunc / e_ij_initial;
            de_pno_[ij] = de_pno_ij;

            // account for symmetry
            if (i < j) {
                K_iajb_[ji] = K_iajb_[ij]->transpose();
                T_iajb_[ji] = T_iajb_[ij]->transpose();
                Tt_iajb_[ji] = Tt_iajb_[ij]->transpose();
                X_pno_[ji] = X_pno_[ij];
                e_pno_[ji] = e_pno_[ij];
                n_pno_[ji] = n_pno_[ij];
                occ_pno_[ji] = occ_pno_[ij];
                trace_pno_[ji] = trace_pno_[ij];
                e_ratio_pno_[ji] = e_ratio_pno_[ij];
                de_pno_[ji] = de_pno_ij;
            }
        }
    } // end for (ij pairs)

    if constexpr (!crude) {
        // Print out PNO domain information
        int pno_count_total = 0, pno_count_min = nbf, pno_count_max = 0;
        double occ_number_total = 0.0, occ_number_min = 2.0, occ_number_max = 0.0;
        double trace_total = 0.0, trace_min = 1.0, trace_max = 0.0;
        double energy_total = 0.0, energy_min = 1.0, energy_max = 0.0;
        de_pno_total_ = 0.0;
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            pno_count_total += n_pno_[ij];
            pno_count_min = std::min(pno_count_min, n_pno_[ij]);
            pno_count_max = std::max(pno_count_max, n_pno_[ij]);

            occ_number_total += occ_pno_[ij];
            occ_number_min = std::min(occ_number_min, occ_pno_[ij]);
            occ_number_max = std::max(occ_number_max, occ_pno_[ij]);

            trace_total += trace_pno_[ij];
            trace_min = std::min(trace_min, trace_pno_[ij]);
            trace_max = std::max(trace_max, trace_pno_[ij]);

            energy_total += e_ratio_pno_[ij];
            energy_min = std::min(energy_min, e_ratio_pno_[ij]);
            energy_max = std::max(energy_max, e_ratio_pno_[ij]);
            
            de_pno_total_ += de_pno_[ij];
        }

        outfile->Printf("  \n");
        outfile->Printf("    Natural Orbitals per Local MO pair:\n");
        outfile->Printf("      Avg: %3d NOs \n", pno_count_total / n_lmo_pairs);
        outfile->Printf("      Min: %3d NOs \n", pno_count_min);
        outfile->Printf("      Max: %3d NOs \n\n", pno_count_max);

        outfile->Printf("      Avg Occ Number Tol: %.3e \n", occ_number_total / n_lmo_pairs);
        outfile->Printf("      Min Occ Number Tol: %.3e \n", occ_number_min);
        outfile->Printf("      Max Occ Number Tol: %.3e \n\n", occ_number_max);

        outfile->Printf("      Avg Trace Sum: %.6f \n", trace_total / n_lmo_pairs);
        outfile->Printf("      Min Trace Sum: %.6f \n", trace_min);
        outfile->Printf("      Max Trace Sum: %.6f \n\n", trace_max);

        outfile->Printf("      Avg Energy Ratio: %.6f \n", energy_total / n_lmo_pairs);
        outfile->Printf("      Min Energy Ratio: %.6f \n", energy_min);
        outfile->Printf("      Max Energy Ratio: %.6f \n\n", energy_max);

        outfile->Printf("    PNO truncation energy = %.12f\n", de_pno_total_);
    }

    return e_ijs;
}

std::vector<double> DLPNOCCSD::pno_lmp2_iterations() {
    /* Code borrowed and adapted from Zach Glick's DLPNO-MP2 */

    int nbf = basisset_->nbf();
    int naocc = i_j_to_ij_.size();
    int n_lmo_pairs = ij_to_i_j_.size();

    outfile->Printf("\n  ==> PNO-LMP2 Memory Estimate <== \n\n");

    double F_CUT_MP2 = options_.get_double("F_CUT");

    size_t pno_transform_memory = 0, amplitude_memory = 0, pno_overlap_memory = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : pno_transform_memory, amplitude_memory, pno_overlap_memory)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];
        int npao_ij = lmopair_to_paos_[ij].size();

        pno_transform_memory += npao_ij * n_pno_[ij];
        amplitude_memory += n_pno_[ij] * n_pno_[ij];

        if (i > j) continue;

        for (int k = 0; k < naocc; ++k) {
            int kj = i_j_to_ij_[k][j];
            int ik = i_j_to_ij_[i][k];

            if (kj != -1 && i != k && fabs(F_lmo_->get(i, k)) > F_CUT_MP2) {
                pno_overlap_memory += n_pno_[ij] * n_pno_[kj];
            }
            if (ik != -1 && j != k && fabs(F_lmo_->get(k, j)) > F_CUT_MP2) {
                pno_overlap_memory += n_pno_[ij] * n_pno_[ik];
            }
        }
    }

    size_t total_memory = qia_memory_ + pno_transform_memory + 3 * amplitude_memory + pno_overlap_memory;

    // 1 GB = 1000^3 = 10^9 Bytes
    const double DOUBLES_TO_GB = pow(10.0, -9) * sizeof(double);
    const double WORDS_TO_GB = pow(10.0, -9);

    outfile->Printf("    (q | i a) [AUX, LMO, PAO]  : %8.3f [GB]\n", qia_memory_ * DOUBLES_TO_GB);
    outfile->Printf("    X^{PNO}_{u_{ij} a_{ij}}    : %8.3f [GB]\n", pno_transform_memory * DOUBLES_TO_GB);
    outfile->Printf("    (i a_{ij} | j b_{ij})      : %8.3f [GB]\n", amplitude_memory * DOUBLES_TO_GB);
    outfile->Printf("    T_{ij}^{a_{ij} b_{ij}}     : %8.3f [GB]\n", amplitude_memory * DOUBLES_TO_GB);
    outfile->Printf("    U_{ij}^{a_{ij} b_{ij}}     : %8.3f [GB]\n", amplitude_memory * DOUBLES_TO_GB);
    outfile->Printf("    PNO overlaps               : %8.3f [GB]\n", pno_overlap_memory * DOUBLES_TO_GB);
    outfile->Printf("    Total Memory Required      : %8.3f [GB]\n", total_memory * DOUBLES_TO_GB);
    outfile->Printf("    Total Memory Given         : %8.3f [GB]\n\n", memory_ * WORDS_TO_GB);
    
    if (total_memory * sizeof(double) > 0.9 * memory_) {
        outfile->Printf("  Total Required Memory is more than 90%% of Available Memory!\n");
        outfile->Printf("    We exhausted all of our options!!! This computation cannot continue...\n");
        outfile->Printf("    Pro Tip: Try adjusting the T_CUT_PNO_MP2, T_CUT_TRACE_MP2, and/or T_CUT_ENERGY_MP2 parameters\n");

        throw PSIEXCEPTION("  Too little memory given for PNO-LMP2 Algorithm!");
    }

    // => Form PNO overlap integrals for preceeding DLPNO-MP2 iterations

    std::vector<std::vector<SharedMatrix>> S_pno_ij_ik(n_lmo_pairs);
    std::vector<std::vector<SharedMatrix>> S_pno_ij_kj(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];

        if (i > j) continue;

        S_pno_ij_ik[ij].resize(naocc);
        S_pno_ij_kj[ij].resize(naocc);

        for (int k = 0; k < naocc; ++k) {
            int kj = i_j_to_ij_[k][j];
            int ik = i_j_to_ij_[i][k];

            if (kj != -1 && i != k && fabs(F_lmo_->get(i, k)) > F_CUT_MP2) {
                S_pno_ij_kj[ij][k] = submatrix_rows_and_cols(*S_pao_, lmopair_to_paos_[ij], lmopair_to_paos_[kj]);
                S_pno_ij_kj[ij][k] = linalg::triplet(X_pno_[ij], S_pno_ij_kj[ij][k], X_pno_[kj], true, false, false);
            }
            if (ik != -1 && j != k && fabs(F_lmo_->get(k, j)) > F_CUT_MP2) {
                S_pno_ij_ik[ij][k] = submatrix_rows_and_cols(*S_pao_, lmopair_to_paos_[ij], lmopair_to_paos_[ik]);
                S_pno_ij_ik[ij][k] = linalg::triplet(X_pno_[ij], S_pno_ij_ik[ij][k], X_pno_[ik], true, false, false);
            }
        } // end k
    } // end ij

    // Store the energy for each pair (used to filter out strong and weak pairs later)
    std::vector<double> e_ijs(n_lmo_pairs);

    // => Computing Truncated LMP2 energies (basically running DLPNO-MP2 here)

    outfile->Printf("\n  ==> Iterative Local MP2 with Pair Natural Orbitals (PNOs) <==\n\n");
    outfile->Printf("    E_CONVERGENCE = %.2e\n", 0.01 * options_.get_double("E_CONVERGENCE"));
    outfile->Printf("    R_CONVERGENCE = %.2e\n\n", 0.01 * options_.get_double("R_CONVERGENCE"));
    outfile->Printf("                         Corr. Energy    Delta E     Max R     Time (s)\n");

    std::vector<SharedMatrix> R_iajb(n_lmo_pairs);

    int iteration = 1, max_iteration = options_.get_int("DLPNO_MAXITER");
    double e_curr = 0.0, e_prev = 0.0, r_curr = 0.0;
    bool e_converged = false, r_converged = false;
    DIISManager diis(options_.get_int("DIIS_MAX_VECS"), "LMP2 DIIS", DIISManager::RemovalPolicy::LargestError, DIISManager::StoragePolicy::InCore);

    while (!(e_converged && r_converged)) {
        // RMS of residual per LMO pair, for assessing convergence
        std::vector<double> R_iajb_rms(n_lmo_pairs, 0.0);

        std::time_t time_start = std::time(nullptr);

        // Calculate residuals from current amplitudes
#pragma omp parallel for schedule(dynamic, 1)
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            int i, j;
            std::tie(i, j) = ij_to_i_j_[ij];

            if (i > j) continue;

            R_iajb[ij] = std::make_shared<Matrix>("Residual", n_pno_[ij], n_pno_[ij]);

            if (n_pno_[ij] == 0) continue;

            for (int a = 0; a < n_pno_[ij]; ++a) {
                for (int b = 0; b < n_pno_[ij]; ++b) {
                    R_iajb[ij]->set(a, b,
                                    K_iajb_[ij]->get(a, b) +
                                        (e_pno_[ij]->get(a) + e_pno_[ij]->get(b) - F_lmo_->get(i, i) - F_lmo_->get(j, j)) *
                                            T_iajb_[ij]->get(a, b));
                }
            }

            for (int k = 0; k < naocc; ++k) {
                int kj = i_j_to_ij_[k][j];
                int ik = i_j_to_ij_[i][k];

                if (kj != -1 && i != k && fabs(F_lmo_->get(i, k)) > options_.get_double("F_CUT") && n_pno_[kj] > 0) {
                    auto S_ij_kj = S_pno_ij_kj[ij][k];
                    auto temp = linalg::triplet(S_ij_kj, T_iajb_[kj], S_ij_kj, false, false, true);
                    temp->scale(-1.0 * F_lmo_->get(i, k));
                    R_iajb[ij]->add(temp);
                }
                if (ik != -1 && j != k && fabs(F_lmo_->get(k, j)) > options_.get_double("F_CUT") && n_pno_[ik] > 0) {
                    auto S_ij_ik = S_pno_ij_ik[ij][k];
                    auto temp = linalg::triplet(S_ij_ik, T_iajb_[ik], S_ij_ik, false, false, true);
                    temp->scale(-1.0 * F_lmo_->get(k, j));
                    R_iajb[ij]->add(temp);
                }
            }

            R_iajb_rms[ij] = R_iajb[ij]->rms();

            if (i < j) {
                int ji = ij_to_ji_[ij];
                R_iajb[ji] = R_iajb[ij]->transpose();
                R_iajb_rms[ji] = R_iajb_rms[ij];
            }
        }

// use residuals to get next amplitudes
#pragma omp parallel for schedule(dynamic, 1)
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            int i, j;
            std::tie(i, j) = ij_to_i_j_[ij];
            for (int a = 0; a < n_pno_[ij]; ++a) {
                for (int b = 0; b < n_pno_[ij]; ++b) {
                    T_iajb_[ij]->add(a, b, -R_iajb[ij]->get(a, b) / ((e_pno_[ij]->get(a) + e_pno_[ij]->get(b)) -
                                                                    (F_lmo_->get(i, i) + F_lmo_->get(j, j))));
                }
            }
        }

        // DIIS extrapolation
        auto T_iajb_flat = flatten_mats(T_iajb_);
        auto R_iajb_flat = flatten_mats(R_iajb);

        if (iteration == 1) {
            diis.set_error_vector_size(R_iajb_flat.get());
            diis.set_vector_size(T_iajb_flat.get());
        }

        diis.add_entry(R_iajb_flat.get(), T_iajb_flat.get());
        diis.extrapolate(T_iajb_flat.get());

        copy_flat_mats(T_iajb_flat, T_iajb_);

#pragma omp parallel for schedule(dynamic, 1)
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            Tt_iajb_[ij]->copy(T_iajb_[ij]);
            Tt_iajb_[ij]->scale(2.0);
            Tt_iajb_[ij]->subtract(T_iajb_[ij]->transpose());
        }

        // evaluate convergence using current amplitudes and residuals
        e_prev = e_curr;
        e_curr = 0.0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : e_curr)
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            int i, j;
            std::tie(i, j) = ij_to_i_j_[ij];

            e_ijs[ij] = K_iajb_[ij]->vector_dot(Tt_iajb_[ij]);
            e_curr += e_ijs[ij];
        }

        r_curr = *max_element(R_iajb_rms.begin(), R_iajb_rms.end());

        r_converged = (fabs(r_curr) < 0.01 * options_.get_double("R_CONVERGENCE"));
        e_converged = (fabs(e_curr - e_prev) < 0.01 * options_.get_double("E_CONVERGENCE"));

        std::time_t time_stop = std::time(nullptr);

        outfile->Printf("  @PNO-LMP2 iter %3d: %16.12f %10.3e %10.3e %8d\n", iteration, e_curr, e_curr - e_prev, r_curr, (int)time_stop - (int)time_start);

        iteration++;

        if (iteration > max_iteration + 1) {
            throw PSIEXCEPTION("Maximum DLPNO iterations exceeded.");
        }
    }

    // Set reference LMP2 reference energy to MP2 energy this iteration
    e_lmp2_ = e_curr;

    return e_ijs;
}

void DLPNOCCSD::recompute_pnos() {
    /* Recompute Pair Natural Orbitals for LCCSD iterations */

    timer_on("Compute PNOs (CCSD)");

    int nbf = basisset_->nbf();
    int naocc = i_j_to_ij_.size();
    int n_lmo_pairs = ij_to_i_j_.size();
    const int MIN_PNOS = options_.get_int("MIN_PNOS");

    outfile->Printf("\n  ==> Forming Pair Natural Orbitals (for LCCSD) <==\n");

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];
        int ji = ij_to_ji_[ij];

        if (i > j) continue;

        auto S_pao_ij = submatrix_rows_and_cols(*S_pao_, lmopair_to_paos_[ij], lmopair_to_paos_[ij]);
        auto F_pao_ij = submatrix_rows_and_cols(*F_pao_, lmopair_to_paos_[ij], lmopair_to_paos_[ij]);

        SharedMatrix X_pao_ij;  // canonical transformation of this domain's PAOs to
        SharedVector e_pao_ij;  // energies of the canonical PAOs
        std::tie(X_pao_ij, e_pao_ij) = orthocanonicalizer(S_pao_ij, F_pao_ij);

        S_pao_ij = linalg::triplet(X_pao_ij, S_pao_ij, X_pao_ij, true, false, false);
        auto F_pno_ij_init = linalg::triplet(X_pno_[ij], F_pao_ij, X_pno_[ij], true, false, false);

        // Construct pair density from amplitudes
        auto D_ij = linalg::doublet(Tt_iajb_[ij], T_iajb_[ij], false, true);
        D_ij->add(linalg::doublet(Tt_iajb_[ij], T_iajb_[ij], true, false));
        if (i == j) D_ij->scale(0.5);

        int nvir_ij = F_pno_ij_init->rowspi(0);

        // Diagonalization of pair density gives PNOs (in basis of the LMO's virtual domain) and PNO occ numbers
        auto X_pno_ij = std::make_shared<Matrix>("eigenvectors", nvir_ij, nvir_ij);
        Vector pno_occ("eigenvalues", nvir_ij);
        D_ij->diagonalize(*X_pno_ij, pno_occ, descending);

        double t_cut_scale = 1.0;
        if (i == j) t_cut_scale *= T_CUT_PNO_DIAG_SCALE_;
        if (i < ncore_ || j < ncore_) t_cut_scale *= T_CUT_PNO_CORE_SCALE_;

        double occ_total = 0.0;
        for (size_t a = 0; a < nvir_ij; ++a) {
            occ_total += pno_occ.get(a);
        }

        double e_pno = 0.0;
        double e_ij_total = K_iajb_[ij]->vector_dot(Tt_iajb_[ij]);
        double occ_pno = 0.0;

        int nvir_ij_final = 0;
        std::vector<int> a_curr;
        auto K_pno_init = linalg::triplet(X_pno_ij, K_iajb_[ij], X_pno_ij, true, false, false);
        auto Tt_pno_init = linalg::triplet(X_pno_ij, Tt_iajb_[ij], X_pno_ij, true, false, false);

        for (size_t a = 0; a < nvir_ij; ++a) {
            if (fabs(pno_occ.get(a)) >= t_cut_scale * T_CUT_PNO_ || occ_pno / occ_total < T_CUT_TRACE_ ||
                    std::fabs(e_pno) < T_CUT_ENERGY_ * std::fabs(e_ij_total) || a < MIN_PNOS) {
                a_curr.push_back(a);

                // Energy criteria
                e_pno = submatrix_rows_and_cols(*K_pno_init, a_curr, a_curr)->vector_dot(submatrix_rows_and_cols(*Tt_pno_init, a_curr, a_curr));

                // Trace criteria
                occ_pno += pno_occ.get(a);

                nvir_ij_final++;
            }
        }

        // Make sure there is at least one PNO per pair :)
        nvir_ij_final = std::max(1, nvir_ij_final);

        Dimension zero(1);
        Dimension dim_final(1);
        dim_final.fill(nvir_ij_final);

        // This transformation gives orbitals that are orthonormal but not canonical
        X_pno_ij = X_pno_ij->get_block({zero, X_pno_ij->rowspi()}, {zero, dim_final});
        pno_occ = pno_occ.get_block({zero, dim_final});

        SharedMatrix pno_canon;
        SharedVector e_pno_ij;
        std::tie(pno_canon, e_pno_ij) = canonicalizer(X_pno_ij, F_pno_ij_init);

        // This transformation gives orbitals that are orthonormal and canonical
        X_pno_ij = linalg::doublet(X_pno_ij, pno_canon, false, false);

        auto K_pno_ij = linalg::triplet(X_pno_ij, K_iajb_[ij], X_pno_ij, true, false, false);
        auto T_pno_ij = linalg::triplet(X_pno_ij, T_iajb_[ij], X_pno_ij, true, false, false);
        auto Tt_pno_ij = linalg::triplet(X_pno_ij, Tt_iajb_[ij], X_pno_ij, true, false, false);

        // (additional) truncation error
        double de_pno_ij = K_iajb_[ij]->vector_dot(Tt_iajb_[ij]) - K_pno_ij->vector_dot(Tt_pno_ij);

        // New PNO transformation matrix
        X_pno_ij = linalg::doublet(X_pno_[ij], X_pno_ij, false, false);

        K_iajb_[ij] = K_pno_ij;
        T_iajb_[ij] = T_pno_ij;
        Tt_iajb_[ij] = Tt_pno_ij;
        X_pno_[ij] = X_pno_ij;
        e_pno_[ij] = e_pno_ij;
        n_pno_[ij] = X_pno_ij->colspi(0);
        occ_pno_[ij] = pno_occ.get(n_pno_[ij] - 1);
        trace_pno_[ij] = occ_pno / occ_total;
        e_ratio_pno_[ij] = e_pno / e_ij_total;
        de_pno_[ij] += de_pno_ij;

        // account for symmetry
        if (i < j) {
            K_iajb_[ji] = K_iajb_[ij]->transpose();
            T_iajb_[ji] = T_iajb_[ij]->transpose();
            Tt_iajb_[ji] = Tt_iajb_[ij]->transpose();
            X_pno_[ji] = X_pno_[ij];
            e_pno_[ji] = e_pno_[ij];
            n_pno_[ji] = n_pno_[ij];
            occ_pno_[ji] = occ_pno_[ij];
            trace_pno_[ji] = trace_pno_[ij];
            e_ratio_pno_[ji] = e_ratio_pno_[ij];
            de_pno_[ji] += de_pno_ij;
        } // end if (i < j)
    }

    // Print out PNO domain information
    int pno_count_total = 0, pno_count_min = nbf, pno_count_max = 0;
    double occ_number_total = 0.0, occ_number_min = 2.0, occ_number_max = 0.0;
    double trace_total = 0.0, trace_min = 1.0, trace_max = 0.0;
    double energy_total = 0.0, energy_min = 1.0, energy_max = 0.0;
    de_weak_ = 0.0, de_pno_total_ = 0.0;

    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];

        if (i_j_to_ij_strong_[i][j] == -1) de_weak_ += K_iajb_[ij]->vector_dot(Tt_iajb_[ij]);

        pno_count_total += n_pno_[ij];
        pno_count_min = std::min(pno_count_min, n_pno_[ij]);
        pno_count_max = std::max(pno_count_max, n_pno_[ij]);

        occ_number_total += occ_pno_[ij];
        occ_number_min = std::min(occ_number_min, occ_pno_[ij]);
        occ_number_max = std::max(occ_number_max, occ_pno_[ij]);

        trace_total += trace_pno_[ij];
        trace_min = std::min(trace_min, trace_pno_[ij]);
        trace_max = std::max(trace_max, trace_pno_[ij]);

        energy_total += e_ratio_pno_[ij];
        energy_min = std::min(energy_min, e_ratio_pno_[ij]);
        energy_max = std::max(energy_max, e_ratio_pno_[ij]);

        de_pno_total_ += de_pno_[ij];
    }

    outfile->Printf("  \n");
    outfile->Printf("    Natural Orbitals per Local MO pair:\n");
    outfile->Printf("      Avg: %3d NOs \n", pno_count_total / n_lmo_pairs);
    outfile->Printf("      Min: %3d NOs \n", pno_count_min);
    outfile->Printf("      Max: %3d NOs \n\n", pno_count_max);

    outfile->Printf("      Avg Occ Number Tol: %.3e \n", occ_number_total / n_lmo_pairs);
    outfile->Printf("      Min Occ Number Tol: %.3e \n", occ_number_min);
    outfile->Printf("      Max Occ Number Tol: %.3e \n\n", occ_number_max);

    outfile->Printf("      Avg Trace Sum: %.6f \n", trace_total / n_lmo_pairs);
    outfile->Printf("      Min Trace Sum: %.6f \n", trace_min);
    outfile->Printf("      Max Trace Sum: %.6f \n\n", trace_max);

    outfile->Printf("      Avg Energy Ratio: %.6f \n", energy_total / n_lmo_pairs);
    outfile->Printf("      Min Energy Ratio: %.6f \n", energy_min);
    outfile->Printf("      Max Energy Ratio: %.6f \n\n", energy_max);

    outfile->Printf("    LMP2 Weak Pair energy = %.12f\n", de_weak_);
    outfile->Printf("    PNO truncation energy = %.12f\n", de_pno_total_);

    timer_off("Compute PNOs (CCSD)");
}

template<bool crude> double DLPNOCCSD::filter_pairs(const std::vector<double>& e_ijs) {
    /* If crude, this function filters out semicanonical MP2 pairs. Otherwise, it filters out the weak pairs from the strong pairs */

    int natom = molecule_->natom();
    int nbf = basisset_->nbf();
    int naocc = i_j_to_ij_.size();
    int n_lmo_pairs = ij_to_i_j_.size();
    int naux = ribasis_->nbf();
    int npao = C_pao_->colspi(0);  // same as nbf

    if constexpr (crude) {
        // Clear information from old maps
        std::vector<std::vector<int>> i_j_to_ij_new;
        std::vector<std::pair<int,int>> ij_to_i_j_new;
        std::vector<int> ij_to_ji_new;

        // Allocate space for new intermediates
        i_j_to_ij_new.resize(naocc);

        for (size_t i = 0; i < naocc; i++) {
            i_j_to_ij_new[i].resize(naocc, -1);
        }

        double delta_e_crude = 0.0;

        int ij_new = 0;
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            auto &[i, j] = ij_to_i_j_[ij];
            if (std::fabs(e_ijs[ij]) >= T_CUT_PAIRS_MP2_) { // If this pair survives, it continues in the computation
                i_j_to_ij_new[i][j] = ij_new;
                ij_to_i_j_new.push_back(std::make_pair(i, j));
                ++ij_new;
            } else { // Otherwise, it dies forever :)
                delta_e_crude += e_ijs[ij];
            }
        }

        for (size_t ij = 0; ij < ij_to_i_j_new.size(); ++ij) {
            auto &[i, j] = ij_to_i_j_new[ij];
            ij_to_ji_new.push_back(i_j_to_ij_new[j][i]);
        }

        i_j_to_ij_ = i_j_to_ij_new;
        ij_to_i_j_ = ij_to_i_j_new;
        ij_to_ji_ = ij_to_ji_new;

        return delta_e_crude;

    } else {
        i_j_to_ij_strong_.resize(naocc);
        i_j_to_ij_weak_.resize(naocc);

        for (size_t i = 0; i < naocc; i++) {
            i_j_to_ij_strong_[i].resize(naocc, -1);
            i_j_to_ij_weak_[i].resize(naocc, -1);
        }

        double delta_e_weak = 0.0;

        int ij_strong = 0, ij_weak = 0;
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            auto &[i, j] = ij_to_i_j_[ij];
            if (std::fabs(e_ijs[ij]) >= T_CUT_PAIRS_) { // Pair is strong pair
                i_j_to_ij_strong_[i][j] = ij_strong;
                ij_to_i_j_strong_.push_back(std::make_pair(i, j));
                ++ij_strong;
            } else { // Pair is weak pair
                i_j_to_ij_weak_[i][j] = ij_weak;
                ij_to_i_j_weak_.push_back(std::make_pair(i, j));
                delta_e_weak += e_ijs[ij];
                ++ij_weak;
            } // end else
        } // end ij

        ij_to_ji_strong_.clear();
        ij_to_ji_weak_.clear();

        for (size_t ij = 0; ij < ij_to_i_j_strong_.size(); ++ij) {
            auto &[i, j] = ij_to_i_j_strong_[ij];
            ij_to_ji_strong_.push_back(i_j_to_ij_strong_[j][i]);
        } // end ij
        
        for (size_t ij = 0; ij < ij_to_i_j_weak_.size(); ++ij) {
            auto &[i, j] = ij_to_i_j_weak_[ij];
            ij_to_ji_weak_.push_back(i_j_to_ij_weak_[j][i]);
        } // end ij

        return delta_e_weak;
    }
}

template<bool crude> void DLPNOCCSD::pair_prescreening() {
    
    int naocc = i_j_to_ij_.size();

    if constexpr (crude) {
        outfile->Printf("\n  ==> Initial semi-canonical MP2 prescreening of pairs <==\n");

        int n_lmo_pairs_init = ij_to_i_j_.size();

        const std::vector<double>& e_ijs_crude = compute_pair_energies<true>();
        de_lmp2_eliminated_ = filter_pairs<true>(e_ijs_crude);

        int n_lmo_pairs_final = ij_to_i_j_.size();
        int n_eliminated_pairs = n_lmo_pairs_init - n_lmo_pairs_final;

        outfile->Printf("    Eliminated Pairs (SC-LMP2)         = %d\n", n_eliminated_pairs);
        outfile->Printf("    Surviving Pairs                    = %d\n", n_lmo_pairs_final);
        outfile->Printf("    Surviving Pairs / Non-dipole Pairs = (%.2f %%)\n", (100.0 * n_lmo_pairs_final) / (n_lmo_pairs_init));
        outfile->Printf("    Eliminated Pair dE                 = %.12f\n\n", de_lmp2_eliminated_);
    } else {
        outfile->Printf("\n  ==> Determining Strong and Weak Pairs (Refined Prescreening Step) <==\n");

        int n_lmo_pairs = ij_to_i_j_.size();

        compute_pair_energies<false>();

        // Compute full LMP2 energies for remaining pairs for a better energy estimate
        timer_on("PNO-LMP2 Iterations");
        const std::vector<double>& e_ijs = pno_lmp2_iterations();
        timer_off("PNO-LMP2 Iterations");

        filter_pairs<false>(e_ijs);

        int n_strong_pairs = ij_to_i_j_strong_.size();
        int n_weak_pairs = ij_to_i_j_weak_.size();

        outfile->Printf("\n  ==> Final Strong and Weak Pairs <==\n\n");
        outfile->Printf("    Weak Pairs                      = %d\n", n_weak_pairs);
        outfile->Printf("    Strong Pairs                    = %d\n", n_strong_pairs);
        outfile->Printf("    Strong Pairs / Total Pairs      = (%.2f %%)\n", (100.0 * n_strong_pairs) / n_lmo_pairs);
    }
}

void DLPNOCCSD::compute_pno_integrals() {
    outfile->Printf("    Computing integrals in the PNO basis from PAO integrals...\n\n");

    int n_lmo_pairs = ij_to_i_j_.size();
    
    // 1 virtual
    K_mibj_.resize(n_lmo_pairs);
    J_ijmb_.resize(n_lmo_pairs);
    L_mibj_.resize(n_lmo_pairs);

    // 2 virtual
    L_iajb_.resize(n_lmo_pairs);

    // 2-virtual non-projected
    J_ikac_non_proj_.resize(n_lmo_pairs);
    K_iakc_non_proj_.resize(n_lmo_pairs);

    // 3 virtual
    K_ivvv_.resize(n_lmo_pairs);

    // DF integrals (only allocate if writing to RAM)
    if (!write_qia_pno_) {
        Qma_ij_.resize(n_lmo_pairs);
    }
    if (!write_qab_pno_) {
        Qab_ij_.resize(n_lmo_pairs);
    }

    i_Qa_ij_.resize(n_lmo_pairs);
    i_Qk_ij_.resize(n_lmo_pairs);

    psio_->open(PSIF_DLPNO_QIA_PNO, PSIO_OPEN_NEW);
    psio_->open(PSIF_DLPNO_QAB_PNO, PSIO_OPEN_NEW);

    std::time_t time_start = std::time(nullptr);
    std::time_t time_lap = std::time(nullptr);

    // Sort pairs by the approximate number of operations (for maximal parallel efficiency)
    std::vector<std::pair<int, size_t>> ij_cost_tuple(n_lmo_pairs);
    
#pragma omp parallel for
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];

        const int npao_ij = lmopair_to_paos_[ij].size();
        const int naux_ij = lmopair_to_ribfs_[ij].size();
        const int nlmo_ij = lmopair_to_lmos_[ij].size();
        
        size_t cost = 0;

        if (i <= j) {
            // Cost of transforming q_vv from PAO to PNO basis
            cost += naux_ij * npao_ij * npao_ij * n_pno_[ij];
            cost += naux_ij * npao_ij * n_pno_[ij] * n_pno_[ij];
            // Cost of transforming q_ov from PAO to PNO basis
            cost += naux_ij * npao_ij * nlmo_ij * n_pno_[ij];

            // For strong pairs, there are the non-projected integrals
            if (i_j_to_ij_strong_[i][j] != -1) {

                std::vector<int> extended_pao_domain = lmopair_to_paos_[ij];

                for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
                    int k = lmopair_to_lmos_[ij][k_ij];
                    extended_pao_domain = merge_lists(extended_pao_domain, lmo_to_paos_[k]);
                }

                const int npao_ext_ij = extended_pao_domain.size();
                cost += naux_ij * npao_ext_ij * npao_ij * n_pno_[ij]; // Most expensive step in forming non-projected integrals
            } // end if
        }

        ij_cost_tuple[ij] = std::make_pair(ij, cost);
    }
    
    std::sort(ij_cost_tuple.begin(), ij_cost_tuple.end(), [&](const std::pair<int, size_t>& a, const std::pair<int, size_t>& b) {
        return (a.second > b.second);
    });

    std::vector<int> ij_sorted_by_cost(n_lmo_pairs);
    
#pragma omp parallel for
    for (int ij_idx = 0; ij_idx < n_lmo_pairs; ++ij_idx) {
        ij_sorted_by_cost[ij_idx] = ij_cost_tuple[ij_idx].first;
    }

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij_idx = 0; ij_idx < n_lmo_pairs; ++ij_idx) {
        int ij = ij_sorted_by_cost[ij_idx];

        auto &[i, j] = ij_to_i_j_[ij];
        const int ji = ij_to_ji_[ij];
        const bool is_strong_pair = i_j_to_ij_strong_[i][j] != -1;

        // number of PNOs in the pair domain
        const int npno_ij = n_pno_[ij];

        if (i > j) continue;

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        // number of LMOs in the pair domain
        const int nlmo_ij = lmopair_to_lmos_[ij].size();
        // number of PAOs in the pair domain (before removing linear dependencies)
        const int npao_ij = lmopair_to_paos_[ij].size();
        // number of auxiliary functions in the pair domain
        const int naux_ij = lmopair_to_ribfs_[ij].size();

        auto q_pair = std::make_shared<Matrix>(naux_ij, 1);

        auto q_io = std::make_shared<Matrix>(naux_ij, nlmo_ij);
        auto q_jo = std::make_shared<Matrix>(naux_ij, nlmo_ij);

        auto q_iv = std::make_shared<Matrix>(naux_ij, npno_ij);
        auto q_jv = std::make_shared<Matrix>(naux_ij, npno_ij);

        auto q_ov = std::make_shared<Matrix>(naux_ij, nlmo_ij * npno_ij);
        auto q_vv = std::make_shared<Matrix>(naux_ij, npno_ij * npno_ij);
        
        J_ikac_non_proj_[ij].resize(nlmo_ij);
        if (i != j) J_ikac_non_proj_[ji].resize(nlmo_ij);

        K_iakc_non_proj_[ij].resize(nlmo_ij);
        if (i != j) K_iakc_non_proj_[ji].resize(nlmo_ij);

        // Extended PAO domain of ij is formed from the union of PAOs
        // from all interacting k_ij (LMOs k such that ik AND kj form a valid pair)
        std::vector<int> extended_pao_domain;
        int npao_ext_ij;
        
        extended_pao_domain = lmopair_to_paos_[ij];
        for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
            int k = lmopair_to_lmos_[ij][k_ij];
            extended_pao_domain = merge_lists(extended_pao_domain, lmo_to_paos_[k]);
        }
        npao_ext_ij = extended_pao_domain.size();

        if (thread == 0) timer_on("DLPNO-CCSD: Setup Integrals");

        for (int q_ij = 0; q_ij < naux_ij; q_ij++) {
            const int q = lmopair_to_ribfs_[ij][q_ij];
            const int centerq = ribasis_->function_to_center(q);

            const int i_sparse = riatom_to_lmos_ext_dense_[centerq][i];
            const int j_sparse = riatom_to_lmos_ext_dense_[centerq][j];
            const std::vector<int> i_slice(1, i_sparse);
            const std::vector<int> j_slice(1, j_sparse);

            const auto sparse_lmo_list = index_list(riatom_to_lmos_ext_[centerq], lmopair_to_lmos_[ij]);
            const auto sparse_pao_list = index_list(riatom_to_paos_ext_[centerq], lmopair_to_paos_[ij]);

            q_pair->set(q_ij, 0, (*qij_[q])(i_sparse, j_sparse));
            
            auto q_io_tmp = submatrix_rows_and_cols(*qij_[q], i_slice, sparse_lmo_list);
            ::memcpy(&(*q_io)(q_ij, 0), &(*q_io_tmp)(0,0), nlmo_ij * sizeof(double));

            auto q_jo_tmp = submatrix_rows_and_cols(*qij_[q], j_slice, sparse_lmo_list);
            ::memcpy(&(*q_jo)(q_ij, 0), &(*q_jo_tmp)(0,0), nlmo_ij * sizeof(double));

            auto q_iv_tmp = submatrix_rows_and_cols(*qia_[q], i_slice, sparse_pao_list);
            q_iv_tmp = linalg::doublet(q_iv_tmp, X_pno_[ij], false, false);
            ::memcpy(&(*q_iv)(q_ij, 0), &(*q_iv_tmp)(0,0), npno_ij * sizeof(double));

            auto q_jv_tmp = submatrix_rows_and_cols(*qia_[q], j_slice, sparse_pao_list);
            q_jv_tmp = linalg::doublet(q_jv_tmp, X_pno_[ij], false, false);
            ::memcpy(&(*q_jv)(q_ij, 0), &(*q_jv_tmp)(0,0), npno_ij * sizeof(double));

            auto q_ov_tmp = submatrix_rows_and_cols(*qia_[q], sparse_lmo_list, sparse_pao_list);
            q_ov_tmp = linalg::doublet(q_ov_tmp, X_pno_[ij], false, false);
            ::memcpy(&(*q_ov)(q_ij, 0), &(*q_ov_tmp)(0,0), nlmo_ij * npno_ij * sizeof(double));

            SharedMatrix q_vv_tmp = std::make_shared<Matrix>(npao_ij, npao_ij);
            
            for (int u_ij = 0; u_ij < npao_ij; ++u_ij) {
                int u = lmopair_to_paos_[ij][u_ij];
                for (int v_ij = 0; v_ij < npao_ij; ++v_ij) {
                    int v = lmopair_to_paos_[ij][v_ij];
                    int uv_idx = riatom_to_pao_pairs_dense_[centerq][u][v];
                    if (uv_idx == -1) continue;
                    q_vv_tmp->set(u_ij, v_ij, qab_[q]->get(uv_idx, 0));
                }
            }
            
            q_vv_tmp = linalg::triplet(X_pno_[ij], q_vv_tmp, X_pno_[ij], true, false, false);
            ::memcpy(&(*q_vv)(q_ij, 0), &(*q_vv_tmp)(0,0), npno_ij * npno_ij * sizeof(double));
        }

        auto A_solve = submatrix_rows_and_cols(*full_metric_, lmopair_to_ribfs_[ij], lmopair_to_ribfs_[ij]);
        SharedMatrix q_io_clone;
        SharedMatrix q_jo_clone;
        SharedMatrix q_iv_clone;
        SharedMatrix q_jv_clone;
        
        // (P_{ij}|Q_{ij})^{-1} (Q_{ij}|m_{ij}i)
        q_io_clone = q_io->clone();
        C_DGESV_wrapper(A_solve->clone(), q_io_clone);
        if (i != j) {
            q_jo_clone = q_jo->clone();
            C_DGESV_wrapper(A_solve->clone(), q_jo_clone);
        }
        
        // (P_{ij}|Q_{ij})^{-1} (Q_{ij}|a_{ij}i)
        q_iv_clone = q_iv->clone();
        C_DGESV_wrapper(A_solve->clone(), q_iv_clone);
        if (i != j) {
            q_jv_clone = q_jv->clone();
            C_DGESV_wrapper(A_solve->clone(), q_jv_clone);
        }
        
        A_solve->power(0.5, 1.0e-14);

        // (P_{ij}|Q_{ij})^{-1/2} (Q_{ij}|p q)
        C_DGESV_wrapper(A_solve->clone(), q_pair);
        C_DGESV_wrapper(A_solve->clone(), q_io);
        C_DGESV_wrapper(A_solve->clone(), q_jo);
        C_DGESV_wrapper(A_solve->clone(), q_iv);
        C_DGESV_wrapper(A_solve->clone(), q_jv);
        C_DGESV_wrapper(A_solve->clone(), q_ov);
        C_DGESV_wrapper(A_solve->clone(), q_vv);

        if (thread == 0) timer_off("DLPNO-CCSD: Setup Integrals");

        if (is_strong_pair) {

            if (thread == 0) timer_on("DLPNO-CCSD: J_ikac integrals");

            // Partially transformed q_vv with first index transformed into PNO virtual space,
            // second in extended PAO space
            SharedMatrix q_vv_partial = std::make_shared<Matrix>("q_vv_partial", naux_ij, npno_ij * npao_ext_ij);

            // Compute (i k_{ij} | a_{ij} b_{kj}) integrals through density-fitting
            for (int q_ij = 0; q_ij < naux_ij; q_ij++) {
                const int q = lmopair_to_ribfs_[ij][q_ij];
                const int centerq = ribasis_->function_to_center(q);
                
                auto q_cd_temp = std::make_shared<Matrix>(npao_ij, npao_ext_ij);
                q_cd_temp->zero();
                for (int u_ij = 0; u_ij < npao_ij; ++u_ij) {
                    int u = lmopair_to_paos_[ij][u_ij];
                    for (int v_ij = 0; v_ij < npao_ext_ij; ++v_ij) {
                        int v = extended_pao_domain[v_ij];
                        int uv_idx = riatom_to_pao_pairs_dense_[centerq][u][v];
                        if (uv_idx == -1) continue;
                        q_cd_temp->set(u_ij, v_ij, qab_[q]->get(uv_idx, 0));
                    }
                }
                q_cd_temp = linalg::doublet(X_pno_[ij], q_cd_temp, true, false);
                ::memcpy(&(*q_vv_partial)(q_ij, 0), q_cd_temp->get_pointer(), npno_ij * npao_ext_ij * sizeof(double));
            }

            SharedMatrix K_iovv_partial = linalg::doublet(q_io_clone, q_vv_partial, true, false);
            SharedMatrix K_jovv_partial;

            if (i != j) {
                K_jovv_partial = linalg::doublet(q_jo_clone, q_vv_partial, true, false);
            }

            for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
                int k = lmopair_to_lmos_[ij][k_ij];
                int ik = i_j_to_ij_[i][k], kj = i_j_to_ij_[k][j];

                auto K_iovv_slice = submatrix_rows(*K_iovv_partial, std::vector<int>(1, k_ij));
                K_iovv_slice->reshape(npno_ij, npao_ext_ij);
                K_iovv_slice = submatrix_cols(*K_iovv_slice, index_list(extended_pao_domain, lmopair_to_paos_[kj]));
                J_ikac_non_proj_[ij][k_ij] = linalg::doublet(K_iovv_slice, X_pno_[kj]);

                if (i != j) {
                    auto K_jovv_slice = submatrix_rows(*K_jovv_partial, std::vector<int>(1, k_ij));
                    K_jovv_slice->reshape(npno_ij, npao_ext_ij);
                    K_jovv_slice = submatrix_cols(*K_jovv_slice, index_list(extended_pao_domain, lmopair_to_paos_[ik]));
                    J_ikac_non_proj_[ji][k_ij] = linalg::doublet(K_jovv_slice, X_pno_[ik]);
                } // end if
            } // end for

            if (thread == 0) timer_off("DLPNO-CCSD: J_ikac integrals");

            if (thread == 0) timer_on("DLPNO-CCSD: K_iakc integrals");

            // Global q_ov_ext containing lmos in ij and paos in extended domain of ij
            SharedMatrix q_ov_ext = std::make_shared<Matrix>("q_ov_ext", naux_ij, nlmo_ij * npao_ext_ij);

            // Compute (i k_{ij} | a_{ij} b_{kj}) integrals through density-fitting

            for (int q_ij = 0; q_ij < naux_ij; q_ij++) {
                const int q = lmopair_to_ribfs_[ij][q_ij];
                const int centerq = ribasis_->function_to_center(q);

                const auto sparse_lmo_list = index_list(riatom_to_lmos_ext_[centerq], lmopair_to_lmos_[ij]);
                const auto sparse_pao_list = index_list(riatom_to_paos_ext_[centerq], extended_pao_domain);

                auto q_ov_ext_tmp = submatrix_rows_and_cols(*qia_[q], sparse_lmo_list, sparse_pao_list);
                ::memcpy(&(*q_ov_ext)(q_ij, 0), q_ov_ext_tmp->get_pointer(), nlmo_ij * npao_ext_ij * sizeof(double));
            }

            SharedMatrix K_oviv = linalg::doublet(q_ov_ext, q_iv_clone, true, false);
            K_oviv->reshape(nlmo_ij, npao_ext_ij * npno_ij);
            SharedMatrix K_ovjv;

            if (i != j) {
                K_ovjv = linalg::doublet(q_ov_ext, q_jv_clone, true, false);
                K_ovjv->reshape(nlmo_ij, npao_ext_ij * npno_ij);
            }

            for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
                int k = lmopair_to_lmos_[ij][k_ij];
                int ik = i_j_to_ij_[i][k], kj = i_j_to_ij_[k][j];

                auto K_oviv_slice = submatrix_rows(*K_oviv, std::vector<int>(1, k_ij));
                K_oviv_slice->reshape(npao_ext_ij, npno_ij);
                K_oviv_slice = submatrix_rows(*K_oviv_slice, index_list(extended_pao_domain, lmopair_to_paos_[kj]));
                K_oviv_slice = linalg::doublet(X_pno_[kj], K_oviv_slice, true, false);
                K_iakc_non_proj_[ij][k_ij] = K_oviv_slice->transpose();

                if (i != j) {
                    auto K_ovjv_slice = submatrix_rows(*K_ovjv, std::vector<int>(1, k_ij));
                    K_ovjv_slice->reshape(npao_ext_ij, npno_ij);
                    K_ovjv_slice = submatrix_rows(*K_ovjv_slice, index_list(extended_pao_domain, lmopair_to_paos_[ik]));
                    K_ovjv_slice = linalg::doublet(X_pno_[ik], K_ovjv_slice, true, false);
                    K_iakc_non_proj_[ji][k_ij] = K_ovjv_slice->transpose();
                } // end if
            }

            if (thread == 0) timer_off("DLPNO-CCSD: K_iakc integrals");

        } // end if

        if (thread == 0) timer_on("DLPNO-CCSD: Contract Integrals");
        
        K_mibj_[ij] = linalg::doublet(q_io, q_jv, true, false);
        J_ijmb_[ij] = linalg::doublet(q_pair, q_ov, true, false);
        J_ijmb_[ij]->reshape(nlmo_ij, npno_ij);
        K_ivvv_[ij] = linalg::doublet(q_iv, q_vv, true, false);

        L_mibj_[ij] = K_mibj_[ij]->clone();

        if (i != j) {
            K_mibj_[ji] = linalg::doublet(q_jo, q_iv, true, false);
            J_ijmb_[ji] = J_ijmb_[ij];
            K_ivvv_[ji] = linalg::doublet(q_jv, q_vv, true, false);

            L_mibj_[ij]->scale(2.0);
            L_mibj_[ij]->subtract(K_mibj_[ji]);

            L_mibj_[ji] = K_mibj_[ji]->clone();
            L_mibj_[ji]->scale(2.0);
            L_mibj_[ji]->subtract(L_mibj_[ij]);
        }

        // Save DF integrals (only for strong pairs)
        if (is_strong_pair) {
            i_Qk_ij_[ij] = q_io;
            i_Qa_ij_[ij] = q_iv;

            if (i != j) {
                i_Qk_ij_[ji] = q_jo;
                i_Qa_ij_[ji] = q_jv;
            }
        }
        
        if (is_strong_pair) {
            if (!write_qia_pno_) {
                Qma_ij_[ij].resize(naux_ij);
                for (int q_ij = 0; q_ij < naux_ij; ++q_ij) {
                    // Save transformed (Q_ij | m_ij a_ij) integrals
                    Qma_ij_[ij][q_ij] = std::make_shared<Matrix>(nlmo_ij, npno_ij);
                    ::memcpy(&(*Qma_ij_[ij][q_ij])(0, 0), &(*q_ov)(q_ij,0), nlmo_ij * npno_ij * sizeof(double));
                }
            } else {
                std::stringstream toc_entry;
                toc_entry << "QIA (PNO) " << ij;
                q_ov->set_name(toc_entry.str());
    #pragma omp critical
                q_ov->save(psio_, PSIF_DLPNO_QIA_PNO, psi::Matrix::SubBlocks);
            }
        }

        if (is_strong_pair) {
            if (!write_qab_pno_) {
                Qab_ij_[ij].resize(naux_ij);
                for (int q_ij = 0; q_ij < naux_ij; ++q_ij) {
                    // Save transformed (Q_ij | a_ij b_ij) integrals
                    Qab_ij_[ij][q_ij] = std::make_shared<Matrix>(npno_ij, npno_ij);
                    ::memcpy(&(*Qab_ij_[ij][q_ij])(0, 0), &(*q_vv)(q_ij,0), npno_ij * npno_ij * sizeof(double));
                }
            } else {
                std::stringstream toc_entry;
                toc_entry << "QAB (PNO) " << ij;
                q_vv->set_name(toc_entry.str());
    #pragma omp critical
                q_vv->save(psio_, PSIF_DLPNO_QAB_PNO, psi::Matrix::ThreeIndexLowerTriangle);
            }
        }

        // L_iajb
        L_iajb_[ij] = K_iajb_[ij]->clone();
        L_iajb_[ij]->scale(2.0);
        L_iajb_[ij]->subtract(K_iajb_[ij]->transpose());

        if (i != j) {
            L_iajb_[ji] = L_iajb_[ij]->transpose();
        }

        if (thread == 0) timer_off("DLPNO-CCSD: Contract Integrals");

        if (thread == 0) {
            std::time_t time_curr = std::time(nullptr);
            int time_elapsed = (int) time_curr - (int) time_lap;
            if (time_elapsed > 60) {
                outfile->Printf("  Time Elapsed from last checkpoint %4d (s), Progress %2d %%, Integrals for (%4d / %4d) Pairs Computed\n", time_elapsed, 
                                    (100 * ij_idx) / n_lmo_pairs, ij_idx, n_lmo_pairs);
                time_lap = std::time(nullptr);
            }
        }
    }

    std::time_t time_stop = std::time(nullptr);
    int time_elapsed = (int) time_stop - (int) time_start;
    outfile->Printf("    Integral Computation Complete!!! Time Elapsed: %4d seconds\n\n", time_elapsed);
}

void DLPNOCCSD::t1_ints() {

    timer_on("DLPNO-CCSD: T1 Ints");

    int naocc = nalpha_ - nfrzc();
    int n_lmo_pairs = ij_to_i_j_.size();

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];

        int nlmo_ij = lmopair_to_lmos_[ij].size();
        int naux_ij = lmopair_to_ribfs_[ij].size();
        int npno_ij = n_pno_[ij];
        int pair_idx = (i > j) ? ij_to_ji_[ij] : ij;
        int i_ij = lmopair_to_lmos_dense_[ij][i], j_ij = lmopair_to_lmos_dense_[ij][j];

        // These are only needed to be computed over strong pairs
        if (i_j_to_ij_strong_[i][j] == -1) continue;

        // => Jiang Eq. 91 <= //
        // \widetilde{B}^{Q}_{ki} = B^{Q}_{ki} + B^{Q}_{ka} T_{i}^{a}
        i_Qk_t1_[ij] = i_Qk_ij_[ij]->clone(); // (Q, k)

        auto qma_ij = QIA_PNO(ij);
        for (int q_ij = 0; q_ij < naux_ij; ++q_ij) {
            for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
                for (int a_ij = 0; a_ij < npno_ij; ++a_ij) {
                    (*i_Qk_t1_[ij])(q_ij, k_ij) += (*qma_ij[q_ij])(k_ij, a_ij) * (*T_n_ij_[ij])(i_ij, a_ij);
                }
            }
        }

        // => Jiang Eq. 92 <= //
        // \widetilde{B}^{Q}_{ai} = B^{Q}_{ai} - T_{k}^{a} B^{Q}_{ki} + B^{Q}_{ab} T_{i}^{b} - T_{k}^{a} B^{Q}_{kb} T_{i}^{b}
        // = B^{Q}_{ai} - T_{k}^{a} \widetilde{B}^{Q}_{ki} + B^{Q}_{ab} T_{i}^{b}
        i_Qa_t1_[ij] = i_Qa_ij_[ij]->clone(); // (Q, a)
        i_Qa_t1_[ij]->subtract(linalg::doublet(i_Qk_t1_[ij], T_n_ij_[ij])); // (Q, k) (k, a)

        auto qab_ij = QAB_PNO(ij);
        for (int q_ij = 0; q_ij < naux_ij; ++q_ij) {
            for (int a_ij = 0; a_ij < npno_ij; ++a_ij) {
                for (int b_ij = 0; b_ij < npno_ij; ++b_ij) {
                    (*i_Qa_t1_[ij])(q_ij, a_ij) += (*qab_ij[q_ij])(a_ij, b_ij) * (*T_n_ij_[ij])(i_ij, b_ij); // (Q, a, b) (i, b)
                } // end b_ij
            } // end a_ij
        } // end q_ij
    }

    timer_off("DLPNO-CCSD: T1 Ints");
}

void DLPNOCCSD::t1_fock() {

    timer_on("DLPNO-CCSD: T1 Fock");

    int naocc = nalpha_ - nfrzc();
    int n_lmo_pairs = ij_to_i_j_.size();

    // => Step 1: Dressing over the contracted indices <= //

    SharedMatrix Fij_bar = F_lmo_->clone(); // (i, j)
    std::vector<SharedMatrix> Fkc_bar(n_lmo_pairs); // (k_{ij}, c_{ij})
    std::vector<SharedMatrix> Fai_bar(naocc); // (a_{ii})
    std::vector<SharedMatrix> Fab_bar(n_lmo_pairs); // (a_{ij}, b_{ij})

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];
        int ji = ij_to_ji_[ij];
        int pair_idx = (i > j) ? ji : ij;

        int naux_ij = lmopair_to_ribfs_[ij].size();
        int nlmo_ij = lmopair_to_lmos_[ij].size();
        int npno_ij = n_pno_[ij];

        // Partially dress Fij (Jiang Eq. 98)
        // \overline{F}_{ij} = f_{ij} + [2(ij|kc) - (ic|kj)] T_{k}^{c} =>
        // \overline{F}_{ij} = f_{ij} + [2J_{ij}^{kc} - K_{ji}^{kc})] T_{k}^{c}
        (*Fij_bar)(i, j) += 2.0 * T_n_ij_[ij]->vector_dot(J_ijmb_[ij]);
        (*Fij_bar)(i, j) -= T_n_ij_[ij]->vector_dot(K_mibj_[ji]);

        if (i > j || i_j_to_ij_strong_[i][j] == -1) continue;

        // Partially dress Fia and Fab (Jiang Eq. 99 and 101)
        // \overline{F}_{kc} = f_{kc} + [2(kc|me) - (ke|mc)] T_{m}^{e}
        // In closed-shell RHF reference, f_{kc} is zero
        Fkc_bar[ij] = std::make_shared<Matrix>(nlmo_ij, npno_ij);

        // \overline{F}_{ab} = f_{ab} + [2(ab|me) - (ae|mb)] T_{m}^{e}
        // In canonical PNO representation, f_{ab} is diagonal
        Fab_bar[ij] = std::make_shared<Matrix>(npno_ij, npno_ij);
        Fab_bar[ij]->set_diagonal(e_pno_[ij]);

        auto qma_ij = QIA_PNO(ij);
        auto qab_ij = QAB_PNO(ij);
        for (int q_ij = 0; q_ij < naux_ij; ++q_ij) {
            auto Qma = qma_ij[q_ij]; // (k, c)
            auto Qab = qab_ij[q_ij]; // (a, b)
            // Common intermediate
            // \Gamma_{Q} = B^{Q}_{me} t_{m}^{e}
            double gamma = Qma->vector_dot(T_n_ij_[ij]);

            // J like contributions
            // \overline{F}_{kc} += 2 B^{Q}_{kc} \Gamma_{Q}
            auto Jcont = Qma->clone();
            Jcont->scale(2.0 * gamma);
            Fkc_bar[ij]->add(Jcont);

            // \overline{F}_{ab} += 2 B^{Q}_{ab} \Gamma_{Q}
            Jcont = Qab->clone();
            Jcont->scale(2.0 * gamma);
            Fab_bar[ij]->add(Jcont);

            // K like contributions
            // \overline{F}_{ab} -= B^{Q}_{ke} T_{m}^{e} B^{Q}_{mc}
            auto Kcont = linalg::triplet(Qma, T_n_ij_[ij], Qma, false, true, false); // (k, e) (m, e) (m, c) -> (k, c)
            Fkc_bar[ij]->subtract(Kcont);

            // \overline{F}_{ab} -= B^{Q}_{ae} T_{m}^{e} B^{Q}_{mb}
            Kcont = linalg::triplet(Qab, T_n_ij_[ij], Qma, false, true, false); // (a, e) (m, e) (m, b) -> (a, b)
            Fab_bar[ij]->subtract(Kcont);
        }

        // Partially dress Fai (Jiang Eq. 100)
        if (i == j) {
            // \overline{F}_{ai} = f_{ai} + [2(ai|me) - (ae|mi)] T_{m}^{e}
            // In closed-shell RHF reference, f_{kc} is zero
            Fai_bar[i] = std::make_shared<Matrix>(npno_ij, 1);

            auto Qia = i_Qa_ij_[ij];
            auto Qik = i_Qk_ij_[ij];

            for (int q_i = 0; q_i < naux_ij; ++q_i) {
                auto Qma = qma_ij[q_i];
                double gamma = Qma->vector_dot(T_n_ij_[ij]); // (m, e) (m, e)

                auto Qab = qab_ij[q_i];
                auto lambda = linalg::doublet(Qab, T_n_ij_[ij], false, true); // (a, e) (m, e) -> (a, m)

                for (int a_i = 0; a_i < npno_ij; ++a_i) {
                    // J like contribution
                    (*Fai_bar[i])(a_i, 0) += 2.0 * gamma * (*Qia)(q_i, a_i); // (a, i)
                    for (int k_i = 0; k_i < nlmo_ij; ++k_i) {
                        // K like contribution
                        (*Fai_bar[i])(a_i, 0) -= (*lambda)(a_i, k_i) * (*Qik)(q_i, k_i); // (a, m) (m, i)
                    } // end k_i
                } // end a_i   
            } // end q_i
        } // end i == j

    } // end ij

    // => Step 2: Dressing over the free/non-contracted indices <= //

    Fkj_ = Fij_bar->clone();
    Fkc_.resize(n_lmo_pairs);
    Fai_.resize(naocc);
    Fab_.resize(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];
        int i_ij = lmopair_to_lmos_dense_[ij][i], j_ij = lmopair_to_lmos_dense_[ij][j];
        int ji = ij_to_ji_[ij], jj = i_j_to_ij_[j][j];
        int pair_idx = (i > j) ? ji : ij;

        int naux_ij = lmopair_to_ribfs_[ij].size();
        int nlmo_ij = lmopair_to_lmos_[ij].size();
        int npno_ij = n_pno_[ij];

        // Fully dress Fkj matrices (Jiang Eq. 94)
        // \widetilde{F}_{ij} = \overline{F}_{ij} (initialized earlier) + \overline{F}_{ic} T_{j}^{c}
        int i_jj = lmopair_to_lmos_dense_[jj][i];
        for (int a_jj = 0; a_jj < n_pno_[jj]; ++a_jj) {
            (*Fkj_)(i, j) += (*Fkc_bar[jj])(i_jj, a_jj) * (*T_ia_[j])(a_jj, 0);
        }

        // Fkc matrices (Jiang Eq. 95)
        // (built separately since \overline{F}_{kc} intermediate is NOT built over weak pairs)
        // \widetilde{F}_{ia} = \overline{F}_{ia} = f_{ia} + [2(ia|kc) - (ic|ka)] T_{k}^{c}
        // => L_{ik}^{ac} T_{k}^{c}
        Fkc_[ij] = std::make_shared<Matrix>(npno_ij, 1);

        for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
            int k = lmopair_to_lmos_[ij][k_ij];
            int ik = i_j_to_ij_[i][k], kk = i_j_to_ij_[k][k];

            // Fkc contributions (projection used to convert from domain of ik to ij)
            auto T_k = linalg::doublet(S_PNO(ik, kk), T_ia_[k]);
            Fkc_[ij]->add(linalg::triplet(S_PNO(ij, ik), L_iajb_[ik], T_k)); // (a, c) (c, 1) -> (a, 1)
        }

        // The rest of these quantities need only to be computed over strong pairs
        // i <= j is only needed due to the equivalence of PNO domains ij and ji
        if (i_j_to_ij_strong_[i][j] == -1 || i > j) continue;

        // Fully dress Fab matrices (Jiang Eq. 97)
        // \widetilde{F}_{ab} = \overline{F}_{ab} - T_{k}^{a} \overline{F}_{kb} 
        Fab_[ij] = Fab_bar[ij]->clone(); // (a, b)
        Fab_[ij]->subtract(linalg::doublet(T_n_ij_[ij], Fkc_bar[ij], true, false)); // (k, a) (k, b) -> (a, b)

        // Fully dress Fai matrices (Jiang Eq. 96)
        // \widetilde{F}_{ai} = \overline{F}_{ai} - T_{k}^{a} \overline{F}_{ki}
        // + \overline{F}_{ab} T_{i}^{b} - T_{k}^{a} \overline{F}_{kb} T_{i}^{b}
        if (i == j) {
            Fai_[i] = Fai_bar[i]->clone(); // (a, 1)
            Fai_[i]->add(linalg::doublet(Fab_bar[ij], T_ia_[i], false, false)); // (a, b) (i, b)
            Fai_[i]->subtract(linalg::triplet(T_n_ij_[ij], Fkc_bar[ij], T_ia_[i], true, false, false)); // (k, a) (k, b) (i, b)

            for (int a_i = 0; a_i < npno_ij; ++a_i) {
                for (int k_i = 0; k_i < nlmo_ij; ++k_i) {
                    int k = lmopair_to_lmos_[ij][k_i];
                    (*Fai_[i])(a_i, 0) -= (*T_n_ij_[ij])(k_i, a_i) * (*Fij_bar)(k, i); // (k, a) (k, i)
                } // end k_i
            } // end a_i
            
        } // end i == j
    }

    timer_off("DLPNO-CCSD: T1 Fock");
}

std::vector<SharedMatrix> DLPNOCCSD::compute_beta() {
    timer_on("DLPNO-CCSD: beta");

    int naocc = nalpha_ - nfrzc();
    int n_lmo_pairs = ij_to_i_j_.size();

    // Jiang Eq. 82
    // \beta_{ij}^{kl} = \widetilde{B}^{Q}_{ki} \widetilde{B}^{Q}_{lj} + t_{ij}^{cd} B^{Q}_{kc} B^{Q}_{ld}
    std::vector<SharedMatrix> beta(n_lmo_pairs); // Stored as n_lmo_pairs * (nlmo_ij, nlmo_ij)

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];
        int ji = ij_to_ji_[ij];

        // This only needs to be computed over strong pairs
        if (i_j_to_ij_strong_[i][j] == -1) continue;

        int naux_ij = lmopair_to_ribfs_[ij].size();
        int nlmo_ij = lmopair_to_lmos_[ij].size();
        int pair_idx = (i > j) ? ji : ij;
        int i_ij = lmopair_to_lmos_dense_[ij][i], j_ij = lmopair_to_lmos_dense_[ij][j];

        // Jiang Eq. 82a
        beta[ij] = linalg::doublet(i_Qk_t1_[ij], i_Qk_t1_[ji], true, false); // (Q, k) (Q, l) -> (k, l)

        // Jiang Eq. 82b
        auto qma_ij = QIA_PNO(ij);
        for (int q_ij = 0; q_ij < naux_ij; ++q_ij) {
            beta[ij]->add(linalg::triplet(qma_ij[q_ij], T_iajb_[ij], qma_ij[q_ij], false, false, true)); // (k, c) (c, d) (l, d) -> (k, l)
        }
    }

    timer_off("DLPNO-CCSD: beta");

    return beta;
}

std::vector<SharedMatrix> DLPNOCCSD::compute_gamma() {
    
    timer_on("DLPNO-CCSD: gamma");

    int naocc = nalpha_ - nfrzc();
    int n_lmo_pairs = ij_to_i_j_.size();

    // Jiang Eq. 83
    // \gamma_{ki}^{ac} = -t_{l}^{a} (ki|lc) + t_{i}^{b} (kb|ac) - t_{i}^{b} (kb|lc) t_{l}^{a} - 0.5 t_{li}^{ad} (kd|lc)
    std::vector<SharedMatrix> gamma(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ki = 0; ki < n_lmo_pairs; ++ki) {
        auto &[k, i] = ij_to_i_j_[ki];
        int ii = i_j_to_ij_[i][i];

        int naux_ki = lmopair_to_ribfs_[ki].size();
        int nlmo_ki = lmopair_to_lmos_[ki].size();
        int npno_ki = n_pno_[ki];
        int k_ki = lmopair_to_lmos_dense_[ki][k], i_ki = lmopair_to_lmos_dense_[ki][i];
        int pair_idx = (k > i) ? ij_to_ji_[ki] : ki;

        gamma[ki] = std::make_shared<Matrix>(npno_ki, npno_ki);
        gamma[ki]->zero();

        // Jiang Eq. 83a \gamma_{ki}^{ac} -= t_{l}^{a} (ki|lc)
        gamma[ki]->subtract(linalg::doublet(T_n_ij_[ki], J_ijmb_[ki], true, false)); // (l, a) (l, c)

        // Jiang Eq. 83b \gamma_{ki}^{ac} += t_{i}^{b} (kb|ac)
        auto T_i = linalg::doublet(S_PNO(ki, ii), T_ia_[i]); // Project from PNO space ii to PNO space ki
        auto K_temp = linalg::doublet(T_i, K_ivvv_[ki], true, false); // (i, b) (b, a * c) -> (a * c)
        K_temp->reshape(n_pno_[ki], n_pno_[ki]); // (a * c) -> (a, c)
        gamma[ki]->add(K_temp);

        // Jiang Eq. 83c \gamma_{ki}^{ac} -= t_{i}^{b} (kb|lc) t_{l}^{a}
        for (int l_ki = 0; l_ki < nlmo_ki; ++l_ki) {
            int l = lmopair_to_lmos_[ki][l_ki];
            int kl = i_j_to_ij_[k][l], ll = i_j_to_ij_[l][l];

            auto T_l = linalg::doublet(S_PNO(ki, ll), T_ia_[l]); // Project from PNO space ll to PNO space ki
            auto T_i_kl = linalg::doublet(S_PNO(kl, ii), T_ia_[i]); // Project from PNO space ii to PNO space kl

            // (b, c) (b, 1) -> (c, 1)... c is also projected to the PNO space of ki
            auto K_kl = linalg::triplet(S_PNO(ki, kl), K_iajb_[kl], T_i_kl, false, true, false);

            // GERRRR! This performs the outer product (a) (c) -> (a, c)
            C_DGER(npno_ki, npno_ki, -1.0, T_l->get_pointer(), 1, K_kl->get_pointer(), 1, gamma[ki]->get_pointer(), npno_ki);
        }
        
        // Jiang Eq. 83d \gamma_{ki}^{ac} -= 0.5 t_{li}^{ad} (kd|lc)
        for (int l_ki = 0; l_ki < nlmo_ki; ++l_ki) {
            int l = lmopair_to_lmos_[ki][l_ki];
            int li = i_j_to_ij_[l][i], kl = i_j_to_ij_[k][l];

            // (a_{li}, d_{li}) (d_{li}, d_{kl}) (d_{kl}, c_{kl}) -> (a_{li}, c_{kl})
            auto gamma_temp = linalg::triplet(T_iajb_[li], S_PNO(li, kl), K_iajb_[kl]);

            // (a_{ki}, a_{li}) (a_{li}, c_{kl}) (c_{kl}, c_{ki}) -> (a_{ki}, c_{ki})
            gamma_temp = linalg::triplet(S_PNO(ki, li), gamma_temp, S_PNO(kl, ki));
            gamma_temp->scale(0.5);

            gamma[ki]->subtract(gamma_temp);
        }
    }

    timer_off("DLPNO-CCSD: gamma");

    return gamma;
}

std::vector<SharedMatrix> DLPNOCCSD::compute_delta() {

    timer_on("DLPNO-CCSD: delta");

    int naocc = nalpha_ - nfrzc();
    int n_lmo_pairs = ij_to_i_j_.size();

    // Jiang Eq. 84 (This term is a doozy)
    // \delta_{ik}^{ac} = -t_{l}^{a} [2(il|kc) - (ik|lc)] + t_{i}^{b} [2(kc|ab) - (kb|ca)]
    // - t_{i}^{b} [2(lb|kc) - (lc|kb)] t_{l}^{a} + 0.5 u_{il}^{ad} [2(kc|ld) - (kd|lc)]
    std::vector<SharedMatrix> delta(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ik = 0; ik < n_lmo_pairs; ++ik) {
        auto &[i, k] = ij_to_i_j_[ik];
        int ii = i_j_to_ij_[i][i], ki = ij_to_ji_[ik];

        int naux_ik = lmopair_to_ribfs_[ik].size();
        int nlmo_ik = lmopair_to_lmos_[ik].size();
        int npno_ik = n_pno_[ik];
        int i_ik = lmopair_to_lmos_dense_[ik][i], k_ik = lmopair_to_lmos_dense_[ik][k];
        int pair_idx = (i > k) ? ki : ik;

        delta[ik] = std::make_shared<Matrix>(npno_ik, npno_ik);
        delta[ik]->zero();

        // Jiang Eq. 84a \delta_{ik}^{ac} -= t_{l}^{a} [2(il|kc) - (ik|lc)]
        auto M_iklc = K_mibj_[ik]->clone(); // (i l_{ik} | k c_{ik})
        M_iklc->scale(2.0);
        M_iklc->subtract(J_ijmb_[ik]); // (i k | l_{ik} c_{ik})
        delta[ik]->subtract(linalg::doublet(T_n_ij_[ik], M_iklc, true, false)); // (l, a) (l, c) -> (a, c)

        // Jiang Eq. 84b \delta_{ik}^{ac} += t_{i}^{b} [2(kc|ab) - (kb|ca)], this is further

        // \delta_{ik}^{ac} += 2.0 t_{i}^{b} (kc|ab)
        auto T_i = linalg::doublet(S_PNO(ik, ii), T_ia_[i]); // (i, b_{ii}) -> (i, b_{ik})
        auto L_temp = K_ivvv_[ki]->clone(); // (k c_{ik} | a_{ik} b_{ik}), stored as (c_{ik}, a_{ik} * b_{ik})
        L_temp->scale(2.0);
        L_temp->reshape(npno_ik * npno_ik, npno_ik); // (c_{ik}, a_{ik} * b_{ik}) -> (c_{ik} * a_{ik}, b_{ik})
        L_temp = linalg::doublet(L_temp, T_i); // (c_{ik} * a_{ik}, b_{ik}) (i, b_{ik}) -> (c_{ik} * a_{ik}, 1)
        L_temp->reshape(npno_ik, npno_ik); // (c_{ik} * a_{ik}, 1) -> (c_{ik}, a_{ik})
        delta[ik]->add(L_temp->transpose()); // (c_{ik}, a_{ik}) -> (a_{ik}, c_{ik})

        // \delta_{ik}^{ac} -= t_{i}^{b} (kb|ca) => t_{i}^{b} (kb|ac) (by symmetry)
        L_temp = K_ivvv_[ki]->clone(); // (k b_{ik} | a_{ik} c_{ik}), stored as (b_{ik}, a_{ik} * c_{ik})
        L_temp = linalg::doublet(T_i, L_temp, true, false); // (b_{ik}, 1) (b_{ik}, a_{ik} * c_{ik}) -> (1, a_{ik} * c_{ik})
        L_temp->reshape(npno_ik, npno_ik); // (1, a_{ik} * c_{ik}) -> (a_{ik}, c_{ik})
        delta[ik]->subtract(L_temp);

        // Jiang Eq. 84c \delta_{ik}^{ac} -= t_{l}^{a} [2(lb|kc) - (lc|kb)] t_{i}^{b} 
        // => t_{i}^{b} L_{lk}^{bc} t_{l}^{a}
        for (int l_ik = 0; l_ik < nlmo_ik; ++l_ik) {
            int l = lmopair_to_lmos_[ik][l_ik];
            int ll = i_j_to_ij_[l][l], lk = i_j_to_ij_[l][k];

            auto T_l = linalg::doublet(S_PNO(ik, ll), T_ia_[l]); // (l, a_{ll}) -> (l, a_{ik})
            auto T_i_lk = linalg::doublet(S_PNO(lk, ii), T_ia_[i]); // (i, b_{ii}) -> (i, b_{lk})

             // (b_{lk}, 1) (b_{lk}, c_{lk}) (c_{lk}, c_{ik}) -> (1, c_{ik})
            auto L_lk = linalg::triplet(T_i_lk, L_iajb_[lk], S_PNO(lk, ik), true, false, false);

            // GERRRR! This performs the outer product (a_{ik}) (c_{ik}) -> (a_{ik}, c_{ik})
            C_DGER(npno_ik, npno_ik, -1.0, T_l->get_pointer(), 1, L_lk->get_pointer(), 1, delta[ik]->get_pointer(), npno_ik);
        }

        // Jiang Eq. 84d \delta_{ik}^{ac} += 0.5 u_{il}^{ad} [2(kc|ld) - (kd|lc)]
        // => 0.5 u_{il}^{ad} L_{kl}^{cd} or 0.5 u_{il}^{ad} L_{lk}^{dc}
        for (int l_ik = 0; l_ik < nlmo_ik; ++l_ik) {
            int l = lmopair_to_lmos_[ik][l_ik];
            int il = i_j_to_ij_[i][l], lk = i_j_to_ij_[l][k];

            // (a_{il}, d_{il}) (d_{il}, d_{lk}) (d_{lk}, c_{lk}) -> (a_{il}, c_{lk})
            auto delta_temp = linalg::triplet(Tt_iajb_[il], S_PNO(il, lk), L_iajb_[lk]);

            // (a_{ik}, a_{il}) (a_{il}, c_{lk}) (c_{lk}, c_{ik}) -> (a_{ik}, c_{ik})
            delta_temp = linalg::triplet(S_PNO(ik, il), delta_temp, S_PNO(lk, ik));
            delta_temp->scale(0.5);

            delta[ik]->add(delta_temp);
        }
    }

    timer_off("DLPNO-CCSD: delta");

    return delta;
}

SharedMatrix DLPNOCCSD::compute_Fkj_double_tilde() {

    timer_on("DLPNO-CCSD: Fkj double tilde");

    int naocc = nalpha_ - nfrzc();
    int n_lmo_pairs = ij_to_i_j_.size();

    // Jiang Eq. 86
    // F_{kj}'' = F_{kj}' + u_{lj}^{cd} (lc|kd) (primes represent tildes)
    // this is equivalent to F_{ij}'' = F_{ij}' + u_{lj}^{cd} (lc|id) by replacing k with i
    SharedMatrix Fkj_double_tilde = Fkj_->clone();

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];
        int nlmo_ij = lmopair_to_lmos_[ij].size();

        for (int l_ij = 0; l_ij < nlmo_ij; ++l_ij) {
            int l = lmopair_to_lmos_[ij][l_ij];

            int li = i_j_to_ij_[l][i], lj = i_j_to_ij_[l][j];
            if (li == -1 || lj == -1) continue;

            // F_{ij}'' += u_{lj}^{cd} K_{li}^{cd}

            // (c_{li}, c_{lj}) (c_{lj}, d_{lj}) (d_{lj}, d_{li}) -> (c_{li}, d_{li})
            auto U_lj = linalg::triplet(S_PNO(li, lj), Tt_iajb_[lj], S_PNO(lj, li));

            // (c_{li}, d_{li}) (c_{li}, d_{li}) -> scalar
            (*Fkj_double_tilde)(i, j) += K_iajb_[li]->vector_dot(U_lj);
        } // end l_ij
    } // end ij

    timer_off("DLPNO-CCSD: Fkj double tilde");

    return Fkj_double_tilde;
}

void DLPNOCCSD::compute_R_ia(std::vector<SharedMatrix>& R_ia, std::vector<std::vector<SharedMatrix>>& R_ia_buffer) {
    
    timer_on("DLPNO-CCSD: Compute R1");

    int n_lmo_pairs = ij_to_i_j_.size();
    int naocc = nalpha_ - nfrzc();

    // Thread and OMP Parallel info
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif

    // Initialize R1 residuals, Jiang Eq. 87a
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < naocc; ++i) {
        int ii = i_j_to_ij_[i][i];
        // Initialize Ria as T1-dressed Fock matrix element Fai
        R_ia[i]->copy(Fai_[i]);
    }

    // Zero out buffers
    for (int thread = 0; thread < nthreads; ++thread) {
        for (int i = 0; i < naocc; ++i) {
            R_ia_buffer[thread][i]->zero();
        }
    }

    // Compute residual for singles amplitude (A and C contributions)
#pragma omp parallel for schedule(dynamic, 1)
    for (int ik = 0; ik < n_lmo_pairs; ++ik) {
        auto &[i, k] = ij_to_i_j_[ik];
        int ki = ij_to_ji_[ik];
        int k_ki = lmopair_to_lmos_dense_[ki][k]; // Grabs the index of k within domain ki
        int pair_idx = (i > k) ? ki : ik;

        int nlmo_ik = lmopair_to_lmos_[ik].size();
        int naux_ik = lmopair_to_ribfs_[ik].size();
        int npno_ik = n_pno_[ik];

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        int i_ik = lmopair_to_lmos_dense_[ik][i], k_ik = lmopair_to_lmos_dense_[ik][k];
        std::vector<int> k_ik_slice = std::vector<int>(1, k_ik);
        int ii = i_j_to_ij_[i][i];
        
        // A_{i}^{a} = u_{ki}^{cd} [(kc|da) - t_{l}^{a}(ld|kc)] (Jiang Eq. 88)

        // A_{i}^{a} += u_{ki}^{cd} (kc|da)
        auto K_kcad = K_ivvv_[ki]->clone(); // (k c_{ki} | d_{ki} a_{ki}), stored as (c_{ki}, d_{ki} * a_{ki})
        K_kcad->reshape(n_pno_[ki] * n_pno_[ki], n_pno_[ki]); // (c_{ki}, d_{ki} * a_{ki}) -> (c_{ki} * d_{ki}, a_{ki})
        auto Uki = Tt_iajb_[ki]->clone();   // (c_{ki}, d_{ki})
        Uki->reshape(npno_ik * npno_ik, 1); // (c_{ki}, d_{ki}) -> (c_{ki} * d_{ki}, 1)

        // (a_{ki}, a_{ii}) (c_{ki} * d_{ki}, a_{ki}) (c_{ki} * d_{ki}, 1) -> (a_{ii})
        R_ia_buffer[thread][i]->add(linalg::triplet(S_PNO(ki, ii), K_kcad, Uki, true, true, false));

        // A_{i}^{a} -= u_{ki}^{cd} (kc|ld) t_{l}^{a}
        for (int l_ik = 0; l_ik < nlmo_ik; ++l_ik) {
            int l = lmopair_to_lmos_[ik][l_ik];
            int kl = i_j_to_ij_[k][l], ll = i_j_to_ij_[l][l];

            auto T_l = linalg::doublet(S_PNO(ii, ll), T_ia_[l]); // (l, a_{ll}) -> (l, a_{ii})

            // (c_{kl}, c_{ki}) (c_{ki}, d_{ki}) (d_{ki}, d_{kl}) -> (c_{kl}, d_{kl})
            auto U_ki = linalg::triplet(S_PNO(kl, ki), Tt_iajb_[ki], S_PNO(ki, kl));

            // (c_{kl}, d_{kl}) * (c_{kl}, d_{kl}) * (a_{ii}) -> (a_{ii})
            T_l->scale(U_ki->vector_dot(K_iajb_[kl]));
            R_ia_buffer[thread][i]->subtract(T_l);
        }

        // C_{i}^{a} = F_{kc}U_{ik}^{ac} (Jiang Eq. 90)
        // (a_{ik}, a_{ii}) (a_{ik}, c_{ik}) (k, c_{ik}) -> (a_{ii})
        R_ia_buffer[thread][i]->add(linalg::triplet(S_PNO(ik, ii), Tt_iajb_[ik], Fkc_[ki], true, false, false));
    } // end ki

    // B contributions
#pragma omp parallel for schedule(dynamic, 1)
    for (int kl = 0; kl < n_lmo_pairs; ++kl) {
        auto &[k, l] = ij_to_i_j_[kl];

        int naux_kl = lmopair_to_ribfs_[kl].size();
        int nlmo_kl = lmopair_to_lmos_[kl].size();
        int npno_kl = n_pno_[kl];
        int pair_idx = (k > l) ? ij_to_ji_[kl] : kl;
        int k_kl = lmopair_to_lmos_dense_[kl][k], l_kl = lmopair_to_lmos_dense_[kl][l];

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        // B_{i}^{a} = -u_{kl}^{ac}[(ki|lc) + t_{i}^{b} (kb|lc)] (Jiang Eq. 89)
        auto K_kilc = K_mibj_[kl]->clone(); // (k i_{kl} | l c_{kl}), stored as (i_{kl}, c_{kl})
        K_kilc->add(linalg::doublet(T_n_ij_[kl], K_iajb_[kl])); // (i_{kl}, b_{kl}) (b_{kl}, c_{kl}) -> (i_{kl}, c_{kl})
        auto B_ia = linalg::doublet(K_kilc, Tt_iajb_[kl], false, true); // (i_{kl}, c_{kl}) (a_{kl}, c_{kl})  -> (i_{kl}, a_{kl})

        // Flush results into buffer
        for (int i_kl = 0; i_kl < nlmo_kl; ++i_kl) {
            int i = lmopair_to_lmos_[kl][i_kl];
            int ii = i_j_to_ij_[i][i];
            std::vector<int> i_kl_slice(1, i_kl);

            // (i_{kl}, a_{kl}) (a_{kl}, a_{ii}) -> (a_{ii})
            auto B_ia_proj = linalg::doublet(submatrix_rows(*B_ia, i_kl_slice), S_PNO(kl, ii), false, false);
            R_ia_buffer[thread][i]->subtract(B_ia_proj->transpose());
        } // end i_kl
    } // end kl

    // Add R_ia buffers to R_ia
    for (int i = 0; i < naocc; ++i) {
        for (int thread = 0; thread < nthreads; ++thread) {
            R_ia[i]->add(R_ia_buffer[thread][i]);
        }
    }

    timer_off("DLPNO-CCSD: Compute R1");
}

void DLPNOCCSD::compute_R_iajb(std::vector<SharedMatrix>& R_iajb, std::vector<SharedMatrix>& Rn_iajb) {

    timer_on("DLPNO-CCSD: Compute R2");

    int n_lmo_pairs = ij_to_i_j_.size();
    int naocc = nalpha_ - nfrzc();

    // Thread and OMP Parallel info
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif

    // Compute intermediates to make our lives easier
    auto beta = compute_beta(); // n_lmo_pairs * (nlmo_ij, nlmo_ij)
    auto gamma = compute_gamma(); // n_lmo_pairs * (nlmo_ij, npno_ij)
    auto delta = compute_delta(); // n_lmo_pairs * (nlmo_ij, npno_ij)
    auto Fkj_double_tilde = compute_Fkj_double_tilde(); // (nlmo, nlmo)

    // Zero out residuals
#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        R_iajb[ij]->zero();
        Rn_iajb[ij]->zero();
    }

    // Sort pairs by the approximate number of operations (for maximal parallel efficiency)
    std::vector<std::pair<int, size_t>> ij_cost_tuple(n_lmo_pairs);
    
#pragma omp parallel for
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];

        const int npao_ij = lmopair_to_paos_[ij].size();
        const int naux_ij = lmopair_to_ribfs_[ij].size();
        const int nlmo_ij = lmopair_to_lmos_[ij].size();
        
        size_t cost = 0;

        if (i <= j) {
            // Cost of PNO projections
            for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
                int k = lmopair_to_lmos_[ij][k_ij];
                for (int l_ij = 0; l_ij < nlmo_ij; ++l_ij) {
                    int l = lmopair_to_lmos_[ij][l_ij];
                    int kl = i_j_to_ij_[k][l];
                    if (kl == -1) continue;

                    cost += n_pno_[ij] * n_pno_[ij] * n_pno_[kl];
                    cost += n_pno_[ij] * n_pno_[kl] * n_pno_[kl];
                } // end l_ij
            } // end k_ij
        } else {
            // Cost of PNO projections
            for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
                int k = lmopair_to_lmos_[ij][k_ij];
                int ik = i_j_to_ij_[i][k], jk = i_j_to_ij_[j][k];

                cost += n_pno_[ij] * n_pno_[ij] * n_pno_[ik];
                cost += n_pno_[ij] * n_pno_[ik] * n_pno_[ik];
                cost += n_pno_[ij] * n_pno_[ij] * n_pno_[jk];
                cost += n_pno_[ij] * n_pno_[jk] * n_pno_[jk];
            } // end k_ij
        } // end else

        ij_cost_tuple[ij] = std::make_pair(ij, cost);
    }
    
    std::sort(ij_cost_tuple.begin(), ij_cost_tuple.end(), [&](const std::pair<int, size_t>& a, const std::pair<int, size_t>& b) {
        return (a.second > b.second);
    });

    std::vector<int> ij_sorted_by_cost(n_lmo_pairs);
    
#pragma omp parallel for
    for (int ij_idx = 0; ij_idx < n_lmo_pairs; ++ij_idx) {
        ij_sorted_by_cost[ij_idx] = ij_cost_tuple[ij_idx].first;
    }

    // Compute residual for doubles amplitude
#pragma omp parallel for schedule(dynamic, 1)
    for (int ij_idx = 0; ij_idx < n_lmo_pairs; ++ij_idx) {
        int ij = ij_sorted_by_cost[ij_idx];

        auto &[i, j] = ij_to_i_j_[ij];
        bool is_weak_pair = (i_j_to_ij_strong_[i][j] == -1);
        int ji = ij_to_ji_[ij];

        int nlmo_ij = lmopair_to_lmos_[ij].size();
        int naux_ij = lmopair_to_ribfs_[ij].size();
        int npno_ij = n_pno_[ij];

        // Skip if this pair is a "weak pair"
        if (is_weak_pair) continue;

        int pair_idx = (i > j) ? ji : ij;

        if (i <= j) {
            // R_{ij}^{ab} += \widetilde{B}^{Q}_{ai} \widetilde{B}^{Q}_{bj} (Jiang Eq. 75)
            auto K_ij = linalg::doublet(i_Qa_t1_[ij], i_Qa_t1_[ji], true, false); // (Q, a) (Q, b) -> (a, b)
            R_iajb[ij]->add(K_ij);
            if (i != j) R_iajb[ji]->add(K_ij->transpose());

            // A_{ij}^{ab} = \widetilde{B}^{Q}_{ac} * t_{ij}^{cd} * \widetilde{B}^{Q}_{bd} (Jiang Eq. 76)
            auto A_ij = std::make_shared<Matrix>(npno_ij, npno_ij);
            A_ij->zero();

            auto qma_ij = QIA_PNO(ij); // naux_ij * (nlmo_ij, npno_ij)
            auto qab_ij = QAB_PNO(ij); // naux_ij * (npno_ij, npno_ij)
            for (int q_ij = 0; q_ij < naux_ij; ++q_ij) {
                // This performs the T1-dressing of Qab on the fly, as this intermeidate is only used once
                // \widetilde{B}^{Q}_{ab} = B^{Q}_{ab} - t_{k}^{a} B^{Q}_{kb} (Jiang Eq. 93)
                auto Qab_t1 = qab_ij[q_ij]->clone(); // (a, b)
                Qab_t1->subtract(linalg::doublet(T_n_ij_[ij], qma_ij[q_ij], true, false)); // (k, a) (k, b) -> (a, b)
                
                // A_{ij}^{ab} = \widetilde{B}^{Q}_{ac} * t_{ij}^{cd} * \widetilde{B}^{Q}_{bd} (Jiang Eq. 76)
                A_ij->add(linalg::triplet(Qab_t1, T_iajb_[ij], Qab_t1, false, false, true)); // (a, c) (c, d) (b, d)
            } // end q_ij
            R_iajb[ij]->add(A_ij);
            if (i != j) R_iajb[ji]->add(A_ij->transpose());

            // This intermediate is computed if the "semi-direct" algorithm is used for the PNO overlap matrices
            // The rows of this intermediate contains an extended PAO domain of ij, while the columns are transformed
            // to the PNO space of ij
            SharedMatrix S_ij;
            std::vector<int> pair_ext_domain;
            if (low_memory_overlap_) {
                for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
                    int k = lmopair_to_lmos_[ij][k_ij];
                    for (int l_ij = 0; l_ij < nlmo_ij; ++l_ij) {
                        int l = lmopair_to_lmos_[ij][l_ij];
                        int kl = i_j_to_ij_[k][l];
                        if (kl == -1 || n_pno_[kl] == 0) continue;
                        pair_ext_domain = merge_lists(pair_ext_domain, lmopair_to_paos_[kl]);
                    } // end k
                } // end l
                S_ij = submatrix_rows_and_cols(*S_pao_, pair_ext_domain, lmopair_to_paos_[ij]);
                S_ij = linalg::doublet(S_ij, X_pno_[ij], false, false);
            } // end if

            // => These two intermediates involve the expensive S(a_{kl}, a_{ij}) PNO projection matrices

            // B_{ij}^{ab} = t_{kl}^{ab} * \beta_{ij}^{kl} (Jiang Eq. 77)
            auto B_ij = std::make_shared<Matrix>(npno_ij, npno_ij);
            B_ij->zero();
            
            // F_{bc}'' = F_{bc}' - u_{kl}^{bd} K_{kl}^{cd} (Jiang Eq. 85)
            auto F_bc_double_tilde = Fab_[ij]->clone();

            for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
                int k = lmopair_to_lmos_[ij][k_ij];
                for (int l_ij = 0; l_ij < nlmo_ij; ++l_ij) {
                    int l = lmopair_to_lmos_[ij][l_ij];
                    int kl = i_j_to_ij_[k][l];
                    if (kl == -1 || n_pno_[kl] == 0) continue;

                    SharedMatrix S_kl_ij = (low_memory_overlap_) ? 
                            linalg::doublet(X_pno_[kl], submatrix_rows(*S_ij, index_list(pair_ext_domain, lmopair_to_paos_[kl])), true, false) : S_PNO(kl, ij);

                    // B contributions
                    // (a_{kl}, b_{kl}) -> (a_{ij}, b_{ij})
                    auto T_kl = linalg::triplet(S_kl_ij, T_iajb_[kl], S_kl_ij, true, false, false);
                    // \beta_{ij}^{kl} T_{kl}^{ab}
                    T_kl->scale((*beta[ij])(k_ij, l_ij));
                    B_ij->add(T_kl);

                    // F_double_tilde contributions
                    auto E_temp = linalg::doublet(Tt_iajb_[kl], K_iajb_[kl], false, true); // (b, d) (c, d) -> (b, c)

                    // (b_{kl}, c_{kl}) -> (b_{ij}, c_{ij})
                    F_bc_double_tilde->subtract(linalg::triplet(S_kl_ij, E_temp, S_kl_ij, true, false, false));
                } // end l_ij
            } // end k_ij
            R_iajb[ij]->add(B_ij);
            if (i != j) R_iajb[ji]->add(B_ij->transpose());

            // E_{ij}^{ab} = t_{ij}^{ac} F_{bc}'' (Jiang Eq. 80)
            // For the residual contribution, this needs to be symmetrized
            // P_{ij}^{ab} t_{ij}^{ac} F_{bc}'' => t_{ij}^{ac} F_{bc}'' + F_{ac}'' t_{ij}^{cb}
            auto E_ij = linalg::doublet(T_iajb_[ij], F_bc_double_tilde, false, true);
            E_ij->add(linalg::doublet(F_bc_double_tilde, T_iajb_[ij], false, false));

            R_iajb[ij]->add(E_ij);
            if (i != j) R_iajb[ji]->add(E_ij->transpose());
        } // end if
        
        // C_{ij}^{ab} = [-\gamma_{ki}^{ac} - (i k | a_{ij} c_{kj})] t_{kj}^{bc} (Jiang Eq. 78)
        auto C_ij = std::make_shared<Matrix>(npno_ij, npno_ij);
        C_ij->zero();

        for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
            int k = lmopair_to_lmos_[ij][k_ij];
            int ki = i_j_to_ij_[k][i], kj = i_j_to_ij_[k][j];
            
            auto gamma_total = J_ikac_non_proj_[ij][k_ij]->clone(); // (i k | a_{ij} c_{kj})
            gamma_total->add(linalg::triplet(S_PNO(ij, ki), gamma[ki], S_PNO(ki, kj))); // (a_{ki}, c_{ki}) -> (a_{ij}, c_{kj})

            // (a_{ij}, c_{kj}) (b_{kj}, c_{kj}) (b_{kj}, b_{ij}) -> (a_{ij}, b_{ij})
            C_ij->subtract(linalg::triplet(gamma_total, T_iajb_[kj], S_PNO(kj, ij), false, true, false));
        }
        // Add all the C terms to the non-symmetrized R buffer
        // P_{ij}^{ab} [0.5 C_{ij}^{ab} + C_{ij}^{ba}] (Jiang Eq. 19)
        auto C_ij_total = C_ij->clone();
        C_ij_total->scale(0.5);
        C_ij_total->add(C_ij->transpose());
        Rn_iajb[ij]->add(C_ij_total);

        // D_{ij}^{ab} = 0.5 [\delta_{ik}^{ac} + 2 (i a_{ij} | k c_{jk}) - (i k | a_{ij} c_{jk})] u_{jk}^{bc}  (Jiang Eq. 79)
        auto D_ij = std::make_shared<Matrix>(npno_ij, npno_ij);
        D_ij->zero();

        for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
            int k = lmopair_to_lmos_[ij][k_ij];
            int ik = i_j_to_ij_[i][k], jk = i_j_to_ij_[j][k];
            int j_ik = lmopair_to_lmos_dense_[ik][j];

            // 2 (i a_{ij} | k c_{jk}) - (i k | a_{ij} c_{jk})
            auto delta_total = K_iakc_non_proj_[ij][k_ij]->clone();
            delta_total->scale(2.0);
            delta_total->subtract(J_ikac_non_proj_[ij][k_ij]);

            // (a_{ik}, c_{ik}) -> (a_{ij}, c_{jk})
            delta_total->add(linalg::triplet(S_PNO(ij, ik), delta[ik], S_PNO(ik, jk)));
            
            // (a_{ij}, c_{jk}) (b_{jk} c_{jk}) (b_{jk}, b_{ij}) -> (a_{ij}, b_{ij})
            D_ij->add(linalg::triplet(delta_total, Tt_iajb_[jk], S_PNO(jk, ij), false, true, false));
        }
        D_ij->scale(0.5);
        Rn_iajb[ij]->add(D_ij);

        // G_{ij}^{ab} = -t_{ik}^{ab} F_{kj}'' (Jiang Eq. 81)
        auto G_ij = std::make_shared<Matrix>(npno_ij, npno_ij);
        G_ij->zero();

        for (int k = 0; k < naocc; ++k) {
            int ik = i_j_to_ij_[i][k];
            if (ik == -1) continue;

            // (a_{ik}, b_{ik}) -> (a_{ij}, b_{ij})
            auto T_ik = linalg::triplet(S_PNO(ij, ik), T_iajb_[ik], S_PNO(ik, ij), false, false, false);
            T_ik->scale((*Fkj_double_tilde)(k, j));
            G_ij->subtract(T_ik);
        }
        Rn_iajb[ij]->add(G_ij);
    } // end ij

    // Add the contributions from projection corrections
#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        int i, j;
        std::tie(i, j) = ij_to_i_j_[ij];
        int ji = ij_to_ji_[ij];

        if (i_j_to_ij_strong_[i][j] != -1) {
            R_iajb[ij]->add(Rn_iajb[ij]);
            R_iajb[ij]->add(Rn_iajb[ji]->transpose());
        } else {
            R_iajb[ij] = std::make_shared<Matrix>(n_pno_[ij], n_pno_[ij]);
            R_iajb[ij]->zero();
        }
    }

    timer_off("DLPNO-CCSD: Compute R2");
}

void DLPNOCCSD::lccsd_iterations() {

    int n_lmo_pairs = ij_to_i_j_.size();
    int naocc = nalpha_ - nfrzc();

    // Thread and OMP Parallel info
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif

    outfile->Printf("\n  ==> Local CCSD <==\n\n");
    outfile->Printf("    E_CONVERGENCE = %.2e\n", options_.get_double("E_CONVERGENCE"));
    outfile->Printf("    R_CONVERGENCE = %.2e\n\n", options_.get_double("R_CONVERGENCE"));
    outfile->Printf("                      Corr. Energy    Delta E     Max R1     Max R2     Time (s)\n");

    // => Initialize Residuals and Amplitudes <= //

    std::vector<SharedMatrix> R_ia(naocc);
    std::vector<SharedMatrix> Rn_iajb(n_lmo_pairs);
    std::vector<SharedMatrix> R_iajb(n_lmo_pairs);

    // => Initialize Singles Residuals and Amplitudes <= //

    T_ia_.resize(naocc);

#pragma omp parallel for
    for (int i = 0; i < naocc; ++i) {
        int ii = i_j_to_ij_[i][i];
        T_ia_[i] = std::make_shared<Matrix>(n_pno_[ii], 1);
        R_ia[i] = std::make_shared<Matrix>(n_pno_[ii], 1);
    }

    // => Initialize Doubles Residuals and Amplitudes <= //
    
#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        R_iajb[ij] = std::make_shared<Matrix>(n_pno_[ij], n_pno_[ij]);
        Rn_iajb[ij] = std::make_shared<Matrix>(n_pno_[ij], n_pno_[ij]);
    }

    // => Thread buffers <= //

    std::vector<std::vector<SharedMatrix>> R_ia_buffer(nthreads);
    for (int thread = 0; thread < nthreads; ++thread) {
        R_ia_buffer[thread].resize(naocc);
        for (int i = 0; i < naocc; ++i) {
            int ii = i_j_to_ij_[i][i];
            R_ia_buffer[thread][i] = std::make_shared<Matrix>(n_pno_[ii], 1);
        }
    }

    int iteration = 1, max_iteration = options_.get_int("DLPNO_MAXITER");
    double e_curr = 0.0, e_prev = 0.0, e_weak = 0.0, r1_curr = 0.0, r2_curr = 0.0;
    bool e_converged = false, r_converged = false;

    DIISManager diis(options_.get_int("DIIS_MAX_VECS"), "LCCSD DIIS", DIISManager::RemovalPolicy::LargestError, DIISManager::StoragePolicy::InCore);

    double F_CUT = options_.get_double("F_CUT");

    i_Qk_t1_.resize(n_lmo_pairs);
    i_Qa_t1_.resize(n_lmo_pairs);

    T_n_ij_.resize(n_lmo_pairs);

    while (!(e_converged && r_converged)) {
        // RMS of residual per single LMO, for assesing convergence
        std::vector<double> R_ia_rms(naocc, 0.0);
        // RMS of residual per LMO pair, for assessing convergence
        std::vector<double> R_iajb_rms(n_lmo_pairs, 0.0);

        std::time_t time_start = std::time(nullptr);

        // Step 1: Create T_n intermediate (Jiang Eq. 70)
        // T_{n_{ij}}^{a_{ij}} = S(a_{ij}, a_{nn}) T_{n}^{a_{nn}}
        // n_{ij} is all n such that in and jn form valid pairs
#pragma omp parallel for schedule(dynamic, 1)
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            auto &[i, j] = ij_to_i_j_[ij];

            int nlmo_ij = lmopair_to_lmos_[ij].size();
            int npno_ij = n_pno_[ij];
            
            T_n_ij_[ij] = std::make_shared<Matrix>(nlmo_ij, npno_ij);

            for (int n_ij = 0; n_ij < nlmo_ij; ++n_ij) {
                int n = lmopair_to_lmos_[ij][n_ij];
                int nn = i_j_to_ij_[n][n];

                // (a_{ij}, a_{nn}) (a_{nn}, 1) -> (a_{ij}, 1)
                auto T_n_temp = linalg::doublet(S_PNO(ij, nn), T_ia_[n], false, false);
                
                for (int a_ij = 0; a_ij < npno_ij; ++a_ij) {
                    (*T_n_ij_[ij])(n_ij, a_ij) = (*T_n_temp)(a_ij, 0);
                } // end a_ij
            } // end n_ij
        } // end ij

        // Step 2: T1-dress integrals and Fock matrices
        t1_ints();
        t1_fock();

        // Step 3: Compute R1 residual
        compute_R_ia(R_ia, R_ia_buffer);

        // Get rms of R_ia
#pragma omp parallel for schedule(dynamic, 1)
        for (int i = 0; i < naocc; ++i) {
            R_ia_rms[i] = R_ia[i]->rms();
        }

        // Step 4: Compute R2 residual
        compute_R_iajb(R_iajb, Rn_iajb);

        // Get rms of R_iajb
#pragma omp parallel for schedule(dynamic, 1)
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            R_iajb_rms[ij] = R_iajb[ij]->rms();
        }

        // Update Singles Amplitude (Jiang Eq. 103)
#pragma omp parallel for
        for (int i = 0; i < naocc; ++i) {
            int ii = i_j_to_ij_[i][i];
            for (int a_ii = 0; a_ii < n_pno_[ii]; ++a_ii) {
                (*T_ia_[i])(a_ii, 0) -= (*R_ia[i])(a_ii, 0) / (e_pno_[ii]->get(a_ii) - F_lmo_->get(i,i));
            }
        }

        // Update Doubles Amplitude (Jiang Eq. 104)
#pragma omp parallel for schedule(dynamic, 1)
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            auto &[i, j] = ij_to_i_j_[ij];

            for (int a_ij = 0; a_ij < n_pno_[ij]; ++a_ij) {
                for (int b_ij = 0; b_ij < n_pno_[ij]; ++b_ij) {
                    (*T_iajb_[ij])(a_ij, b_ij) -= (*R_iajb[ij])(a_ij, b_ij) / 
                                    (e_pno_[ij]->get(a_ij) + e_pno_[ij]->get(b_ij) - F_lmo_->get(i,i) - F_lmo_->get(j,j));
                }
            }
            
        }

        // DIIS Extrapolation
        std::vector<SharedMatrix> T_vecs;
        T_vecs.reserve(T_ia_.size() + T_iajb_.size());
        T_vecs.insert(T_vecs.end(), T_ia_.begin(), T_ia_.end());
        T_vecs.insert(T_vecs.end(), T_iajb_.begin(), T_iajb_.end());

        std::vector<SharedMatrix> R_vecs;
        R_vecs.reserve(R_ia.size() + R_iajb.size());
        R_vecs.insert(R_vecs.end(), R_ia.begin(), R_ia.end());
        R_vecs.insert(R_vecs.end(), R_iajb.begin(), R_iajb.end());

        auto T_vecs_flat = flatten_mats(T_vecs);
        auto R_vecs_flat = flatten_mats(R_vecs);

        if (iteration == 1) {
            diis.set_error_vector_size(R_vecs_flat);
            diis.set_vector_size(T_vecs_flat);
        }

        diis.add_entry(R_vecs_flat.get(), T_vecs_flat.get());
        diis.extrapolate(T_vecs_flat.get());

        copy_flat_mats(T_vecs_flat, T_vecs);

        // Update symmetrized doubles amplitude
        // Tt_iajb (or U_iajb)
        // u_{ij}^{ab} = 2t_{ij}^{ab} - t_{ij}^{ba}
#pragma omp parallel for schedule(dynamic, 1)
        for (int ij = 0; ij < n_lmo_pairs; ij++) {
            Tt_iajb_[ij] = T_iajb_[ij]->clone();
            Tt_iajb_[ij]->scale(2.0);
            Tt_iajb_[ij]->subtract(T_iajb_[ij]->transpose());
        }

        // evaluate convergence using current amplitudes and residuals
        e_prev = e_curr;

        // Compute LCCSD energy (Jiang Eq. 45)
        // E_{CCSD} = (t_{ij}^{ab} + t_{i}^{a}t_{j}^{b}) L_{ij}^{ab}
        // Note, the singles contribution to weak pairs is still computed!!!
        e_curr = 0.0, e_weak = 0.0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : e_curr, e_weak)
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            auto &[i, j] = ij_to_i_j_[ij];
            int ii = i_j_to_ij_[i][i], jj = i_j_to_ij_[j][j];
            
            auto T_i = linalg::doublet(S_PNO(ij, ii), T_ia_[i]);
            auto T_j = linalg::doublet(S_PNO(ij, jj), T_ia_[j]);

            // \tau_{ij}^{ab} = (t_{ij}^{ab} + t_{i}^{a}t_{j}^{b})
            auto tau = T_iajb_[ij]->clone();

            for (int a_ij = 0; a_ij < n_pno_[ij]; ++a_ij) {
                for (int b_ij = 0; b_ij < n_pno_[ij]; ++b_ij) {
                    (*tau)(a_ij, b_ij) += (*T_i)(a_ij, 0) * (*T_j)(b_ij, 0);
                } // end b_ij
            } // end a_ij

            double e_ij = tau->vector_dot(L_iajb_[ij]);
            
            e_curr += e_ij;
            if (i_j_to_ij_strong_[i][j] == -1) e_weak += e_ij;
        }
        double r_curr1 = *max_element(R_ia_rms.begin(), R_ia_rms.end());
        double r_curr2 = *max_element(R_iajb_rms.begin(), R_iajb_rms.end());

        r_converged = (fabs(r_curr1) < options_.get_double("R_CONVERGENCE"));
        r_converged &= (fabs(r_curr2) < options_.get_double("R_CONVERGENCE"));
        e_converged = (fabs(e_curr - e_prev) < options_.get_double("E_CONVERGENCE"));

        std::time_t time_stop = std::time(nullptr);

        outfile->Printf("  @LCCSD iter %3d: %16.12f %10.3e %10.3e %10.3e %8d\n", iteration, e_curr, e_curr - e_prev, r_curr1, r_curr2, (int)time_stop - (int)time_start);

        iteration++;

        if (iteration > max_iteration + 1) {
            throw PSIEXCEPTION("Maximum DLPNO iterations exceeded.");
        }
    }

    e_lccsd_ = e_curr - e_weak;
    de_weak_ = e_weak;
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

    timer_on("Compute Metric");
    compute_metric();
    timer_off("Compute Metric");

    // Adjust parameters for "crude" prescreening
    T_CUT_MKN_ *= 100;
    T_CUT_DO_ *= 2;

    outfile->Printf("  ==> Crude Prescreening (Determines Semicanonical-MP2 Pairs) <==\n\n");
    outfile->Printf("    T_CUT_MKN set to %6.3e\n", T_CUT_MKN_);
    outfile->Printf("    T_CUT_DO  set to %6.3e\n\n", T_CUT_DO_);

    timer_on("Sparsity");
    prep_sparsity(true, false);
    timer_off("Sparsity");

    timer_on("Crude DF Ints");
    compute_qia();
    timer_off("Crude DF Ints");

    timer_on("Initial Pair Prescreening");
    pair_prescreening<true>();
    timer_off("Initial Pair Prescreening");

    // Reset Sparsity After
    T_CUT_MKN_ *= 0.01;
    T_CUT_DO_ *= 0.5;

    outfile->Printf("  ==> Refined Prescreening (Determines Strong and Weak Pairs) <==\n\n");
    outfile->Printf("    T_CUT_MKN reset to %6.3e\n", T_CUT_MKN_);
    outfile->Printf("    T_CUT_DO  reset to %6.3e\n\n", T_CUT_DO_);

    timer_on("Sparsity");
    prep_sparsity(false, false);
    timer_off("Sparsity");

    timer_on("Refined DF Ints");
    print_integral_sparsity();
    compute_qia();
    timer_off("Refined DF Ints");

    timer_on("Refined Pair Prescreening");
    pair_prescreening<false>();
    timer_off("Refined Pair Prescreening");

    // Set variables from LMP2
    double e_scf = reference_wavefunction_->energy();
    double e_lmp2_corr = e_lmp2_ + de_lmp2_eliminated_ + de_dipole_ + de_pno_total_;
    double e_lmp2_total = e_scf + e_lmp2_corr;

    set_scalar_variable("MP2 CORRELATION ENERGY", e_lmp2_corr);
    set_scalar_variable("CURRENT CORRELATION ENERGY", e_lmp2_corr);
    set_scalar_variable("MP2 TOTAL ENERGY", e_lmp2_total);
    set_scalar_variable("CURRENT ENERGY", e_lmp2_total);

    outfile->Printf("  \n");
    outfile->Printf("  Total DLPNO-MP2 Correlation Energy: %16.12f \n", e_lmp2_ + de_lmp2_eliminated_ + de_pno_total_ + de_dipole_);
    outfile->Printf("    MP2 Correlation Energy:           %16.12f \n", e_lmp2_);
    outfile->Printf("    Semicanonical MP2 Correction:     %16.12f \n", de_lmp2_eliminated_);
    outfile->Printf("    Dipole Correction:                %16.12f \n", de_dipole_);
    outfile->Printf("    PNO Truncation Correction:        %16.12f \n", de_pno_total_);
    outfile->Printf("\n\n  @Total DLPNO-MP2 Energy: %16.12f \n", variables_["SCF TOTAL ENERGY"] + e_lmp2_ + de_lmp2_eliminated_ + de_pno_total_ + de_dipole_);
    outfile->Printf("\n   * WARNING: This answer will likely vary from one obtained by a energy('dlpno-mp2') call");
    outfile->Printf("\n                due to lack of a semi-canonical MP2 prescreening step in DLPNO-MP2, as well");
    outfile->Printf("\n                as slightly tighter cutoffs utilized to increase accuracy in the context of CC!!!\n\n");

    // Now we do the hard stuff (CCSD)
    timer_on("Sparsity");
    prep_sparsity(false, true);
    timer_off("Sparsity");

    recompute_pnos();

    timer_on("DF Ints");
    compute_qij();
    compute_qab();
    timer_off("DF Ints");

    timer_on("PNO Integrals");
    estimate_memory();
    compute_pno_integrals();
    timer_off("PNO Integrals");

    timer_on("PNO Overlaps");
    compute_pno_overlaps();
    timer_off("PNO Overlaps");

    timer_on("LCCSD");
    lccsd_iterations();
    timer_off("LCCSD");

    // Bye bye (Q_ij | m_ij a_ij) integrals. You won't be missed
    psio_->close(PSIF_DLPNO_QIA_PNO, 0);
    // Bye bye (Q_ij | a_ij b_ij) integrals. You won't be missed
    psio_->close(PSIF_DLPNO_QAB_PNO, 0);

    print_results();

    timer_off("DLPNO-CCSD");
    
    double e_ccsd_corr = e_lccsd_ + de_weak_ + de_lmp2_eliminated_ + de_dipole_ + de_pno_total_;
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
    outfile->Printf("              DOI: 10.1063/5.0219963           \n");
    outfile->Printf("   --------------------------------------------\n\n");
    outfile->Printf("  DLPNO convergence set to %s.\n\n", options_.get_str("PNO_CONVERGENCE").c_str());
    outfile->Printf("  Detailed DLPNO thresholds and cutoffs:\n");
    outfile->Printf("    T_CUT_PNO        = %6.4e \n", T_CUT_PNO_);
    outfile->Printf("    T_DIAG_SCALE     = %6.4e \n", T_CUT_PNO_DIAG_SCALE_);
    outfile->Printf("    T_CORE_SCALE     = %6.4e \n", T_CUT_PNO_CORE_SCALE_);
    outfile->Printf("    T_CUT_TRACE      = %6.4e \n", T_CUT_TRACE_);
    outfile->Printf("    T_CUT_ENERGY     = %6.4e \n", T_CUT_ENERGY_);
    outfile->Printf("    T_CUT_PAIRS      = %6.4e \n", T_CUT_PAIRS_);
    outfile->Printf("    T_CUT_PAIRS_MP2  = %6.4e \n", T_CUT_PAIRS_MP2_);
    outfile->Printf("    T_CUT_PRE        = %6.4e \n", T_CUT_PRE_);
    outfile->Printf("    T_CUT_DO_PRE     = %6.4e \n", options_.get_double("T_CUT_DO_PRE"));
    outfile->Printf("    T_CUT_MKN        = %6.4e \n", T_CUT_MKN_);
    outfile->Printf("    T_CUT_PNO_MP2    = %6.4e \n", T_CUT_PNO_MP2_);
    outfile->Printf("    T_CUT_TRACE_MP2  = %6.4e \n", T_CUT_TRACE_MP2_);
    outfile->Printf("    T_CUT_ENERGY_MP2 = %6.4e \n", T_CUT_ENERGY_MP2_);
    outfile->Printf("    T_CUT_DO         = %6.4e \n", T_CUT_DO_);
    outfile->Printf("    T_CUT_DO_ij      = %6.4e \n", options_.get_double("T_CUT_DO_ij"));
    outfile->Printf("    T_CUT_DO_uv      = %6.4e \n", options_.get_double("T_CUT_DO_uv"));
    outfile->Printf("    T_CUT_CLMO       = %6.4e \n", options_.get_double("T_CUT_CLMO"));
    outfile->Printf("    T_CUT_CPAO       = %6.4e \n", options_.get_double("T_CUT_CPAO"));
    outfile->Printf("    S_CUT            = %6.4e \n", options_.get_double("S_CUT"));
    outfile->Printf("    F_CUT            = %6.4e \n", options_.get_double("F_CUT"));
    outfile->Printf("    INTS_TOL (AO)    = %6.4e \n", options_.get_double("DLPNO_AO_INTS_TOL"));
    outfile->Printf("    MIN_PNOS         = %6d   \n", options_.get_int("MIN_PNOS"));
    outfile->Printf("\n\n");

    outfile->Printf("  ==> Basis Set Information <==\n\n");
    outfile->Printf("   ----------------------------------------------\n");
    outfile->Printf("      NBF   NFRZC   NACT   NDOCC   NVIR   NAUX   \n");
    outfile->Printf("   ----------------------------------------------\n");
    outfile->Printf("    %5d  %5d  %5d   %5d   %5d  %5d\n", basisset_->nbf(), nfrzc(), nalpha_ - nfrzc(), nalpha_,
                                                                                basisset_->nbf() - nalpha_, ribasis_->nbf());
    outfile->Printf("   ----------------------------------------------\n\n");

}

void DLPNOCCSD::print_results() {
    int naocc = i_j_to_ij_.size();
    double t1diag = 0.0;
#pragma omp parallel for reduction(+ : t1diag)
    for (int i = 0; i < naocc; ++i) {
        t1diag += T_ia_[i]->vector_dot(T_ia_[i]);
    }
    t1diag = std::sqrt(t1diag / (2.0 * naocc));
    outfile->Printf("\n  T1 Diagnostic: %8.8f \n", t1diag);
    if (t1diag > 0.02) {
        outfile->Printf("    WARNING: T1 Diagnostic is greater than 0.02, CCSD results may be unreliable!\n");
    }
    set_scalar_variable("CC T1 DIAGNOSTIC", t1diag);

    outfile->Printf("  \n");
    outfile->Printf("  Total DLPNO-CCSD Correlation Energy: %16.12f \n", e_lccsd_ + de_weak_ + de_lmp2_eliminated_ + de_pno_total_ + de_dipole_);
    outfile->Printf("    CCSD Correlation Energy:           %16.12f \n", e_lccsd_);
    outfile->Printf("    Weak Pair Contribution:            %16.12f \n", de_weak_);
    outfile->Printf("    Semicanonical MP2 Correction:      %16.12f \n", de_lmp2_eliminated_);
    outfile->Printf("    Dipole Pair Correction:            %16.12f \n", de_dipole_);
    outfile->Printf("    PNO Truncation Correction:         %16.12f \n", de_pno_total_);
    outfile->Printf("\n\n  @Total DLPNO-CCSD Energy: %16.12f \n", variables_["SCF TOTAL ENERGY"] + e_lccsd_ + de_lmp2_eliminated_ + de_weak_ + de_pno_total_ + de_dipole_);
    outfile->Printf("    *** Wow, that was fast! A thousand hallelujahs!!! \n\n");
}

}  // namespace dlpno
}  // namespace psi