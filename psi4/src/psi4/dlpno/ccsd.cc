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

inline Tensor<double, 3> DLPNOCCSD::QIA_PNO_EINSUMS(const int ij) {
    auto &[i, j] = ij_to_i_j_[ij];
    const int pair_idx = (i > j) ? ij_to_ji_[ij] : ij;

    int naux_ij = lmopair_to_ribfs_[ij].size();
    int nlmo_ij = lmopair_to_lmos_[ij].size();
    int npno_ij = n_pno_[ij];

    if (write_qia_pno_) {
        std::stringstream toc_entry;
        toc_entry << "QIA (PNO) EINSUMS " << pair_idx;

        auto q_ov = std::make_shared<Matrix>(toc_entry.str(), naux_ij, nlmo_ij * npno_ij);
    #pragma omp critical
        q_ov->load(psio_, PSIF_DLPNO_QIA_PNO, psi::Matrix::SubBlocks);

        // TODO: An efficient disk tensor I/O would eliminate the need for this back and forth copying
        Tensor<double, 3> qia_pno("qia_pno", naux_ij, nlmo_ij, npno_ij);
        ::memcpy(qia_pno.data(), q_ov->get_pointer(), naux_ij * nlmo_ij * npno_ij * sizeof(double));

        return qia_pno;
    } else {
        return q_ov_pno_[pair_idx];
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

inline Tensor<double, 3> DLPNOCCSD::QAB_PNO_EINSUMS(const int ij) {
    auto &[i, j] = ij_to_i_j_[ij];
    const int pair_idx = (i > j) ? ij_to_ji_[ij] : ij;

    int naux_ij = lmopair_to_ribfs_[ij].size();
    int npno_ij = n_pno_[ij];

    if (write_qab_pno_) {
        std::stringstream toc_entry;
        toc_entry << "QAB (PNO) EINSUMS " << pair_idx;

        auto q_vv = std::make_shared<Matrix>(toc_entry.str(), naux_ij, npno_ij * npno_ij);
    #pragma omp critical
        q_vv->load(psio_, PSIF_DLPNO_QAB_PNO, psi::Matrix::ThreeIndexLowerTriangle);

        // TODO: An efficient disk tensor I/O would eliminate the need for this back and forth copying
        Tensor<double, 3> qab_pno("qab_pno", naux_ij, npno_ij, npno_ij);
        ::memcpy(qab_pno.data(), q_vv->get_pointer(), naux_ij * npno_ij * npno_ij * sizeof(double));

        return qab_pno;
    } else {
        return q_vv_pno_[pair_idx];
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

    size_t low_overlap_memory = 0, pno_overlap_memory = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : pno_overlap_memory)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        int i, j;
        std::tie(i, j) = ij_to_i_j_[ij];

        const int npno_ij = n_pno_[ij];
        const int nlmo_ij = lmopair_to_lmos_[ij].size();

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

        if (i <= j && !low_memory_overlap_) {
            for (int mn_ij = 0; mn_ij < nlmo_ij * nlmo_ij; mn_ij++) {
                const int m_ij = mn_ij / nlmo_ij, n_ij = mn_ij % nlmo_ij;
                const int m = lmopair_to_lmos_[ij][m_ij], n = lmopair_to_lmos_[ij][n_ij];
                const int mn = i_j_to_ij_[m][n];
                // If these are true, the quantity has already been formed in a previous intermediate
                if (i == m || i == n || j == m || j == n || m == n) continue;
                if (mn == -1 || m_ij > n_ij) continue;

                pno_overlap_memory += n_pno_[ij] * n_pno_[mn];
            }
        }
    }

    size_t oooo = 0;
    size_t ooov = 0;
    size_t oovv = 0;
    size_t oovv_large = 0;
    size_t ovvv = 0;
    size_t qov = 0;
    size_t qvv = 0;

#pragma omp parallel for schedule(dynamic) reduction(+ : oooo, ooov, oovv, oovv_large, ovvv, qov, qvv)
    for (int ij = 0; ij < n_lmo_pairs; ij++) {
        int i, j;
        std::tie(i, j) = ij_to_i_j_[ij];

        const int naux_ij = lmopair_to_ribfs_[ij].size();
        const int nlmo_ij = lmopair_to_lmos_[ij].size();
        const int npno_ij = n_pno_[ij];

        int dij = (i == j) ? 1.0 : 0.0;

        oooo += nlmo_ij * nlmo_ij;
        ooov += (3 - dij) * nlmo_ij * npno_ij;
        oovv += (4 - dij) * npno_ij * npno_ij;
        ovvv += npno_ij * npno_ij * npno_ij;

        bool is_strong_pair = i_j_to_ij_strong_[i][j] != -1;

        if (is_strong_pair) {
            for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
                int k = lmopair_to_lmos_[ij][k_ij];
                int ik = i_j_to_ij_[i][k], jk = i_j_to_ij_[j][k];
                oovv_large += n_pno_[ij] * n_pno_[jk];
                oovv_large += n_pno_[ij] * n_pno_[jk];
            }
        }

        if (i <= j && is_strong_pair) {
            qov += naux_ij * nlmo_ij * npno_ij;
        }

        if (i <= j && is_strong_pair) {
            qvv += naux_ij * npno_ij * npno_ij;
        }
    }

    low_memory_overlap_ = options_.get_bool("LOW_MEMORY_OVERLAP");
    if (low_memory_overlap_) pno_overlap_memory = low_overlap_memory;

    write_qia_pno_ = options_.get_bool("WRITE_QIA_PNO");
    if (write_qia_pno_) qov = 0;

    write_qab_pno_ = options_.get_bool("WRITE_QAB_PNO");
    if (write_qab_pno_) qvv = 0;

    const size_t total_df_memory = qij_memory_ + qia_memory_ + qab_memory_;
    const size_t total_pno_int_memory = oooo + ooov + oovv + oovv_large + ovvv + qov + qvv;
    size_t total_memory = total_df_memory + pno_overlap_memory + total_pno_int_memory;

    // 2^30 bytes per GiB
    outfile->Printf("    (q | i j) [AUX, LMO]       : %.3f [GiB]\n", qij_memory_ * pow(2.0, -30) * sizeof(double));
    outfile->Printf("    (q | i a) [AUX, LMO, PAO]  : %.3f [GiB]\n", qia_memory_ * pow(2.0, -30) * sizeof(double));
    outfile->Printf("    (q | a b) [AUX, PAO]       : %.3f [GiB]\n", qab_memory_ * pow(2.0, -30) * sizeof(double));
    outfile->Printf("    0-virtual LMO ERI          : %.3f [GiB]\n", oooo * pow(2.0, -30) * sizeof(double));
    outfile->Printf("    1-virtual LMO/PNO ERI      : %.3f [GiB]\n", ooov * pow(2.0, -30) * sizeof(double));
    outfile->Printf("    2-virtual LMO/PNO ERI      : %.3f [GiB]\n", oovv * pow(2.0, -30) * sizeof(double));
    outfile->Printf("    2-virtual non-projected    : %.3f [GiB]\n", oovv_large * pow(2.0, -30) * sizeof(double));
    outfile->Printf("    3-virtual LMO/PNO ERI      : %.3f [GiB]\n", ovvv * pow(2.0, -30) * sizeof(double));
    outfile->Printf("    (Q_ij | m_ij a_ij)         : %.3f [GiB]\n", qov * pow(2.0, -30) * sizeof(double));
    outfile->Printf("    (Q_ij | a_ij b_ij)         : %.3f [GiB]\n", qvv * pow(2.0, -30) * sizeof(double));
    outfile->Printf("    PNO overlaps               : %.3f [GiB]\n\n", pno_overlap_memory * pow(2.0, -30) * sizeof(double));
    outfile->Printf("    Total Memory Given         : %.3f [GiB]\n", memory_ * pow(2.0, -30));
    outfile->Printf("    Total Memory Required      : %.3f [GiB]\n\n", total_memory * pow(2.0, -30) * sizeof(double));

    // Memory checks!!!
    bool memory_changed = false;

    if (total_memory * sizeof(double) > 0.75 * memory_) {
        outfile->Printf("  Total Required Memory is more than 75%% of Available Memory!\n");
        outfile->Printf("    Attempting to switch to semi-direct low memory PNO overlap algorithm...\n");

        total_memory += (low_overlap_memory - pno_overlap_memory);
        low_memory_overlap_ = true;
        memory_changed = true;
        pno_overlap_memory = low_overlap_memory;
        outfile->Printf("    Required Memory Reduced to %.3f [GiB]\n\n", total_memory * pow(2.0, -30) * sizeof(double));
    }

    if (total_memory * sizeof(double) > 0.75 * memory_) {
        outfile->Printf("  Total Required Memory is (still) more than 75%% of Available Memory!\n");
        outfile->Printf("    Attempting to switch to disk IO for (Q_{ij}|m_{ij} a_{ij}) integrals...\n");

        total_memory -= qov;
        write_qia_pno_ = true;
        memory_changed = true;
        qov = 0L;
        outfile->Printf("    Required Memory Reduced to %.3f [GiB]\n\n", total_memory * pow(2.0, -30) * sizeof(double));
    }

    if (total_memory * sizeof(double) > 0.75 * memory_) {
        outfile->Printf("  Total Required Memory is (still) more than 75%% of Available Memory!\n");
        outfile->Printf("    Attempting to switch to disk IO for (Q_{ij}|a_{ij} b_{ij}) integrals...\n");

        total_memory -= qvv;
        write_qab_pno_ = true;
        memory_changed = true;
        qvv = 0L;
        outfile->Printf("    Required Memory Reduced to %.3f [GiB]\n\n", total_memory * pow(2.0, -30) * sizeof(double));
    }

    if (total_memory * sizeof(double) > 0.75 * memory_) {
        outfile->Printf("  Total Required Memory is (still) more than 75%% of Available Memory!\n");
        outfile->Printf("    We exhausted all of our options!!! This computation cannot continue...\n");

        throw PSIEXCEPTION("   Too little memory given for DLPNO-CCSD Algorithm!");
    }

    if (memory_changed) {
        outfile->Printf("  ==> (Updated) DLPNO-CCSD Memory Requirements <== \n\n");

        outfile->Printf("    (q | i j) [AUX, LMO]       : %.3f [GiB]\n", qij_memory_ * pow(2.0, -30) * sizeof(double));
        outfile->Printf("    (q | i a) [AUX, LMO, PAO]  : %.3f [GiB]\n", qia_memory_ * pow(2.0, -30) * sizeof(double));
        outfile->Printf("    (q | a b) [AUX, PAO]       : %.3f [GiB]\n", qab_memory_ * pow(2.0, -30) * sizeof(double));
        outfile->Printf("    0-virtual LMO ERI          : %.3f [GiB]\n", oooo * pow(2.0, -30) * sizeof(double));
        outfile->Printf("    1-virtual LMO/PNO ERI      : %.3f [GiB]\n", ooov * pow(2.0, -30) * sizeof(double));
        outfile->Printf("    2-virtual LMO/PNO ERI      : %.3f [GiB]\n", oovv * pow(2.0, -30) * sizeof(double));
        outfile->Printf("    2-virtual non-projected    : %.3f [GiB]\n", oovv_large * pow(2.0, -30) * sizeof(double));
        outfile->Printf("    3-virtual LMO/PNO ERI      : %.3f [GiB]\n", ovvv * pow(2.0, -30) * sizeof(double));
        outfile->Printf("    (Q_ij | m_ij a_ij)         : %.3f [GiB]\n", qov * pow(2.0, -30) * sizeof(double));
        outfile->Printf("    (Q_ij | a_ij b_ij)         : %.3f [GiB]\n", qvv * pow(2.0, -30) * sizeof(double));
        outfile->Printf("    PNO overlaps               : %.3f [GiB]\n\n", pno_overlap_memory * pow(2.0, -30) * sizeof(double));
        outfile->Printf("    Total Memory Given         : %.3f [GiB]\n", memory_ * pow(2.0, -30));
        outfile->Printf("    Total Memory Required      : %.3f [GiB]\n\n", total_memory * pow(2.0, -30) * sizeof(double));
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

            double t_cut_scale = (i == j) ? T_CUT_PNO_DIAG_SCALE_ : 1.0;

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

void DLPNOCCSD::pno_lmp2_iterations() {
    /* Code borrowed and adapted from Zach Glick's DLPNO-MP2 */

    int nbf = basisset_->nbf();
    int naocc = i_j_to_ij_.size();
    int n_lmo_pairs = ij_to_i_j_.size();

    outfile->Printf("\n  ==> PNO-LMP2 Memory Estimate <== \n\n");

    double F_CUT_MP2 = options_.get_double("F_CUT");

    size_t pno_overlap_memory = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : pno_overlap_memory)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];

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

    // => Form PNO overlap integrals for preceeding DLPNO-MP2 iterations
        
    // 2^30 bytes per GiB
    outfile->Printf("    PNO overlaps           : %.3f [GiB]\n", pno_overlap_memory * pow(2.0, -30) * sizeof(double));
    outfile->Printf("    Total Memory Given     : %.3f [GiB]\n\n", memory_ * pow(2.0, -30));

    if (pno_overlap_memory * sizeof(double) > 0.9 * memory_) {
        outfile->Printf("  Total Required Memory is more than 90%% of Available Memory!\n");
        outfile->Printf("    We exhausted all of our options!!! This computation cannot continue...\n");
        outfile->Printf("    Pro Tip: Try adjusting the T_CUT_PNO_MP2, T_CUT_TRACE_MP2, and/or T_CUT_ENERGY_MP2 parameters\n");

        throw PSIEXCEPTION("  Too little memory given for PNO-LMP2 Algorithm!");
    }

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

            e_curr += K_iajb_[ij]->vector_dot(Tt_iajb_[ij]);
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

        double t_cut_scale = (i == j) ? T_CUT_PNO_DIAG_SCALE_ : 1.0;

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
        outfile->Printf("\n  ==> Determining Strong and Weak Pairs (Refined Prescreening Step) <==\n\n");

        int n_lmo_pairs = ij_to_i_j_.size();

        const std::vector<double>& e_ijs = compute_pair_energies<false>();
        filter_pairs<false>(e_ijs);

        int n_strong_pairs = ij_to_i_j_strong_.size();
        int n_weak_pairs = ij_to_i_j_weak_.size();

        outfile->Printf("    Weak Pairs                      = %d\n", n_weak_pairs);
        outfile->Printf("    Strong Pairs                    = %d\n", n_strong_pairs);
        outfile->Printf("    Strong Pairs / Total Pairs      = (%.2f %%)\n", (100.0 * n_strong_pairs) / n_lmo_pairs);
    }
}

void DLPNOCCSD::compute_cc_integrals() {
    outfile->Printf("    Computing CC integrals...\n\n");

    int n_lmo_pairs = ij_to_i_j_.size();

    // 0 virtual
    K_mnij_.resize(n_lmo_pairs);
    // 1 virtual
    K_bar_chem_.resize(n_lmo_pairs);
    K_bar_.resize(n_lmo_pairs);
    L_bar_.resize(n_lmo_pairs);
    // 2 virtual
    lmopair_to_paos_ext_.resize(n_lmo_pairs);
    J_ijab_.resize(n_lmo_pairs);
    L_iajb_.resize(n_lmo_pairs);
    M_iajb_.resize(n_lmo_pairs);
    // 2-virtual non-projected
    J_ij_kj_.resize(n_lmo_pairs);
    K_ij_kj_.resize(n_lmo_pairs);
    // 3 virtual
    K_tilde_chem_.resize(n_lmo_pairs);

    // Einsums integrals
    q_io_list_.resize(n_lmo_pairs);
    q_iv_list_.resize(n_lmo_pairs);

    // DF integrals (used in DLPNO-CCSD with T1 Transformed Hamiltonian)
    if (!write_qia_pno_) {
        Qma_ij_.resize(n_lmo_pairs);
        q_ov_pno_.resize(n_lmo_pairs);
    }
    if (!write_qab_pno_) {
        Qab_ij_.resize(n_lmo_pairs);
        q_vv_pno_.resize(n_lmo_pairs);
    }

    i_Qa_ij_.resize(n_lmo_pairs);
    i_Qk_ij_.resize(n_lmo_pairs);

    size_t qvv_memory = 0;

    psio_->open(PSIF_DLPNO_QIA_PNO, PSIO_OPEN_NEW);
    psio_->open(PSIF_DLPNO_QAB_PNO, PSIO_OPEN_NEW);

    std::time_t time_start = std::time(nullptr);
    std::time_t time_lap = std::time(nullptr);

#pragma omp parallel for schedule(dynamic, 1) reduction(+ : qvv_memory)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        int i, j;
        std::tie(i, j) = ij_to_i_j_[ij];
        const int ji = ij_to_ji_[ij];

        // number of PNOs in the pair domain
        const int npno_ij = n_pno_[ij];
        if (i > j || npno_ij == 0) continue;
        bool is_strong_pair = i_j_to_ij_strong_[i][j] != -1;

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

        auto q_oo = std::make_shared<Matrix>(naux_ij, nlmo_ij * nlmo_ij);
        auto q_ov = std::make_shared<Matrix>(naux_ij, nlmo_ij * npno_ij);
        auto q_vv = std::make_shared<Matrix>(naux_ij, npno_ij * npno_ij);
        
        J_ij_kj_[ij].resize(nlmo_ij);
        if (i != j) J_ij_kj_[ji].resize(nlmo_ij);

        K_ij_kj_[ij].resize(nlmo_ij);
        if (i != j) K_ij_kj_[ji].resize(nlmo_ij);

        std::vector<int> extended_pao_domain;
        int npao_ext_ij;
        
        extended_pao_domain = lmopair_to_paos_[ij];
        for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
            int k = lmopair_to_lmos_[ij][k_ij];
            extended_pao_domain = merge_lists(extended_pao_domain, lmo_to_paos_[k]);
        }
        npao_ext_ij = extended_pao_domain.size();

        lmopair_to_paos_ext_[ij] = extended_pao_domain;
        if (i != j) lmopair_to_paos_ext_[ji] = extended_pao_domain;

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

            auto q_oo_tmp = submatrix_rows_and_cols(*qij_[q], sparse_lmo_list, sparse_lmo_list);
            ::memcpy(&(*q_oo)(q_ij, 0), &(*q_oo_tmp)(0,0), nlmo_ij * nlmo_ij * sizeof(double));

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
        
        q_io_clone = q_io->clone();
        C_DGESV_wrapper(A_solve->clone(), q_io_clone);
        if (i != j) {
            q_jo_clone = q_jo->clone();
            C_DGESV_wrapper(A_solve->clone(), q_jo_clone);
        }

        q_iv_clone = q_iv->clone();
        C_DGESV_wrapper(A_solve->clone(), q_iv_clone);
        if (i != j) {
            q_jv_clone = q_jv->clone();
            C_DGESV_wrapper(A_solve->clone(), q_jv_clone);
        }
        
        A_solve->power(0.5, 1.0e-14);

        C_DGESV_wrapper(A_solve->clone(), q_pair);
        C_DGESV_wrapper(A_solve->clone(), q_io);
        C_DGESV_wrapper(A_solve->clone(), q_jo);
        C_DGESV_wrapper(A_solve->clone(), q_iv);
        C_DGESV_wrapper(A_solve->clone(), q_jv);
        C_DGESV_wrapper(A_solve->clone(), q_oo);
        C_DGESV_wrapper(A_solve->clone(), q_ov);
        C_DGESV_wrapper(A_solve->clone(), q_vv);

        // Copy necessary integrals to Einsums
        q_io_list_[ij] = Tensor<double, 2>("q_io", naux_ij, nlmo_ij);
        ::memcpy(q_io_list_[ij].data(), q_io->get_pointer(), naux_ij * nlmo_ij * sizeof(double));
        q_iv_list_[ij] = Tensor<double, 2>("q_iv", naux_ij, npno_ij);
        ::memcpy(q_iv_list_[ij].data(), q_iv->get_pointer(), naux_ij * npno_ij * sizeof(double));
        if (i != j) {
            q_io_list_[ji] = Tensor<double, 2>("q_jo", naux_ij, nlmo_ij);
            ::memcpy(q_io_list_[ji].data(), q_jo->get_pointer(), naux_ij * nlmo_ij * sizeof(double));
            q_iv_list_[ji] = Tensor<double, 2>("q_jv", naux_ij, npno_ij);
            ::memcpy(q_iv_list_[ji].data(), q_jv->get_pointer(), naux_ij * npno_ij * sizeof(double));
        }

        // Write (Q_{ij} | m_{ij} a_{ij}) integrals to disk
        if (write_qia_pno_) {
            std::stringstream qov_entry;
            qov_entry << "QIA (PNO) EINSUMS " << ij;
            q_ov->set_name(qov_entry.str());
    #pragma omp critical
            q_ov->save(psio_, PSIF_DLPNO_QIA_PNO, psi::Matrix::SubBlocks);
        } else {
            q_ov_pno_[ij] = Tensor<double, 3>("QIA (PNO)", naux_ij, nlmo_ij, npno_ij);
            ::memcpy(q_ov_pno_[ij].data(), q_ov->get_pointer(), naux_ij * nlmo_ij * npno_ij * sizeof(double));
        }

        // Write (Q_{ij} | a_{ij} b_{ij}) integrals to disk
        if (write_qab_pno_) {
            std::stringstream qvv_entry;
            qvv_entry << "QAB (PNO) EINSUMS " << ij;
            q_vv->set_name(qvv_entry.str());
    #pragma omp critical
            q_vv->save(psio_, PSIF_DLPNO_QAB_PNO, psi::Matrix::ThreeIndexLowerTriangle);
        } else {
            q_vv_pno_[ij] = Tensor<double, 3>("QAB (PNO)", naux_ij, npno_ij, npno_ij);
            ::memcpy(q_vv_pno_[ij].data(), q_vv->get_pointer(), naux_ij * npno_ij * npno_ij * sizeof(double));
        }

        if (thread == 0) timer_off("DLPNO-CCSD: Setup Integrals");

        if (is_strong_pair) {

            if (thread == 0) timer_on("DLPNO-CCSD: J_ij_kj Integrals");
            
            for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
                int k = lmopair_to_lmos_[ij][k_ij];
                int ik = i_j_to_ij_[i][k], kj = i_j_to_ij_[k][j];
                J_ij_kj_[ij][k_ij] = std::make_shared<Matrix>(n_pno_[ij], n_pno_[kj]);
                if (i != j) J_ij_kj_[ji][k_ij] = std::make_shared<Matrix>(n_pno_[ij], n_pno_[ik]);
            }

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

                /*
                for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
                    int k = lmopair_to_lmos_[ij][k_ij];
                    int ik = i_j_to_ij_[i][k], kj = i_j_to_ij_[k][j];
                    auto kj_idx = index_list(extended_pao_domain, lmopair_to_paos_[kj]);

                    auto J_ij_kj_temp = linalg::doublet(submatrix_cols(*q_cd_temp, kj_idx), X_pno_[kj], false, false);
                    J_ij_kj_temp->scale((*q_io_clone)(q_ij, k_ij));
                    J_ij_kj_[ij][k_ij]->add(J_ij_kj_temp);

                    if (i != j) {
                        auto ik_idx = index_list(extended_pao_domain, lmopair_to_paos_[ik]);

                        auto J_ij_ik_temp = linalg::doublet(submatrix_cols(*q_cd_temp, ik_idx), X_pno_[ik], false, false);
                        J_ij_ik_temp->scale((*q_jo_clone)(q_ij, k_ij));
                        J_ij_kj_[ji][k_ij]->add(J_ij_ik_temp);
                    }
                }
                */
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
                auto K_iovv_kj = submatrix_cols(*K_iovv_slice, index_list(extended_pao_domain, lmopair_to_paos_[kj]));
                J_ij_kj_[ij][k_ij] = linalg::doublet(K_iovv_kj, X_pno_[kj]);

                if (i != j) {
                    auto K_jovv_slice = submatrix_rows(*K_jovv_partial, std::vector<int>(1, k_ij));
                    K_jovv_slice->reshape(npno_ij, npao_ext_ij);
                    auto K_jovv_ik = submatrix_cols(*K_jovv_slice, index_list(extended_pao_domain, lmopair_to_paos_[ik]));
                    J_ij_kj_[ji][k_ij] = linalg::doublet(K_jovv_ik, X_pno_[ik]);
                } // end if
            } // end for

            if (thread == 0) timer_off("DLPNO-CCSD: J_ij_kj Integrals");

            if (thread == 0) timer_on("DLPNO-CCSD: K_ij_kj Integrals");
            
            for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
                int k = lmopair_to_lmos_[ij][k_ij];
                int ik = i_j_to_ij_[i][k], kj = i_j_to_ij_[k][j];
                K_ij_kj_[ij][k_ij] = std::make_shared<Matrix>(n_pno_[ij], n_pno_[kj]);
                if (i != j) K_ij_kj_[ji][k_ij] = std::make_shared<Matrix>(n_pno_[ij], n_pno_[ik]);
            }

            for (int q_ij = 0; q_ij < naux_ij; q_ij++) {
                const int q = lmopair_to_ribfs_[ij][q_ij];
                const int centerq = ribasis_->function_to_center(q);

                for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
                    int k = lmopair_to_lmos_[ij][k_ij];
                    int ik = i_j_to_ij_[i][k], kj = i_j_to_ij_[k][j];
                    int k_sparse = riatom_to_lmos_ext_dense_[centerq][k];
                    auto kj_idx = index_list(riatom_to_paos_ext_[centerq], lmopair_to_paos_[kj]);
                    
                    auto K_ij_kj_temp = linalg::doublet(submatrix_rows_and_cols(*qia_[q], std::vector<int>(1, k_sparse), kj_idx), X_pno_[kj]);

                    for (int a_ij = 0; a_ij < n_pno_[ij]; ++a_ij) {
                        for (int c_kj = 0; c_kj < n_pno_[kj]; ++c_kj) {
                            (*K_ij_kj_[ij][k_ij])(a_ij, c_kj) += (*q_iv_clone)(q_ij, a_ij) * (*K_ij_kj_temp)(0, c_kj);
                        }
                    }
                    
                    if (i != j) {
                        auto ik_idx = index_list(riatom_to_paos_ext_[centerq], lmopair_to_paos_[ik]);
                        auto K_ij_ik_temp = linalg::doublet(submatrix_rows_and_cols(*qia_[q], std::vector<int>(1, k_sparse), ik_idx), X_pno_[ik]);

                        for (int a_ij = 0; a_ij < n_pno_[ij]; ++a_ij) {
                            for (int c_ik = 0; c_ik < n_pno_[ik]; ++c_ik) {
                                (*K_ij_kj_[ji][k_ij])(a_ij, c_ik) += (*q_jv_clone)(q_ij, a_ij) * (*K_ij_ik_temp)(0, c_ik);
                            }
                        }
                    } // end (i != j)

                }
            }

            if (thread == 0) timer_off("DLPNO-CCSD: K_ij_kj Integrals");

            // Take the difference with respect to projected integrals
            /*
            for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
                int k = lmopair_to_lmos_[ij][k_ij];
                int kj = i_j_to_ij_[k][j], ik = i_j_to_ij_[i][k];
                
                auto S_ij_kj = submatrix_rows_and_cols(*S_pao_, lmopair_to_paos_[ij], lmopair_to_paos_[kj]);
                S_ij_kj = linalg::triplet(X_pno_[ij], S_ij_kj, X_pno_[kj], true, false, false);

                auto J_ij_kj = std::make_shared<Matrix>(npno_ij, npno_ij);
                J_ij_kj->zero();

                auto K_ij_kj = std::make_shared<Matrix>(npno_ij, npno_ij);
                K_ij_kj->zero();

                for (int a_ij = 0; a_ij < npno_ij; ++a_ij) {
                    for (int b_ij = 0; b_ij < npno_ij; ++b_ij) {
                        for (int q_ij = 0; q_ij < naux_ij; ++q_ij) {
                            (*J_ij_kj)(a_ij, b_ij) += (*q_io)(q_ij, k_ij) * (*q_vv)(q_ij, a_ij * npno_ij + b_ij);
                            (*K_ij_kj)(a_ij, b_ij) += (*q_iv)(q_ij, a_ij) * (*q_ov)(q_ij, k_ij * npno_ij + b_ij);
                        } // end q_ij
                    } // end b_ij
                } // end a_ij

                // Take projection difference
                J_ij_kj_[ij][k_ij]->subtract(linalg::doublet(J_ij_kj, S_ij_kj, false, false));
                K_ij_kj_[ij][k_ij]->subtract(linalg::doublet(K_ij_kj, S_ij_kj, false, false));

                if (i != j) {
                    auto S_ij_ik = submatrix_rows_and_cols(*S_pao_, lmopair_to_paos_[ij], lmopair_to_paos_[ik]);
                    S_ij_ik = linalg::triplet(X_pno_[ij], S_ij_ik, X_pno_[ik], true, false, false);

                    auto J_ij_ik = std::make_shared<Matrix>(npno_ij, npno_ij);
                    J_ij_ik->zero();

                    auto K_ij_ik = std::make_shared<Matrix>(npno_ij, npno_ij);
                    K_ij_ik->zero();

                    for (int a_ij = 0; a_ij < npno_ij; ++a_ij) {
                        for (int b_ij = 0; b_ij < npno_ij; ++b_ij) {
                            for (int q_ij = 0; q_ij < naux_ij; ++q_ij) {
                                (*J_ij_ik)(a_ij, b_ij) += (*q_jo)(q_ij, k_ij) * (*q_vv)(q_ij, a_ij * npno_ij + b_ij);
                                (*K_ij_ik)(a_ij, b_ij) += (*q_jv)(q_ij, a_ij) * (*q_ov)(q_ij, k_ij * npno_ij + b_ij);
                            } // end q_ij
                        } // end b_ij
                    } // end a_ij

                    // Take projection difference
                    J_ij_kj_[ji][k_ij]->subtract(linalg::doublet(J_ij_ik, S_ij_ik, false, false));
                    K_ij_kj_[ji][k_ij]->subtract(linalg::doublet(K_ij_ik, S_ij_ik, false, false));
                } // end if
            }
            */
        
        } // end if

        if (thread == 0) timer_on("DLPNO-CCSD: Contract Integrals");

        K_mnij_[ij] = linalg::doublet(q_io, q_jo, true, false);
        K_bar_[ij] = linalg::doublet(q_io, q_jv, true, false);
        J_ijab_[ij] = linalg::doublet(q_pair, q_vv, true, false);
        J_ijab_[ij]->reshape(npno_ij, npno_ij);

        K_bar_chem_[ij] = linalg::doublet(q_pair, q_ov, true, false);
        K_bar_chem_[ij]->reshape(nlmo_ij, npno_ij);
        K_tilde_chem_[ij] = linalg::doublet(q_iv, q_vv, true, false);

        if (i != j) {
            K_mnij_[ji] = K_mnij_[ij]->transpose();
            K_bar_[ji] = linalg::doublet(q_jo, q_iv, true, false);
            J_ijab_[ji] = J_ijab_[ij];
            K_bar_chem_[ji] = K_bar_chem_[ij];
            K_tilde_chem_[ji] = linalg::doublet(q_jv, q_vv, true, false);
        }

        if (is_strong_pair) {
            i_Qk_ij_[ij] = q_io;
            i_Qa_ij_[ij] = q_iv;

            if (i != j) {
                i_Qk_ij_[ji] = q_jo;
                i_Qa_ij_[ji] = q_jv;
            }
        }

        // Save DF integrals (only for strong pairs)
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

        // Lt_iajb
        M_iajb_[ij] = K_iajb_[ij]->clone();
        M_iajb_[ij]->scale(2.0);
        M_iajb_[ij]->subtract(J_ijab_[ij]);

        if (i != j) {
            M_iajb_[ji] = M_iajb_[ij]->transpose();
        }

        if (thread == 0) timer_off("DLPNO-CCSD: Contract Integrals");

        if (thread == 0) {
            std::time_t time_curr = std::time(nullptr);
            int time_elapsed = (int) time_curr - (int) time_lap;
            if (time_elapsed > 60) {
                outfile->Printf("  Time Elapsed from last checkpoint %4d (s), Progress %2d %%, Integrals for (%4d / %4d) Pairs Computed\n", time_elapsed, 
                                    (100 * ij) / n_lmo_pairs, ij, n_lmo_pairs);
                time_lap = std::time(nullptr);
            }
        }
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

        if (i_j_to_ij_strong_[i][j] == -1) continue;

        // => Jiang Eq. 91 <= //

        i_Qk_t1_[ij] = i_Qk_ij_[ij]->clone();

        auto qma_ij = QIA_PNO(ij);
        for (int q_ij = 0; q_ij < naux_ij; ++q_ij) {
            for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
                for (int a_ij = 0; a_ij < npno_ij; ++a_ij) {
                    (*i_Qk_t1_[ij])(q_ij, k_ij) += (*qma_ij[q_ij])(k_ij, a_ij) * (*T_n_ij_[ij])(i_ij, a_ij);
                }
            }
        }

        // => Jiang Eq. 92 <= //

        i_Qa_t1_[ij] = i_Qa_ij_[ij]->clone();
        i_Qa_t1_[ij]->subtract(linalg::doublet(i_Qk_ij_[ij], T_n_ij_[ij]));

        auto qab_ij = QAB_PNO(ij);
        for (int q_ij = 0; q_ij < naux_ij; ++q_ij) {
            auto Qtemp = qab_ij[q_ij]->clone();
            Qtemp->subtract(linalg::doublet(T_n_ij_[ij], qma_ij[q_ij], true, false));

            for (int a_ij = 0; a_ij < npno_ij; ++a_ij) {
                for (int b_ij = 0; b_ij < npno_ij; ++b_ij) {
                    (*i_Qa_t1_[ij])(q_ij, a_ij) += (*Qtemp)(a_ij, b_ij) * (*T_n_ij_[ij])(i_ij, b_ij);
                } // end b_ij
            } // end a_ij
        } // end q_ij
    }

    // Einsums part
    q_io_t1_.resize(n_lmo_pairs);
    q_iv_t1_.resize(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];
        int i_ij = lmopair_to_lmos_dense_[ij][i];

        Tensor<double, 3> Qma_ij = QIA_PNO_EINSUMS(ij); // Read in (Q_{ij} | m_{ij} a_{ij}) integrals from disk
        Tensor<double, 3> Qab_ij = QAB_PNO_EINSUMS(ij); // Read in (Q_{ij} | a_{ij} b_{ij}) integrals from disk
        Tensor<double, 1> T_i = T_n_ij_sums_[ij](i_ij, All); // Slice i out of T_n_ij

        // Jiang Eq. 36
        q_io_t1_[ij] = q_io_list_[ij];
        einsum(1.0, Indices{index::Q, index::k}, &q_io_t1_[ij], 1.0, Indices{index::Q, index::k, index::a}, Qma_ij, 
                Indices{index::a}, T_i);

        // Jiang Eq. 38
        q_iv_t1_[ij] = q_iv_list_[ij];
        einsum(1.0, Indices{index::Q, index::a}, &q_iv_t1_[ij], -1.0, Indices{index::Q, index::k}, q_io_t1_[ij], 
                Indices{index::k, index::a}, T_n_ij_sums_[ij]);
        einsum(1.0, Indices{index::Q, index::a}, &q_iv_t1_[ij], 1.0, Indices{index::Q, index::a, index::b}, Qab_ij, 
                Indices{index::b}, T_i);
    }

    timer_off("DLPNO-CCSD: T1 Ints");
}

void DLPNOCCSD::t1_fock() {

    timer_on("DLPNO-CCSD: T1 Fock");

    int naocc = nalpha_ - nfrzc();
    int n_lmo_pairs = ij_to_i_j_.size();

    // => Step 1: Dressing over the contracted indices <= //

    SharedMatrix Fij_bar = F_lmo_->clone();
    std::vector<SharedMatrix> Fia_bar(n_lmo_pairs);
    std::vector<SharedMatrix> Fai_bar(naocc);
    std::vector<SharedMatrix> Fab_bar(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];
        int ji = ij_to_ji_[ij];
        int pair_idx = (i > j) ? ji : ij;

        int naux_ij = lmopair_to_ribfs_[ij].size();
        int nlmo_ij = lmopair_to_lmos_[ij].size();
        int npno_ij = n_pno_[ij];

        // Partially dress Fij (Jiang Eq. 98)
        (*Fij_bar)(i, j) += 2.0 * T_n_ij_[ij]->vector_dot(K_bar_chem_[ij]);
        (*Fij_bar)(i, j) -= T_n_ij_[ij]->vector_dot(K_bar_[ji]);

        if (i > j || i_j_to_ij_strong_[i][j] == -1) continue;

        // Partially dress Fia and Fab (Jiang Eq. 99 and 101)
        Fia_bar[ij] = std::make_shared<Matrix>(nlmo_ij, npno_ij);

        Fab_bar[ij] = std::make_shared<Matrix>(npno_ij, npno_ij);
        Fab_bar[ij]->set_diagonal(e_pno_[ij]);

        auto qma_ij = QIA_PNO(ij);
        auto qab_ij = QAB_PNO(ij);
        for (int q_ij = 0; q_ij < naux_ij; ++q_ij) {
            auto Qma = qma_ij[q_ij];
            auto Qab = qab_ij[q_ij];
            double gamma = Qma->vector_dot(T_n_ij_[ij]);

            // J like contributions
            auto Jcont = Qma->clone();
            Jcont->scale(2.0 * gamma);
            Fia_bar[ij]->add(Jcont);

            Jcont = Qab->clone();
            Jcont->scale(2.0 * gamma);
            Fab_bar[ij]->add(Jcont);

            // K like contributions
            auto Kcont = linalg::triplet(Qma, T_n_ij_[ij], Qma, false, true, false);
            Fia_bar[ij]->subtract(Kcont);

            Kcont = linalg::triplet(Qab, T_n_ij_[ij], Qma, false, true, false);
            Fab_bar[ij]->subtract(Kcont);
        }

        // Partially dress Fai (Jiang Eq. 100)
        if (i == j) {
            Fai_bar[i] = std::make_shared<Matrix>(npno_ij, 1);

            auto Qia = i_Qa_ij_[ij]->clone();
            auto Qik = i_Qk_ij_[ij]->clone();

            for (int q_i = 0; q_i < naux_ij; ++q_i) {
                auto Qma = qma_ij[q_i];
                double gamma = Qma->vector_dot(T_n_ij_[ij]);

                auto Qab = qab_ij[q_i];
                auto lambda = linalg::doublet(Qab, T_n_ij_[ij], false, true);

                for (int a_i = 0; a_i < npno_ij; ++a_i) {
                    // J like contribution
                    (*Fai_bar[i])(a_i, 0) += 2.0 * gamma * (*Qia)(q_i, a_i);
                    for (int k_i = 0; k_i < nlmo_ij; ++k_i) {
                        // K like contribution
                        (*Fai_bar[i])(a_i, 0) -= (*lambda)(a_i, k_i) * (*Qik)(q_i, k_i);
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
        int i_jj = lmopair_to_lmos_dense_[jj][i];
        for (int a_jj = 0; a_jj < n_pno_[jj]; ++a_jj) {
            (*Fkj_)(i, j) += (*Fia_bar[jj])(i_jj, a_jj) * (*T_ia_[j])(a_jj, 0);
        }

        // Fkc matrices (built separately since Fia intermediate not built for weak pairs)
        // Jiang Eq. 95
        Fkc_[ij] = std::make_shared<Matrix>(1, npno_ij);

        for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
            int k = lmopair_to_lmos_[ij][k_ij];
            int ik = i_j_to_ij_[i][k], kk = i_j_to_ij_[k][k];

            // Fkc contributions
            auto T_k = linalg::doublet(S_PNO(ik, kk), T_ia_[k]);
            Fkc_[ij]->add(linalg::triplet(S_PNO(ij, ik), L_iajb_[ik], T_k)->transpose());
        }

        if (i_j_to_ij_strong_[i][j] == -1 || i > j) continue;

        // Fully dress Fab matrices (Jiang Eq. 97)
        Fab_[ij] = Fab_bar[ij]->clone();
        Fab_[ij]->subtract(linalg::doublet(T_n_ij_[ij], Fia_bar[ij], true, false));

        // Fully dress Fai matrices (Jiang Eq. 96)
        if (i == j) {
            Fai_[i] = Fai_bar[i]->clone();
            Fai_[i]->add(linalg::doublet(Fab_bar[ij], T_ia_[i], false, false));
            Fai_[i]->subtract(linalg::triplet(T_n_ij_[ij], Fia_bar[ij], T_ia_[i], true, false, false));

            for (int a_i = 0; a_i < npno_ij; ++a_i) {
                for (int k_i = 0; k_i < nlmo_ij; ++k_i) {
                    int k = lmopair_to_lmos_[ij][k_i];
                    (*Fai_[i])(a_i, 0) -= (*T_n_ij_[ij])(k_i, a_i) * (*Fij_bar)(k, i);
                } // end k_i
            } // end a_i
            
        } // end i == j
    }

    // => Einsums part... old stuff goes away once we confirm this works :)

    // T1-dressed Fock matrix elements (Jiang Eq. 94 - 97)
    F_ki_t1_.resize(n_lmo_pairs);
    F_kc_t1_.resize(n_lmo_pairs);
    F_ai_t1_.resize(naocc); // These only need to be computed for singles
    F_ab_t1_.resize(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];
        // lmopair_to_lmos_dense_ tracks all orbitals m in the molecule for a domain ij.
        // If m is in the domain of ij, it gives the corresponding m_ij, else, returns -1
        int i_ij = lmopair_to_lmos_dense_[ij][i], j_ij = lmopair_to_lmos_dense_[ij][j];
        int ji = ij_to_ji_[ij], jj = i_j_to_ij_[j][j];
        int pair_idx = (i > j) ? ji : ij;

        if (i_j_to_ij_strong_[i][j] == -1 || i > j) continue;

        int naux_ij = lmopair_to_ribfs_[ij].size(); // Local auxiliary basis (defined by Mulliken population)
        int nlmo_ij = lmopair_to_lmos_[ij].size();  // All occ orbitals k, such that ik AND jk are valid pairs
        int npno_ij = n_pno_[ij]; // Pair natural orbitals per ij

        // Read in 3-center integrals from disk (or core... depending on option)
        Tensor<double, 3> Qma_ij = QIA_PNO_EINSUMS(ij); // (naux_ij * nlmo_ij * npno_ij)
        Tensor<double, 3> Qab_ij = QAB_PNO_EINSUMS(ij); // (naux_ij * npno_ij * npno_ij)

        // Common intermediate is Gamma_{Q_{ij}}, which is used in every Fock matrix
        // \Gamma_{Q} = B^{Q}_{me} t_{m}^{e}
        // T_n_ij_sums = t_{n_{ij}}^{a_{nn}} S_{a_{nn}}^{a_{ij}} (Jiang Eq. 70)
        // Singles in domain of ij with relavant amplitudes projected onto PNO space of ij
        Tensor<double, 1> Gamma_Q("Gamma_Q", naux_ij);
        einsum(0.0, Indices{index::Q}, &Gamma_Q, 1.0, Indices{index::Q, index::m, index::e}, Qma_ij,
                Indices{index::m, index::e}, T_n_ij_sums_[ij]);

        // => Part 1: J-like contributions <= //

        // \overline{F}_{ki} = f_{ki} + [2B^{Q}_{ki}B^{Q}_{me} - B^{Q}_{ke}B^{Q}_{mi}]t_{m}^{e}
        // = f_{ki} + 2B^{Q}_{ki} \Gamma_{Q} - B^{Q}_{ke}B^{Q}_{mi}t_{m}^{e}
        // Reindex of Eq. 98, using Eq. 44 as baseline
        Tensor<double, 1> F_ki_bar("F_ki_bar", nlmo_ij);

        // Copy LMO Fock matrix elements into F_ki_bar
        for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
            int k = lmopair_to_lmos_[ij][k_ij];
            F_ki_bar(k_ij) = (*F_lmo_)(k, i);
        }

        // J-like contribution to F_ki_bar (2B^{Q}_{ki} \Gamma_{Q})
        einsum(1.0, Indices{index::k}, &F_ki_bar, 2.0, Indices{index::Q, index::k}, q_io_list_[ij],
                Indices{index::Q}, Gamma_Q);
        // Einsums sequence
        // F_ki = A (1.0) * F_{ki} + B (2.0) * (Q, k) (q_io_list_[ij]) * (Q) (Gamma_Q)
        // (Q, k) (Q) -> GEMV (matrix-vector product)

        // \overline{F}_{kc} = f_{kc} + [2B^{Q}_{kc}B^{Q}_{me} - B^{Q}_{ke}B^{Q}_{mc}]t_{m}^{e} 
        // Reindex of (Jiang Eq. 99), using Eq. 44 as baseline
        // f_{kc} is zero in restricted closed-shell case
        Tensor<double, 2> F_kc_bar("F_kc_bar", nlmo_ij, npno_ij);

        // J-like contribution to F_kc_bar \overline{F}_{kc} = 2 B^{Q}_{kc} \Gamma_{Q}
        einsum(0.0, Indices{index::k, index::c}, &F_kc_bar, 2.0, Indices{index::Q, index::k, index::c},
                Qma_ij, Indices{index::Q}, Gamma_Q);

        // \overline{F}_{ab} = f_{ab} + [2B^{Q}_{ab}B^{Q}_{me} - B^{Q}_{ae}B^{Q}_{mb}]t_{m}^{e} 
        // Reindex of (Jiang Eq. 101), using Eq. 44 as baseline
        // f_{a_{ij}b_{ij}} is diagonal in each PNO domain
        Tensor<double, 2> F_ab_bar("F_ab_bar", npno_ij, npno_ij);
        F_ab_bar.zero();

        for (int a_ij = 0; a_ij < npno_ij; ++a_ij) {
            (F_ab_bar)(a_ij, a_ij) = (*e_pno_[ij])(a_ij);
        }

        // J-like contribution to F_ab_bar \overline{F}_{ab} = 2 B^{Q}_{ab} \Gamma_{Q}
        einsum(1.0, Indices{index::a, index::b}, &F_ab_bar, 2.0, Indices{index::Q, index::a, index::b},
                Qab_ij, Indices{index::Q}, Gamma_Q);

        Tensor<double, 1> F_ai_bar;
        if (i == j) {
            // \overline{F}_{ai} = f_{ai} + [2B^{Q}_{ai}B^{Q}_{me} - B^{Q}_{ae}B^{Q}_{mi}]t_{m}^{e}
            // Reindex of (Jiang Eq. 100), using Eq. 44 as baseline
            // f_{a_{ii}i} is zero in restricted closed-shell reference
            F_ai_bar = Tensor<double, 1>("F_ai_bar", npno_ij);
            F_ai_bar.zero();

            // J-like contribution to F_ai_bar \overline{F}_{ai} = 2 B^{Q}_{ai} \Gamma_{Q}
            einsum(1.0, Indices{index::a}, &F_ai_bar, 2.0, Indices{index::Q, index::a},
                    q_iv_list_[ij], Indices{index::Q}, Gamma_Q);
        }

        // => Part 2: K-like contributions <= //
        
        for (int q_ij = 0; q_ij < naux_ij; ++q_ij) {
            // Integral slices
            Tensor<double, 2> Qov_slice = (Qma_ij)(q_ij, All, All);
            Tensor<double, 2> Qvv_slice = (Qab_ij)(q_ij, All, All);
            Tensor<double, 1> Qmi_slice = (q_io_list_[ij])(q_ij, All);

            // Common intermediates

            // This intermediate is used for F_ki and F_ai
            // F_xi_int (e) = Qmi_slice (m) * T_n_ij_sums (m, e)
            Tensor<double, 1> F_xi_int("F_xi_int", npno_ij);
            einsum(0.0, Indices{index::e}, &F_xi_int, 1.0, Indices{index::m, index::e}, T_n_ij_sums_[ij],
                    Indices{index::m}, Qmi_slice);

            // K-like contribution to F_ki_bar
            // \overline{F}_{ki} -= B^{Q}_{ke} B^{Q}_{mi} t_{m}^{e} (k, e) (m) (m, e)
            // or \overline{F}_{ki} -= B^{Q}_{ke} (k, e) F_xi_int (e)
            einsum(1.0, Indices{index::k}, &F_ki_bar, -1.0, Indices{index::k, index::e}, Qov_slice, 
                    Indices{index::e}, F_xi_int);

            // K-like contribution to F_kc_bar
            // \overline{F}_{kc} -= B^{Q}_{ke} B^{Q}_{mc} t_{m}^{e} (k, e) (m, c) (m, e)
            Tensor<double, 2> F_kc_int("F_kc_int", nlmo_ij, nlmo_ij);
            einsum(0.0, Indices{index::k, index::m}, &F_kc_int, 1.0, Indices{index::k, index::e},
                    Qov_slice, Indices{index::m, index::e}, T_n_ij_sums_[ij]);
            einsum(1.0, Indices{index::k, index::c}, &F_kc_bar, -1.0, Indices{index::k, index::m}, 
                    F_kc_int, Indices{index::m, index::c}, Qov_slice);

            // K-like contribution to F_ab_bar
            // \overline{F}_{ab} -= B^{Q}_{ae} B^{Q}_{mb} t_{m}^{e} (a, e) (m, b) (m, e)
            Tensor<double, 2> F_ab_int("F_ab_int", npno_ij, nlmo_ij);
            einsum(0.0, Indices{index::a, index::m}, &F_ab_int, 1.0, Indices{index::a, index::e},
                    Qvv_slice, Indices{index::m, index::e}, T_n_ij_sums_[ij]);
            einsum(1.0, Indices{index::a, index::b}, &F_ab_bar, -1.0, Indices{index::a, index::m}, 
                    F_ab_int, Indices{index::m, index::b}, Qov_slice);

            if (i == j) {
                // K-like contribution to F_ai_bar
                // \overline{F}_{ai} -= B^{Q}_{ae} B^{Q}_{mi} t_{m}^{e} (a, e) (m) (m, e)
                // \overline{F}_{ai} -= B^{Q}_{ae} (a, e) F_xi_int (e)
                einsum(1.0, Indices{index::a}, &F_ai_bar, -1.0, Indices{index::a, index::e}, 
                        Qvv_slice, Indices{index::e}, F_xi_int);
            }
        } // end q_ij

        // Build Fully T1-Dressed Fock Matrix Intermediates

        // \widetilde{F}_{ki} = \overline{F}_{ki} + \overline{F}_{kc}t_{i}^{c}
        // (Jiang Eq. 94)... small modification... using more robust virtual space of ij
        // compared to jj
        F_ki_t1_[ij] = F_ki_bar; // in einsums, this is copy
        Tensor<double, 1> T_i = (T_n_ij_sums_[ij])(i_ij, All); // Grab projected i from all projected m_ij amplitudes
        einsum(1.0, Indices{index::k}, &F_ki_t1_[ij], 1.0, Indices{index::k, index::c}, F_kc_bar,
                Indices{index::c}, T_i);

        // We do not need to make modifications for fully T1-dressed Fkc
        // (Jiang Eq. 95 assumed that the index was restricted over strong pairs only,
        // which is no longer the case due to strong/weak pair coupling
        F_kc_t1_[ij] = F_kc_bar; // in einsums, this is copy

        // \widetilde{F}_{ab} = \overline{F}_{ab} - t_{k}^{a} \overline{F}_{kb}
        // (Jiang Eq. 96)... no mods
        F_ab_t1_[ij] = F_ab_bar;
        einsum(1.0, Indices{index::a, index::b}, &F_ab_t1_[ij], -1.0, Indices{index::k, index::a},
                T_n_ij_sums_[ij], Indices{index::k, index::b}, F_kc_bar);

        // Only do if i == j, since only needed for singles
        if (i == j) {
            // \widetilde{F}_{ai} = \overline{F}_{ai} - t_{k}^{a} \overline{F}_{ki}
            // + \overline{F}_{ab} t_{i}^{b} - t_{k}^{a} \overline{F}_{kb} t_{i}^{b}
            // => factor out - t_{k}^{a}
            // \overline{F}_{ai} + \overline{F}_{ab} t_{i}^{b} - t_{k}^{a} 
            // (\overline{F}_{ki} + \overline{F}_{kb} t_{i}^{b}) <- Eq. 94 (Jaden is genius)
            // \widetilde{F}_{ai} = \overline{F}_{ai} + \overline{F}_{ab} t_{i}^{b}
            // - t_{k}^{a} \widetilde{F}_{ki}
            F_ai_t1_[i] = F_ai_bar;
            einsum(1.0, Indices{index::a}, &F_ai_t1_[i], 1.0, Indices{index::a, index::b}, F_ab_bar,
                    Indices{index::b}, T_i);
            einsum(1.0, Indices{index::a}, &F_ai_t1_[i], -1.0, Indices{index::k, index::a}, T_n_ij_sums_[ij],
                    Indices{index::k}, F_ki_t1_[ij]);
        } // end if

    } // end for

    timer_off("DLPNO-CCSD: T1 Fock");
}

std::vector<SharedMatrix> DLPNOCCSD::compute_beta() {
    timer_on("DLPNO-CCSD: beta");

    int naocc = nalpha_ - nfrzc();
    int n_lmo_pairs = ij_to_i_j_.size();

    std::vector<SharedMatrix> beta(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];
        int ji = ij_to_ji_[ij];

        if (i_j_to_ij_strong_[i][j] == -1) continue;

        int naux_ij = lmopair_to_ribfs_[ij].size();
        int nlmo_ij = lmopair_to_lmos_[ij].size();
        int pair_idx = (i > j) ? ji : ij;
        int i_ij = lmopair_to_lmos_dense_[ij][i], j_ij = lmopair_to_lmos_dense_[ij][j];

        // Jiang Eq. 82a
        beta[ij] = linalg::doublet(i_Qk_t1_[ij], i_Qk_t1_[ji], true, false);

        // Jiang Eq. 82b
        auto qma_ij = QIA_PNO(ij);
        for (int q_ij = 0; q_ij < naux_ij; ++q_ij) {
            beta[ij]->add(linalg::triplet(qma_ij[q_ij], T_iajb_[ij], qma_ij[q_ij], false, false, true));
        }
    }

    timer_off("DLPNO-CCSD: beta");

    return beta;
}

std::vector<SharedMatrix> DLPNOCCSD::compute_gamma() {
    
    timer_on("DLPNO-CCSD: gamma");

    int naocc = nalpha_ - nfrzc();
    int n_lmo_pairs = ij_to_i_j_.size();

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

        // Jiang Eq. 83a
        gamma[ki]->subtract(linalg::doublet(T_n_ij_[ki], K_bar_chem_[ki], true, false));

        // Jiang Eq. 83b
        auto T_i = linalg::doublet(S_PNO(ki, ii), T_ia_[i]);
        auto K_temp = linalg::doublet(T_i, K_tilde_chem_[ki], true, false);
        K_temp->reshape(n_pno_[ki], n_pno_[ki]);
        gamma[ki]->add(K_temp);

        // Jiang Eq. 83c
        for (int l_ki = 0; l_ki < nlmo_ki; ++l_ki) {
            int l = lmopair_to_lmos_[ki][l_ki];
            int kl = i_j_to_ij_[k][l], ll = i_j_to_ij_[l][l];

            auto T_l = linalg::doublet(S_PNO(ki, ll), T_ia_[l]);
            auto T_i_kl = linalg::doublet(S_PNO(kl, ii), T_ia_[i]);
            auto K_kl = linalg::triplet(S_PNO(ki, kl), K_iajb_[kl], T_i_kl, false, true, false);

            C_DGER(npno_ki, npno_ki, -1.0, T_l->get_pointer(), 1, K_kl->get_pointer(), 1, gamma[ki]->get_pointer(), npno_ki);
        }
        
        // Jiang Eq. 83d
        for (int l_ki = 0; l_ki < nlmo_ki; ++l_ki) {
            int l = lmopair_to_lmos_[ki][l_ki];
            int li = i_j_to_ij_[l][i], kl = i_j_to_ij_[k][l];

            auto gamma_temp = linalg::triplet(T_iajb_[li], S_PNO(li, kl), K_iajb_[kl]);
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

        // Jiang Eq. 84a
        auto L_bar_temp = K_bar_[ik]->clone();
        L_bar_temp->scale(2.0);
        L_bar_temp->subtract(K_bar_chem_[ik]);
        delta[ik]->subtract(linalg::doublet(T_n_ij_[ik], L_bar_temp, true, false));

        // Jiang Eq. 84b
        auto T_i = linalg::doublet(S_PNO(ik, ii), T_ia_[i]);
        auto L_temp = K_tilde_chem_[ki]->clone();
        L_temp->scale(2.0);
        L_temp->reshape(npno_ik * npno_ik, npno_ik);
        L_temp = linalg::doublet(L_temp, T_i);
        L_temp->reshape(npno_ik, npno_ik);
        delta[ik]->add(L_temp->transpose());        
        L_temp = K_tilde_chem_[ki]->clone();
        L_temp = linalg::doublet(T_i, L_temp, true, false);
        L_temp->reshape(npno_ik, npno_ik);
        delta[ik]->subtract(L_temp->transpose());

        // Jiang Eq. 84c
        for (int l_ik = 0; l_ik < nlmo_ik; ++l_ik) {
            int l = lmopair_to_lmos_[ik][l_ik];
            int ll = i_j_to_ij_[l][l], lk = i_j_to_ij_[l][k];

            auto T_l = linalg::doublet(S_PNO(ik, ll), T_ia_[l]);
            auto T_i_lk = linalg::doublet(S_PNO(lk, ii), T_ia_[i]);
            auto L_lk = linalg::triplet(T_i_lk, L_iajb_[lk], S_PNO(lk, ik), true, false, false);

            C_DGER(npno_ik, npno_ik, -1.0, T_l->get_pointer(), 1, L_lk->get_pointer(), 1, delta[ik]->get_pointer(), npno_ik);
        }

        // Jiang Eq. 84d
        for (int l_ik = 0; l_ik < nlmo_ik; ++l_ik) {
            int l = lmopair_to_lmos_[ik][l_ik];
            int il = i_j_to_ij_[i][l], lk = i_j_to_ij_[l][k];

            auto delta_temp = linalg::triplet(Tt_iajb_[il], S_PNO(il, lk), L_iajb_[lk]);
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

    SharedMatrix Fkj_double_tilde = Fkj_->clone();

    // Jiang Eq. 86
#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < naocc * naocc; ++ij) {
        int i = ij / naocc, j = ij % naocc;
        for (int l = 0; l < naocc; ++l) {
            int il = i_j_to_ij_[i][l], lj = i_j_to_ij_[l][j];
            if (il == -1 || lj == -1) continue;

            auto U_lj = linalg::triplet(S_PNO(il, lj), Tt_iajb_[lj], S_PNO(lj, il));
            (*Fkj_double_tilde)(i, j) += K_iajb_[il]->vector_dot(U_lj->transpose());
        } // end l
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

    // Initialize R1 residuals, DePrince 2013 Eq. 19a, Jiang Eq. 87a
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < naocc; ++i) {
        int ii = i_j_to_ij_[i][i];
        // R_ia[i]->copy(Fai_[i]);
        // Initialize Ria as T1-dressed Fock matrix element Fai
        ::memcpy(R_ia[i]->get_pointer(), F_ai_t1_[i].data(), n_pno_[ii] * sizeof(double));
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
        
        // A_{i}^{a} = u_{ik}^{dc} [(ad|kc) - t_{l}^{a}(ld|kc)] (Jiang Eq. 88)
        // (d_{ik}, c_{ik}) * (Q_{ik}, a_{ik}, d_{ik}) * (Q_{ik}, k, c_{ik})
        // (d_{ik}, c_{ik}) * (l_{ik}, a_{ik}) * (Q_{ik}, l_{ik}, d_{ik}) * (Q_{ik}, k, c_{ik})
        auto K_kcad = K_tilde_chem_[ki]->clone();
        K_kcad->reshape(n_pno_[ki] * n_pno_[ki], n_pno_[ki]);
        auto Uki = Tt_iajb_[ki]->clone();
        Uki->reshape(npno_ik * npno_ik, 1);
        R_ia_buffer[thread][i]->add(linalg::triplet(S_PNO(ki, ii), K_kcad, Uki, true, true, false));

        // C_{i}^{a} = F_{kc}U_{ik}^{ac} (DePrince 2013 Eq. 22, Jiang Eq. 90)
        R_ia_buffer[thread][i]->add(linalg::triplet(S_PNO(ik, ii), Tt_iajb_[ik], Fkc_[ki], true, false, true));
    } // end ki

    // A2 and B contributions
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

        // B_{i}^{a} = -u_{kl}^{ac}B^{Q}_{ki}B^{Q}_{lc} (DePrince 2013 Eq. 21, Jiang Eq. 89a)
        auto K_kilc = K_bar_[kl]->clone();
        K_kilc->add(linalg::doublet(T_n_ij_[kl], K_iajb_[kl]));
        auto B_ia = linalg::doublet(Tt_iajb_[kl], K_kilc, false, true);

        for (int i_kl = 0; i_kl < nlmo_kl; ++i_kl) {
            int i = lmopair_to_lmos_[kl][i_kl];
            int ii = i_j_to_ij_[i][i], ki = i_j_to_ij_[k][i];

            // A2 contribution (Jiang Eq. 88b)
            std::vector<int> l_ii_slice(1, lmopair_to_lmos_dense_[ii][l]);
            auto U_ki = linalg::triplet(S_PNO(kl, ki), Tt_iajb_[ki], S_PNO(ki, kl));
            auto T_l = submatrix_rows(*T_n_ij_[ii], l_ii_slice);
            T_l->scale(K_iajb_[kl]->vector_dot(U_ki));
            R_ia_buffer[thread][i]->subtract(T_l->transpose());

            // B contribution (Jiang Eq. 89b)
            std::vector<int> i_kl_slice(1, i_kl);
            R_ia_buffer[thread][i]->subtract(linalg::doublet(S_PNO(kl, ii), submatrix_cols(*B_ia, i_kl_slice), true, false));
        }
    }

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

    auto beta = compute_beta();
    auto gamma = compute_gamma();
    auto delta = compute_delta();
    auto Fkj_double_tilde = compute_Fkj_double_tilde();

    // Zero out residuals
#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        R_iajb[ij]->zero();
        Rn_iajb[ij]->zero();
    }

    // Compute residual for doubles amplitude
#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];
        bool is_weak_pair = (i_j_to_ij_strong_[i][j] == -1);
        int ji = ij_to_ji_[ij];

        int nlmo_ij = lmopair_to_lmos_[ij].size();
        int naux_ij = lmopair_to_ribfs_[ij].size();
        int npno_ij = n_pno_[ij];

        // Skip if pair is weak or contains no pair natural orbitals
        if (is_weak_pair || npno_ij == 0) continue;

        int pair_idx = (i > j) ? ji : ij;

        if (i <= j) {
            // => SECTION 0: SETUP (buffers, indices, etc.) <= //

            // Useful information for integral slices
            int i_ij = lmopair_to_lmos_dense_[ij][i], j_ij = lmopair_to_lmos_dense_[ij][j];

            // Buffers for contractions and holding data
            Tensor<double, 2> R_ij("R_ij", npno_ij, npno_ij);
            R_ij.zero();

            Tensor<double, 2> pno_buffer_a("pno_buffer_a", npno_ij, npno_ij);
            Tensor<double, 2> pno_buffer_b("pno_buffer_b", npno_ij, npno_ij);

            // => SECTION 1: INTEGRALS <= //

            // Expensive DF integrals are read from disk once to minimize Disk I/O overhead
            Tensor<double, 3> Qma_ij = QIA_PNO_EINSUMS(ij);
            Tensor<double, 3> Qab_ij = QAB_PNO_EINSUMS(ij);

            // Form T1-dressed version of Qab_ij (Jiang Eq. 39)
            Tensor<double, 3> Qab_t1("Qab_t1", naux_ij, npno_ij, npno_ij);
            for (int q_ij = 0; q_ij < naux_ij; ++q_ij) {
                Tensor<double, 2> Qma_slice = Qma_ij(q_ij, All, All);
                Tensor<double, 2> Qab_slice = Qab_ij(q_ij, All, All);
                einsum(1.0, Indices{index::a, index::b}, &Qab_slice, -1.0, Indices{index::k, index::a}, T_n_ij_sums_[ij], 
                        Indices{index::k, index::b}, Qma_slice);
                
                Qab_t1(q_ij, All, All) = Qab_slice;
            }

            // => SECTION 2: PROJECTED AMPLITUDES <= //

            // Father of all overlap integrals (all overlap integrals are built using semi-direct algorithm)
            // (using the extended PAO domain, which combines PAOs with all pairs ij with PAOs from all possible pairs
            // ik, jk, il, jl, as well as kl
            SharedMatrix S_ij;
            std::vector<int> pair_ext_domain;
            for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
                int k = lmopair_to_lmos_[ij][k_ij];
                for (int l_ij = 0; l_ij < nlmo_ij; ++l_ij) {
                    int l = lmopair_to_lmos_[ij][l_ij];
                    pair_ext_domain = merge_lists(pair_ext_domain, lmo_to_paos_[k]);
                    pair_ext_domain = merge_lists(pair_ext_domain, lmo_to_paos_[l]);
                } // end k_ij
            } // end l_ij
            S_ij = submatrix_rows_and_cols(*S_pao_, lmopair_to_paos_[ij], pair_ext_domain);
            S_ij = linalg::doublet(X_pno_[ij], S_ij, true, false);

            // Amplitude intermediates
            Tensor<double, 2> T_ij_ein("T_ij_ein", npno_ij, npno_ij);
            ::memcpy(T_ij_ein.data(), T_iajb_[ij]->get_pointer(), npno_ij * npno_ij * sizeof(double));

            // => SECTION 3: Contract Amplitudes and Integrals into R2 Residual! <= //

            // R_{ij}^{ab} += (i a_ij | q_ij)' * (q_ij | j b_ij)' (DePrince 2013 Eq. 10a, Jiang Eq. 75)
            einsum(0.0, Indices{index::a, index::b}, &R_ij, 1.0, Indices{index::Q, index::a}, q_iv_t1_[ij],
                    Indices{index::Q, index::b}, q_iv_t1_[ji]);

            // R_{ij}^{ab} += B^{Q}_{ac} * t_{ij}^{cd} * B^{Q}_{bd} (DePrince 2013 Eq. 11, Jiang Eq. 76)
            for (int q_ij = 0; q_ij < naux_ij; ++q_ij) {
                Tensor<double, 2> Qab_slice = Qab_t1(q_ij, All, All);

                einsum(0.0, Indices{index::a, index::d}, &pno_buffer_a, 1.0, Indices{index::a, index::c}, Qab_slice,
                        Indices{index::c, index::d}, T_ij_ein);
                einsum(0.0, Indices{index::a, index::b}, &pno_buffer_b, 1.0, Indices{index::a, index::d}, pno_buffer_a,
                        Indices{index::b, index::d}, Qab_slice);
                R_ij += pno_buffer_b;
            }

            // beta_{ij}^{kl} = B^{Q}_{ki} B^{Q}_{lj} + B^{Q}_{kc} B^{Q}_{ld} t_{ij}^{cd} (Jiang Eq. 27)
            Tensor<double, 2> beta_ij("beta_ij", nlmo_ij, nlmo_ij);
            einsum(0.0, Indices{index::k, index::l}, &beta_ij, 1.0, Indices{index::Q, index::k}, q_io_t1_[ij], 
                    Indices{index::Q, index::l}, q_io_t1_[ji]);

            for (int q_ij = 0; q_ij < naux_ij; ++q_ij) {
                Tensor<double, 2> Qma_slice = Qma_ij(q_ij, All, All);

                // beta_{ij}^{kl} += B^{Q}_{kc} B^{Q}_{ld} t_{ij}^{cd}
                Tensor<double, 2> beta_temp("beta_temp", npno_ij, nlmo_ij);
                einsum(0.0, Indices{index::c, index::l}, &beta_temp, 1.0, Indices{index::c, index::d}, T_ij_ein,
                        Indices{index::l, index::d}, Qma_slice);

                einsum(1.0, Indices{index::k, index::l}, &beta_ij, 1.0, Indices{index::k, index::c}, Qma_slice,
                        Indices{index::c, index::l}, beta_temp);
            }
            
            // Make B and E intermediates

            // F_double_tilde_{bc} = Fbc - U_{kl}^{bd}[B^{Q}_{ld}B^{Q}_{kc}]) (Jiang Eq. 30)
            Tensor<double, 2> F_bc_double_tilde = F_ab_t1_[ij];
            auto F_bc_buffer = std::make_shared<Matrix>(n_pno_[ij], n_pno_[ij]);
            F_bc_buffer->zero();
            
            for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
                int k = lmopair_to_lmos_[ij][k_ij];
                for (int l_ij = 0; l_ij < nlmo_ij; ++l_ij) {
                    int l = lmopair_to_lmos_[ij][l_ij];
                    int kl = i_j_to_ij_[k][l];
                    if (kl == -1) continue;
                    
                    SharedMatrix S_ij_kl;
                    if (low_memory_overlap_) { // Compute S_ij_kl using the "semi-direct" algorithm on the fly
                        S_ij_kl = submatrix_cols(*S_ij, index_list(pair_ext_domain, lmopair_to_paos_[kl]));
                        S_ij_kl = linalg::doublet(S_ij_kl, X_pno_[kl], false, false);
                    } else { // Big chunk of memory
                        S_ij_kl = S_PNO(ij, kl);
                    }
                    
                    // Perform PNO projection and copy data
                    auto T_kl_proj = linalg::triplet(S_ij_kl, T_iajb_[kl], S_ij_kl, false, false, true);

                    // R_{ij}^{ab} += t_{kl}^{ab} * beta_{ij}^{kl} (Jiang Eq. 22)
                    T_kl_proj->scale((beta_ij)(k_ij, l_ij));
                    R_iajb[ij]->add(T_kl_proj);
                    if (i != j) R_iajb[ji]->add(T_kl_proj->transpose());

                    // F_double_tilde_{bc} -= U_{kl}^{bd}K_{kl}^{cd} (one can form this intermediate outside loop)
                    auto F_bc_double_tilde_cont = linalg::doublet(Tt_iajb_[kl], K_iajb_[kl], false, true);
                    F_bc_buffer->subtract(linalg::triplet(S_ij_kl, F_bc_double_tilde_cont, S_ij_kl, false, false, true));
                } // end l_ij
            } // end k_ij
            ::memcpy(pno_buffer_a.data(), F_bc_buffer->get_pointer(), n_pno_[ij] * n_pno_[ij] * sizeof(double));
            F_bc_double_tilde += pno_buffer_a;
            
            // R_{ij}^{ab} += P_{ij}^{ab} [t_{ij}^{ac} (F_double_tilde_{bc})] (Jiang Eq. 25)
            // or t_{ij}^{ac} F_double_tilde_{bc} + F_double_tilde_{ac} t_{ij}^{cb}
            einsum(1.0, Indices{index::a, index::b}, &R_ij, 1.0, Indices{index::a, index::c}, T_ij_ein,
                    Indices{index::b, index::c}, F_bc_double_tilde);
            einsum(1.0, Indices{index::a, index::b}, &R_ij, 1.0, Indices{index::a, index::c}, F_bc_double_tilde,
                    Indices{index::c, index::b}, T_ij_ein);
            
            // Flush buffers
            // TODO: Get rid of this copying in the future
            auto R_ij_psi = std::make_shared<Matrix>(npno_ij, npno_ij);
            ::memcpy(R_ij_psi->get_pointer(), R_ij.data(), npno_ij * npno_ij * sizeof(double));
            R_iajb[ij]->add(R_ij_psi);
            if (i != j) R_iajb[ji]->add(R_ij_psi->transpose());
        } // end if
        /*
        // (Projection Error Correction Terms)
        // Projection Error correction terms from C_{ij}^{ab} and D_{ij}^{ab}
        auto dC_ij = std::make_shared<Matrix>(npno_ij, npno_ij);
        auto dD_ij = std::make_shared<Matrix>(npno_ij, npno_ij);
        dC_ij->zero();
        dD_ij->zero();
        for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
            int k = lmopair_to_lmos_[ij][k_ij];
            int jk = i_j_to_ij_[j][k];

            auto S_ij_jk = submatrix_rows_and_cols(*S_pao_, lmopair_to_paos_[ij], lmopair_to_paos_[jk]);
            S_ij_jk = linalg::triplet(X_pno_[ij], S_ij_jk, X_pno_[jk], true, false, false);

            // (ki|ac) correction (Adapted from Jiang Eq. 78)
            // dC_ij->subtract(linalg::triplet(dJ_ij_kj_[ij][k_ij], T_iajb_[jk], S_ij_jk, false, false, true));

            // 2.0 (ai|kc) - (ki|ac) correction (Adapted from Jiang Eq. 79)
            // auto L_aikc = dK_ij_kj_[ij][k_ij]->clone();
            // L_aikc->scale(2.0);
            // L_aikc->subtract(dJ_ij_kj_[ij][k_ij]);
            // auto D_temp = linalg::triplet(L_aikc, Tt_iajb_[jk], S_ij_jk, false, true, true);
            // D_temp->scale(0.5);
            // dD_ij->add(D_temp);
            dD_ij->add(linalg::triplet(dK_ij_kj_[ij][k_ij], Tt_iajb_[jk], S_ij_jk, false, true, true));
        }
        // Add all the C terms to the non-symmetrized R buffer
        auto dC_ij_total = dC_ij->clone();
        dC_ij_total->scale(0.5);
        dC_ij_total->add(dC_ij->transpose());
        Rn_iajb[ij]->add(dC_ij_total);
        // Add the D contribution to the non-symmetrized R buffer
        Rn_iajb[ij]->add(dD_ij);
        */
        
        // C_{ij}^{ab} = -t_{kj}^{bc}[B^{Q}_{ki}B^{Q}_{ac} - 0.5t_{li}^{ad}(B^{Q}_{kd}B^{Q}_{lc})] 
        // (DePrince 2013 Eq. 13, Jiang Eq. 78)
        auto C_ij = std::make_shared<Matrix>(npno_ij, npno_ij);
        for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
            int k = lmopair_to_lmos_[ij][k_ij];
            int ik = i_j_to_ij_[i][k], ki = i_j_to_ij_[k][i], kj = i_j_to_ij_[k][j];

            auto gamma_total = J_ij_kj_[ij][k_ij]->clone();
            gamma_total->add(linalg::triplet(S_PNO(ij, ik), gamma[ki], S_PNO(ik, kj)));
            C_ij->subtract(linalg::triplet(gamma_total, T_iajb_[kj], S_PNO(kj, ij), false, true, false));
        }
        // Add all the C terms to the non-symmetrized R buffer
        auto C_ij_total = C_ij->clone();
        C_ij_total->scale(0.5);
        C_ij_total->add(C_ij->transpose());
        Rn_iajb[ij]->add(C_ij_total);

        // D_{ij}^{ab} = u_{jk}^{bc}(L_{aikc} + 0.5[u_{il}^{ad}L_{ldkc}]) 
        // (DePrince 2013 Eq. 14, Jiang Eq. 79)
        auto D_ij = R_iajb[ij]->clone();
        D_ij->zero();
        for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
            int k = lmopair_to_lmos_[ij][k_ij];
            int ik = i_j_to_ij_[i][k], jk = i_j_to_ij_[j][k];
            int j_ik = lmopair_to_lmos_dense_[ik][j];

            auto U_jk = linalg::triplet(S_PNO(ij, jk), Tt_iajb_[jk], S_PNO(jk, ik));
            auto D_temp = linalg::triplet(S_PNO(ij, ik), delta[ik], U_jk, false, false, true);

            auto L_aikc = K_ij_kj_[ij][k_ij]->clone();
            L_aikc->scale(2.0);
            L_aikc->subtract(J_ij_kj_[ij][k_ij]);
            
            D_temp->add(linalg::triplet(L_aikc, Tt_iajb_[jk], S_PNO(jk, ij), false, true, false));
            D_temp->scale(0.5);
            D_ij->add(D_temp);
        }
        Rn_iajb[ij]->add(D_ij);

        // G_{ij}^{ab} = -t_{ik}^{ab} (Fkj + U_{lj}^{cd}[B^{Q}_{kd}B^{Q}_{lc}]) 
        // (DePrince 2013 Eq. 16, Jiang Eq. 81)
        auto G_ij = R_iajb[ij]->clone();
        G_ij->zero();

        for (int k = 0; k < naocc; ++k) {
            int ik = i_j_to_ij_[i][k];
            if (ik == -1 || n_pno_[ik] == 0) continue;

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

    std::vector<SharedMatrix> tau(n_lmo_pairs);
#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        if (n_pno_[ij] == 0) continue;

        tau[ij] = T_iajb_[ij]->clone();
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
    T_n_ij_sums_.resize(n_lmo_pairs);

    while (!(e_converged && r_converged)) {
        // RMS of residual per single LMO, for assesing convergence
        std::vector<double> R_ia_rms(naocc, 0.0);
        // RMS of residual per LMO pair, for assessing convergence
        std::vector<double> R_iajb_rms(n_lmo_pairs, 0.0);

        std::time_t time_start = std::time(nullptr);

        // Step 1: Create T_n intermediate (Jiang Eq. 70)
#pragma omp parallel for schedule(dynamic, 1)
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            auto &[i, j] = ij_to_i_j_[ij];

            int nlmo_ij = lmopair_to_lmos_[ij].size();
            int npno_ij = n_pno_[ij];
            
            T_n_ij_[ij] = std::make_shared<Matrix>(nlmo_ij, npno_ij);

            for (int n_ij = 0; n_ij < nlmo_ij; ++n_ij) {
                int n = lmopair_to_lmos_[ij][n_ij];
                int nn = i_j_to_ij_[n][n];
                auto T_n_temp = linalg::doublet(S_PNO(ij, nn), T_ia_[n], false, false);
                
                for (int a_ij = 0; a_ij < npno_ij; ++a_ij) {
                    (*T_n_ij_[ij])(n_ij, a_ij) = (*T_n_temp)(a_ij, 0);
                } // end a_ij
            } // end n_ij

            T_n_ij_sums_[ij] = Tensor<double, 2>("T_n_ij", nlmo_ij, npno_ij);
            ::memcpy(T_n_ij_sums_[ij].data(), T_n_ij_[ij]->get_pointer(), nlmo_ij * npno_ij * sizeof(double));
        }

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

            auto S_ij_ii = S_PNO(ij, ii);
            auto S_ij_jj = S_PNO(ij, jj);
            auto tia_temp = linalg::doublet(S_ij_ii, T_ia_[i]);
            auto tjb_temp = linalg::doublet(S_ij_jj, T_ia_[j]);

            for (int a_ij = 0; a_ij < n_pno_[ij]; ++a_ij) {
                for (int b_ij = 0; b_ij < n_pno_[ij]; ++b_ij) {
                    double t1_cont = tia_temp->get(a_ij, 0) * tjb_temp->get(b_ij, 0);
                    double t2_cont = T_iajb_[ij]->get(a_ij, b_ij);

                    tau[ij]->set(a_ij, b_ij, t2_cont + t1_cont);
                }
            }
        }

        // evaluate convergence using current amplitudes and residuals
        e_prev = e_curr;
        // Compute LCCSD energy (Jiang Eq. 45)
        e_curr = 0.0;
        e_weak = 0.0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : e_curr, e_weak)
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            int i, j;
            std::tie(i, j) = ij_to_i_j_[ij];
            double e_ij = tau[ij]->vector_dot(L_iajb_[ij]);
            
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

    outfile->Printf("  Starting Crude Prescreening...\n");
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

    outfile->Printf("  Starting Refined Prescreening...\n");
    outfile->Printf("    T_CUT_MKN reset to %6.3e\n", T_CUT_MKN_);
    outfile->Printf("    T_CUT_DO  reset to %6.3e\n\n", T_CUT_DO_);

    timer_on("Sparsity");
    prep_sparsity(false, false);
    timer_off("Sparsity");

    timer_on("Refined DF Ints");
    compute_qia();
    timer_off("Refined DF Ints");

    timer_on("Refined Pair Prescreening");
    pair_prescreening<false>();
    timer_off("Refined Pair Prescreening");

    timer_on("Sparsity");
    prep_sparsity(false, true);
    timer_off("Sparsity");

    timer_on("PNO-LMP2 Iterations");
    pno_lmp2_iterations();
    timer_off("PNO-LMP2 Iterations");

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

    // Now we do the hard stuff (CCSD)
    recompute_pnos();

    timer_on("DF Ints");
    print_integral_sparsity();
    compute_qij();
    compute_qab();
    timer_off("DF Ints");

    timer_on("CC Integrals");
    estimate_memory();
    compute_cc_integrals();
    timer_off("CC Integrals");

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
    outfile->Printf("   --------------------------------------------\n\n");
    outfile->Printf("  DLPNO convergence set to %s.\n\n", options_.get_str("PNO_CONVERGENCE").c_str());
    outfile->Printf("  Detailed DLPNO thresholds and cutoffs:\n");
    outfile->Printf("    T_CUT_PNO        = %6.4e \n", T_CUT_PNO_);
    outfile->Printf("    T_DIAG_SCALE     = %6.4e \n", T_CUT_PNO_DIAG_SCALE_);
    outfile->Printf("    T_CUT_TRACE      = %6.4e \n", T_CUT_TRACE_);
    outfile->Printf("    T_CUT_ENERGY     = %6.4e \n", T_CUT_ENERGY_);
    outfile->Printf("    T_CUT_PAIRS      = %6.4e \n", T_CUT_PAIRS_);
    outfile->Printf("    MIN_PNOS         = %6d   \n", options_.get_int("MIN_PNOS"));
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
}

}  // namespace dlpno
}  // namespace psi