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

DLPNOCCSD_Lambda::DLPNOCCSD_Lambda(SharedWavefunction ref_wfn, Options& options) : DLPNOCCSD(ref_wfn, options) {}
DLPNOCCSD_Lambda::~DLPNOCCSD_Lambda() {}

void DLPNOCCSD_Lambda::estimate_memory() {

    int n_lmo_pairs = ij_to_i_j_.size();

    outfile->Printf(" ==> CCSD_Lambda Memory Estimate <== \n\n");

    size_t delta_imae_size = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : delta_imae_size)
    for (int im = 0; im < n_lmo_pairs; ++im) {
        auto &[i, m] = ij_to_i_j_[im];
        int ii = i_j_to_ij_[i][i], mm = i_j_to_ij_[m][m];

        delta_imae_size += n_pno_[ii] * n_pno_[mm];
    } // end im

    // Memory Estimate for K_{ma_{ii}}^{e_{mi} f_{mi}} intermediate
    size_t K_maef_dt_size = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : K_maef_dt_size)
    for (int mi = 0; mi < n_lmo_pairs; ++mi) {
        auto &[m, i] = ij_to_i_j_[mi];
        int ii = i_j_to_ij_[i][i];

        int nlmo_mi = lmopair_to_lmos_[mi].size();
        K_maef_dt_size += n_pno_[ii] * n_pno_[mi] * n_pno_[mi];
    } // end mi

    // Memory Estimate for K_{e_{mn} i}^{m n} intermediate
    size_t K_eimn_dt_size = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : K_eimn_dt_size)
    for (int mn = 0; mn < n_lmo_pairs; ++mn) {
        auto &[m, n] = ij_to_i_j_[mn];
        
        int nlmo_mn = lmopair_to_lmos_[mn].size();
        K_eimn_dt_size += nlmo_mn * n_pno_[mn];
    } // end mn

    size_t M_kace_bar_size = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : M_kace_bar_size)
    for (int ki = 0; ki < n_lmo_pairs; ++ki) {
        auto &[k, i] = ij_to_i_j_[ki];
        int ii = i_j_to_ij_[i][i];
        
        M_kace_bar_size += n_pno_[ki] * n_pno_[ki] * n_pno_[ii];
    } // end ki

    size_t M_mkic_bar_size = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : M_mkic_bar_size)
    for (int mk = 0; mk < n_lmo_pairs; ++mk) {
        auto &[m, k] = ij_to_i_j_[mk];
        
        int nlmo_mk = lmopair_to_lmos_[mk].size();
        M_mkic_bar_size += nlmo_mk * n_pno_[mk];
    } // end mk

    size_t F_knia_hat_size = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : F_knia_hat_size)
    for (int kn = 0; kn < n_lmo_pairs; ++kn) {
        auto &[k, n] = ij_to_i_j_[kn];
        
        int nlmo_kn = lmopair_to_lmos_[kn].size();

        for (int i_kn = 0; i_kn < nlmo_kn; ++i_kn) {
            int i = lmopair_to_lmos_[kn][i_kn];
            int ii = i_j_to_ij_[i][i];

            F_knia_hat_size += n_pno_[ii];
        } // end i_kn
    } // end kn

    size_t L_ieab_bar_size = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : L_ieab_bar_size)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];
        int jj = i_j_to_ij_[j][j];

        L_ieab_bar_size += n_pno_[jj] * n_pno_[ij] * n_pno_[ij];
    } // end ij
    
    size_t K_mbij_bar_size = 0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : K_mbij_bar_size)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];

        int nlmo_ij = lmopair_to_lmos_[ij].size();

        K_mbij_bar_size += nlmo_ij * n_pno_[ij];
    } // end ij

    // 1 GB = 1000^3 = 10^9 Bytes
    const double DOUBLES_TO_GB = pow(10.0, -9) * sizeof(double);
    size_t total_size = delta_imae_size + K_maef_dt_size + K_eimn_dt_size + 2 * M_kace_bar_size + 2 * M_mkic_bar_size + F_knia_hat_size + L_ieab_bar_size + K_mbij_bar_size;

    outfile->Printf("     delta_{im}^{a_{ii} e_{mm}}   : %8.3f [GB]\n", delta_imae_size * DOUBLES_TO_GB);
    outfile->Printf("    (a_{ii}, b_{ij}, c_{ij})-like : %8.3f [GB]\n", (K_maef_dt_size + 2 * M_kace_bar_size + L_ieab_bar_size) * DOUBLES_TO_GB);
    outfile->Printf("    (k_{ij}, a_{ij})-like         : %8.3f [GB]\n", (K_eimn_dt_size + 2 * M_mkic_bar_size + K_mbij_bar_size) * DOUBLES_TO_GB);
    outfile->Printf("    F_knia_hat                    : %8.3f [GB]\n", F_knia_hat_size * DOUBLES_TO_GB);
    outfile->Printf("    Total Memory Required         : %8.3f [GB]\n\n", total_size * DOUBLES_TO_GB);

} // end function

void DLPNOCCSD_Lambda::form_goo() {
    // Number of active occupied orbitals
    int naocc = nalpha_ - nfrzc();
    // Number of surviving pairs after DLPNO screening
    int n_lmo_pairs = ij_to_i_j_.size();

    // \rho^{OO}_{nk} = \sum_{m, e, f} \lambda_{mn}^{e_{mn} f_{mn}} [S(e_{mn}, e_{mk}) T_{mk}^{e_{mk} f_{mk}} S(f_{mk}, f_{mn})]
    rho_oo_ = std::make_shared<Matrix>("rho_oo", naocc, naocc);

#pragma omp parallel for schedule(dynamic, 1)
    for (int nk = 0; nk < n_lmo_pairs; ++nk) {
        auto &[n, k] = ij_to_i_j_[nk];

        for (int m_nk = 0; m_nk < lmopair_to_lmos_[nk].size(); ++m_nk) {
            int m = lmopair_to_lmos_[nk][m_nk];
            int mn = i_j_to_ij_[m][n], mk = i_j_to_ij_[m][k];

            // Transform amplitude from domain of mk to mn [S(e_{mn}, e_{mk}) T_{mk}^{e_{mk} f_{mk}} S(f_{mk}, f_{mn})]
            auto T_mk_to_mn = linalg::triplet(S_PNO(mn, mk), T_iajb_[mk], S_PNO(mk, mn)); // (e_{mn}, f_{mn})

            // \rho^{OO}_{nk} += \lambda_{mn}^{e_{mn} f_{mn}} T_mk_to_mn(e_{mn} f_{mn})
            double val = lambda_iajb_[mn]->vector_dot(T_mk_to_mn) + rho_oo_->get(n, k);
            rho_oo_->set(n, k, val);
        } // end m_nk
    } // end nk

    // \rho^{VV}_{f_{mn} c{mn}} = \sum_{e_{mn}} \lambda_{mn}^{e_{mn} f_{mn}} T_{mn}^{e_{mn} c_{mn}}
    rho_vv_.resize(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int mn = 0; mn < n_lmo_pairs; ++mn) {
        auto &[m, n] = ij_to_i_j_[mn];

        rho_vv_[mn] = linalg::doublet(lambda_iajb_[mn], T_iajb_[mn], true, false);
    } // end
}

void DLPNOCCSD_Lambda::compute_lambda_intermediates() {

    outfile->Printf("   ==> Computing Lambda Intermediates <== \n\n");

    // Number of active occupied orbitals
    int naocc = nalpha_ - nfrzc();
    // Number of surviving pairs after DLPNO screening
    int n_lmo_pairs = ij_to_i_j_.size();

    outfile->Printf("   T1-dressing integrals and Fock matrices from converged T1...");

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

    std::time_t time_stop = std::time(nullptr);
    outfile->Printf("   %3d seconds\n\n", (int) time_stop - (int) time_start);

    outfile->Printf("   Computing beta, gamma, and delta from converged T2...");

    time_start = std::time(nullptr);
    
    beta_ = compute_beta();
    gamma_ = compute_gamma();
    delta_ = compute_delta();

    time_stop = std::time(nullptr);
    outfile->Printf("   %3d seconds\n\n", (int) time_stop - (int) time_start);

    outfile->Printf("   Forming K_maef_dt...");

    time_start = std::time(nullptr);
    
    // Toth Eq. 47 (\widetilde{\widetilde{K}}_{m a_{ii}}^{e_{mi} f_{mi}})
    K_maef_dt_.resize(n_lmo_pairs);
    
#pragma omp parallel for
    for (int mi = 0; mi < n_lmo_pairs; ++mi) {
        auto &[m, i] = ij_to_i_j_[mi];
        int im = ij_to_ji_[mi];
        int ii = i_j_to_ij_[i][i];

        int naux_mi = lmopair_to_ribfs_[mi].size();
        int nlmo_mi = lmopair_to_lmos_[mi].size();
        int npno_mi = n_pno_[mi];

        auto q_vo_t1 = i_Qa_t1_[mi];
        auto q_oo_t1 = i_Qk_t1_[mi];
        auto q_ov = QIA_PNO(mi);
        auto q_vv = QAB_PNO(mi);

        K_maef_dt_[mi] = std::make_shared<Matrix>(npno_mi, npno_mi * npno_mi); // (a_{mi}, e_{mi} f_{mi}) -> (a_{ii}, e_{mi} f_{mi}) later
        K_maef_dt_[mi]->zero();

        // (Toth Eq. 47a) +1.0 \widetilde{B}^{Q_{mi}}_{e_{mi} m} \widetilde{B}^{Q_{mi}}_{f_{mi} a_{mi}} S^{a_{mi}}_{a_{ii}}

        for (int q_mi = 0; q_mi < naux_mi; ++q_mi) {
            auto q_vv_t1 = q_vv[q_mi]->clone();
            q_vv_t1->subtract(linalg::doublet(T_n_ij_[mi], q_ov[q_mi], true, false)); // (k_{mi}, f_{mi}) (k_{mi}, a_{mi})

            for (int a_mi = 0; a_mi < n_pno_[mi]; ++a_mi) {
                for (int e_mi = 0; e_mi < n_pno_[mi]; ++e_mi) {
                    for (int f_mi = 0; f_mi < n_pno_[mi]; ++f_mi) {
                        double val = K_maef_dt_[mi]->get(a_mi, e_mi * n_pno_[mi] + f_mi) + q_vv_t1->get(f_mi, a_mi) * q_vo_t1->get(q_mi, e_mi);
                        K_maef_dt_[mi]->set(a_mi, e_mi * n_pno_[mi] + f_mi, val);
                        // (*K_maef_dt_[mi])(a_mi, e_mi * n_pno_[mi] + f_mi) += (*q_vv_t1)(f_mi, a_mi) * (*q_vo_t1)(q_mi, e_mi);
                    } // end f_mi
                } // end e_mi
            } // end a_ii
        } // end q_mi

        K_maef_dt_[mi] = linalg::doublet(S_PNO(ii, mi), K_maef_dt_[mi]); // (a_{mi}, e_{mi} f_{mi}) -> (a_{ii}, e_{mi} f_{mi})

        for (int k_mi = 0; k_mi < nlmo_mi; ++k_mi) {
            int k = lmopair_to_lmos_[mi][k_mi];
            int k_ii = lmopair_to_lmos_dense_[ii][k];
            int mk = i_j_to_ij_[m][k];

            // (Toth Eq. 47c) -1.0 (S^{e_{mi}}_{e_{mk}} T_{mk}^{e_{mk} f_{mk}} S^{f_{mi}}_{f_{mk}}) \overline{F}_{k_{ii} a_{ii}}

            auto glizzy_sticker = linalg::triplet(S_PNO(mi, mk), T_iajb_[mk], S_PNO(mk, mi));

            for (int a_ii = 0; a_ii < n_pno_[ii]; ++a_ii) {
                for (int e_mi = 0; e_mi < n_pno_[mi]; ++e_mi) {
                    for (int f_mi = 0; f_mi < n_pno_[mi]; ++f_mi) {
                        double val = K_maef_dt_[mi]->get(a_ii, e_mi * npno_mi + f_mi) - Fkc_bar_[ii]->get(k_ii, a_ii) * glizzy_sticker->get(e_mi, f_mi);
                        K_maef_dt_[mi]->set(a_ii, e_mi * npno_mi + f_mi, val);
                        // (*K_maef_dt_[mi])(a_ii, e_mi * npno_mi + f_mi) -= (*Fkc_bar_[ii])(k_ii, a_ii) * (*glizzy_sticker)(e_mi, f_mi) * (*Fkc_bar_[ii])(k_ii, a_ii);
                    } // end f_mi
                } // end e_mi
            } // end a_ii

            for (int l_mi = 0; l_mi < nlmo_mi; ++l_mi) {
                int l = lmopair_to_lmos_[mi][l_mi];
                int kl = i_j_to_ij_[k][l];
                if (kl == -1) continue; // checks to make sure kl is a pair

                // (Toth Eq. 47b) +1.0 (S^{e_{mi}}_{e_{kl}} T_{kl}^{e_{kl} f_{kl}} S^{f_{mi}}_{f_{kl}}) \widetilde{B}^{Q_{mi}}_{k_{mi} m} B^{Q_{mi}}_{l_{mi} a_{mi}} S^{a_{mi}}_{a_{ii}}
                auto ender_dragon = linalg::triplet(S_PNO(mi, kl), T_iajb_[kl], S_PNO(kl, mi));

                // \widetilde{B}^{Q_{mi}}_{k_{mi} m} B^{Q_{mi}}_{l_{mi} a_{mi}}
                auto lo_mein = std::make_shared<Matrix>(n_pno_[mi], 1);
                lo_mein->zero();

                for (int q_mi = 0; q_mi < naux_mi; ++q_mi) {
                    // (S^{a_{ii}}_{a_{mi}}) * B^{Q_{mi}}_{l_{mi}a_{mi}} -> (a_{ii}, 1)
                    for (int a_mi = 0; a_mi < n_pno_[mi]; ++a_mi) {
                        double val = lo_mein->get(a_mi, 0) + q_ov[q_mi]->get(l_mi, a_mi) * q_oo_t1->get(q_mi, k_mi);
                        lo_mein->set(a_mi, 0, val);
                    } // end a_mi
                } // end q_mi

                lo_mein = linalg::doublet(S_PNO(ii, mi), lo_mein); // (a_{mi}, 1) -> (a_{ii}, 1)

                for (int a_ii = 0; a_ii < n_pno_[ii]; ++a_ii) {
                    for (int e_mi = 0; e_mi < n_pno_[mi]; ++e_mi) {
                        for (int f_mi = 0; f_mi < n_pno_[mi]; ++f_mi) {
                                double val = K_maef_dt_[mi]->get(a_ii, e_mi * npno_mi + f_mi) + ender_dragon->get(e_mi, f_mi) * lo_mein->get(a_ii, 0);
                                K_maef_dt_[mi]->set(a_ii, e_mi * npno_mi + f_mi, val);
                        } // end f_mi
                    } // end e_mi
                } // end a_ii

            } // end l_mi
        } // end k_mi
    } // end mi

    time_stop = std::time(nullptr);
    outfile->Printf("   %3d seconds\n\n", (int) time_stop - (int) time_start);

    outfile->Printf("   Forming K_eimn_dt...");

    time_start = std::time(nullptr);

    // Toth Eq. 48 \widetilde{\widetilde{K}}_{e_{mn} i}^{m n}

    K_eimn_dt_.resize(n_lmo_pairs);
#pragma omp parallel for schedule(dynamic, 1)
    for (int mn = 0; mn < n_lmo_pairs; ++mn) {
        auto &[m, n] = ij_to_i_j_[mn];
        int nm = ij_to_ji_[mn];

        int naux_mn = lmopair_to_ribfs_[mn].size();
        int nlmo_mn = lmopair_to_lmos_[mn].size();
        int npno_mn = n_pno_[mn];

        auto q_vv = QAB_PNO(mn);
        auto q_ov = QIA_PNO(mn);
        auto q_vo_t1 = i_Qa_t1_[mn];
        auto q_oo_t1 = i_Qk_t1_[nm];

        // (Toth Eq. 48a) +1.0 \widetilde{B}^{Q_{mn}}_{e_{mn} m} \widetilde{B}^{Q_{mn}}_{i_{mn} n}
        K_eimn_dt_[mn] = linalg::doublet(q_vo_t1, q_oo_t1, true, false); // (Q, e) (Q, i) -> (e, i)

        for (int q_mn = 0; q_mn < naux_mn; ++q_mn) {
            auto q_vv_t1 = q_vv[q_mn]->clone();
            q_vv_t1->subtract(linalg::doublet(T_n_ij_[mn], q_ov[q_mn], true, false)); // (k_{mi}, f_{mi}) (k_{mi}, a_{mi})

            // (Toth Eq. 48b) +1.0 \widetilde{B}^{Q_{mn}}_{e_{mn} c_{mn}} T_{mn}^{c_{mn} d_{mn}} B^{Q_{mn}}_{i_{mn} d_{mn}}
            K_eimn_dt_[mn]->add(linalg::triplet(q_vv_t1, T_iajb_[mn], q_ov[q_mn], false, false, true)); // (e, c) (c, d) (i, d) -> (e, i)
        } // end q_mn

        // (Toth Eq. 48c) +1.0 T_{mn}^{e_{mn} c_{mn}} \overline{F}_{i_{mn} c_{mn}}
        int mn_idx = (m < n) ? mn : nm;
        K_eimn_dt_[mn]->add(linalg::doublet(T_iajb_[mn], Fkc_bar_[mn_idx], false, true)); // (e, c) (i, c) -> (e, i)
    } // end mn

    time_stop = std::time(nullptr);
    outfile->Printf("   %3d seconds\n\n", (int) time_stop - (int) time_start);

    outfile->Printf("   Forming M_kace_bar_...");

    time_start = std::time(nullptr);

    // Toth Eq. 29
    M_kace_bar_.resize(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ki = 0; ki < n_lmo_pairs; ++ki) {
        auto &[k, i] = ij_to_i_j_[ki];
        int ii = i_j_to_ij_[i][i];

        int nlmo_ki = lmopair_to_lmos_[ki].size();
        
        M_kace_bar_[ki] = std::make_shared<Matrix>(n_pno_[ki], n_pno_[ki] * n_pno_[ki]);
        
        for (int a_ki = 0; a_ki < n_pno_[ki]; ++a_ki) {
            for (int c_ki = 0; c_ki < n_pno_[ki]; ++c_ki) {
                for (int e_ki = 0; e_ki < n_pno_[ki]; ++e_ki) {
                    double val = 2.0 * K_ivvv_[ki]->get(c_ki, a_ki * n_pno_[ki] + e_ki) - K_ivvv_[ki]->get(a_ki, c_ki * n_pno_[ki] + e_ki);
                    M_kace_bar_[ki]->set(a_ki, c_ki * n_pno_[ki] + e_ki, val);
                    // (*M_kace_bar_[ki])(a_ki, c_ki * n_pno_[ki] + e_ki) = 2.0 * (*K_ivvv_[ki])(c_ki, a_ki * n_pno_[ki] + e_ki)
                            // - (*K_ivvv_[ki])(a_ki, c_ki * n_pno_[ki] + e_ki);
                } // end e_ki
            } // end c_i
        } // end a_ki

        M_kace_bar_[ki] = linalg::doublet(S_PNO(ii, ki), M_kace_bar_[ki]);

        for (int l_ki = 0; l_ki < nlmo_ki; ++l_ki) {
            int l = lmopair_to_lmos_[ki][l_ki];
            int lk = i_j_to_ij_[l][k], ll = i_j_to_ij_[l][l];

            auto forg = linalg::triplet(S_PNO(ii, lk), L_iajb_[lk], S_PNO(lk, ki));
            auto T_l_ki = linalg::doublet(S_PNO(ki, ll), T_ia_[l]);

            for (int a_ii = 0; a_ii < n_pno_[ii]; ++a_ii) {
                for (int c_ki = 0; c_ki < n_pno_[ki]; ++c_ki) {
                    for (int e_ki = 0; e_ki < n_pno_[ki]; ++e_ki) {
                        double val = -forg->get(a_ii, c_ki) * T_l_ki->get(e_ki, 0)
                            + M_kace_bar_[ki]->get(a_ii, c_ki * n_pno_[ki] + e_ki);
                        M_kace_bar_[ki]->set(a_ii, c_ki * n_pno_[ki] + e_ki, val);
                        // (*M_kace_bar_[ki])(a_ii, c_ki * n_pno_[ki] + e_ki) -= (*forg)(a_ii, c_ki) * (*T_l_ki)(e_ki, 0);
                    } // end e_ki
                } // end c_ki
            } // end a_ki
        } // end l_ki

    } // end ki

    time_stop = std::time(nullptr);
    outfile->Printf("   %3d seconds\n\n", (int) time_stop - (int) time_start);

    outfile->Printf("   Forming M_mkic_bar...");

    time_start = std::time(nullptr);

    // Toth Eq. 30
    M_mkic_bar_.resize(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int mk = 0; mk < n_lmo_pairs; ++mk) {
        auto &[m, k] = ij_to_i_j_[mk];
        int mm = i_j_to_ij_[m][m];

        int nlmo_mk = lmopair_to_lmos_[mk].size();
        int npno_mk = n_pno_[mk];

        M_mkic_bar_[mk] = K_mibj_[mk]->clone();
        M_mkic_bar_[mk]->scale(2.0);
        M_mkic_bar_[mk]->subtract(J_ijmb_[mk]);

        for (int i_mk = 0; i_mk < nlmo_mk; ++i_mk) {
            int i = lmopair_to_lmos_[mk][i_mk];
            int ik = i_j_to_ij_[i][k];

            auto T_m_ik = linalg::doublet(S_PNO(ik, mm), T_ia_[m]); // (f_ik, 1)

            auto temp_m = linalg::triplet(T_m_ik, L_iajb_[ik], S_PNO(ik, mk), true, false, false);

            for (int c_mk = 0; c_mk < n_pno_[mk]; ++c_mk) {
                double val = temp_m->get(0, c_mk) + M_mkic_bar_[mk]->get(i_mk, c_mk);
                M_mkic_bar_[mk]->set(i_mk, c_mk, val);
                // (*M_mkic_bar_[mk])(i_mk, c_mk) += (*temp_m)(0, c_mk);
            }
        } // end i_mk
    } // end mk

    time_stop = std::time(nullptr);
    outfile->Printf("   %3d seconds\n\n", (int) time_stop - (int) time_start);

    outfile->Printf("   Forming J_kmic_bar...");

    time_start = std::time(nullptr);

    // Toth Eq. 31
    // \overline{J}_{km}^{ic} = (km | i c_{km}) + \widetilde{T}_{m}^{f_{ki}} (k f_{ki} | i c_{ki})
    //  S_{c_{ki}}^{c_{km}}
    J_kmic_bar_.resize(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int km = 0; km < n_lmo_pairs; ++km) {
        auto &[k, m] = ij_to_i_j_[km];
        int mm = i_j_to_ij_[m][m];

        J_kmic_bar_[km] = J_ijmb_[km]->clone();

        for (int i_km = 0; i_km < lmopair_to_lmos_[km].size(); ++i_km) {
            int i = lmopair_to_lmos_[km][i_km];
            int ki = i_j_to_ij_[k][i];

            auto T_m_ki = linalg::doublet(S_PNO(ki, mm), T_ia_[m]);
            auto f_steak_university = linalg::triplet(T_m_ki, K_iajb_[ki], S_PNO(ki, km), true, false, false); // a column vector of dimension (1, n_pno_[km])

            for (int c_km = 0; c_km < n_pno_[km]; ++c_km) {
                double val = f_steak_university->get(0, c_km) + J_kmic_bar_[km]->get(i_km, c_km);
                J_kmic_bar_[km]->set(i_km, c_km, val);
                // (*J_kmic_bar_[km])(i_km, c_km) += (*f_steak_university)(0, c_km);
            } // end c_km
        } // end i_km
    } // end km

    time_stop = std::time(nullptr);
    outfile->Printf("   %3d seconds\n\n", (int) time_stop - (int) time_start);

    outfile->Printf("   Forming J_kaec_bar...");
    time_start = std::time(nullptr);

    // Toth Eq. 32
    // \overline{J}_{ka_{ii}}^{e_{ki}c_{ki}} = S^{a_{ii}}_{a_{ki}} (ka_{ki}|e_{ki}c_{ki}) -
    // S_{a_{ii}}^{a_{kl}} (k a_{kl} | l c_{kl}) S_{c_{kl}}^{c_{ki}} \widetilde{T}_{l}^{e_{ki}}
    J_kaec_bar_.resize(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ki = 0; ki < n_lmo_pairs; ++ki) {
        auto &[k, i] = ij_to_i_j_[ki];
        int ii = i_j_to_ij_[i][i];

        J_kaec_bar_[ki] = linalg::doublet(S_PNO(ii, ki), K_ivvv_[ki]); // (a, e * c)

        for (int l_ki = 0; l_ki < lmopair_to_lmos_[ki].size(); ++l_ki) {
            int l = lmopair_to_lmos_[ki][l_ki];
            int kl = i_j_to_ij_[k][l], ll = i_j_to_ij_[l][l];

            // S_{a_{ii}}^{a_{kl}} (k a_{kl} | l c_{kl}) S_{c_{kl}}^{c_{ki}}
            auto south_ohio = linalg::triplet(S_PNO(ii, kl), K_iajb_[kl], S_PNO(kl, ki));
            auto T_l_ki = linalg::doublet(S_PNO(ki, ll), T_ia_[l]); // \widetilde{T}_{l}^{e_{ki}} (Jiang Eq. 70?)

            for (int a_ii = 0; a_ii < n_pno_[ii]; ++a_ii) {
                for (int e_ki = 0; e_ki < n_pno_[ki]; ++e_ki) {
                    for (int c_ki = 0; c_ki < n_pno_[ki]; ++c_ki) {
                        double val = J_kaec_bar_[ki]->get(a_ii, e_ki * n_pno_[ki] + c_ki) - south_ohio->get(a_ii, c_ki) * T_l_ki->get(e_ki, 0);
                        J_kaec_bar_[ki]->set(a_ii, e_ki * n_pno_[ki] + c_ki, val);
                    } // end c_ki
                } // end e_ki
            } // end a_ii

        } // end l_ki
    } // end for ki

    time_stop = std::time(nullptr);
    outfile->Printf("   %3d seconds\n\n", (int) time_stop - (int) time_start);

    // Toth Eq. 33
    /*
    F_fcia_hat_.resize(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int mn = 0; mn < n_lmo_pairs; ++mn) {
        auto &[m, n] = ij_to_i_j_[mn];

        int naux_mn = lmopair_to_ribfs_[mn].size();
        int nlmo_mn = lmopair_to_lmos_[mn].size();
        int npno_mn = n_pno_[mn];

        F_fcia_hat_[mn].resize(nlmo_mn);

        auto qia_mn = QIA_PNO(mn); // naux_mn * (nlmo_mn, npno_mn)
        auto qab_mn = QAB_PNO(mn); // naux_mn * (npno_mn, npno_mn)

        SharedMatrix F_fcia_temp = std::make_shared<Matrix>("F_fcia_temp", nlmo_mn * npno_mn, npno_mn * npno_mn);

        // Eq 33a and b
        // TODO: Add tilde to B^{Q}_{fc} term
        for (int q_mn = 0; q_mn < naux_mn; ++q_mn) {
            auto q_vv = qab_mn[q_mn]->clone(); // (npno_mn, npno_mn)
            auto q_ov = qia_mn[q_mn]->clone(); // (nlmo_mn, npno_mn)

            // This performs a "T1-dressing" on Qvv, 
            // B^{Q_{mn}}_{f_{mn}c_{mn}} -= \widetilde{T}_{k_{mn}}^{f_{mn}} B^{Q_{mn}}_{k_{mn}c_{mn}}
            q_vv->subtract(linalg::doublet(T_n_ij_[mn], q_ov, true, false));

            for (int i_mn = 0; i_mn < nlmo_mn; ++i_mn) {
                for (int a_mn = 0; a_mn < npno_mn; ++a_mn) {
                    for (int f_mn = 0; f_mn < npno_mn; ++f_mn) {
                        for (int c_mn = 0; c_mn < npno_mn; ++c_mn) {
                            double val = 2.0 * q_vv->get(f_mn, c_mn) * q_ov->get(i_mn, a_mn) - q_vv->get(f_mn, a_mn) * q_ov->get(i_mn, c_mn);
                            F_fcia_temp->set(i_mn * npno_mn + a_mn, f_mn * npno_mn + c_mn, 
                                                F_fcia_temp->get(i_mn * npno_mn + a_mn, f_mn * npno_mn + c_mn) + val);
                        } // end c_mn
                    } // end f_mn
                } // end a_mn
            } // end i_mn
        }

        // Change the dimensions of F_fcia_temp to be (i_mn, a_mn * f_mn * c_mn)
        F_fcia_temp->reshape(nlmo_mn, npno_mn * npno_mn * npno_mn); // (i_mn, a_mn * f_mn * c_mn)

        // Package it up into F_fcia_hat_[mn][i_mn]
        for (int i_mn = 0; i_mn < nlmo_mn; ++i_mn) {
            int i = lmopair_to_lmos_[mn][i_mn];
            int ii = i_j_to_ij_[i][i];

            // This is the slice of F_fcia_temp corresponding to the current i_mn
            F_fcia_hat_[mn][i_mn] = submatrix_rows(*F_fcia_temp, std::vector<int>(1, i_mn)); // (1, a_mn * f_mn * c_mn)
            F_fcia_hat_[mn][i_mn]->reshape(npno_mn, npno_mn * npno_mn); // (a_mn, f_mn * c_mn)
            F_fcia_hat_[mn][i_mn] = linalg::doublet(S_PNO(ii, mn), F_fcia_hat_[mn][i_mn]); // (a_ii, f_mn * c_mn)
        } // end i_mn
    } // end mn
     */

    outfile->Printf("   Forming F_knia_hat...");
    time_start = std::time(nullptr);

    // Toth Eq. 34a
    F_knia_hat_.resize(n_lmo_pairs);
#pragma omp parallel for schedule(dynamic, 1)
    for (int kn = 0; kn < n_lmo_pairs; ++kn) {
        auto &[k, n] = ij_to_i_j_[kn];
        int nk = ij_to_ji_[kn];

        int naux_kn = lmopair_to_ribfs_[kn].size();
        int nlmo_kn = lmopair_to_lmos_[kn].size();
        int npno_kn = n_pno_[kn];

        F_knia_hat_[kn].resize(nlmo_kn);

        auto F_knia_temp = J_ijmb_[kn]->clone(); // (i_kn, a_kn)
        F_knia_temp->scale(2.0);
        F_knia_temp->subtract(K_mibj_[nk]); // (i_kn, a_kn)

        for (int i_kn = 0; i_kn < nlmo_kn; ++i_kn) {
            int i = lmopair_to_lmos_[kn][i_kn];
            int ii = i_j_to_ij_[i][i];

            F_knia_hat_[kn][i_kn] = submatrix_rows(*F_knia_temp, std::vector<int>(1, i_kn)); // (1, a_kn)
            F_knia_hat_[kn][i_kn]->reshape(npno_kn, 1); // (a_kn, 1)
            F_knia_hat_[kn][i_kn] = linalg::doublet(S_PNO(ii, kn), F_knia_hat_[kn][i_kn]); // (a_ii, 1)
        } // end i_kn
    } // end kn

    // Toth Eq. 34b
#pragma omp parallel for schedule(dynamic, 1)
    for (int n = 0; n < naocc; ++n) {
        int nn = i_j_to_ij_[n][n];

        int naux_nn = lmopair_to_ribfs_[nn].size(); // Number of auxiliary functions in domain of nn
        int nlmo_nn = lmopair_to_lmos_[nn].size(); // Number of LMOs in domain of pair nn
        int npno_nn = n_pno_[nn];

        auto qov_nn = QIA_PNO(nn); // naux_nn * (nlmo_nn, npno_nn)

        for (int q_nn = 0; q_nn < naux_nn; ++q_nn) {
            auto qov = qov_nn[q_nn]->clone(); // (nlmo_nn, npno_nn)
            auto qov_contracted = linalg::doublet(qov, T_ia_[n], false, false); // (nlmo_nn, npno_nn) * (npno_nn, 1) -> (nlmo_nn, 1)

            for (int k_nn = 0; k_nn < nlmo_nn; ++k_nn) {
                int k = lmopair_to_lmos_[nn][k_nn];
                int kn = i_j_to_ij_[k][n];
                for (int i_nn = 0; i_nn < nlmo_nn; ++i_nn) {
                    int i = lmopair_to_lmos_[nn][i_nn];
                    int ii = i_j_to_ij_[i][i];
                    int i_kn = lmopair_to_lmos_dense_[kn][i]; // Index of i in the domain of kn
                    if (i_kn == -1) continue; // If i is not in the domain of kn, skip

                    auto F_knia_temp = std::make_shared<Matrix>("F_knia_temp", npno_nn, 1);

                    for (int a_nn = 0; a_nn < npno_nn; ++a_nn) {
                        double val = 2.0 * qov_contracted->get(k_nn, 0) * qov->get(i_nn, a_nn) - 
                            qov_contracted->get(i_nn, 0) * qov->get(k_nn, a_nn);
                        F_knia_temp->set(a_nn, 0, F_knia_temp->get(a_nn, 0) + val);
                    }

                    F_knia_hat_[kn][i_kn]->add(linalg::doublet(S_PNO(ii, nn), F_knia_temp)); // (a_ii, 1)
                } // end i_nn
            } // end k_nn
        } // q_nn

    } // end n

    time_stop = std::time(nullptr);
    outfile->Printf("   %3d seconds\n\n", (int) time_stop - (int) time_start);

    outfile->Printf("   Forming L_ieab_bar_ and K_ijmb_bar_...");
    time_start = std::time(nullptr);

    L_ieab_bar_.resize(n_lmo_pairs);
    K_ijmb_bar_.resize(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];
        int jj = i_j_to_ij_[j][j];

        int naux_ij = lmopair_to_ribfs_[ij].size();
        int nlmo_ij = lmopair_to_lmos_[ij].size();
        int npno_ij = n_pno_[ij];

        // Toth Eq. 41 (S_{e_{ij}}^{e_{jj}} L_{ie_{ij}}^{a_{ij}b_{ij}}) 
        // - \widetilde{T}_{l}^{e_{jj}} [S_{a_{il}}^{a_{ij}}L_{il}^{a_{il}b_{il}}]S_{b_{il}}^{b_{ij}}
        L_ieab_bar_[ij] = std::make_shared<Matrix>("L_ieab_bar", npno_ij, npno_ij * npno_ij); // (e_ij, a_ij * b_ij)

        for (int a_ij = 0; a_ij < npno_ij; ++a_ij) {
            for (int e_ij = 0; e_ij < npno_ij; ++e_ij) {
                for (int b_ij = 0; b_ij < npno_ij; ++b_ij) {
                    // L_{ie_{ij}}^{a_{ij} b_{ij}} = 2 (i a_{ij} | e_{ij} b_{ij}) - (i b_{ij} | e_{ij} a_{ij})
                    double val = 2.0 * K_ivvv_[ij]->get(a_ij, e_ij * n_pno_[ij] + b_ij) - K_ivvv_[ij]->get(b_ij, e_ij * n_pno_[ij] + a_ij);
                    L_ieab_bar_[ij]->set(e_ij, a_ij * npno_ij + b_ij, val);
                } // end b_ij
            } // end e_ij
        } // end a_ij

        L_ieab_bar_[ij] = linalg::doublet(S_PNO(jj, ij), L_ieab_bar_[ij]); // (e_ij, a_ij * b_ij) -> (e_jj, a_ij * b_ij)

        for (int l_ij = 0; l_ij < nlmo_ij; ++l_ij) {
            int l = lmopair_to_lmos_[ij][l_ij];
            int il = i_j_to_ij_[i][l], ll = i_j_to_ij_[l][l];

            auto T_l_j = linalg::doublet(S_PNO(jj, ll), T_ia_[l]); // (e_ll, 1) -> (e_jj, 1)
            auto jon_arbuckle = linalg::triplet(S_PNO(ij, il), L_iajb_[il], S_PNO(il, ij)); // (a_ij, b_ij)
            
            for (int e_jj = 0; e_jj < n_pno_[jj]; ++e_jj) {
                for (int a_ij = 0; a_ij < n_pno_[ij]; ++a_ij) {
                    for (int b_ij = 0; b_ij < n_pno_[ij]; ++b_ij) {
                        double val = (L_ieab_bar_[ij])->get(e_jj, a_ij * n_pno_[ij] + b_ij) - T_l_j->get(e_jj, 0) * jon_arbuckle->get(a_ij, b_ij);
                        L_ieab_bar_[ij]->set(e_jj, a_ij * n_pno_[ij] + b_ij, val);
                    } // end b_ij
                } // end a_ij
            } // end e_jj
        } // end l_ij

        // Toth Eq. 42 K_{ij}^{mb_{ij}} + \widetilde{T}_m^{e_{ij}} K_{ij}^{e_{ij}b_{ij}}
        K_ijmb_bar_[ij] = K_mibj_[ij]->clone();
        K_ijmb_bar_[ij]->add(linalg::doublet(T_n_ij_[ij], K_iajb_[ij])); // (m, e) (e, b) -> (m, b)
    } // end ij

    time_stop = std::time(nullptr);
    outfile->Printf("   %3d seconds\n\n", (int) time_stop - (int) time_start);

    outfile->Printf("   Forming delta_imae_tilde ...");
    time_start = std::time(nullptr);

    // Toth Eq. 44 and 45
    delta_imae_tilde_.resize(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int im = 0; im < n_lmo_pairs; ++im) {
        auto &[i, m] = ij_to_i_j_[im];
        int ii = i_j_to_ij_[i][i], mm = i_j_to_ij_[m][m];
        int i_mm = lmopair_to_lmos_dense_[mm][i]; // Index of LMO (occupied) i in the domain of mm

        int nlmo_im = lmopair_to_lmos_[im].size();

        // The relavent density-fitted two-electron integrals needed
        auto Qov = QIA_PNO(mm); // (Q_{mm} | i_{mm} a_{mm})
        auto Qvv = QAB_PNO(mm); // (Q_{mm} | a_{mm} e_{mm})
        auto Qme = i_Qa_ij_[mm]; // \tilde(Q_{mm} | m e_{mm})
        auto Qmi = i_Qk_ij_[mm]; // \tilde(Q_{mm} | m i_{mm})

        // => Toth Eq. 44a S(a_{ii}, a_{mm}) M_{im}^{a_{mm}e_{mm}} <= //
        
        auto M_imae_tilde = std::make_shared<Matrix>(n_pno_[mm], n_pno_[mm]);
        M_imae_tilde->zero();
        
        // M_{im}^{a_{mm}e_{mm}} = 2(ia_{mm}|me_{mm}) - (im|a_{mm}e_{mm})
        for (int q_mm = 0; q_mm < lmopair_to_ribfs_[mm].size(); ++q_mm) {
            for (int a_mm = 0; a_mm < n_pno_[mm]; ++a_mm) {
                for (int e_mm = 0; e_mm < n_pno_[mm]; ++e_mm) {
                    double val = 2.0 * Qov[q_mm]->get(i_mm, a_mm) * Qme->get(q_mm, e_mm)
                                    - Qmi->get(q_mm, i_mm) * Qvv[q_mm]->get(a_mm, e_mm);
                    M_imae_tilde->set(a_mm, e_mm, M_imae_tilde->get(a_mm, e_mm) + val);
                } // end e_mm
            } // end a_mm
        } // end q_mm

        M_imae_tilde = linalg::doublet(S_PNO(ii, mm), M_imae_tilde);

        // => Toth Eq. 44c S(a_{ii}, a_{mm}) M_{if_{mm}}^{a_{mm}e_{mm}} T_{m}^{f_{mm}} <= //

        // M_{if_{mm}}^{a_{mm}e_{mm}} = 2 (ia_{mm} | f_{mm} e_{mm}) - (if_{mm} | a_{mm} e_{mm}) (24c)
        auto M_ifae = std::make_shared<Matrix>(n_pno_[mm], n_pno_[mm] * n_pno_[mm]);

        for (int q_mm = 0; q_mm < lmopair_to_ribfs_[mm].size(); ++q_mm) {
            for (int f_mm = 0; f_mm < n_pno_[mm]; ++f_mm) {
                for (int a_mm = 0; a_mm < n_pno_[mm]; ++a_mm) {
                    for (int e_mm = 0; e_mm < n_pno_[mm]; ++e_mm) {
                        double val = 2.0 * Qov[q_mm]->get(i_mm, a_mm) * Qvv[q_mm]->get(f_mm, e_mm)
                                        - Qov[q_mm]->get(i_mm, f_mm) * Qvv[q_mm]->get(a_mm, e_mm);
                        M_ifae->set(f_mm, a_mm * n_pno_[mm] + e_mm, M_ifae->get(f_mm, a_mm * n_pno_[mm] + e_mm) + val);
                    } // end e_mm
                } // end a_mm
            } // end f_mm
        } // end q_mm
        M_ifae = linalg::doublet(T_ia_[m], M_ifae, true, false); // (f_mm, 1) x (f_mm, a_mm * e_mm) -> (1, a_mm * e_mm) 
        M_ifae->reshape(n_pno_[mm], n_pno_[mm]); // (1, a_mm * e_mm) -> (a_mm, e_mm)
        M_imae_tilde->add(linalg::doublet(S_PNO(ii, mm), M_ifae)); // (a_ii, a_mm) x (a_mm, e_mm) -> (a_ii, e_mm)

        // Toth Eq. 44d -\widetilde{T}^{e_{mm}}_k L^{a_{mm}f_{mm}}_{ik} S^{a_{mm}}_{a_{ii}} T^{f_{mm}}_m

        // L_{ik}^{a_{mm}f_{mm}} = 2(ia_{mm} | kf_{mm})  - (if_{mm} | ka_{mm})
        for (int k_im = 0; k_im < nlmo_im; ++k_im) {
            int k = lmopair_to_lmos_[im][k_im];
            int mk = i_j_to_ij_[m][k], kk = i_j_to_ij_[k][k];
            int i_mk = lmopair_to_lmos_dense_[mk][i];
            int k_mm = lmopair_to_lmos_dense_[mm][k];

            // \widetilde{T}^{e_{mm}}_k
            auto T_k_m = submatrix_rows(*T_n_ij_[mm], std::vector<int>(1, k_mm)); // (1, e_mm)

            auto L_ikaf = std::make_shared<Matrix>(n_pno_[mm], n_pno_[mm]);
            L_ikaf->zero();

            for (int q_mm = 0; q_mm < lmopair_to_ribfs_[mm].size(); ++q_mm) {
                for (int a_mm = 0; a_mm < n_pno_[mm]; ++a_mm) {
                    for (int f_mm = 0; f_mm < n_pno_[mm]; ++f_mm) {
                        double val = 2.0 * Qov[q_mm]->get(i_mm, a_mm) * Qov[q_mm]->get(k_mm, f_mm)
                            - Qov[q_mm]->get(i_mm, f_mm) * Qov[q_mm]->get(k_mm, a_mm);
                        
                        L_ikaf->set(a_mm, f_mm, L_ikaf->get(a_mm, f_mm) + val);
                    } // end f_mm
                } // end a_mm
            } // end q_mm

            // S_{a_{ii}}^{a_{mm}} L^{a_{mm}f_{mm}}_{ik}  T^{f_{mm}}_m
            L_ikaf = linalg::triplet(S_PNO(ii, mm), L_ikaf, T_ia_[m]); // (a_ii, a_mm) x (a_mm, f_mm) x (f_mm, 1) -> (a_ii, 1)

            for (int a_ii = 0; a_ii < n_pno_[ii]; ++a_ii) {
                for (int e_mm = 0; e_mm < n_pno_[mm]; ++e_mm) {
                    double val = L_ikaf->get(a_ii, 0) * T_k_m->get(0, e_mm);
                    M_imae_tilde->set(a_ii, e_mm, M_imae_tilde->get(a_ii, e_mm) - val);
                } // end e_mm
            } // end a_ii
        } // end k_im

        // => Toth Eq. 45b S^{a_{ii}}_{a_{ik}} L_{ik}^{a_{ik}c_{ik}} S^{c_{ik}}_{c_{km}} U_{km}^{c_{km}e_{km}} S^{e_{km}}_{e_{mm}}

        // Compute delta_imae_tilde
        delta_imae_tilde_[im] = M_imae_tilde->clone();

        for (int k_im = 0; k_im < lmopair_to_lmos_[im].size(); ++k_im) {
            int k = lmopair_to_lmos_[im][k_im];
            int ik = i_j_to_ij_[i][k], km = i_j_to_ij_[k][m];
            
            auto bear = linalg::triplet(L_iajb_[ik], S_PNO(ik, km), Tt_iajb_[km]);
            delta_imae_tilde_[im]->add(linalg::triplet(S_PNO(ii, ik), bear, S_PNO(km, mm)));
        }
        
    } // end im

#pragma omp parallel for schedule(dynamic, 1)
    for (int mk = 0; mk < n_lmo_pairs; ++mk) {
        auto &[m, k] = ij_to_i_j_[mk];
        int mm = i_j_to_ij_[m][m], kk = i_j_to_ij_[k][k];
        int nlmo_mk = lmopair_to_lmos_[mk].size();

        // => Toth Eq. 44b -\widetilde{T}_{k}^{e_{mm}} N_{mk}^{i a_{mk}} S(a_{mk}, a_{ii}) <= //

        // T_k_m = \widetilde{T}_{k}^{e_{mm}}
        auto T_k_m = linalg::doublet(S_PNO(mm, kk), T_ia_[k]); // (e_{mm}, 1)

        auto N_mkia = J_ijmb_[mk]->clone(); // (i_{mk}, a_{mk})
        N_mkia->scale(2.0);
        N_mkia->subtract(K_mibj_[mk]); // (i_{mk}, a_{mk})

        // N_{mk}^{ia_{mk}} = 2 (mk | ia_{mk}) - (mi | ka_{mk})
        for (int i_mk = 0; i_mk < nlmo_mk; ++i_mk) {
            int i = lmopair_to_lmos_[mk][i_mk];
            int ii = i_j_to_ij_[i][i], im = i_j_to_ij_[i][m];
            
            auto N_mkia_slice = submatrix_rows(*N_mkia, std::vector<int>(1, i_mk))->transpose();
            N_mkia_slice = linalg::doublet(S_PNO(ii, mk), N_mkia_slice);

            for (int a_ii = 0; a_ii < n_pno_[ii]; ++a_ii) {
                for (int e_mm = 0; e_mm < n_pno_[mm]; ++e_mm) {
#pragma omp atomic
                    (*delta_imae_tilde_[im])(a_ii, e_mm) -= N_mkia_slice->get(a_ii, 0) * T_k_m->get(e_mm, 0);
                } // end e_mm
            } // end a_ii
        } // end k_im
    }

    time_stop = std::time(nullptr);
    outfile->Printf("   %3d seconds\n\n", (int) time_stop - (int) time_start);

    outfile->Printf("   Forming F_vv_double_tilde...");
    time_start = std::time(nullptr);

    F_vv_double_tilde_.resize(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];
        int ji = ij_to_ji_[ij];

        int nlmo_ij = lmopair_to_lmos_[ij].size();
        int ij_idx = (i > j) ? ji : ij;

        F_vv_double_tilde_[ij] = Fab_[ij_idx]->clone();

        for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
            int k = lmopair_to_lmos_[ij][k_ij];
    
            for (int l_ij = 0; l_ij < nlmo_ij; ++l_ij) {
                int l = lmopair_to_lmos_[ij][l_ij];
                int kl = i_j_to_ij_[k][l];
                if (kl == -1) continue;

                auto chilly_glizzy = linalg::doublet(Tt_iajb_[kl], K_iajb_[kl], false, true);
                F_vv_double_tilde_[ij]->subtract(linalg::triplet(S_PNO(ij, kl), chilly_glizzy, S_PNO(kl, ij)));

            } // end l_ij
        } // end k_ij
    } // end ij

    time_stop = std::time(nullptr);
    outfile->Printf("   %3d seconds\n\n", (int) time_stop - (int) time_start);

    outfile->Printf("   Forming F_im_double_tilde...");
    time_start = std::time(nullptr);
    
    F_im_double_tilde_ = Fkj_->clone();

#pragma omp parallel for schedule(dynamic, 1)
    for (int im = 0; im < n_lmo_pairs; ++im) {
        auto &[i, m] = ij_to_i_j_[im];

        int nlmo_im = lmopair_to_lmos_[im].size();
        
        for (int l_im = 0; l_im < nlmo_im; ++l_im) {
            int l = lmopair_to_lmos_[im][l_im];
            int li = i_j_to_ij_[l][i], lm = i_j_to_ij_[l][m];

            auto sussy_baka = linalg::triplet(S_PNO(li, lm), Tt_iajb_[lm], S_PNO(li, lm), false, false, true);
            double val = F_im_double_tilde_->get(i, m) + sussy_baka->vector_dot(K_iajb_[li]);
            F_im_double_tilde_->set(i, m, val);
        } // end l_im
    } // end im

    time_stop = std::time(nullptr);
    outfile->Printf("   %3d seconds\n\n", (int) time_stop - (int) time_start);

    outfile->Printf("   Forming gamma_double_tilde and delta_double_tilde...");
    time_start = std::time(nullptr);

    gamma_double_tilde_.resize(n_lmo_pairs);
    delta_double_tilde_.resize(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];

        int nlmo_ij = lmopair_to_lmos_[ij].size();

        gamma_double_tilde_[ij].resize(nlmo_ij);
        delta_double_tilde_[ij].resize(nlmo_ij);

        for (int n_ij = 0; n_ij < nlmo_ij; ++n_ij) {
            int n = lmopair_to_lmos_[ij][n_ij];
            int in = i_j_to_ij_[i][n], jn = i_j_to_ij_[j][n], nj = i_j_to_ij_[n][j];
            int i_nj = lmopair_to_lmos_dense_[nj][i];

            gamma_double_tilde_[ij][n_ij] = std::make_shared<Matrix>(n_pno_[jn], n_pno_[ij]);
            gamma_double_tilde_[ij][n_ij]->zero();

            for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
                int k = lmopair_to_lmos_[ij][k_ij];
                int kn = i_j_to_ij_[k][n], ik = i_j_to_ij_[i][k], kj = i_j_to_ij_[k][j];
                if (kn == -1) continue;

                auto T_kn = linalg::triplet(S_PNO(jn, kn), T_iajb_[kn], S_PNO(kn, ik));
                auto K_ik = linalg::doublet(K_iajb_[ik], S_PNO(ik, ij));

                gamma_double_tilde_[ij][n_ij]->add(linalg::doublet(T_kn, K_ik));
                
            } // end k_ij

        } // end n_ij

        for (int n_ij = 0; n_ij < nlmo_ij; ++n_ij) {
            int n = lmopair_to_lmos_[ij][n_ij];
            int nj = i_j_to_ij_[n][j], in = i_j_to_ij_[i][n], ni = i_j_to_ij_[n][i];
            int j_ni = lmopair_to_lmos_dense_[ni][j];

            delta_double_tilde_[ij][n_ij] = std::make_shared<Matrix>(n_pno_[ni], n_pno_[ij]);
            delta_double_tilde_[ij][n_ij]->zero();
            
            for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
                int k = lmopair_to_lmos_[ij][k_ij];
                int nk = i_j_to_ij_[n][k], ik = i_j_to_ij_[i][k], kj = i_j_to_ij_[k][j];
                if (nk == -1) continue;

                auto U_nk = linalg::triplet(S_PNO(ni, nk), Tt_iajb_[nk], S_PNO(nk, kj));
                auto L_kj = linalg::doublet(L_iajb_[kj], S_PNO(kj, ij));

                delta_double_tilde_[ij][n_ij]->add(linalg::doublet(U_nk, L_kj));
            } // end k_ij
        }
    } // end ij

    time_stop = std::time(nullptr);
    outfile->Printf("   %3d seconds\n\n", (int) time_stop - (int) time_start);
}

void DLPNOCCSD_Lambda::compute_L_ia(std::vector<SharedMatrix>& L_ia, std::vector<std::vector<SharedMatrix>> &L_ia_buffer) {

    timer_on("DLPNO-CCSD Lambda : Compute L1");

    int naocc = nalpha_ - nfrzc();
    int n_lmo_pairs = ij_to_i_j_.size();

    // Thread and OMP Parallel Info
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif

    // Zero out buffers
    for (int thread = 0; thread < nthreads; ++thread) {
        for (int i = 0; i < naocc; ++i) {
            L_ia_buffer[thread][i]->zero();
        }
    }

    // \lambda_{i}^{e_{ii}} \widetilde{\widetilde{F}}_(e_{ii}, a_{ii}) (Toth Eq. 46c)
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < naocc; ++i) {
        int ii = i_j_to_ij_[i][i];

        L_ia[i] = linalg::doublet(F_vv_double_tilde_[ii], lambda_ia_[i], true, false); // (e, a) x (e, 1) -> (a, 1)
    } // end for

    // + 2.0 * (S^{a_{ii}}_{a_{in}} L_{in}^{a_{in}f_{in}} S^{f_{in}}_{f_{nn}}) T_{n}^{f_{nn}}
#pragma omp parallel for schedule(dynamic, 1)
    for (int in = 0; in < n_lmo_pairs; ++in) {
        auto &[i, n] = ij_to_i_j_[in];
        int ii = i_j_to_ij_[i][i], nn = i_j_to_ij_[n][n];

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif
        /*
        auto chicken = linalg::triplet(S_PNO(ii, in), L_iajb_[in], S_PNO(nn, in), false, false, true);
        auto sandwich = linalg::doublet(chicken, T_ia_[n]);
        sandwich->scale(2.0);
        */

        auto chicken = linalg::triplet(L_iajb_[in], S_PNO(in, nn), T_ia_[n]);
        auto sandwich = linalg::doublet(S_PNO(ii, in), chicken);
        sandwich->scale(2.0);

        L_ia_buffer[thread][i]->add(sandwich);
    } // end in

#pragma omp parallel for schedule(dynamic, 1)
    for (int im = 0; im < n_lmo_pairs; ++im) {
        auto &[i, m] = ij_to_i_j_[im];
        int ii = i_j_to_ij_[i][i], mm = i_j_to_ij_[m][m], mi = i_j_to_ij_[m][i];

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        // \widetilde{\delta}_{im}^{a_{ii} e_{mm}} \lambda_{m}^{e_{mm}} (Toth Eq. 46a)
        L_ia_buffer[thread][i]->add(linalg::doublet(delta_imae_tilde_[im], lambda_ia_[m]));

        // - \widetilde{widetilde{F}}_{im} \lambda_{m}^{a_{mm}} S_{a_{mm}}^{a_{ii}} (Toth Eq. 46b)
        auto john_big_back_buffer = linalg::doublet(S_PNO(ii, mm), lambda_ia_[m]);
        john_big_back_buffer->scale(F_im_double_tilde_->get(i, m));
        L_ia_buffer[thread][i]->subtract(john_big_back_buffer);
        
        /* l_{i}^{a_{ii}} += \widetilde{\widetilde{K}}_{ma_{ii}}^{e_{mi}f_{mi}} \widetilde{\lambda}_{mi}^{e_{mi}f_{mi}} (Toth Eq. 55a) */
        auto lambda_mi_slice = lambda_iajb_[mi]->clone();
        lambda_mi_slice->reshape(n_pno_[mi] * n_pno_[mi], 1);
        L_ia_buffer[thread][i]->add(linalg::doublet(K_maef_dt_[mi], lambda_mi_slice)); // (a_{ii}, e_{mi} * f_{mi}) (e_{mi} * f_{mi}, 1)
    } // end im

    // rip John Pork :( bigback = ibs (irritable bigback syndrome)
    //    ^     ^
    // (  *     *  )
    // (  ( * * )  )
    // (           )
    // OINK!!! John Pork is calling...
#pragma omp parallel for schedule(dynamic, 1)
    for (int mn = 0; mn < n_lmo_pairs; ++mn) {
        auto &[m, n] = ij_to_i_j_[mn];

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        int nlmo_mn = lmopair_to_lmos_[mn].size();
	    int naux_mn = lmopair_to_ribfs_[mn].size();

        auto qia_mn = QIA_PNO(mn); // naux_mn * (nlmo_mn, npno_mn)
        auto qab_mn = QAB_PNO(mn); // naux_mn * (npno_mn, npno_mn)

        /* l_{i}^{a_{ii}} -= \widetilde{\widetilde{K}}_{e_{mn}i}^{mn} \widetilde{\lambda}_{mn}^{e_{mn}a_{mn}}S_{a_{mn}}^{a_{ii}} (Toth Eq. 55b) */
        auto bruvver = linalg::doublet(K_eimn_dt_[mn], lambda_iajb_[mn], true, false); // (e_{mn}, i_{mn}) x (e_{mn}, a_{mn}) -> (i_{mn}, a_{mn})

        for (int i_mn = 0; i_mn < lmopair_to_lmos_[mn].size(); ++i_mn) {
            int i = lmopair_to_lmos_[mn][i_mn];
            int ii = i_j_to_ij_[i][i];
            auto antonios_slice = submatrix_rows(*bruvver, std::vector<int>(1, i_mn));
            L_ia_buffer[thread][i]->subtract(linalg::doublet(S_PNO(ii, mn), antonios_slice, false, true)); // (a_{ii}, a_{mn}) x (1, a_{mn}) -> (a_{ii}, 1)
        }

        // l_{i}^{a_{ii}} \mathrel{+}= \rho^{\mathrm{VV}}_{f_{mn}c_{mn}}\hat{F}^{ia_{ii}}_{f_{mn}c_{mn}} - \rho^{\mathrm{OO}}_{nm} \hat{F}_{mn}^{ia_{ii}} (Toth Eq. 56a)
        for (int q_mn = 0; q_mn < naux_mn; ++q_mn) {
            auto q_vv = qab_mn[q_mn]->clone(); // (npno_mn, npno_mn)
            auto q_ov = qia_mn[q_mn]->clone(); // (nlmo_mn, npno_mn)

            // This performs a "T1-dressing" on Qvv, 
            // B^{Q_{mn}}_{f_{mn}c_{mn}} -= \widetilde{T}_{k_{mn}}^{f_{mn}} B^{Q_{mn}}_{k_{mn}c_{mn}}
            q_vv->subtract(linalg::doublet(T_n_ij_[mn], q_ov, true, false));

            auto Gvv_temp = q_ov->clone();
            Gvv_temp->scale(2.0 * rho_vv_[mn]->vector_dot(q_vv));
            Gvv_temp->subtract(linalg::triplet(q_ov, rho_vv_[mn], q_vv, false, true, false)); // (i, c) (f, c) (f, a) -> (i, a)

            for (int i_mn = 0; i_mn < lmopair_to_lmos_[mn].size(); ++i_mn) {
                int i = lmopair_to_lmos_[mn][i_mn];
                int ii = i_j_to_ij_[i][i];

                auto Gvv_slice = submatrix_rows(*Gvv_temp, std::vector<int>(1, i_mn)); // (1, a)
                L_ia_buffer[thread][i]->add(linalg::doublet(S_PNO(ii, mn), Gvv_slice, false, true));
            } // end i_mn
        } // end q_mn

        // l_{i}^{a_{ii}} \mathrel{-}= \rho^{OO}_{nm} \hat{F}_{mn}^{ia_{ii}} (Toth Eq. 56b)
        for (int i_mn = 0; i_mn < lmopair_to_lmos_[mn].size(); ++i_mn) {
            int i = lmopair_to_lmos_[mn][i_mn];

            /*
            auto Gvv_temp = linalg::doublet(F_fcia_hat_[mn][i_mn], Gvv_slice, false, false); // (a_ii, f_mn * c_mn) (f_mn * c_mn, 1) -> (a_ii, 1)
            L_ia_buffer[thread][i]->add(Gvv_temp);
            */

            auto F_mnia_slice = F_knia_hat_[mn][i_mn]->clone(); // (a_ii, 1)
            F_mnia_slice->scale(rho_oo_->get(n, m)); // Needs to be (n, m) not mn
            L_ia_buffer[thread][i]->subtract(F_mnia_slice);
        }
    } // end mn

    /* \textcolor{blue}{\begin{equation}
        l^{a_{ii}}_{i} \mathrel{+}= [S^{a_{ii}}_{a_{km}}S^{a_{km}}_{a_{mn}}\overline{\lambda}^{a_{mn}f_{mn}}_{mn}S^{f_{kn}}_{f_{mn}}]T^{f_{kn}c_{kn}}_{kn}S^{c_{km}}_{c_{kn}}\overline{J}^{ic_{km}}_{km} (Toth Eq. 57)
        \end{equation}}
    */
#pragma omp parallel for schedule(dynamic, 1)
    for (int mn = 0; mn < n_lmo_pairs; ++mn) {
        auto &[m, n] = ij_to_i_j_[mn];

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        for (int k_mn = 0; k_mn < lmopair_to_lmos_[mn].size(); ++k_mn) {
            int k = lmopair_to_lmos_[mn][k_mn];
            int km = i_j_to_ij_[k][m], kn = i_j_to_ij_[k][n];
            auto gus = linalg::doublet(lambda_iajb_bar_[mn], S_PNO(mn, kn));
            auto charlie = linalg::triplet(gus, T_iajb_[kn], S_PNO(kn, km));

            // Done! (From Toth Eq. 31)
            auto airbuds = linalg::doublet(charlie, J_kmic_bar_[km], false, true); // (a, c) * (i, c)
            
            for (int i_km = 0; i_km < lmopair_to_lmos_[km].size(); ++i_km) {
                int i = lmopair_to_lmos_[km][i_km];
                int ii = i_j_to_ij_[i][i];
                auto ryan = submatrix_cols(*airbuds, std::vector<int>(1, i_km));
                L_ia_buffer[thread][i]->add(linalg::doublet(S_PNO(ii, mn), ryan)); // Acts like a (7)5-year old
            } // end for

        } // end for
    } // end for

    // l_{i}^{a_{ii}} -= (S^{e_{ki}}_{e_{in}}\overline{\lambda}^{e_{in}f_{in}}_{in}S^{f_{kn}}_{f_{in}})T^{f_{kn}c_{kn}}_{kn}S^{c_{ki}}_{c_{kn}}\overline{J}^{e_{ki}c_{ki}}_{ka_{ii}} (Toth Eq. 38)
#pragma omp parallel for schedule(dynamic, 1) 
    for (int in = 0; in < n_lmo_pairs; ++in) {
        auto &[i, n] = ij_to_i_j_[in];

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        for (int k_in = 0; k_in < lmopair_to_lmos_[in].size(); ++k_in) {
            int k = lmopair_to_lmos_[in][k_in];
            int kn = i_j_to_ij_[k][n], in = i_j_to_ij_[i][n], ki = i_j_to_ij_[k][i];
            auto garfield = linalg::triplet(S_PNO(ki, in), lambda_iajb_bar_[in], S_PNO(in, kn));
            auto lasagne = linalg::triplet(garfield, T_iajb_[kn], S_PNO(kn, ki));

            lasagne->reshape(n_pno_[ki] * n_pno_[ki], 1);
            auto steak = linalg::doublet(J_kaec_bar_[ki], lasagne); // (a, e * c) (e * c, 1) => (a, 1) (steak sauce)
            L_ia_buffer[thread][i]->subtract(steak);
        } // end k_in
    } // end in
        
    /*
    \textcolor{blue}{\begin{equation}
        l^{a_{ii}}_i \mathrel{+}= \frac{1}{2} \widetilde{\lambda}_{in}^{e_{in}f_{in}} [S^{e_{in}}_{e_{ik}}(S^{c_{nk}}_{c_{ik}}u^{f_{nk}c_{nk}}_{nk}S^{f_{in}}_{f_{nk}})]\overline{M}^{c_{ki}e_{ki}}_{ka_{ii}} (Toth Eq. 59)
    \end{equation}}
    */
#pragma omp parallel for schedule(dynamic, 1)
    for (int in = 0; in < n_lmo_pairs; ++in) {
        auto &[i, n] = ij_to_i_j_[in];

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        for (int k_in = 0; k_in < lmopair_to_lmos_[in].size(); ++k_in) {
            int k = lmopair_to_lmos_[in][k_in];
            int nk = i_j_to_ij_[n][k], in = i_j_to_ij_[i][n], ki = i_j_to_ij_[k][i];
            auto costco = linalg::triplet(S_PNO(ki, in), lambda_iajb_[in], S_PNO(in, nk)); // (e_{ki}, e_{in}) (e_{in}, f_{in}) (f_{in}, f_{nk}) -> (e_{ki}, f_{nk})
            auto pizza = linalg::triplet(costco, Tt_iajb_[nk], S_PNO(nk, ki));

            pizza->transpose_this(); // (e, c) -> (c, e)
            pizza->reshape(n_pno_[ki] * n_pno_[ki], 1);
            auto breadsticks = linalg::doublet(M_kace_bar_[ki], pizza); // (a, c * e) (c * e, 1) => (a, 1) (steak sauce)
            breadsticks->scale(0.5);
            L_ia_buffer[thread][i]->add(breadsticks);
        } // end k_in
    } // end in

    /*
        \textcolor{blue}{\begin{equation}
            l^{a_{ii}}_i \mathrel{-}= \frac{1}{2}[S^{a_{ii}}_{a_{mk}}(S^{a_{mk}}_{a_{mn}}\widetilde{\lambda}^{a_{mn}f_{mn}}_{mn}S^{f_{mn}}_{f_{nk}})]u^{f_{nk}c_{nk}}_{nk}S^{c_{nk}}_{c_{mk}}\overline{M}^{ic_{mk}}_{mk} (Toth Eq. 40)
        \end{equation}}
    */
#pragma omp parallel for schedule(dynamic, 1)
    for (int mn = 0; mn < n_lmo_pairs; ++mn) {
        auto &[m, n] = ij_to_i_j_[mn];

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        for (int k_mn = 0; k_mn < lmopair_to_lmos_[mn].size(); ++k_mn) {
            int k = lmopair_to_lmos_[mn][k_mn];
            int nk = i_j_to_ij_[n][k], mk = i_j_to_ij_[m][k];
            auto johnpork = linalg::doublet(lambda_iajb_[mn], S_PNO(mn, nk)); // (a, f)
            auto squiddy = linalg::triplet(johnpork, Tt_iajb_[nk], S_PNO(nk, mk)); // (a, f) (f, c) -> (a, c)

            for (int i_mk = 0; i_mk < lmopair_to_lmos_[mk].size(); ++i_mk) {
                int i = lmopair_to_lmos_[mk][i_mk];
                int ii = i_j_to_ij_[i][i];
                auto M_slice = submatrix_rows(*M_mkic_bar_[mk], std::vector<int>(1, i_mk));
                auto calimari = linalg::triplet(S_PNO(ii, mn), squiddy, M_slice, false, false, true);
                calimari->scale(-0.5);
                L_ia_buffer[thread][i]->add(calimari);
            } // end i_mk
        } // end k_mn
    } // end mn

    // Add L_ia_buffers to L_ia
    for (int i = 0; i < naocc; ++i) {
        for (int thread = 0; thread < nthreads; ++thread) {
            L_ia[i]->add(L_ia_buffer[thread][i]);
        } // end thread
    } // end int i

    timer_off("DLPNO-CCSD Lambda : Compute L1");
}

std::vector<SharedMatrix> DLPNOCCSD_Lambda::compute_alpha_ijkl() {
    timer_on("DLPNO-CCSD Lambda : alpha ijkl");

    int n_lmo_pairs = ij_to_i_j_.size();
    int naocc = nalpha_ - nfrzc();

    std::vector<SharedMatrix> alpha_ijkl(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];

        if (i > j) continue;

        int nlmo_ij = lmopair_to_lmos_[ij].size();
        alpha_ijkl[ij] = std::make_shared<Matrix>("alpha_ijkl", nlmo_ij, nlmo_ij);

        // alpha_{ij}^{kl} = \lambda_{ij}^{e_{ij} f_{ij}} (S_{e_{kl}}^{e_{ij}} T_{kl}^{e_{kl}f_{kl}} S_{f_{kl}}^{f_{ij}}) (Toth Eq. 64)
        for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
            int k = lmopair_to_lmos_[ij][k_ij];
            for (int l_ij = 0; l_ij < nlmo_ij; ++l_ij) {
                int l = lmopair_to_lmos_[ij][l_ij];
                int kl = i_j_to_ij_[k][l];
                if (kl == -1) continue;
                auto gremlin = linalg::triplet(S_PNO(ij, kl), T_iajb_[kl], S_PNO(kl, ij));
                double val = gremlin->vector_dot(lambda_iajb_[ij]);
                alpha_ijkl[ij]->set(k_ij, l_ij, val);
            } // end l_ij
        } // end k_ij
    }

    timer_off("DLPNO-CCSD Lambda : alpha ijkl");

    return alpha_ijkl;
}

void DLPNOCCSD_Lambda::compute_L_iajb(std::vector<SharedMatrix>& L_iajb, std::vector<SharedMatrix>& Ln_iajb) {

    timer_on("DLPNO-CCSD Lambda : Compute L2");

    int n_lmo_pairs = ij_to_i_j_.size();
    int naocc = nalpha_ - nfrzc();

    // Thread and OMP Parallel info
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif

    // Zero out residuals
#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        L_iajb[ij]->zero();
        Ln_iajb[ij]->zero();
    } // end for

    auto alpha_ijkl = compute_alpha_ijkl();

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];
        int ji = ij_to_ji_[ij], ii = i_j_to_ij_[i][i];

        const int npao_ij = lmopair_to_paos_[ij].size();
        const int naux_ij = lmopair_to_ribfs_[ij].size();
        const int nlmo_ij = lmopair_to_lmos_[ij].size();
        const int npno_ij = n_pno_[ij];

        Ln_iajb[ij] = L_iajb_[ij]->clone();

        // These are the slow delinquent terms we apply to make our code faster
        if (i <= j) {
            // Necessary three-center integrals
            auto qma_ij = QIA_PNO(ij); // naux_ij * (nlmo_ij, npno_ij)
            auto qab_ij = QAB_PNO(ij); // naux_ij * (npno_ij, npno_ij)

            // Toth Eq. 70
            for (int q_ij = 0; q_ij < naux_ij; ++q_ij) {
                // This performs the T1-dressing of Qab on the fly, as this intermeidate is only used once
                // \widetilde{B}^{Q}_{ab} = B^{Q}_{ab} - t_{k}^{a} B^{Q}_{kb} (Jiang Eq. 93)
                auto Qab_t1 = qab_ij[q_ij]->clone(); // (a, b)
                Qab_t1->subtract(linalg::doublet(T_n_ij_[ij], qma_ij[q_ij], true, false)); // (k, a) (k, b) -> (a, b)

                auto L_temp = std::make_shared<Matrix>(n_pno_[ij], n_pno_[ij]);
                L_temp->zero();
                // l^{a_{ij}b_{ij}}_{ij} += 0.5 * \widetilde{\lambda}^{e_{ij}f_{ij}}_{ij}[\widetilde{B}^{Q_{ij}}_{e_{ij}a_{ij}}\widetilde{B}^{Q_{ij}}_{f_{ij}b_{ij}} (Toth Eq. 50a)
                L_temp->add(linalg::triplet(Qab_t1, lambda_iajb_[ij], Qab_t1, true, false, false)); // (e, a) (e, f) (f, b)
                // l^{a_{ij}b_{ij}}_{ij} += 0.5 * B^{Q_{ij}}_{k_{ij}a_{ij}} B^{Q_{ij}}_{l_{ij}b_{ij}} \alpha_{ij}^{k_{ij}l_{ij}} (Toth Eq. 50b)
                L_temp->add(linalg::triplet(qma_ij[q_ij], alpha_ijkl[ij], qma_ij[q_ij], true, false, false)); // (k, a) (k, l) (l, b)

                L_iajb[ij]->add(L_temp);
                if (i != j) L_iajb[ji]->add(L_temp->transpose());
            } // end q_ij

            // l_{ij}^{a_{ij}b_{ij}} \mathrel{+}= \frac{1}{2} (S_{a_{mn}}^{a_{ij}} \widetilde{\lambda}_{mn}^{a_{mn}b_{mn}}S_{b_{mn}}^{b_{ij}})\beta_{mn}^{ij} (Toth Eq. 51)
            for (int m_ij = 0; m_ij < nlmo_ij; ++m_ij) {
                int m = lmopair_to_lmos_[ij][m_ij];
                for (int n_ij = 0; n_ij < nlmo_ij; ++n_ij) {
                    int n = lmopair_to_lmos_[ij][n_ij];
                    int mn = i_j_to_ij_[m][n];
                    if (mn == -1) continue;
                    int i_mn = lmopair_to_lmos_dense_[mn][i], j_mn = lmopair_to_lmos_dense_[mn][j];

                    auto ethan = linalg::triplet(S_PNO(ij, mn), lambda_iajb_[mn], S_PNO(mn, ij));
                    ethan->scale(beta_[mn]->get(i_mn, j_mn));
                    L_iajb[ij]->add(ethan);
                    if (i != j) L_iajb[ji]->add(ethan->transpose());
                } // end n_ij
            } // end m_ij

            // l^{a_{ij}b_{ij}}_{ij} += \widetilde{\lambda}^{a_{ij}f_{ij}}_{ij}\widetilde{\widetilde{F}}_{f_{ij}b_{ij}} - (2 - P_{ab}) \rho^{\mathrm{VV}}_{a_{mn}c_{mn}}
            // S^{a_{mn}}_{a_{ij}} K^{c_{ij}b_{ij}}_{ij}S^{c_{mn}}_{c_{ij}} (Toth Eq. 54)
            auto E_temp = linalg::doublet(lambda_iajb_[ij], F_vv_double_tilde_[ij], false, false); // (a, f) (f, b) -> (a, b)
            E_temp->add(linalg::doublet(F_vv_double_tilde_[ij], lambda_iajb_[ij], true, false)); // (f, a) (f, b) -> (a, b)

            auto big_poob = std::make_shared<Matrix>(n_pno_[ij], n_pno_[ij]);
            big_poob->zero();
            
            for (int m_ij = 0; m_ij < nlmo_ij; ++m_ij) {
                int m = lmopair_to_lmos_[ij][m_ij];
                for (int n_ij = 0; n_ij < nlmo_ij; ++n_ij) {
                    int n = lmopair_to_lmos_[ij][n_ij];
                    int mn = i_j_to_ij_[m][n];
                    if (mn == -1) continue;

                    big_poob->add(linalg::triplet(S_PNO(ij, mn), rho_vv_[mn], S_PNO(mn, ij)));
                } // end n_ij
            } // end m_ij
            E_temp->subtract(linalg::doublet(big_poob, L_iajb_[ij])); // (a, c) (c, b) -> (a, b)
            E_temp->subtract(linalg::doublet(L_iajb_[ij], big_poob, false, true)); // (a, c) (b, c) -> (a, b)

            L_iajb[ij]->add(E_temp);
            if (i != j) L_iajb[ji]->add(E_temp->transpose());
        } // end i <= j

        // l^{a_{ij}b_{ij}}_{ij} += \lambda^{e_{jj}}_j\overline{L}^{a_{ij}b_{ij}}_{ie_{jj}} (Toth Eq. 43a)
        auto L_temp = linalg::doublet(lambda_ia_[j], L_ieab_bar_[ij], true, false); // (e, 1) (e, a * b) -> (1, a * b)
        L_temp->reshape(n_pno_[ij], n_pno_[ij]);
        Ln_iajb[ij]->add(L_temp);

        // l^{a_{ij}b_{ij}}_{ij} -= (2 - P_{ab})[S^{a_{ij}}_{a_{mm}} \lambda^{a_{mm}}_{m} \overline{K}^{mb_{ij}}_{ij}] (Toth Eq. 43b)
        for (int m_ij = 0; m_ij < nlmo_ij; ++m_ij) {
            int m = lmopair_to_lmos_[ij][m_ij];
            int mm = i_j_to_ij_[m][m];

            auto lambda_m_ij = linalg::doublet(S_PNO(ij, mm), lambda_ia_[m]);
            for (int a_ij = 0; a_ij < npno_ij; ++a_ij) {
                for (int b_ij = 0; b_ij < npno_ij; ++b_ij) {
                    double val = Ln_iajb[ij]->get(a_ij, b_ij) - 2.0 * lambda_m_ij->get(a_ij, 0) * K_ijmb_bar_[ij]->get(m_ij, b_ij)
                                    + lambda_m_ij->get(b_ij, 0) * K_ijmb_bar_[ij]->get(m_ij, a_ij);
                    Ln_iajb[ij]->set(a_ij, b_ij, val);
                } // end b_ij
            } // end a_ij
        } // end m_ij

        // l^{a_{ij}b_{ij}}_{ij} += (2 - P_{ab}) \lambda^{a_{ii}}_{i} S^{a_{ij}}_{a_{ii}} \widetilde{F}_{jb_{ij}} (Toth Eq. 43c)
        auto lambda_i = linalg::doublet(S_PNO(ij, ii), lambda_ia_[i]);
        for (int a_ij = 0; a_ij < npno_ij; ++a_ij) {
            for (int b_ij = 0; b_ij < npno_ij; ++b_ij) {
                double val = Ln_iajb[ij]->get(a_ij, b_ij) + 2.0 * lambda_i->get(a_ij, 0) * Fkc_[ji]->get(b_ij, 0)
                                - lambda_i->get(b_ij, 0) * Fkc_[ji]->get(a_ij, 0);
                Ln_iajb[ij]->set(a_ij, b_ij, val);
            } // end b_ij
        } // end a_ij

        for (int n_ij = 0; n_ij < nlmo_ij; ++n_ij) {
            int n = lmopair_to_lmos_[ij][n_ij];
            int in = i_j_to_ij_[i][n], jn = i_j_to_ij_[j][n], nj = i_j_to_ij_[n][j];
            int i_nj = lmopair_to_lmos_dense_[nj][i];

            // l_{ij}^{a_{ij}b_{ij}} -= \overline{\lambda}_{jn}^{a_{jn}f_{jn}}[S_{a_{nj}}^{a_{ij}}\widetilde{\gamma}_{in}^{f_{nj}b_{ij}} (Toth Eq. 52a)
            auto ron = linalg::doublet(S_PNO(ij, jn), lambda_iajb_bar_[jn]);
            auto uncle_andy = linalg::triplet(S_PNO(nj, in), gamma_[in], S_PNO(in, ij));
            uncle_andy->add(J_ikac_non_proj_[nj][i_nj]); // Uncle Andy stays delinquent
            auto weasley = linalg::doublet(ron, uncle_andy);
            Ln_iajb[ij]->subtract(weasley);

            // l_{ij}^{a_{ij}b_{ij}} -= 0.5 * S_{a_{kj}}^{a_{jn}}(S_{a_{ij}}^{a_{kj}}S_{b_{ki}}^{b_{ij}} K_{ki}^{b_{ki}c_{ki}}S_{c_{ki}}^{c_{kj}})S_{c_{kj}}^{c_{kn}}T_{kn}^{f_{kn}c_{kn}}S_{f_{kn}}^{f_{jn}}] (Toth Eq. 52b)
            auto gamma_dos = gamma_double_tilde_[ij][n_ij]->clone();
            gamma_dos->scale(0.5);
            Ln_iajb[ij]->add(linalg::triplet(S_PNO(ij, jn), lambda_iajb_bar_[jn], gamma_dos));
        } // end n_ij
        
        for (int n_ij = 0; n_ij < nlmo_ij; ++n_ij) {
            int n = lmopair_to_lmos_[ij][n_ij];
            int nj = i_j_to_ij_[n][j], in = i_j_to_ij_[i][n], ni = i_j_to_ij_[n][i];
            int j_ni = lmopair_to_lmos_dense_[ni][j];

            // l^{a_{ij}b_{ij}}_{ij} += (2 - P_{ab}) [\frac{1}{2}\widetilde{\lambda}^{f_{ni}a_{ni}}_{ni}[\widetilde{\delta}^{f_{ni}b_{ij}}_{nj}S^{a_{ni}}_{a_{ij}} (Toth Eq. 53a) 
            auto forg = linalg::doublet(S_PNO(ij, in), lambda_iajb_[in]);
            /// [ij][k_ij] => (i a_{ij} | k c_{kj}) 
            /// [ij][k_ij] => (i k | a_{ij} c_{kj})

            /// [ni][j_ni] => (n a_{ni} | j c_{ji})

            /// (n f_{ni} | j b_{ij})
            auto glizzy = K_iakc_non_proj_[ni][j_ni]->clone(); // I forgor a clone :)
            glizzy->scale(2.0);
            /// (n j | f_{ni} b_{ij})
            glizzy->subtract(J_ikac_non_proj_[ni][j_ni]);
            glizzy->add(linalg::triplet(S_PNO(ni, nj), delta_[nj], S_PNO(nj, ij))); // (f_{ni}, b_{ij})
            auto parts = linalg::doublet(forg, glizzy);
            
            Ln_iajb[ij]->add(parts);
            parts->scale(-0.5);
            Ln_iajb[ij]->add(parts->transpose());

            auto delta_dos = delta_double_tilde_[ij][n_ij]->clone();
            delta_dos->scale(0.5);
            
            auto delta_dos_temp = linalg::triplet(S_PNO(ij, in), lambda_iajb_[in], delta_dos);
            Ln_iajb[ij]->add(delta_dos_temp);
            delta_dos_temp->scale(-0.5);
            Ln_iajb[ij]->add(delta_dos_temp->transpose());
        }

        // l^{a_{ij}b_{ij}}_{ij} -= (S^{a_{ij}}_{a_{in}}\widetilde{\lambda}^{a_{in}b_{in}}_{in}S^{b_{ij}}_{b_{in}})\widetilde{\widetilde{F}}_{jn} + 
        // \rho^{\mathrm{OO}}_{jk}(S^{a_{ik}}_{a_{ij}}L^{a_{ik}b_{ik}}_{ik}S^{b_{ik}}_{b_{ij}}) (Toth Eq. 55)
        for (int n_ij = 0; n_ij < nlmo_ij; ++n_ij) {
            int n = lmopair_to_lmos_[ij][n_ij];
            int in = i_j_to_ij_[i][n];
            auto lambda_temp = linalg::triplet(S_PNO(ij, in), lambda_iajb_[in], S_PNO(in, ij));
            lambda_temp->scale(F_im_double_tilde_->get(j, n));
            Ln_iajb[ij]->subtract(lambda_temp);

            auto ohio = linalg::triplet(S_PNO(ij, in), L_iajb_[in], S_PNO(in, ij));
            ohio->scale(rho_oo_->get(j, n));
            Ln_iajb[ij]->subtract(ohio);
        } // end n_ij
    } // end ij

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        int i, j;
        std::tie(i, j) = ij_to_i_j_[ij];
        int ji = ij_to_ji_[ij];
        
        L_iajb[ij]->add(Ln_iajb[ij]);
        L_iajb[ij]->add(Ln_iajb[ji]->transpose());
    }

    timer_off("DLPNO-CCSD Lambda : Compute L2");
}

void DLPNOCCSD_Lambda::lambda_ccsd_iterations() {

    int n_lmo_pairs = ij_to_i_j_.size();
    int naocc = nalpha_ - nfrzc();

    // Thread and OMP Parallel info
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif

    outfile->Printf("\n  ==> Lambda DLPNO-CCSD <==\n\n");
    outfile->Printf("    E_CONVERGENCE = %.2e\n", options_.get_double("E_CONVERGENCE"));
    outfile->Printf("    R_CONVERGENCE = %.2e\n\n", options_.get_double("R_CONVERGENCE"));
    outfile->Printf("                            Corr. Energy    Delta E     Max L1     Max L2     Time (s)\n");

    // => Initialize Residuals and Amplitudes <= //

    std::vector<SharedMatrix> L_ia(naocc);
    std::vector<SharedMatrix> L_ia_prev(naocc);
    std::vector<SharedMatrix> Ln_iajb(n_lmo_pairs);
    std::vector<SharedMatrix> L_iajb(n_lmo_pairs);
    std::vector<SharedMatrix> L_iajb_prev(n_lmo_pairs);

    // => Initialize Singles and Doubles Residuals and Amplitudes <= //
    lambda_ia_.resize(naocc);
    lambda_iajb_.resize(n_lmo_pairs);
    lambda_iajb_bar_.resize(n_lmo_pairs);

#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < naocc; ++i) {
        int ii = i_j_to_ij_[i][i];
        lambda_ia_[i] = std::make_shared<Matrix>(n_pno_[ii], 1);
        L_ia[i] = std::make_shared<Matrix>(n_pno_[ii], 1);
        L_ia_prev[i] = std::make_shared<Matrix>(n_pno_[ii], 1);
    }

    // => Initialize Doubles Residuals and Amplitudes <= //
    
#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        lambda_iajb_[ij] = std::make_shared<Matrix>(n_pno_[ij], n_pno_[ij]);
        lambda_iajb_bar_[ij] = std::make_shared<Matrix>(n_pno_[ij], n_pno_[ij]);
        L_iajb[ij] = std::make_shared<Matrix>(n_pno_[ij], n_pno_[ij]);
        Ln_iajb[ij] = std::make_shared<Matrix>(n_pno_[ij], n_pno_[ij]);
        L_iajb_prev[ij] = std::make_shared<Matrix>(n_pno_[ij], n_pno_[ij]);
    }

    // => Thread buffers <= //

    std::vector<std::vector<SharedMatrix>> L_ia_buffer(nthreads);
    for (int thread = 0; thread < nthreads; ++thread) {
        L_ia_buffer[thread].resize(naocc);
        for (int i = 0; i < naocc; ++i) {
            int ii = i_j_to_ij_[i][i];
            L_ia_buffer[thread][i] = std::make_shared<Matrix>(n_pno_[ii], 1);
        } // end i
    } // end thread

    int iteration = 1, max_iteration = options_.get_int("DLPNO_MAXITER");
    double damping = options_.get_double("DLPNO_LAMBDA_DAMPING");
    double e_curr = 0.0, e_prev = 0.0, l1_curr = 0.0, l2_curr = 0.0;
    bool e_converged = false, l_converged = false;

    DIISManager diis(options_.get_int("DIIS_MAX_VECS"), "LCCSD DIIS", DIISManager::RemovalPolicy::LargestError, DIISManager::StoragePolicy::InCore);

    while (!(e_converged && l_converged)) {
        // RMS of residual per single LMO, for assesing convergence
        std::vector<double> L_ia_rms(naocc, 0.0);
        // RMS of residual per LMO pair, for assessing convergence
        std::vector<double> L_iajb_rms(n_lmo_pairs, 0.0);

        std::time_t time_start = std::time(nullptr);

        // Form intermediates
        form_goo();

        // Step 4: Compute R2 residual
        compute_L_iajb(L_iajb, Ln_iajb);

        // Get rms of R_iajb
#pragma omp parallel for schedule(dynamic, 1)
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            L_iajb_rms[ij] = L_iajb[ij]->rms();
        }

        // Update Doubles Amplitude (Jiang Eq. 104)
#pragma omp parallel for schedule(dynamic, 1)
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            auto &[i, j] = ij_to_i_j_[ij];

            // Dynamic Damping
            /*
            double m = (iteration > 10) ? -L_iajb[ij]->vector_dot(K_iajb_[ij]) / L_iajb_prev[ij]->vector_dot(K_iajb_[ij]) : -0.3;
            double alpha = (m > 0.0) ? 1.0 : 1.0 / (1.0 - m);
            */

            for (int a_ij = 0; a_ij < n_pno_[ij]; ++a_ij) {
                for (int b_ij = 0; b_ij < n_pno_[ij]; ++b_ij) {
                    double val = lambda_iajb_[ij]->get(a_ij, b_ij) - (1.0 - damping) * L_iajb[ij]->get(a_ij, b_ij) /
                                    (e_pno_[ij]->get(a_ij) + e_pno_[ij]->get(b_ij) - F_lmo_->get(i,i) - F_lmo_->get(j,j));
                    lambda_iajb_[ij]->set(a_ij, b_ij, val);
                }
            }
            L_iajb_prev[ij] = L_iajb[ij]->clone();
            // L_iajb[ij]->scale(alpha);
        }

        if (!brueckner_orbs_) {
            // Form Goo a second time (using updated lambda)
            form_goo();

            // Step 1: Compute R1 residual
            compute_L_ia(L_ia, L_ia_buffer);

            // Get rms of L_ia
#pragma omp parallel for schedule(dynamic, 1)
            for (int i = 0; i < naocc; ++i) {
                L_ia_rms[i] = L_ia[i]->rms();
            }

            // Update Singles Amplitude (Jiang Eq. 103)
#pragma omp parallel for schedule(dynamic, 1)
            for (int i = 0; i < naocc; ++i) {
                int ii = i_j_to_ij_[i][i];

                // Dynamic Damping
                /*
                double m = (iteration > 10) ? -L_ia[i]->vector_dot(linalg::doublet(K_iajb_[ii], L_ia[i])) 
                                            / L_ia_prev[i]->vector_dot(linalg::doublet(K_iajb_[ii], L_ia_prev[i])) : -0.3;
                double alpha = (m > 0.0) ? 1.0 : 1.0 / (1.0 - m);
                */

                for (int a_ii = 0; a_ii < n_pno_[ii]; ++a_ii) {
                    double val = lambda_ia_[i]->get(a_ii, 0) - (1.0 - damping) * L_ia[i]->get(a_ii, 0) / (e_pno_[ii]->get(a_ii) - F_lmo_->get(i,i));
                    lambda_ia_[i]->set(a_ii, 0, val);
                }
                L_ia_prev[i] = L_ia[i]->clone();
                // L_ia[i]->scale(alpha);
            }
        }

        // DIIS Extrapolation
        std::vector<SharedMatrix> lambda_vecs;
        lambda_vecs.reserve(lambda_ia_.size() + lambda_iajb_.size());
        lambda_vecs.insert(lambda_vecs.end(), lambda_ia_.begin(), lambda_ia_.end());
        lambda_vecs.insert(lambda_vecs.end(), lambda_iajb_.begin(), lambda_iajb_.end());

        std::vector<SharedMatrix> L_vecs;
        L_vecs.reserve(L_ia.size() + L_iajb.size());
        L_vecs.insert(L_vecs.end(), L_ia.begin(), L_ia.end());
        L_vecs.insert(L_vecs.end(), L_iajb.begin(), L_iajb.end());

        auto lambda_vecs_flat = flatten_mats(lambda_vecs);
        auto L_vecs_flat = flatten_mats(L_vecs);

        if (iteration == 1) {
            diis.set_error_vector_size(lambda_vecs_flat);
            diis.set_vector_size(L_vecs_flat);
        }
        
        diis.add_entry(lambda_vecs_flat.get(), lambda_vecs_flat.get());
        diis.extrapolate(L_vecs_flat.get());

        copy_flat_mats(lambda_vecs_flat, lambda_vecs);
        
        // Compute lambda CCSD pseudoenergy (and remake lambda_ijab_bar)
        e_curr = 0.0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : e_curr)
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            auto &[i, j] = ij_to_i_j_[ij];

            // TODO: HOMEWORK FOR ZHUDNEY
            lambda_iajb_bar_[ij] = lambda_iajb_[ij]->clone();
            lambda_iajb_bar_[ij]->scale(0.5);
            lambda_iajb_bar_[ij]->add(lambda_iajb_[ij]->transpose());

            e_curr += 0.5 * K_iajb_[ij]->vector_dot(lambda_iajb_[ij]);
        }
        
        double l_curr1 = *max_element(L_ia_rms.begin(), L_ia_rms.end());
        double l_curr2 = *max_element(L_iajb_rms.begin(), L_iajb_rms.end());

        l_converged = (fabs(l_curr1) < options_.get_double("R_CONVERGENCE"));
        l_converged &= (fabs(l_curr2) < options_.get_double("R_CONVERGENCE"));
        e_converged = (fabs(e_curr - e_prev) < options_.get_double("E_CONVERGENCE"));

        std::time_t time_stop = std::time(nullptr);

        outfile->Printf("  @Lambda-CCSD iter %3d: %16.12f %10.3e %10.3e %10.3e %8d\n", iteration, e_curr, e_curr - e_prev, l_curr1, l_curr2, (int)time_stop - (int)time_start);

        iteration++;

        if (iteration > max_iteration + 1) {
            throw PSIEXCEPTION("Maximum DLPNO iterations exceeded.");
        }

        // For next iteration
        e_prev = e_curr;
    } // end iter
}

void DLPNOCCSD_Lambda::compute_opdm() {

    int naocc = i_j_to_ij_.size();
    int n_lmo_pairs = ij_to_i_j_.size();

    form_goo();

    // Toth Eq. 65
    Doo_ = rho_oo_->clone();
    Doo_->scale(-1.0);
#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];
        int ii = i_j_to_ij_[i][i], jj = i_j_to_ij_[j][j];

        auto baby_shark = linalg::triplet(lambda_ia_[i], S_PNO(ii, jj), T_ia_[j], true, false, false);
        (*Doo_)(i, j) -= baby_shark->get(0, 0);
    } // end ij

    // Zero out buffers

    // Thread and OMP Parallel Info
    int nthreads = 1;
#ifdef _OPENMP
    nthreads = Process::environment.get_n_threads();
#endif

    std::vector<std::vector<SharedMatrix>> D_ov_buffer(nthreads);
    for (int thread = 0; thread < nthreads; ++thread) {
        D_ov_buffer[thread].resize(naocc);
        for (int i = 0; i < naocc; ++i) {
            int ii = i_j_to_ij_[i][i];
            D_ov_buffer[thread][i] = std::make_shared<Matrix>(n_pno_[ii], 1);
        }
    }

    // Evil mf (Toth Eq. 66)
    Dov_.resize(naocc);
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < naocc; ++i) {
        int ii = i_j_to_ij_[i][i];
        // 66a
        Dov_[i] = T_ia_[i]->clone();
        Dov_[i]->scale(2.0);
    }

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];
        int ii = i_j_to_ij_[i][i], jj = i_j_to_ij_[j][j];

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif
            
        // 66b
        auto sus = linalg::triplet(S_PNO(ii, ij), Tt_iajb_[ij], S_PNO(ij, jj));
        D_ov_buffer[thread][i]->add(linalg::doublet(sus, lambda_ia_[j]));

        // 66c
        auto T_i_to_j = linalg::doublet(S_PNO(ii, jj), T_ia_[j]);
        T_i_to_j->scale(rho_oo_->get(j, i));
        D_ov_buffer[thread][i]->subtract(T_i_to_j);
    }

#pragma omp parallel for schedule(dynamic, 1)
    for (int mn = 0; mn < n_lmo_pairs; ++mn) {
        auto &[m, n] = ij_to_i_j_[mn];

        int nlmo_mn = lmopair_to_lmos_[mn].size();

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        for (int i_mn = 0; i_mn < nlmo_mn; ++i_mn) {
            int i = lmopair_to_lmos_[mn][i_mn];
            int ii = i_j_to_ij_[i][i];

            // (e_{ii}, 1), (e_{ii}, e_{mn}), (e_{mn}, a_{mn})
            auto bean = linalg::triplet(T_ia_[i], S_PNO(ii, mn), rho_vv_[mn], true, false, false); 
            D_ov_buffer[thread][i]->subtract(linalg::doublet(S_PNO(ii, mn), bean, false, true));
        } // end i_mn
    } // end mn

    for (int i = 0; i < naocc; ++i) {
        for (int thread = 0; thread < nthreads; ++thread) {
            Dov_[i]->add(D_ov_buffer[thread][i]);
        } // end thread
    } // end int i

    Dvv_pair_.resize(n_lmo_pairs);
#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        Dvv_pair_[ij] = rho_vv_[ij]->clone();
    }

    Dvv_singles_.resize(naocc);
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < naocc; ++i) {
        int ii = i_j_to_ij_[i][i];

        Dvv_singles_[i] = std::make_shared<Matrix>(n_pno_[ii], n_pno_[ii]);

        for (int a = 0; a < n_pno_[ii]; ++a) {
            for (int b = 0; b < n_pno_[ii]; ++b) {
                Dvv_singles_[i]->set(a, b, T_ia_[i]->get(a, 0) * lambda_ia_[i]->get(b, 0));
            } // end b
        } // end a
    } // end i
    
} // end function

Vector3 DLPNOCCSD_Lambda::compute_dipole_moment() {

    int naocc = i_j_to_ij_.size();
    int n_lmo_pairs = ij_to_i_j_.size();

    // Nuclear contribution to the dipole moment
    Vector3 nuclear_contribution(0.0, 0.0, 0.0);

    for (int A = 0; A < molecule_->natom(); ++A) {
        Vector3 R_A = molecule_->xyz(A);
        nuclear_contribution[0] += molecule_->Z(A) * R_A[0];
        nuclear_contribution[1] += molecule_->Z(A) * R_A[1];
        nuclear_contribution[2] += molecule_->Z(A) * R_A[2];
    } // end A

    const auto ao_dipole = MintsHelper(basisset_, options_).ao_dipole();

    // Compute AO density matrix (for Hartree-Fock electronic contribution)
    auto C_occ = reference_wavefunction_->Ca_subset("AO", "OCC");
    auto D_ao = linalg::doublet(C_occ, C_occ, false, true);

    Vector3 hf_elec_contribution(0.0, 0.0, 0.0);
    for (int soup = 0; soup < 3; ++soup) {
        hf_elec_contribution[soup] += 2.0 * ao_dipole[soup]->vector_dot(D_ao);
    }

    // Correlated contribution
    Vector3 ccsd_contribution(0.0, 0.0, 0.0);

    for (int soup = 0; soup < 3; ++soup) {
        auto mu_oo = linalg::triplet(C_lmo_, ao_dipole[soup], C_lmo_, true, false, false);
        auto mu_ov = linalg::triplet(C_lmo_, ao_dipole[soup], C_pao_, true, false, false);
        auto mu_vv = linalg::triplet(C_pao_, ao_dipole[soup], C_pao_, true, false, false);

        double dipole_cont = 0.0;
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : dipole_cont)
        for (int ij = 0; ij < n_lmo_pairs; ++ij) {
            auto &[i, j] = ij_to_i_j_[ij];

            // Doo contributions (Toth Eq. 64a)
            dipole_cont += Doo_->get(i, j) * mu_oo->get(i, j);

            // Dvv (doubles) contributions (Toth Eq. 64c)
            auto mu_vv_ij = submatrix_rows_and_cols(*mu_vv, lmopair_to_paos_[ij], lmopair_to_paos_[ij]);
            mu_vv_ij = linalg::triplet(X_pno_[ij], mu_vv_ij, X_pno_[ij], true, false, false);
            dipole_cont += Dvv_pair_[ij]->vector_dot(mu_vv_ij);
        }
        
#pragma omp parallel for schedule(dynamic, 1) reduction(+ : dipole_cont)
        for (int i = 0; i < naocc; ++i) {
            int ii = i_j_to_ij_[i][i];

            // Dov contributions (Toth Eq. 64b)
            auto mu_ov_slice = submatrix_rows_and_cols(*mu_ov, std::vector<int>(1, i), lmopair_to_paos_[ii]);
            auto mu_ov_ii = linalg::doublet(mu_ov_slice, X_pno_[ii], false, false); // <i|x|a_{ii}>
            auto Dov_total = Dov_[i]->clone(); // D_{i}^{a_{ii}}
            Dov_total->add(lambda_ia_[i]);

            dipole_cont += Dov_total->vector_dot(mu_ov_ii->transpose());

            // Dvv (singles) contributions (Toth Eq. 64d)
            auto mu_vv_ii = submatrix_rows_and_cols(*mu_vv, lmopair_to_paos_[ii], lmopair_to_paos_[ii]);
            mu_vv_ii = linalg::triplet(X_pno_[ii], mu_vv_ii, X_pno_[ii], true, false, false);
            dipole_cont += Dvv_singles_[i]->vector_dot(mu_vv_ii);
        }

        ccsd_contribution[soup] += dipole_cont;
    } // end soup

    outfile->Printf("    Nuclear Contribution: \n");
    outfile->Printf("    X: %.6f, Y: %.6f, Z: %.6f\n\n", nuclear_contribution[0], nuclear_contribution[1], nuclear_contribution[2]);
    outfile->Printf("    SCF Electronic Contribution: \n");
    outfile->Printf("    X: %.6f, Y: %.6f, Z: %.6f\n\n", hf_elec_contribution[0], hf_elec_contribution[1], hf_elec_contribution[2]);
    outfile->Printf("    CCSD Correlated Contribution: \n");
    outfile->Printf("    X: %.6f, Y: %.6f, Z: %.6f\n\n", ccsd_contribution[0], ccsd_contribution[1], ccsd_contribution[2]);

    Vector3 total_dipole = nuclear_contribution;
    total_dipole += hf_elec_contribution;
    total_dipole += ccsd_contribution;

    outfile->Printf("    ==> Total Dipole Moment <== \n");
    outfile->Printf("    X: %.6f, Y: %.6f, Z: %.6f\n\n", total_dipole[0], total_dipole[1], total_dipole[2]);

    return total_dipole;
}

void DLPNOCCSD_Lambda::print_header() {}
void DLPNOCCSD_Lambda::print_results() {}

double DLPNOCCSD_Lambda::compute_energy() {
    // Run DLPNO-CCSD
    double e_dlpno_ccsd = DLPNOCCSD::compute_energy();

    estimate_memory();
    compute_lambda_intermediates();

    lambda_ccsd_iterations();

    compute_opdm();
    compute_dipole_moment();

    return e_dlpno_ccsd;
}

} // end dlpno
} // end psi
