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

void DLPNOCCSD_Lambda::compute_lambda_intermediates() {

    // Number of active occupied orbitals
    int naocc = nalpha_ - nfrzc();
    // Number of surviving pairs after DLPNO screening
    int n_lmo_pairs = ij_to_i_j_.size();

    // \rho^{OO}_{nk} = \sum_{m, e, f} \lambda_{mn}^{e_{mn} f_{mn}} [S(e_{mn}, e_{mk}) T_{mk}^{e_{mk} f_{mk}} S(f_{mk}, f_{mn})]
    rho_oo_ = std::make_shared<Matrix>("rho_oo", naocc, naocc);

#pragma omp parallel for
    for (int nk = 0; nk < n_lmo_pairs; ++nk) {
        auto &[n, k] = ij_to_i_j_[nk];

        for (int m_nk = 0; m_nk < lmopair_to_lmos_[nk].size(); ++m_nk) {
            int m = lmopair_to_lmos_[nk][m_nk];
            int mn = i_j_to_ij_[m][n], mk = i_j_to_ij_[m][k];

            // Transform amplitude from domain of mk to mn [S(e_{mn}, e_{mk}) T_{mk}^{e_{mk} f_{mk}} S(f_{mk}, f_{mn})]
            auto T_mk_to_mn = linalg::triplet(S_PNO(mn, mk), T_iajb_[mk], S_PNO(mk, mn)); // (e_{mn}, f_{mn})

            // \rho^{OO}_{nk} += \lambda_{mn}^{e_{mn} f_{mn}} T_mk_to_mn(e_{mn} f_{mn})
            rho_oo_ += lambda_iajb_[mn]->vector_dot(T_mk_to_mn);

        } // end m_nk
    } // end nk

    // \rho^{VV}_{f_{mn} c{mn}} = \sum_{e, m, n} \lambda_{mn}^{e_{mn} f_{mn}} T_{mn}^{e_{mn} c_{mn}}
    rho_vv_.resize(n_lmo_pairs);

#pragma omp parallel for
    for (int mn = 0; mn < n_lmo_pairs; ++mn) {
        auto &[m, n] = ij_to_i_j_[mn];

        rho_vv_[mn] = linalg::doublet(lambda_iajb_[mn], T_iajb_[mn], true, false);
    } // end

    /*
    // Form the \widetilde{M}_imae (Eq. 24)
    M_imae_tilde_.resize(n_lmo_pairs);

#pragma omp parallel for
    for (int im = 0; im < n_lmo_pairs; ++im) {
        auto &[i, m] = ij_to_i_j_[im];
        int ii = i_j_to_ij_[i][i], mm = i_j_to_ij_[m][m];

        // \widetilde{M}_{im}^{a_{ii} e_{mm}} = S^{a_{ii}}_{a_{mm}} M_{im}^{a_{mm} e_{mm}} (Eq. 24a)
        M_imae_tilde_[im] = linalg::doublet(S_PNO(ii, mm), M_imae_[im]);
        // Note M_{im}^{a_{mm} e_{mm}} = 2K_{im}^{a_{mm} e_{mm}} - J_{im}^{a_{mm} e_{mm}} 
        // TODO: Form this later

        for (int k_im = 0; k_im < lmopair_to_lmos_[im].size(); ++k_im) {
            int k = lmopair_to_lmos_[im][k_im];
            int mk = ij_to_i_j_[m][k], ik = ij_to_i_j_[i][k], kk = i_j_to_ij_[k][k];

            // - \widetilde{T}_{k}^{e_{mm}} N_{mk}^{ia_{mk}} S_{a_{mk}}^{a_{ii}} (Eq. 24b)
            auto N_mkia = J_ijmb_[mk]->clone();
            N_mkia->scale(2.0);
            N_mkia->subtract(K_mibj_[mk]);

            // \widetilde{T}_{k}^{e_{mm}} 
            T_k_to_m = linalg::doublet(S_PNO(mm, kk), T_ia_[k]);
            // This is a temporary variable
            N_mkia_ii = linalg::doublet(N_mkia, S_PNO(mk, ii));

            for (int e_mm = 0; e_mm < n_pno_[mm]; ++e_mm) {

            }

        } // end k_im

    } // end im
    */

    // Toth Eq. 31
    // \overline{J}_{km}^{ic} = (km | i c_{km}) + \widetilde{T}_{m}^{f_{ki}} (k f_{ki} | i c_{ki})
    //  S_{c_{ki}}^{c_{km}}
#pragma omp parallel for
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
                (*J_kmic_bar_[km])(i_km, c_km) += (*f_steak_university)(0, c_km);
            } // end c_km
        } // end i_km
    } // end km

    // Toth Eq. 32
    // \overline{J}_{ka_{ii}}^{e_{ki}c_{ki}} = S^{a_{ii}}_{a_{ki}} (ka_{ki}|e_{ki}c_{ki}) -
    // S_{a_{ii}}^{a_{kl}} (k a_{kl} | l c_{kl}) S_{c_{kl}}^{c_{ki}} \widetilde{T}_{l}^{e_{ki}}
#pragma omp parallel for
    for (int ki = 0; ki < n_lmo_pairs; ++ki) {
        auto &[k, i] = ij_to_i_j_[ki];
        int ii = i_j_to_ij_[i][i];

        auto J_kaec_bar_[ki] = linalg::doublet(S_PNO(ii, ki), K_ivvv_[ki]); // (a, e * c)

        for (int l_ki = 0; l_ki < lmopair_to_lmos_[ki].size(); ++l_ki) {
            int l = lmopair_to_lmos_[ki][l_ki];
            int kl = i_j_to_ij_[k][l];

            // S_{a_{ii}}^{a_{kl}} (k a_{kl} | l c_{kl}) S_{c_{kl}}^{c_{ki}}
            auto south_ohio = linalg::triplet(S_PNO(ii, kl), K_iajb_[kl], S_PNO(kl, ki));
            auto T_l_ki = linalg::doublet(S_PNO(ki, ll), T_ia_[l]); // \widetilde{T}_{l}^{e_{ki}} (Jiang Eq. 70?)

            for (int a_ki = 0; a_ki < n_pno_[ki]; ++a_ki) {
                for (int e_ki = 0; e_ki < n_pno_[ki]; ++e_ki) {
                    for (int c_ki = 0; c_ki < n_pno_[ki]; ++c_ki) {
                        (*J_kaec_bar_[ki])(a_ki, e_ki * n_pno_[ki] + c_ki) -= (*south_ohio)(a_ki, c_ki) * (*T_l_ki)(e_ki, 0);
                    } // end c_ki
                } // end e_ki
            } // end a_ki

        } // end l_ki
    } // end for ki
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

    // TODO: Add Crawford Line 471 from GitHub

    // \lambda_{i}^{e_{ii}} \widetilde{\widetilde{F}}_(e_{ii}, a_{ii})
#pragma omp parallel for schedule(dynamic, 1)
    for (int i = 0; i < naocc; ++i) {
        int ii = i_j_to_ij_[i][i];

        L_ia[i] = linalg::doublet(F_vv_double_tilde_[ii], lambda_ia_[i], true, false); // TODO: Make this variable
    } // end for

#pragma omp parallel for schedule(dynamic, 1)
    for (int im = 0; im < n_lmo_pairs; ++im) {
        auto &[i, m] = ij_to_i_j_[im];
        int ii = i_j_to_ij_[i][i], mm = i_j_to_ij_[m][m];

        int thread = 0;
#ifdef _OPENMP
        thread = omp_get_thread_num();
#endif

        // \widetilde{\delta}_{im}^{a_{ii} e_{mm}} \lambda_{m}^{e_{mm}}
        L_ia_buffer[thread][i]->add(linalg::doublet(delta_imae_tilde_[im], lambda_ia_[m]));

        // - \widetilde{widetilde{F}}_{im} \lambda_{m}^{a_{mm}} S_{a_{mm}}^{a_{ii}}
        auto john_big_back_buffer = linalg::doublet(S_PNO(ii, mm), lambda_ia_[m]);
        john_big_back_buffer->scale((*F_im_double_tilde_)(i, m));
        L_ia_buffer[thread][i]->subtract(john_big_back_buffer);
        
        /* l_{i}^{a_{ii}} += \widetilde{\widetilde{K}}_{ma_{ii}}^{e_{mi}f_{mi}} \widetilde{\lambda}_{mi}^{e_{mi}f_{mi}} (Toth Eq. 35a) */
        // TODO: Define K_maef_dt_
        auto lambda_mn_slice = lambda_iajb_[mn]->clone();
        lambda_mn_slice->reshape(n_pno_[mn] * n_pno_[mn], 1);
        L_ia_buffer[thread][i]->add(linalg::doublet(K_maef_dt_[mi], lambda_mn_slice)); // (a_{ii}, e_{mi} * f_{mi}) (e_{mi} * f_{mi}, 1)

    } // end im

    // rip John Pork :( bigback = ibs (irritable bigback syndrome)
    //    ^     ^
    // (  *     *  )
    // (  ( * * )  )
    // (           )
    // OINK!!! John Pork is calling...
#pragma omp parallel for schedule(dynamic, 1)
    for (int mn = 0; mn < n_lmo_pairs; ++mn) {

        int nlmo_mn = lmopair_to_lmos_[mn].size();

        /* l_{i}^{a_{ii}} -= \widetilde{\widetilde{K}}_{e_{mn}i}^{mn} \widetilde{\lambda}_{mn}^{e_{mn}a_{mn}}S_{a_{mn}}^{a_{ii}} (Toth Eq. 35b) */
        // TODO: K_eimn_dt has not been defined yet ()
        auto bruvver = linalg::triplet(K_eimn_dt_[mn], lambda_iajb_[mn], S_PNO(mn, ii), true, false, false);

        for (int i_mn = 0; i_mn < lmopair_to_lmos_[mn].size(); ++i_mn) {
            int i = lmopair_to_lmos_[mn][i_mn];
            auto antonios_slice = submatrix_rows(std::vector<int>(1, i_mn), *bruvver);
            L_ia_buffer[thread][i]->subtract(antonios_slice->transpose());
        }

        // l_{i}^{a_{ii}} \mathrel{+}= \rho^{\mathrm{VV}}_{f_{mn}c_{mn}}\hat{F}^{ia_{ii}}_{f_{mn}c_{mn}} - \rho^{\mathrm{OO}}_{nm} \hat{F}_{mn}^{ia_{ii}} (Toth Eq. 36)
        // TODO: Define F_iafc_hat_ and F_mnia_hat_
        auto Gvv_slice = rho_vv_[mn]->clone();
        Gvv_slice->reshape(n_pno_[mn] * n_pno_[mn], 1);
        for (int i_mn = 0; i_mn < lmopair_to_lmos_[mn].size(); ++i_mn) {
            int i = lmopair_to_lmos_[mn][i_mn];

            auto Gvv_temp = linalg::doublet(F_iafc_hat_[mn][i_mn], Gvv_slice, false, false); // (a_ii, f_mn * c_mn) (f_mn * c_mn, 1) -> (a_ii, 1)
            L_ia_buffer[thread][i]->add(Gvv_temp);

            auto F_mnia_slice = F_mnia_hat_[mn][i_mn]->clone(); // (a_ii, 1)
            F_mnia_slice->scale((*rho_oo_)(m, n));
            L_ia_buffer[thread][i]->subtract(F_mnia_slice);
        }

        /* \textcolor{blue}{\begin{equation}
            l^{a_{ii}}_{i} \mathrel{+}= [S^{a_{ii}}_{a_{km}}S^{a_{km}}_{a_{mn}}\overline{\lambda}^{a_{mn}f_{mn}}_{mn}S^{f_{kn}}_{f_{mn}}]T^{f_{kn}c_{kn}}_{kn}S^{c_{km}}_{c_{kn}}\overline{J}^{ic_{km}}_{km} (Toth Eq. 37)
            \end{equation}}
        */
        for (int mn = 0; mn < n_lmo_pairs; ++mn) {
            auto &[m, n] = ij_to_i_j_[mn];

            for (int k_mn = 0; k_mn < lmopair_to_lmos_[mn].size(); ++k_mn) {
                int k = lmopair_to_lmos_[mn][k_mn];
                int km = i_j_to_ij_[k][m], kn = i_j_to_ij_[k][n];
                auto gus = linalg::triplet(S_PNO(km, mn), lambda_iajb_[mn], S_PNO(mn, kn));
                auto charlie = linalg::triplet(gus, T_iajb_[kn], S(kn, km));

                // Done! (From Toth Eq. 31)
                auto airbuds = linalg::doublet(charlie, J_kmic_bar_[km], false, true); // (a, c) * (i, c)
                
                for (int i_km = 0; i_km < lmopair_to_lmos_[km].size(); ++i_km) {
                    int i = lmopair_to_lmos_[km][i_km];
                    int ii = i_j_to_ij_[i][i];
                    auto ryan = submatrix_cols(std::vector<int>(1, i_km), *airbuds)->transpose();
                    L_ia_buffer[thread][i]->add(linalg::doublet(S_PNO(ii, km), ryan)); // Acts like a (7)5-year old
                } // end for

            } // end for
        } // end for


        // l_{i}^{a_{ii}} -= (S^{e_{ki}}_{e_{in}}\overline{\lambda}^{e_{in}f_{in}}_{in}S^{f_{kn}}_{f_{in}})T^{f_{kn}c_{kn}}_{kn}S^{c_{ki}}_{c_{kn}}\overline{J}^{e_{ki}c_{ki}}_{ka_{ii}} (Toth Eq. 38)
        for (int in = 0; in < n_lmo_pairs; ++in) {
            auto &[i, n] = ij_to_i_j_[in];

            for (int k_in = 0; k_in < lmopair_to_lmos_[in].size(); ++k_in) {
                int k = lmopair_to_lmos_[in][k_in];
                int kn = i_j_to_ij_[k][n], in = i_j_to_ij_[i][n];
                auto garfield = linalg::triplet(S_PNO(ki, in), lambda_iajb_bar_[in], S_PNO(in, kn));
                auto lasagne = linalg::triplet(garfield, T_iajb_[kn], S_PNO(kn, ki));

                lasagne->reshape(n_pno_[ki] * n_pno_[ki], 1);
                auto steak = linalg::doublet(J_kaec_bar_[ki], lasagne); // (a, e * c) (e * c, 1) => (a, 1) (steak sauce)
                L_ia_buffer[thread][i]->subtract(steak);
            } // end k_in
        } // end in
    }
        
    /*
    \textcolor{blue}{\begin{equation}
        l^{a_{ii}}_i \mathrel{+}= \frac{1}{2} \widetilde{\lambda}_{in}^{e_{in}f_{in}} [S^{e_{in}}_{e_{ik}}(S^{c_{nk}}_{c_{ik}}u^{f_{nk}c_{nk}}_{nk}S^{f_{in}}_{f_{nk}})]\overline{M}^{c_{ki}e_{ki}}_{ka_{ii}} (Toth Eq. 39)
    \end{equation}}
    */

    for (int in = 0; in < n_lmo_pairs; ++in) {
        auto &[i, n] = ij_to_i_j_[in];

        for (int k_in = 0; k_in < lmopair_to_lmos_[in].size(); ++k_in) {
            int k = lmopair_to_lmos_[in][k_in];
            int nk = i_j_to_ij_[n][k], in = i_j_to_ij_[i][n];
            auto homeslice = linalg::triplet(S_PNO(ki, in), lambda_iajb_[in], S_PNO(in, kn));
            auto pizza = linalg::triplet(homeslice, Tt_iajb_[nk], S_PNO(nk, ki));

            // TODO: Define M_kace_bar_ (a, c * e)
            homeslice->transpose_this(); // (e, c) -> (c, e)
            homeslice->reshape(n_pno_[ki] * n_pno_[ki], 1);
            auto breadsticks = linalg::doublet(M_kace_bar_[ki], homeslice); // (a, c * e) (c * e, 1) => (a, 1) (steak sauce)
            breadsticks->scale(0.5);
            L_ia_buffer[thread][i]->add(breadsticks);
        } // end k_in
    } // end in

    /*
        \textcolor{blue}{\begin{equation}
            l^{a_{ii}}_i \mathrel{-}= \frac{1}{2}[S^{a_{ii}}_{a_{mk}}(S^{a_{mk}}_{a_{mn}}\widetilde{\lambda}^{a_{mn}f_{mn}}_{mn}S^{f_{mn}}_{f_{nk}})]u^{f_{nk}c_{nk}}_{nk}S^{c_{nk}}_{c_{mk}}\overline{M}^{ic_{mk}}_{mk} (Toth Eq. 40)
        \end{equation}}
    */

    for (int mn = 0; mn < n_lmo_pairs; ++mn) {
        auto &[m, n] = ij_to_i_j_[mn];

        for (int k_mn = 0; k_mn < lmopair_to_lmos_[mn].size(); ++k_mn) {
            int k = lmopair_to_lmos_[mn][k_mn];
            int nk = i_j_to_ij_[n][k], mk = i_j_to_ij_[m][k];
            auto johnpork = linalg::triplet(S_PNO(mk, mn), lambda_iajb_[mn], S_PNO(mn, nk));
            auto squiddy = linalg::triplet(johnpork, Tt_iajb_[nk], S_PNO(nk, mk));

            for (int i_mk = 0; i_mk < lmopair_to_lmos_[mk].size(); ++i_mk) {
                auto M_slice = submatrix_rows(std::vector<int>(1, i_mk), *M_mkic_bar_[mk]); // TODO: Need to define M_mkic_bar
                auto calimari = linalg::triplet(S_PNO(ii, mk), squiddy, M_slice, false, false, true);
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

void DLPNOCCSD_Lambda::compute_alpha_ijkl() {
    timer_on("DLPNO-CCSD Lambda : alpha ijkl");

    int n_lmo_pairs = ij_to_i_j_.size();
    int naocc = nalpha_ - nfrzc();

    std::vector<SharedMatrix> alpha_ijkl(n_lmo_pairs);

#pragma omp parallel for
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];

        int nlmo_ij = lmopair_to_lmos_[ij].size();
        alpha_ijkl[ij] = std::make_shared<Matrix>("alpha_ijkl", nlmo_ij, nlmo_ij);

        // alpha_{ij}^{kl} = \lambda_{ij}^{e_{ij} f_{ij}} (S_{e_{kl}}^{e_{ij}} T_{kl}^{e_{kl}f_{kl}} S_{f_{kl}}^{f_{ij}})
        for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
            int k = lmopair_to_lmos_[ij][k_ij];
            for (int l_ij = 0; l_ij < nlmo_ij; ++l_ij) {
                int l = lmopair_to_lmos_[ij][l_ij];
                int kl = i_j_to_ij_[k][l];
                auto gremlin = linalg::triplet(S_PNO(ij, kl), T_iajb_[kl], S_PNO(kl, ij));
                
                (*alpha_ijkl[ij])(k_ij, l_ij) = gremlin->vector_dot(lambda_iajb_[ij]);

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
    auto gamma = compute_gamma();
    auto delta = compute_delta();

    // TODO: Add Crawford Line 516 from GitHub

#pragma omp parallel for
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        auto &[i, j] = ij_to_i_j_[ij];
        int ji = ij_to_ji_[ij], ii = i_j_to_ij_[i][i];

        const int npao_ij = lmopair_to_paos_[ij].size();
        const int naux_ij = lmopair_to_ribfs_[ij].size();
        const int nlmo_ij = lmopair_to_lmos_[ij].size();
        const int npno_ij = n_pno_[ij];

        // l^{a_{ij}b_{ij}}_{ij} += \lambda^{e_{jj}}_j\overline{L}^{a_{ij}b_{ij}}_{ie_{jj}} (Toth Eq. 43a)
        L_temp = linalg::doublet(lambda_ia_[j], L_ijab_bar_[ij], false, false); // TODO: Define L_ijab_bar_
        L_temp->reshape(n_pno_[ij], n_pno_[ij]);
        Ln_iajb[ij]->add(L_temp);

        // l^{a_{ij}b_{ij}}_{ij} -= (2 - P_{ab})[S^{a_{ij}}_{a_{mm}} \lambda^{a_{mm}}_{m} \overline{K}^{mb_{ij}}_{ij}] (Toth Eq. 43b)
        for (int m_ij = 0; m_ij < nlmo_ij; ++m_ij) {
            int m = lmopair_to_lmos_[ij][m_ij];
            int mm = i_j_to_ij_[m][m];

            auto lambda_m_ij = linalg::doublet(S_PNO(ij, mm), lambda_ia_[m]);
            for (int a_ij = 0; a_ij < npno_ij; ++a_ij) {
                for (int b_ij = 0; b_ij < npno_ij; ++b_ij) {
                    (*Ln_iajb[ij])(a_ij, b_ij) -= 2.0 * (*lambda_m_ij)(a_ij, 0) * (*K_ijmb_bar_[ij])(m_ij, b_ij)
                        - (*lambda_m_ij)(b_ij, 0) * (*K_ijmb_bar_[ij])(m_ij, a_ij); // TODO: Define K_ijmb_bar_
                } // end b_ij
            } // end a_ij
        } // end m_ij

        // l^{a_{ij}b_{ij}}_{ij} += (2 - P_{ab}) \lambda^{a_{ii}}_{i} S^{a_{ij}}_{a_{ii}} \widetilde{F}_{jb_{ij}} (Toth Eq. 43c)
        auto lambda_i = linalg::doublet(S_PNO(ij, ii), lambda_ia_[i]);
        for (int a_ij = 0; a_ij < npno_ij; ++a_ij) {
            for (int b_ij = 0; b_ij < npno_ij; ++b_ij) {
                (*Ln_iajb[ij])(a_ij, b_ij) += 2.0 * (*lambda_i)(a_ij, 0) * (F_jb_tilde_[ji])(b_ij, 0) 
                        - (*lambda_i)(b_ij, 0) * (F_jb_tilde_[ji])(a_ij, 0); // TODO: define F_ij_tilde_
            } // end b_ij
        } // end a_ij

        // Necessary three-center integrals
        auto qma_ij = QIA_PNO(ij); // naux_ij * (nlmo_ij, npno_ij)
        auto qab_ij = QAB_PNO(ij); // naux_ij * (npno_ij, npno_ij)
        
        // Toth Eq. 50
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
            L_temp->add(linalg::triplet(qma_ij, alpha_ijkl[ij], qma_ij, true, false, false)); // (k, a) (k, l) (l, b)
            
            L_temp->scale(0.5);
            Ln_iajb[ij]->add(L_temp);
        } // end q_ij

        // l_{ij}^{a_{ij}b_{ij}} \mathrel{+}= \frac{1}{2} (S_{a_{mn}}^{a_{ij}} \widetilde{\lambda}_{mn}^{a_{mn}b_{mn}}S_{b_{mn}}^{b_{ij}})\beta_{mn}^{ij} (Toth Eq. 51)
        for (int m_ij = 0; m_ij < nlmo_ij; ++m_ij) {
            int m = lmopair_to_lmos_[ij][m_ij];
            for (int n_ij = 0; n_ij < nlmo_ij; ++n_ij) {
                int n = lmopair_to_lmos_[ij][n_ij];
                int mn = i_j_to_ij_[m][n];

                auto ethan = linalg::triplet(S_PNO(ij, mn), lambda_iajb_[mn], S_PNO(mn, ij));
                ethan->scale(0.5 * (*beta_mnij_[ij])(m_ij, n_ij)); // TODO: Define beta_mnij_
                Ln_iajb[ij]->add(ethan);
            } // end n_ij
        } // end m_ij

        for (int n_ij = 0; n_ij < nlmo_ij; ++n_ij) {
            int n = lmopair_to_lmos_[ij][n_ij];
            int in = i_j_to_ij_[i][n], jn = i_j_to_ij_[j][n];

            // l_{ij}^{a_{ij}b_{ij}} -= \overline{\lambda}_{jn}^{a_{jn}f_{jn}}[S_{a_{nj}}^{a_{ij}}\widetilde{\gamma}_{in}^{f_{nj}b_{ij}} (Toth Eq. 52a) (TODO: Add J_ikac_non_proj_ next time)
            auto ron = linalg::triplet(S_PNO(ij, jn), lambda_iajb_bar_[jn], S_PNO(jn, in));
            auto weasley = linalg::triplet(ron, gamma[in], S_PNO(in, ij));
            Ln_iajb[ij]->subtract(weasley);

            // l_{ij}^{a_{ij}b_{ij}} -= 0.5 * S_{a_{kj}}^{a_{jn}}(S_{a_{ij}}^{a_{kj}}S_{b_{ki}}^{b_{ij}} K_{ki}^{b_{ki}c_{ki}}S_{c_{ki}}^{c_{kj}})S_{c_{kj}}^{c_{kn}}T_{kn}^{f_{kn}c_{kn}}S_{f_{kn}}^{f_{jn}}] (Toth Eq. 52b)
            for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
                int k = lmopair_to_lmos_[ij][k_ij];
                int kn = i_j_to_ij_[k][n], ik = i_j_to_ij_[i][k], kj = i_j_to_ij_[k][j];
                auto chunk = linalg::triplet(S_PNO(kj, jn), lambda_iajb_bar_[jn], S_PNO(jn, kn));
                auto punk = linalg::triplet(S_PNO(kj, ik), K_iajb_[ik], S_PNO(ik, ij));
                auto munk = linalg::triplet(S_PNO(ij, kj), chunk, T_iajb_[kn]);
                auto kerplunk = linalg::triplet(munk, S_PNO(kn, kj), punk);
                kerplunk->scale(0.5);
                Ln_iajb[ij]->add(kerplunk);
            } // end k_ij

        } // end n_ij
        
        for (int n_ij = 0; n_ij < nlmo_ij; ++n_ij) {
            int n = lmopair_to_lmos_[ij][n_ij];
            int nj = i_j_to_ij_[n][j], in = i_j_to_ij_[i][n];

            // l^{a_{ij}b_{ij}}_{ij} += (2 - P_{ab}) [\frac{1}{2}\widetilde{\lambda}^{f_{ni}a_{ni}}_{ni}[\widetilde{\delta}^{f_{ni}b_{ij}}_{nj}S^{a_{ni}}_{a_{ij}} (Toth Eq. 53a) 
            // (TODO: Add J_ikac_non_proj_ and K_ikac_non_proj_ next time)
            auto forg = linalg::triplet(S_PNO(ij, in), lambda_iajb_[in], S_PNO(in, nj));
            auto parts = linalg::triplet(forg, delta[nj], S_PNO(nj, ij));
            
            Ln_iajb[ij]->add(parts);
            parts->scale(-0.5);
            Ln_iajb[ij]->add(parts->transpose());

            // 0.5 * u^{f_{nk}c_{nk}}_{nk}S^{f_{in}}_{f_{nk}}S^{a_{in}}_{a_{ik}}(S^{a_{ik}}_{a_{ij}}S^{b_{ij}}_{b_{kj}}L^{c_{kj}b_{kj}}_{kj}S^{c_{ik}}_{c_{kj}})S^{c_{nk}}_{c_{ik}}]] (Toth Eq. 53b)
            for (int k_ij = 0; k_ij < nlmo_ij; ++k_ij) {
                int nk = i_j_to_ij_[n][k], ik = i_j_to_ij_[i][k], kj = i_j_to_ij_[k][j];
                auto jon = linalg::triplet(S_PNO(ik, in), lambda_iajb_[in], S_PNO(in, nk));
                auto arbuckle = linalg::triplet(S_PNO(ik, kj), L_iajb_[kj], S_PNO(kj, ij));
                auto john = linalg::triplet(S_PNO(ij, ik), jon, Tt_iajb_[nk]);
                auto pork = linalg::triplet(john, S_PNO(nk, ik), arbuckle);

                pork->scale(0.5);
                Ln_iajb[ij]->add(pork);
                pork->scale(-0.5);
                Ln_iajb[ij]->add(pork->transpose());
            } // end k_ij
        }

        // l^{a_{ij}b_{ij}}_{ij} += \widetilde{\lambda}^{a_{ij}f_{ij}}_{ij}\widetilde{\widetilde{F}}_{f_{ij}b_{ij}} - (2 - P_{ab}) \rho^{\mathrm{VV}}_{a_{mn}c_{mn}}
        // S^{a_{mn}}_{a_{ij}} K^{c_{ij}b_{ij}}_{ij}S^{c_{mn}}_{c_{ij}} (Toth Eq. 54)
        Ln_iajb->add(linalg::doublet(lambda_iajb_[ij], F_vv_double_tilde_[ij])); // TODO: Make F_vv_double_tilde
        for (int m_ij = 0; m_ij < nlmo_ij; ++m_ij) {
            int m = lmopair_to_lmos_[ij][m_ij];
            for (int n_ij = 0; n_ij < nlmo_ij; ++n_ij) {
                int n = lmopair_to_lmos_[ij][n_ij];
                int mn = i_j_to_ij_[m][n];

                auto poob = linalg::triplet(S_PNO(ij, mn), rho_vv_[mn], S_PNO(mn, ij));
                auto missouri = linalg::doublet(poob, K_iajb_[ij]);
                missouri->scale(-2.0);
                Ln_iajb[ij]->add(missouri);
                missouri->scale(-0.5);
                Ln_iajb[ij]->add(missouri->transpose());
            } // end n_ij
        } // end m_ij

        // l^{a_{ij}b_{ij}}_{ij} -= (S^{a_{ij}}_{a_{in}}\widetilde{\lambda}^{a_{in}b_{in}}_{in}S^{b_{ij}}_{b_{in}})\widetilde{\widetilde{F}}_{jn} + 
        // \rho^{\mathrm{OO}}_{jk}(S^{a_{ik}}_{a_{ij}}L^{a_{ik}b_{ik}}_{ik}S^{b_{ik}}_{b_{ij}}) (Toth Eq. 55)
        for (int n_ij = 0; n_ij < nlmo_ij; ++n_ij) {
            int n = lmopair_to_lmos_[ij][n_ij];
            int in = i_j_to_ij_[i][n];
            auto lambda_temp = linalg::triplet(S_PNO(ij, in), lambda_iajb_[in], S_PNO(in, ij));
            lambda_temp->scale((*F_im_double_tilde_)(j, n)); // TODO: Define F_im_double_tilde_
            Ln_iajb[ij]->subtract(lambda_temp);

            auto ohio = linalg::triplet(S_PNO(ij, in), L_iajb_[in], S_PNO(in, ij));
            ohio->scale((*rho_oo_)(j, n));
            Ln_iajb[ij]->subtract(ohio);
        } // end n_ij

#pragma omp parallel for schedule(dynamic, 1)
    for (int ij = 0; ij < n_lmo_pairs; ++ij) {
        int i, j;
        std::tie(i, j) = ij_to_i_j_[ij];
        int ji = ij_to_ji_[ij];

        L_iajb[ij]->zero();

        if (i_j_to_ij_strong_[i][j] != -1) {
            L_iajb[ij]->add(Ln_iajb[ij]);
            L_iajb[ij]->add(Ln_iajb[ji]->transpose());
        } else {
            L_iajb[ij] = std::make_shared<Matrix>(n_pno_[ij], n_pno_[ij]);
            L_iajb[ij]->zero();
        }
    }

    }

    timer_off("DLPNO-CCSD Lambda : Compute L2");
}

void DLPNOCCSD_Lambda::compute_lambda_pno_integrals() {
    throw PSIEXCEPTION("YOU DELINQUENT! WE HAVE NOT IMPLEMENTED THE CODE YET! RAHHHHHHHH!");
}

double DLPNOCCSD_Lambda::compute_energy() {
    return 0.0;
}

}
}