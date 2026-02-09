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

                // TODO: define J_kmic_bar_
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

                // TODO: Define J_kaec_bar_ (a, e * c)
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

    timer_off("DLPNO-CCSD Lambda : Compute L1");
}


void DLPNOCCSD_Lambda::compute_lambda_pno_integrals() {
    throw PSIEXCEPTION("YOU DELINQUENT! WE HAVE NOT IMPLEMENTED THE CODE YET! RAHHHHHHHH!");
}

double DLPNOCCSD_Lambda::compute_energy() {
    return 0.0;
}

}
}