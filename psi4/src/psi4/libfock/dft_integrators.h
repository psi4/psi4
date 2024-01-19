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

#ifndef LIBFOCK_DFT_INTEGRATORS
#define LIBFOCK_DFT_INTEGRATORS

#include "cubature.h"
#include "points.h"

#include "psi4/libfunctional/LibXCfunctional.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libqt/qt.h"
#include "psi4/psi4-dec.h"

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"

#include <cstdlib>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi {

namespace dft_integrators {

inline std::vector<double> rks_quadrature_integrate(std::shared_ptr<BlockOPoints> block,
                                                    std::shared_ptr<SuperFunctional> fworker,
                                                    std::shared_ptr<PointFunctions> pworker) {
    // Block data
    int npoints = block->npoints();
    auto x = block->x();
    auto y = block->y();
    auto z = block->z();
    auto w = block->w();

    // Superfunctional data
    auto zk = fworker->value("V")->pointer();
    auto QTp = fworker->value("Q_TMP")->pointer();

    // Points data
    auto rho_a = pworker->point_value("RHO_A")->pointer();

    // Build quadrature
    std::vector<double> ret(5);
    ret[0] = C_DDOT(npoints, w, 1, zk, 1);
    for (int P = 0; P < npoints; P++) {
        QTp[P] = w[P] * rho_a[P];
    }
    ret[1] = C_DDOT(npoints, w, 1, rho_a, 1);
    ret[2] = C_DDOT(npoints, QTp, 1, x, 1);
    ret[3] = C_DDOT(npoints, QTp, 1, y, 1);
    ret[4] = C_DDOT(npoints, QTp, 1, z, 1);

    return ret;
}

inline void sap_integrator(std::shared_ptr<BlockOPoints> block, const SharedVector& sap_potential,
                           std::shared_ptr<PointFunctions> pworker, SharedMatrix V) {
    // Block data
    const std::vector<int>& function_map = block->functions_local_to_global();
    int nlocal = function_map.size();
    int npoints = block->npoints();
    double* w = block->w();

    // Scratch is updated
    double** Tp = pworker->scratch()[0]->pointer();

    // Points data
    double** phi = pworker->basis_value("PHI")->pointer();
    size_t coll_funcs = pworker->basis_value("PHI")->ncol();

    // V2 Temporary
    int max_functions = V->ncol();
    double** V2p = V->pointer();
    auto& potential = *sap_potential;

    // => LSDA contribution (symmetrized) <= //
    for (int P = 0; P < npoints; P++) {
        std::fill(Tp[P], Tp[P] + nlocal, 0.0);
        C_DAXPY(nlocal, 0.5 * potential[P] * w[P], phi[P], 1, Tp[P], 1);
    }
    // parallel_timer_off("LSDA Phi_tmp", rank);

    // Collect V terms
    C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, phi[0], coll_funcs, Tp[0], max_functions, 0.0, V2p[0],
            max_functions);

    for (int m = 0; m < nlocal; m++) {
        for (int n = 0; n <= m; n++) {
            V2p[m][n] = V2p[n][m] = V2p[m][n] + V2p[n][m];
        }
    }
}

inline void rks_integrator(std::shared_ptr<BlockOPoints> block, std::shared_ptr<SuperFunctional> fworker,
                           std::shared_ptr<PointFunctions> pworker, SharedMatrix V, int ansatz = -1) {
    ansatz = (ansatz == -1 ? fworker->ansatz() : ansatz);
    // printf("Ansatz %d\n", ansatz);

    // Block data
    const auto& function_map = block->functions_local_to_global();
    auto nlocal = function_map.size();
    auto npoints = block->npoints();
    auto w = block->w();

    // Scratch is updated
    auto Tp = pworker->scratch()[0]->pointer();

    // Points data
    auto phi = pworker->basis_value("PHI")->pointer();
    auto rho_a = pworker->point_value("RHO_A")->pointer();
    auto coll_funcs = pworker->basis_value("PHI")->ncol();

    // V2 Temporary
    auto max_functions = V->ncol();
    auto V2p = V->pointer();

    // => LSDA contribution <= //
    //                                         ∂
    // T := 1/2 einsum("p, p, pn -> pn", w, φ, -- f)
    //                                         ∂ρ
    auto v_rho_a = fworker->value("V_RHO_A")->pointer();
    for (int P = 0; P < npoints; P++) {
        std::fill(Tp[P], Tp[P] + nlocal, 0.0);
        C_DAXPY(nlocal, 0.5 * v_rho_a[P] * w[P], phi[P], 1, Tp[P], 1);
    }
    // parallel_timer_off("LSDA Phi_tmp", rank);

    // => GGA contribution <= //
    if (ansatz >= 1) {
        //                                        ∂
        // T += einsum("p, p, xp, xpn -> pnσ", w, -- f, ∇ρ, ∇φ)
        //                                        ∂Γ
        // parallel_timer_on("GGA Phi_tmp", rank);
        auto phix = pworker->basis_value("PHI_X")->pointer();
        auto phiy = pworker->basis_value("PHI_Y")->pointer();
        auto phiz = pworker->basis_value("PHI_Z")->pointer();
        auto rho_ax = pworker->point_value("RHO_AX")->pointer();
        auto rho_ay = pworker->point_value("RHO_AY")->pointer();
        auto rho_az = pworker->point_value("RHO_AZ")->pointer();
        auto v_gamma_aa = fworker->value("V_GAMMA_AA")->pointer();

        for (int P = 0; P < npoints; P++) {
            C_DAXPY(nlocal, w[P] * (2.0 * v_gamma_aa[P] * rho_ax[P]), phix[P], 1, Tp[P], 1);
            C_DAXPY(nlocal, w[P] * (2.0 * v_gamma_aa[P] * rho_ay[P]), phiy[P], 1, Tp[P], 1);
            C_DAXPY(nlocal, w[P] * (2.0 * v_gamma_aa[P] * rho_az[P]), phiz[P], 1, Tp[P], 1);
        }
        // parallel_timer_off("GGA Phi_tmp", rank);
    }

    // ==> Contract T aginst φ, replacing a point index with  an AO index <==
    C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, phi[0], coll_funcs, Tp[0], max_functions, 0.0, V2p[0],
            max_functions);

    // ==> Add the adjoint to complete the LDA and GGA contributions  <==
    for (int m = 0; m < nlocal; m++) {
        for (int n = 0; n <= m; n++) {
            V2p[m][n] = V2p[n][m] = V2p[m][n] + V2p[n][m];
        }
    }

    // => Meta contribution <= //
    if (ansatz >= 2) {
        // parallel_timer_on("Meta", rank);
        auto phix = pworker->basis_value("PHI_X")->pointer();
        auto phiy = pworker->basis_value("PHI_Y")->pointer();
        auto phiz = pworker->basis_value("PHI_Z")->pointer();
        auto v_tau_a = fworker->value("V_TAU_A")->pointer();

        double** phi_w[3];
        phi_w[0] = phix;
        phi_w[1] = phiy;
        phi_w[2] = phiz;

        for (int i = 0; i < 3; i++) {
            auto phiw = phi_w[i];
            for (int P = 0; P < npoints; P++) {
                std::fill(Tp[P], Tp[P] + nlocal, 0.0);
                C_DAXPY(nlocal, v_tau_a[P] * w[P], phiw[P], 1, Tp[P], 1);
            }
            C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, phiw[0], coll_funcs, Tp[0], max_functions, 1.0, V2p[0],
                    max_functions);
        }
        // parallel_timer_off("Meta", rank);
    }
}

inline void rks_gradient_integrator(std::shared_ptr<BasisSet> primary, std::shared_ptr<BlockOPoints> block,
                                    std::shared_ptr<SuperFunctional> fworker, std::shared_ptr<PointFunctions> pworker,
                                    SharedMatrix G, SharedMatrix U, int ansatz = -1) {
    ansatz = (ansatz == -1 ? fworker->ansatz() : ansatz);

    // => Setup scratch pointers, and associated variables <= //
    auto Gp = G->pointer();

    auto Up = U->pointer();
    auto Tp = pworker->scratch()[0]->pointer();
    auto Dp = pworker->D_scratch()[0]->pointer();

    // Fine for now, but not true once we start caching
    auto max_functions = U->ncol();

    // => Per-block setup <= //
    auto npoints = block->npoints();
    auto w = block->w();
    const auto& function_map = block->functions_local_to_global();
    auto nlocal = function_map.size();

    // => Setup accessors to computed values <= //
    auto phi = pworker->basis_value("PHI")->pointer();
    auto phi_x = pworker->basis_value("PHI_X")->pointer();
    auto phi_y = pworker->basis_value("PHI_Y")->pointer();
    auto phi_z = pworker->basis_value("PHI_Z")->pointer();
    auto rho_a = pworker->point_value("RHO_A")->pointer();
    auto coll_funcs = pworker->basis_value("PHI")->ncol();

    // => phi_x type contributions <= //
    // ==> LSDA Contribution <== //
    //                                      ∂        ∂
    // T:= -2 * einsum("p, p, pm -> pm", w, -- f, φ, -- φ, φ, D, δ)
    //                                      ∂ρ       ∂x
    auto v_rho_a = fworker->value("V_RHO_A")->pointer();
    for (int P = 0; P < npoints; P++) {
        std::fill(Tp[P], Tp[P] + nlocal, 0.0);
        C_DAXPY(nlocal, -2.0 * w[P] * v_rho_a[P], phi[P], 1, Tp[P], 1);
    }

    // ==> GGA Contribution (Term 1) <== //
    if (fworker->is_gga()) {
        auto rho_ax = pworker->point_value("RHO_AX")->pointer();
        auto rho_ay = pworker->point_value("RHO_AY")->pointer();
        auto rho_az = pworker->point_value("RHO_AZ")->pointer();
        auto v_gamma_aa = fworker->value("V_GAMMA_AA")->pointer();

        for (int P = 0; P < npoints; P++) {
            C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_ax[P]), phi_x[P], 1, Tp[P], 1);
            C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_ay[P]), phi_y[P], 1, Tp[P], 1);
            C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_az[P]), phi_z[P], 1, Tp[P], 1);
        }
    }

    // ==> Complete Terms <== //
    // U := einsum("pm, mn -> pn")
    C_DGEMM('N', 'N', npoints, nlocal, nlocal, 1.0, Tp[0], max_functions, Dp[0], max_functions, 0.0, Up[0],
            max_functions);

    //                                      ∂
    // dE += einsum("pn, pnx, ni -> ix", U, -- φ, δ)
    //                                      ∂x
    for (int ml = 0; ml < nlocal; ml++) {
        auto A = primary->function_to_center(function_map[ml]);
        Gp[A][0] += C_DDOT(npoints, &Up[0][ml], max_functions, &phi_x[0][ml], coll_funcs);
        Gp[A][1] += C_DDOT(npoints, &Up[0][ml], max_functions, &phi_y[0][ml], coll_funcs);
        Gp[A][2] += C_DDOT(npoints, &Up[0][ml], max_functions, &phi_z[0][ml], coll_funcs);
    }

    // => GGA Contribution (Term 2) <= //
    if (fworker->is_gga()) {
        double** phi_xx = pworker->basis_value("PHI_XX")->pointer();
        double** phi_xy = pworker->basis_value("PHI_XY")->pointer();
        double** phi_xz = pworker->basis_value("PHI_XZ")->pointer();
        double** phi_yy = pworker->basis_value("PHI_YY")->pointer();
        double** phi_yz = pworker->basis_value("PHI_YZ")->pointer();
        double** phi_zz = pworker->basis_value("PHI_ZZ")->pointer();
        double* rho_ax = pworker->point_value("RHO_AX")->pointer();
        double* rho_ay = pworker->point_value("RHO_AY")->pointer();
        double* rho_az = pworker->point_value("RHO_AZ")->pointer();
        double* v_gamma_aa = fworker->value("V_GAMMA_AA")->pointer();

        C_DGEMM('N', 'N', npoints, nlocal, nlocal, 1.0, phi[0], coll_funcs, Dp[0], max_functions, 0.0, Up[0],
                max_functions);

        // x
        for (int P = 0; P < npoints; P++) {
            std::fill(Tp[P], Tp[P] + nlocal, 0.0);
            C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_ax[P]), Up[P], 1, Tp[P], 1);
        }
        for (int ml = 0; ml < nlocal; ml++) {
            int A = primary->function_to_center(function_map[ml]);
            Gp[A][0] += C_DDOT(npoints, &Tp[0][ml], max_functions, &phi_xx[0][ml], coll_funcs);
            Gp[A][1] += C_DDOT(npoints, &Tp[0][ml], max_functions, &phi_xy[0][ml], coll_funcs);
            Gp[A][2] += C_DDOT(npoints, &Tp[0][ml], max_functions, &phi_xz[0][ml], coll_funcs);
        }

        // y
        for (int P = 0; P < npoints; P++) {
            std::fill(Tp[P], Tp[P] + nlocal, 0.0);
            C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_ay[P]), Up[P], 1, Tp[P], 1);
        }
        for (int ml = 0; ml < nlocal; ml++) {
            int A = primary->function_to_center(function_map[ml]);
            Gp[A][0] += C_DDOT(npoints, &Tp[0][ml], max_functions, &phi_xy[0][ml], coll_funcs);
            Gp[A][1] += C_DDOT(npoints, &Tp[0][ml], max_functions, &phi_yy[0][ml], coll_funcs);
            Gp[A][2] += C_DDOT(npoints, &Tp[0][ml], max_functions, &phi_yz[0][ml], coll_funcs);
        }

        // z
        for (int P = 0; P < npoints; P++) {
            std::fill(Tp[P], Tp[P] + nlocal, 0.0);
            C_DAXPY(nlocal, -2.0 * w[P] * (2.0 * v_gamma_aa[P] * rho_az[P]), Up[P], 1, Tp[P], 1);
        }
        for (int ml = 0; ml < nlocal; ml++) {
            int A = primary->function_to_center(function_map[ml]);
            Gp[A][0] += C_DDOT(npoints, &Tp[0][ml], max_functions, &phi_xz[0][ml], coll_funcs);
            Gp[A][1] += C_DDOT(npoints, &Tp[0][ml], max_functions, &phi_yz[0][ml], coll_funcs);
            Gp[A][2] += C_DDOT(npoints, &Tp[0][ml], max_functions, &phi_zz[0][ml], coll_funcs);
        }
    }

    // => Meta Contribution <= //
    if (fworker->is_meta()) {
        double** phi_xx = pworker->basis_value("PHI_XX")->pointer();
        double** phi_xy = pworker->basis_value("PHI_XY")->pointer();
        double** phi_xz = pworker->basis_value("PHI_XZ")->pointer();
        double** phi_yy = pworker->basis_value("PHI_YY")->pointer();
        double** phi_yz = pworker->basis_value("PHI_YZ")->pointer();
        double** phi_zz = pworker->basis_value("PHI_ZZ")->pointer();
        double* v_tau_a = fworker->value("V_TAU_A")->pointer();

        double** phi_i[3];
        phi_i[0] = phi_x;
        phi_i[1] = phi_y;
        phi_i[2] = phi_z;

        double** phi_ij[3][3];
        phi_ij[0][0] = phi_xx;
        phi_ij[0][1] = phi_xy;
        phi_ij[0][2] = phi_xz;
        phi_ij[1][0] = phi_xy;
        phi_ij[1][1] = phi_yy;
        phi_ij[1][2] = phi_yz;
        phi_ij[2][0] = phi_xz;
        phi_ij[2][1] = phi_yz;
        phi_ij[2][2] = phi_zz;

        for (int i = 0; i < 3; i++) {
            double*** phi_j = phi_ij[i];
            C_DGEMM('N', 'N', npoints, nlocal, nlocal, 1.0, phi_i[i][0], coll_funcs, Dp[0], max_functions, 0.0, Up[0],
                    max_functions);
            for (int P = 0; P < npoints; P++) {
                std::fill(Tp[P], Tp[P] + nlocal, 0.0);
                C_DAXPY(nlocal, -2.0 * w[P] * (v_tau_a[P]), Up[P], 1, Tp[P], 1);
            }
            for (int ml = 0; ml < nlocal; ml++) {
                int A = primary->function_to_center(function_map[ml]);
                Gp[A][0] += C_DDOT(npoints, &Tp[0][ml], max_functions, &phi_j[0][0][ml], coll_funcs);
                Gp[A][1] += C_DDOT(npoints, &Tp[0][ml], max_functions, &phi_j[1][0][ml], coll_funcs);
                Gp[A][2] += C_DDOT(npoints, &Tp[0][ml], max_functions, &phi_j[2][0][ml], coll_funcs);
            }
        }
    }
}

}  // namespace dft_integrators
}  // namespace psi
#endif
