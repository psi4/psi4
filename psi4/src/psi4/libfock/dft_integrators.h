/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
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
#include "psi4/libfunctional/functional.h"
#include "psi4/libfunctional/superfunctional.h"
#include "psi4/libqt/qt.h"
#include "psi4/psi4-dec.h"

#include "psi4/libmints/basisset.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/molecule.h"
#include "psi4/libmints/petitelist.h"
#include "psi4/libmints/vector.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/process.h"

#include <cstdlib>
#include <numeric>
#include <sstream>
#include <string>

#ifdef _OPENMP
#include <omp.h>
#endif


namespace psi {

namespace dft_integrators {

inline std::vector<double> rks_quadrature_integrate(std::shared_ptr<BlockOPoints> block,
                                                    std::shared_ptr<SuperFunctional> fworker,
                                                    std::shared_ptr<PointFunctions> pworker){
        // Block data
    int npoints = block->npoints();
    double* x = block->x();
    double* y = block->y();
    double* z = block->z();
    double* w = block->w();

    // Superfunctional data
    double* zk = fworker->value("V")->pointer();
    double* QTp = fworker->value("Q_TMP")->pointer();

    // Points data
    double* rho_a = pworker->point_value("RHO_A")->pointer();

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

inline void rks_integrator(std::shared_ptr<BlockOPoints> block, std::shared_ptr<SuperFunctional> fworker,
                                    std::shared_ptr<PointFunctions> pworker, SharedMatrix V, int ansatz=-1){


    ansatz = (ansatz == -1 ? fworker->ansatz() : ansatz);
    // printf("Ansatz %d\n", ansatz);

    // Block data
    const std::vector<int>& function_map = block->functions_local_to_global();
    int nlocal = function_map.size();
    int npoints = block->npoints();
    double* w = block->w();

    // Scratch is updated
    double** Tp = pworker->scratch()[0]->pointer();

    // Points data
    double** phi = pworker->basis_value("PHI")->pointer();
    double* rho_a = pworker->point_value("RHO_A")->pointer();

    // Fine for now, but not true once we start caching
    int max_functions = V->ncol();
    double** V2p = V->pointer();

    // => LSDA contribution (symmetrized) <= //
    double* v_rho_a = fworker->value("V_RHO_A")->pointer();
    for (int P = 0; P < npoints; P++) {
        std::fill(Tp[P], Tp[P] + nlocal, 0.0);
        C_DAXPY(nlocal, 0.5 * v_rho_a[P] * w[P], phi[P], 1, Tp[P], 1);
    }
    // parallel_timer_off("LSDA Phi_tmp", rank);

    // => GGA contribution (symmetrized) <= //
    if (ansatz >= 1) {
        // parallel_timer_on("GGA Phi_tmp", rank);
        double** phix = pworker->basis_value("PHI_X")->pointer();
        double** phiy = pworker->basis_value("PHI_Y")->pointer();
        double** phiz = pworker->basis_value("PHI_Z")->pointer();
        double* rho_ax = pworker->point_value("RHO_AX")->pointer();
        double* rho_ay = pworker->point_value("RHO_AY")->pointer();
        double* rho_az = pworker->point_value("RHO_AZ")->pointer();
        double* v_sigma_aa = fworker->value("V_GAMMA_AA")->pointer();

        for (int P = 0; P < npoints; P++) {
            C_DAXPY(nlocal, w[P] * (2.0 * v_sigma_aa[P] * rho_ax[P]), phix[P], 1, Tp[P], 1);
            C_DAXPY(nlocal, w[P] * (2.0 * v_sigma_aa[P] * rho_ay[P]), phiy[P], 1, Tp[P], 1);
            C_DAXPY(nlocal, w[P] * (2.0 * v_sigma_aa[P] * rho_az[P]), phiz[P], 1, Tp[P], 1);
        }
        // parallel_timer_off("GGA Phi_tmp", rank);
    }

    // Collect V terms
    C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, phi[0], max_functions, Tp[0], max_functions, 0.0, V2p[0],
        max_functions);

    for (int m = 0; m < nlocal; m++) {
        for (int n = 0; n <= m; n++) {
            V2p[m][n] = V2p[n][m] = V2p[m][n] + V2p[n][m];
        }
    }

    // => Meta contribution <= //
    if (ansatz >= 2) {
        // parallel_timer_on("Meta", rank);
        double** phix = pworker->basis_value("PHI_X")->pointer();
        double** phiy = pworker->basis_value("PHI_Y")->pointer();
        double** phiz = pworker->basis_value("PHI_Z")->pointer();
        double* v_tau_a = fworker->value("V_TAU_A")->pointer();

        double** phi_w[3];
        phi_w[0] = phix;
        phi_w[1] = phiy;
        phi_w[2] = phiz;

        for (int i = 0; i < 3; i++) {
            double** phiw = phi_w[i];
            for (int P = 0; P < npoints; P++) {
                std::fill(Tp[P], Tp[P] + nlocal, 0.0);
                C_DAXPY(nlocal, v_tau_a[P] * w[P], phiw[P], 1, Tp[P], 1);
            }
            C_DGEMM('T', 'N', nlocal, nlocal, npoints, 1.0, phiw[0], max_functions, Tp[0], max_functions, 1.0,
                    V2p[0], max_functions);
        }
        // parallel_timer_off("Meta", rank);
    }
}


}  // End dft_integrator namespace
}  // End psi namespace
#endif
