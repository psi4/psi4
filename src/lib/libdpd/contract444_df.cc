/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

/*! \file
    \ingroup DPD
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cmath>
#include <libqt/qt.h>
#include <libpsio/psio.h>
#include "dpd.h"
#include <psi4-dec.h>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace psi {
/**
 * @brief contract444_df Forms integrals on the fly from density fitted integrals and contracts them against T amplitudes
 *
 * The integral tensor should arrive packed, in the order (a>=b|Q), where Q is the auxilliary basis index
 *
 * @param B The DF integral tensor, arranged with the auxilliary index as the fast-running index
 * @param tau1_AO The incoming amplitudes
 * @param tau2_AO The outgoing amplitudes
 * @param alpha The prefactor for the contraction
 * @param beta The prefactor for the existing tensor
 */
#define PermSym 1
int DPD::contract444_df(dpdbuf4 *B, dpdbuf4 *tau_in, dpdbuf4 *tau_out, double alpha, double beta)
{
    // Zero out the target buffer
    buf4_scm(tau_out, beta);

    // Open up the entire tau matrix
    for(int h = 0; h < tau_in->params->nirreps; ++h){
        buf4_mat_irrep_init(tau_out, h);
        buf4_mat_irrep_rd(tau_out, h);
        buf4_mat_irrep_init(tau_in, h);
        buf4_mat_irrep_rd(tau_in, h);
    }

    int nthreadz = 1;
#ifdef _OPENMP
    nthreadz = omp_get_max_threads();
#endif

    for(int Gpr = 0; Gpr < B->params->nirreps; ++Gpr){
        buf4_mat_irrep_init(B, Gpr);
        buf4_mat_irrep_rd(B, Gpr);

        int **orbs = B->params->roworb[Gpr];
        for(int pr = 0; pr < B->params->rowtot[Gpr]; ++pr){
            int p = orbs[pr][0];
            int r = orbs[pr][1];
            int psym = B->params->psym[p];
            int rsym = B->params->qsym[r];

            // Scale the diagonals by 0.5, to avoid "if" statements below
            if(p==r)
                for(int Q = 0; Q < B->params->coltot[Gpr]; ++Q)
                    B->matrix[Gpr][pr][Q] *= 0.5;

            for(int qs = 0; qs < pr; ++qs){
                int q = orbs[qs][0];
                int s = orbs[qs][1];
                int qsym = B->params->psym[q];
                int ssym = B->params->qsym[s];

                int len = B->params->coltot[Gpr];
                if(len == 0)
                    continue;

                // Build the integral
                double prqs = alpha * C_DDOT(len, B->matrix[Gpr][pr], 1, B->matrix[Gpr][qs], 1);

                int Gpq = psym ^ qsym;
                int Grs = rsym ^ ssym;
                int Gps = psym ^ ssym;
                int Grq = rsym ^ qsym;
                int pq = tau_out->params->rowidx[p][q];
                int qp = tau_out->params->rowidx[q][p];
                int rs = tau_out->params->rowidx[r][s];
                int sr = tau_out->params->rowidx[s][r];
                int rq = tau_out->params->rowidx[r][q];
                int qr = tau_out->params->rowidx[q][r];
                int ps = tau_out->params->rowidx[p][s];
                int sp = tau_out->params->rowidx[s][p];

                len = tau_in->params->coltot[Gpq];
                if(len){
                    // T_pq_ij <- (pr|qs) T_rs_ij
                    C_DAXPY(len, prqs, tau_in->matrix[Gpq][rs], 1, tau_out->matrix[Gpq][pq], 1);
                    // T_qp_ij <- (qs|pr) T_sr_ij
                    C_DAXPY(len, prqs, tau_in->matrix[Gpq][sr], 1, tau_out->matrix[Gpq][qp], 1);
                    // T_rs_ij <- (rp|sq) T_pq_ij
                    C_DAXPY(len, prqs, tau_in->matrix[Grs][pq], 1, tau_out->matrix[Grs][rs], 1);
                    // T_sr_ij <- (sq|rp) T_qp_ij
                    C_DAXPY(len, prqs, tau_in->matrix[Grs][qp], 1, tau_out->matrix[Grs][sr], 1);
                }
                len = tau_in->params->coltot[Grq];
                if(len){
                    // T_rq_ij <- (rp|qs) T_ps_ij
                    C_DAXPY(len, prqs, tau_in->matrix[Grq][ps], 1, tau_out->matrix[Grq][rq], 1);
                    // T_qr_ij <- (qs|rp) T_sp_ij
                    C_DAXPY(len, prqs, tau_in->matrix[Grq][sp], 1, tau_out->matrix[Grq][qr], 1);
                    // T_ps_ij <- (pr|sq) T_rq_ij
                    C_DAXPY(len, prqs, tau_in->matrix[Gps][rq], 1, tau_out->matrix[Gps][ps], 1);
                    // T_sp_ij <- (sq|pr) T_qr_ij
                    C_DAXPY(len, prqs, tau_in->matrix[Gps][qr], 1, tau_out->matrix[Gps][sp], 1);
                }

            } // End qs loop

            // Now, add in the diagonal elements, pr == qs
            int Gpq = 0;
            int Grq = rsym ^ psym;
            int Grs = 0;
            int Gps = Grq;
            int pq = tau_out->params->rowidx[p][p];
            int rs = tau_out->params->rowidx[r][r];
            int rq = tau_out->params->rowidx[r][p];
            int ps = tau_out->params->rowidx[p][r];

            int len = B->params->coltot[Gpr];
            if(len == 0)
                continue;

            // Build the integral
            double prqs = alpha * C_DDOT(len, B->matrix[Gpr][pr], 1, B->matrix[Gpr][pr], 1);

            len = tau_in->params->coltot[Gpq];
            if(len){
                // T_pq_ij <- (pr|qs) T_rs_ij
                C_DAXPY(len, prqs, tau_in->matrix[Gpq][rs], 1, tau_out->matrix[Gpq][pq], 1);
                // T_rs_ij <- (rp|sq) T_pq_ij
                C_DAXPY(len, prqs, tau_in->matrix[Grs][pq], 1, tau_out->matrix[Grs][rs], 1);
            }
            len = tau_in->params->coltot[Grq];
            if(len){
                // T_rq_ij <- (rp|qs) T_ps_ij
                C_DAXPY(len, prqs, tau_in->matrix[Grq][ps], 1, tau_out->matrix[Grq][rq], 1);
                // T_ps_ij <- (pr|sq) T_rq_ij
                C_DAXPY(len, prqs, tau_in->matrix[Gps][rq], 1, tau_out->matrix[Gps][ps], 1);
            }

        } // End pr loop
        buf4_mat_irrep_close(B, Gpr);
    } // End Gpr loop

    // Close the B matrix
    for(int h = 0; h < B->params->nirreps; ++h){
        buf4_mat_irrep_wrt(tau_out, h);
        buf4_mat_irrep_close(tau_out, h);
        buf4_mat_irrep_close(tau_in, h);
    }

    return 0;
}


}
