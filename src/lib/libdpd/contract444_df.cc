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
namespace psi {
/**
 * @brief contract444_df Forms integrals on the fly from density fitted integrals and contracts them against T amplitudes
 * @param B The DF integral tensor, arranged with the auxilliary index as the fast-running index
 * @param tau1_AO The incoming amplitudes
 * @param tau2_AO The outgoing amplitudes
 * @param alpha The prefactor for the contraction
 * @param beta The prefactor for the existing tensor
 */
int DPD::contract444_df(dpdbuf4 *B, dpdbuf4 *tau_in, dpdbuf4 *tau_out, double alpha, double beta)
{
    // Zero out the target buffer
    buf4_scm(tau_out, beta);

    // Open up the entire B matrix
    for(int h = 0; h < B->params->nirreps; ++h){
        buf4_mat_irrep_init(tau_out, h);
        buf4_mat_irrep_rd(tau_out, h);
        buf4_mat_irrep_init(tau_in, h);
        buf4_mat_irrep_rd(tau_in, h);
    }

    for(int Gpr = 0; Gpr < B->params->nirreps; ++Gpr){
        buf4_mat_irrep_init(B, Gpr);
        buf4_mat_irrep_rd(B, Gpr);

        for(int pr = 0; pr < B->params->rowtot[Gpr]; ++pr){
            int p = tau_out->params->roworb[Gpr][pr][0];
            int r = tau_out->params->roworb[Gpr][pr][1];
            int psym = tau_out->params->psym[p];
            int rsym = tau_out->params->qsym[r];
            for(int qs = 0; qs < B->params->rowtot[Gpr]; ++qs){
                int q = tau_in->params->roworb[Gpr][qs][0];
                int s = tau_in->params->roworb[Gpr][qs][1];
                int qsym = tau_in->params->psym[q];
                int ssym = tau_in->params->qsym[s];

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
                int rq = tau_out->params->rowidx[r][q];
                int ps = tau_out->params->rowidx[p][s];

                // T_pq_ij <- (pr|qs) T_rs_ij
                len = tau_in->params->coltot[Gpq];
                if(len){
                    C_DAXPY(len, prqs, tau_in->matrix[Gpq][rs], 1, tau_out->matrix[Gpq][pq], 1);
                    // T_qp_ij <- -(pr|qs) T_rs_ij
//                    C_DAXPY(len, -prqs, tau_in->matrix[Gpq][rs], 1, tau_out->matrix[Gpq][qp], 1);
                }

#if 0
                // T_rq_ij <- (rp|qs) T_ps_ij
                if(p != r){
                    len = tau1_AO->params->coltot[Grq];
                    if(len)
                        C_DAXPY(len, prqs, tau1_AO->matrix[Grq][ps], 1, tau2_AO->matrix[Grq][rq], 1);
                }

                // T_ps_ij <- (pr|sq) T_rq_ij
                if(q != s){
                    len = tau1_AO->params->coltot[Gps];
                    if(len)
                        C_DAXPY(len, prqs, tau1_AO->matrix[Gps][rq], 1, tau2_AO->matrix[Gps][ps], 1);
                }

                // T_rs_ij <- (rp|sq) T_pq_ij
                if(p != r && q != s){
                    len = tau1_AO->params->coltot[Grs];
                    if(len)
                        C_DAXPY(len, prqs, tau1_AO->matrix[Grs][pq], 1, tau2_AO->matrix[Grs][rs], 1);
                }
#endif
            }
        }

        buf4_mat_irrep_close(B, Gpr);
    }

    // Close the B matrix
    for(int h = 0; h < B->params->nirreps; ++h){
        buf4_mat_irrep_wrt(tau_out, h);
        buf4_mat_irrep_close(tau_out, h);
        buf4_mat_irrep_close(tau_in, h);
    }

    return 0;
}


}
