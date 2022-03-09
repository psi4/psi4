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

/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/

/*! \defgroup CCENERGY ccenergy: Compute the Coupled-Cluster Energy */

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libiwl/iwl.h"
#include "psi4/libqt/qt.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/cc/ccwave.h"

namespace psi {
namespace ccenergy {

int CCEnergyWavefunction::AO_contribute(struct iwlbuf *InBuf, dpdbuf4 *tau1_AO, dpdbuf4 *tau2_AO) {
    int p, q, r, s;
    double value = 0.0;
    int Gp, Gq, Gr, Gs, Gpr, Gps, Gqr, Gqs, Grp, Gsp, Grq, Gsq;
    int pr, ps, qr, qs, rp, rq, sp, sq, pq, rs;
    int count = 0;

    auto lblptr = InBuf->labels;
    auto valptr = InBuf->values;

    for (int idx = 4 * InBuf->idx; InBuf->idx < InBuf->inbuf; InBuf->idx++) {
        p = std::abs((int)lblptr[idx++]);
        q = (int)lblptr[idx++];
        r = (int)lblptr[idx++];
        s = (int)lblptr[idx++];

        value = (double)valptr[InBuf->idx];
        count++;

        Gp = tau1_AO->params->psym[p];
        Gq = tau1_AO->params->psym[q];
        Gr = tau1_AO->params->psym[r];
        Gs = tau1_AO->params->psym[s];

        Gpr = Grp = Gp ^ Gr;
        Gps = Gsp = Gp ^ Gs;
        Gqr = Grq = Gq ^ Gr;
        Gqs = Gsq = Gq ^ Gs;

        pq = tau1_AO->params->rowidx[p][q];
        rs = tau1_AO->params->rowidx[r][s];

        pr = tau1_AO->params->rowidx[p][r];
        rp = tau1_AO->params->rowidx[r][p];
        ps = tau1_AO->params->rowidx[p][s];
        sp = tau1_AO->params->rowidx[s][p];
        qr = tau1_AO->params->rowidx[q][r];
        rq = tau1_AO->params->rowidx[r][q];
        qs = tau1_AO->params->rowidx[q][s];
        sq = tau1_AO->params->rowidx[s][q];

        /* (pq|rs) */
        if (tau1_AO->params->coltot[Gpr])
            C_DAXPY(tau1_AO->params->coltot[Gpr], value, tau1_AO->matrix[Gpr][qs], 1, tau2_AO->matrix[Gpr][pr], 1);

        if (p != q && r != s && pq != rs) {
            /* (pq|sr) */
            if (tau1_AO->params->coltot[Gps])
                C_DAXPY(tau1_AO->params->coltot[Gps], value, tau1_AO->matrix[Gps][qr], 1, tau2_AO->matrix[Gps][ps], 1);

            /* (qp|rs) */
            if (tau1_AO->params->coltot[Gqr])
                C_DAXPY(tau1_AO->params->coltot[Gqr], value, tau1_AO->matrix[Gqr][ps], 1, tau2_AO->matrix[Gqr][qr], 1);

            /* (qp|sr) */
            if (tau1_AO->params->coltot[Gqs])
                C_DAXPY(tau1_AO->params->coltot[Gqs], value, tau1_AO->matrix[Gqs][pr], 1, tau2_AO->matrix[Gqs][qs], 1);

            /* (rs|pq) */
            if (tau1_AO->params->coltot[Grp])
                C_DAXPY(tau1_AO->params->coltot[Grp], value, tau1_AO->matrix[Grp][sq], 1, tau2_AO->matrix[Grp][rp], 1);

            /* (sr|pq) */
            if (tau1_AO->params->coltot[Gsp])
                C_DAXPY(tau1_AO->params->coltot[Gsp], value, tau1_AO->matrix[Gsp][rq], 1, tau2_AO->matrix[Gsp][sp], 1);

            /* (rs|qp) */
            if (tau1_AO->params->coltot[Grq])
                C_DAXPY(tau1_AO->params->coltot[Grq], value, tau1_AO->matrix[Grq][sp], 1, tau2_AO->matrix[Grq][rq], 1);

            /* (sr|qp) */
            if (tau1_AO->params->coltot[Gsq])
                C_DAXPY(tau1_AO->params->coltot[Gsq], value, tau1_AO->matrix[Gsq][rp], 1, tau2_AO->matrix[Gsq][sq], 1);

        } else if (p != q && r != s && pq == rs) {
            /* (pq|sr) */
            if (tau1_AO->params->coltot[Gps])
                C_DAXPY(tau1_AO->params->coltot[Gps], value, tau1_AO->matrix[Gps][qr], 1, tau2_AO->matrix[Gps][ps], 1);

            /* (qp|rs) */
            if (tau1_AO->params->coltot[Gqr])
                C_DAXPY(tau1_AO->params->coltot[Gqr], value, tau1_AO->matrix[Gqr][ps], 1, tau2_AO->matrix[Gqr][qr], 1);

            /* (qp|sr) */
            if (tau1_AO->params->coltot[Gqs])
                C_DAXPY(tau1_AO->params->coltot[Gqs], value, tau1_AO->matrix[Gqs][pr], 1, tau2_AO->matrix[Gqs][qs], 1);

        } else if (p != q && r == s) {
            /* (qp|rs) */
            if (tau1_AO->params->coltot[Gqr])
                C_DAXPY(tau1_AO->params->coltot[Gqr], value, tau1_AO->matrix[Gqr][ps], 1, tau2_AO->matrix[Gqr][qr], 1);

            /* (rs|pq) */
            if (tau1_AO->params->coltot[Grp])
                C_DAXPY(tau1_AO->params->coltot[Grp], value, tau1_AO->matrix[Grp][sq], 1, tau2_AO->matrix[Grp][rp], 1);

            /* (rs|qp) */
            if (tau1_AO->params->coltot[Grq])
                C_DAXPY(tau1_AO->params->coltot[Grq], value, tau1_AO->matrix[Grq][sp], 1, tau2_AO->matrix[Grq][rq], 1);

        }

        else if (p == q && r != s) {
            /* (pq|sr) */
            if (tau1_AO->params->coltot[Gps])
                C_DAXPY(tau1_AO->params->coltot[Gps], value, tau1_AO->matrix[Gps][qr], 1, tau2_AO->matrix[Gps][ps], 1);

            /* (rs|pq) */
            if (tau1_AO->params->coltot[Grp])
                C_DAXPY(tau1_AO->params->coltot[Grp], value, tau1_AO->matrix[Grp][sq], 1, tau2_AO->matrix[Grp][rp], 1);

            /* (sr|pq) */
            if (tau1_AO->params->coltot[Gsp])
                C_DAXPY(tau1_AO->params->coltot[Gsp], value, tau1_AO->matrix[Gsp][rq], 1, tau2_AO->matrix[Gsp][sp], 1);

        }

        else if (p == q && r == s && pq != rs) {
            /* (rs|pq) */
            if (tau1_AO->params->coltot[Grp])
                C_DAXPY(tau1_AO->params->coltot[Grp], value, tau1_AO->matrix[Grp][sq], 1, tau2_AO->matrix[Grp][rp], 1);
        }
    }

    return count;
}

}  // namespace ccenergy
}  // namespace psi
