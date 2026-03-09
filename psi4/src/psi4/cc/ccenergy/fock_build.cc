/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libiwl/iwl.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"
#include "MOInfo.h"
#include "psi4/cc/ccwave.h"

#include "psi4/libmints/wavefunction.h"
#include "psi4/libtrans/mospace.h"
#include "psi4/libmints/matrix.h"

namespace psi {
namespace ccenergy {

#define INDEX(i, j) ((i) >= (j) ? EXPLICIT_IOFF(i) + (j) : EXPLICIT_IOFF(j) + (i))

void CCEnergyWavefunction::rhf_fock_build(double **fock, double **D) {
    int lastbuf, p, q, r, s, pq, rs;
    double value;
    Value *valptr;
    Label *lblptr;
    struct iwlbuf InBuf;

    auto nso = moinfo_.nso;

    auto H = H_->to_block_matrix();

    for (int i = 0; i < nso; i++) {
        for (int j = 0; j <= i; j++) {
            fock[i][j] = fock[j][i] = H[i][j];
        }
    }

    /* two-electron contributions */
    iwl_buf_init(&InBuf, PSIF_SO_TEI, 0.0, 1, 1);
    do {
        lastbuf = InBuf.lastbuf;
        lblptr = InBuf.labels;
        valptr = InBuf.values;

        for (int idx = 4 * InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
            p = std::abs((int)lblptr[idx++]);
            q = (int)lblptr[idx++];
            r = (int)lblptr[idx++];
            s = (int)lblptr[idx++];
            value = (double)valptr[InBuf.idx];

            pq = INDEX(p, q);
            rs = INDEX(r, s);

            /*
            outfile->Printf( "%d %d %d %d [%d] [%d] %20.15f\n", p, q, r, s, pq, rs, value);
            */

            /* (pq|rs) */
            fock[p][q] += 2.0 * D[r][s] * value;
            fock[p][r] -= D[q][s] * value;

            if (p != q && r != s && pq != rs) {
                /* (pq|sr) */
                fock[p][q] += 2.0 * D[s][r] * value;
                fock[p][s] -= D[q][r] * value;

                /* (qp|rs) */
                fock[q][p] += 2.0 * D[r][s] * value;
                fock[q][r] -= D[p][s] * value;

                /* (qp|sr) */
                fock[q][p] += 2.0 * D[s][r] * value;
                fock[q][s] -= D[p][r] * value;

                /* (rs|pq) */
                fock[r][s] += 2.0 * D[p][q] * value;
                fock[r][p] -= D[s][q] * value;

                /* (rs|qp) */
                fock[r][s] += 2.0 * D[q][p] * value;
                fock[r][q] -= D[s][p] * value;

                /* (sr|pq) */
                fock[s][r] += 2.0 * D[p][q] * value;
                fock[s][p] -= D[r][q] * value;

                /* (sr|qp) */
                fock[s][r] += 2.0 * D[q][p] * value;
                fock[s][q] -= D[r][p] * value;
            } else if (p != q && r != s && pq == rs) {
                /* (pq|sr) */
                fock[p][q] += 2.0 * D[s][r] * value;
                fock[p][s] -= D[q][r] * value;

                /* (qp|rs) */
                fock[q][p] += 2.0 * D[r][s] * value;
                fock[q][r] -= D[p][s] * value;

                /* (qp|sr) */
                fock[q][p] += 2.0 * D[s][r] * value;
                fock[q][s] -= D[p][r] * value;

            } else if (p != q && r == s) {
                /* (qp|rs) */
                fock[q][p] += 2.0 * D[r][s] * value;
                fock[q][r] -= D[p][s] * value;

                /* (rs|pq) */
                fock[r][s] += 2.0 * D[p][q] * value;
                fock[r][p] -= D[s][q] * value;

                /* (rs|qp) */
                fock[r][s] += 2.0 * D[q][p] * value;
                fock[r][q] -= D[s][p] * value;

            } else if (p == q && r != s) {
                /* (pq|sr) */
                fock[p][q] += 2.0 * D[s][r] * value;
                fock[p][s] -= D[q][r] * value;

                /* (rs|pq) */
                fock[r][s] += 2.0 * D[p][q] * value;
                fock[r][p] -= D[s][q] * value;

                /* (sr|pq) */
                fock[s][r] += 2.0 * D[p][q] * value;
                fock[s][p] -= D[r][q] * value;

            } else if (p == q && r == s && pq != rs) {
                /* (rs|pq) */
                fock[r][s] += 2.0 * D[p][q] * value;
                fock[r][p] -= D[s][q] * value;
            }
        }

        if (!lastbuf) iwl_buf_fetch(&InBuf);

    } while (!lastbuf);
    iwl_buf_close(&InBuf, 1);
}

void CCEnergyWavefunction::uhf_fock_build(double **fock_a, double **fock_b, double **D_a, double **D_b) {
    int lastbuf, p, q, r, s, pq, rs;
    double value;
    Value *valptr;
    Label *lblptr;
    struct iwlbuf InBuf;

    auto nso = moinfo_.nso;

    auto Dt = block_matrix(nso, nso);
    for (int p = 0; p < nso; p++)
        for (int q = 0; q < nso; q++) Dt[p][q] = D_a[p][q] + D_b[p][q];

    /* one-electron contributions */

    auto H = H_->to_block_matrix();

    for (int i = 0; i < nso; i++) {
        for (int j = 0; j <= i; j++) {
            fock_a[i][j] = fock_a[j][i] = H[i][j];
            fock_b[i][j] = fock_b[j][i] = H[i][j];
        }
    }

    iwl_buf_init(&InBuf, PSIF_SO_TEI, 0.0, 1, 1);
    do {
        lastbuf = InBuf.lastbuf;
        lblptr = InBuf.labels;
        valptr = InBuf.values;

        for (int idx = 4 * InBuf.idx; InBuf.idx < InBuf.inbuf; InBuf.idx++) {
            p = std::abs((int)lblptr[idx++]);
            q = (int)lblptr[idx++];
            r = (int)lblptr[idx++];
            s = (int)lblptr[idx++];
            value = (double)valptr[InBuf.idx];

            pq = INDEX(p, q);
            rs = INDEX(r, s);

            /*
            outfile->Printf( "%d %d %d %d [%d] [%d] %20.15f\n", p, q, r, s, pq, rs, value);
            */

            /* (pq|rs) */
            fock_a[p][q] += Dt[r][s] * value;
            fock_a[p][r] -= D_a[q][s] * value;
            fock_b[p][q] += Dt[r][s] * value;
            fock_b[p][r] -= D_b[q][s] * value;

            if (p != q && r != s && pq != rs) {
                /* (pq|sr) */
                fock_a[p][q] += Dt[s][r] * value;
                fock_a[p][s] -= D_a[q][r] * value;
                fock_b[p][q] += Dt[s][r] * value;
                fock_b[p][s] -= D_b[q][r] * value;

                /* (qp|rs) */
                fock_a[q][p] += Dt[r][s] * value;
                fock_a[q][r] -= D_a[p][s] * value;
                fock_b[q][p] += Dt[r][s] * value;
                fock_b[q][r] -= D_b[p][s] * value;

                /* (qp|sr) */
                fock_a[q][p] += Dt[s][r] * value;
                fock_a[q][s] -= D_a[p][r] * value;
                fock_b[q][p] += Dt[s][r] * value;
                fock_b[q][s] -= D_b[p][r] * value;

                /* (rs|pq) */
                fock_a[r][s] += Dt[p][q] * value;
                fock_a[r][p] -= D_a[s][q] * value;
                fock_b[r][s] += Dt[p][q] * value;
                fock_b[r][p] -= D_b[s][q] * value;

                /* (rs|qp) */
                fock_a[r][s] += Dt[q][p] * value;
                fock_a[r][q] -= D_a[s][p] * value;
                fock_b[r][s] += Dt[q][p] * value;
                fock_b[r][q] -= D_b[s][p] * value;

                /* (sr|pq) */
                fock_a[s][r] += Dt[p][q] * value;
                fock_a[s][p] -= D_a[r][q] * value;
                fock_b[s][r] += Dt[p][q] * value;
                fock_b[s][p] -= D_b[r][q] * value;

                /* (sr|qp) */
                fock_a[s][r] += Dt[q][p] * value;
                fock_a[s][q] -= D_a[r][p] * value;
                fock_b[s][r] += Dt[q][p] * value;
                fock_b[s][q] -= D_b[r][p] * value;
            } else if (p != q && r != s && pq == rs) {
                /* (pq|sr) */
                fock_a[p][q] += Dt[s][r] * value;
                fock_a[p][s] -= D_a[q][r] * value;
                fock_b[p][q] += Dt[s][r] * value;
                fock_b[p][s] -= D_b[q][r] * value;

                /* (qp|rs) */
                fock_a[q][p] += Dt[r][s] * value;
                fock_a[q][r] -= D_a[p][s] * value;
                fock_b[q][p] += Dt[r][s] * value;
                fock_b[q][r] -= D_b[p][s] * value;

                /* (qp|sr) */
                fock_a[q][p] += Dt[s][r] * value;
                fock_a[q][s] -= D_a[p][r] * value;
                fock_b[q][p] += Dt[s][r] * value;
                fock_b[q][s] -= D_b[p][r] * value;

            } else if (p != q && r == s) {
                /* (qp|rs) */
                fock_a[q][p] += Dt[r][s] * value;
                fock_a[q][r] -= D_a[p][s] * value;
                fock_b[q][p] += Dt[r][s] * value;
                fock_b[q][r] -= D_b[p][s] * value;

                /* (rs|pq) */
                fock_a[r][s] += Dt[p][q] * value;
                fock_a[r][p] -= D_a[s][q] * value;
                fock_b[r][s] += Dt[p][q] * value;
                fock_b[r][p] -= D_b[s][q] * value;

                /* (rs|qp) */
                fock_a[r][s] += Dt[q][p] * value;
                fock_a[r][q] -= D_a[s][p] * value;
                fock_b[r][s] += Dt[q][p] * value;
                fock_b[r][q] -= D_b[s][p] * value;

            } else if (p == q && r != s) {
                /* (pq|sr) */
                fock_a[p][q] += Dt[s][r] * value;
                fock_a[p][s] -= D_a[q][r] * value;
                fock_b[p][q] += Dt[s][r] * value;
                fock_b[p][s] -= D_b[q][r] * value;

                /* (rs|pq) */
                fock_a[r][s] += Dt[p][q] * value;
                fock_a[r][p] -= D_a[s][q] * value;
                fock_b[r][s] += Dt[p][q] * value;
                fock_b[r][p] -= D_b[s][q] * value;

                /* (sr|pq) */
                fock_a[s][r] += Dt[p][q] * value;
                fock_a[s][p] -= D_a[r][q] * value;
                fock_b[s][r] += Dt[p][q] * value;
                fock_b[s][p] -= D_b[r][q] * value;

            } else if (p == q && r == s && pq != rs) {
                /* (rs|pq) */
                fock_a[r][s] += Dt[p][q] * value;
                fock_a[r][p] -= D_a[s][q] * value;
                fock_b[r][s] += Dt[p][q] * value;
                fock_b[r][p] -= D_b[s][q] * value;
            }
        }

        if (!lastbuf) iwl_buf_fetch(&InBuf);

    } while (!lastbuf);
    iwl_buf_close(&InBuf, 1);

    free_block(Dt);
}
}
}  // namespace psi
