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

/*! \file
    \ingroup CCEOM
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "psi4/libqt/qt.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace cceom {

struct onestack {
    int i;
    int a;
    double value;
};

void stack_insert(struct onestack *stack, double value, int i, int a, int level, int stacklen);

void local_guess() {
    int nso, nocc, nvir, nroot, i, ii, ij, a, m;
    double *T1bar, *T1tilde;
    double fii, value, norm;
    psio_address next;

    struct onestack *stack;

    char lbl[32];
    dpdfile2 CME;

    nso = local.nso;
    nocc = local.nocc;
    nvir = local.nvir;

    local.pairdom_len = init_int_array(nocc * nocc);
    local.pairdom_nrlen = init_int_array(nocc * nocc);
    local.eps_occ = init_array(nocc);
    psio_read_entry(PSIF_CC_INFO, "Local Pair Domain Length", (char *)local.pairdom_len, nocc * nocc * sizeof(int));
    psio_read_entry(PSIF_CC_INFO, "Local Pair Domain Length (Non-redundant basis)", (char *)local.pairdom_nrlen,
                    nocc * nocc * sizeof(int));
    psio_read_entry(PSIF_CC_INFO, "Local Occupied Orbital Energies", (char *)local.eps_occ, nocc * sizeof(double));

    local.W = (double ***)malloc(nocc * nocc * sizeof(double **));
    local.V = (double ***)malloc(nocc * nocc * sizeof(double **));
    local.eps_vir = (double **)malloc(nocc * nocc * sizeof(double *));
    next = PSIO_ZERO;
    for (ij = 0; ij < nocc * nocc; ij++) {
        local.eps_vir[ij] = init_array(local.pairdom_nrlen[ij]);
        psio_read(PSIF_CC_INFO, "Local Virtual Orbital Energies", (char *)local.eps_vir[ij],
                  local.pairdom_nrlen[ij] * sizeof(double), next, &next);
    }
    next = PSIO_ZERO;
    for (ij = 0; ij < nocc * nocc; ij++) {
        local.V[ij] = block_matrix(nvir, local.pairdom_len[ij]);
        psio_read(PSIF_CC_INFO, "Local Residual Vector (V)", (char *)local.V[ij][0],
                  nvir * local.pairdom_len[ij] * sizeof(double), next, &next);
    }
    next = PSIO_ZERO;
    for (ij = 0; ij < nocc * nocc; ij++) {
        local.W[ij] = block_matrix(local.pairdom_len[ij], local.pairdom_nrlen[ij]);
        psio_read(PSIF_CC_INFO, "Local Transformation Matrix (W)", (char *)local.W[ij][0],
                  local.pairdom_len[ij] * local.pairdom_nrlen[ij] * sizeof(double), next, &next);
    }

    nroot = eom_params.states_per_irrep[0]; /* only C1 allowed */

    stack = (struct onestack *)malloc(nroot * sizeof(struct onestack));
    for (m = 0; m < nroot; m++) {
        stack[m].i = -1;
        stack[m].a = -1;
        stack[m].value = 1e12;
    }

    /* find the nroot lowest excitations in the non-redunant, orthogonal (bar) space */
    for (i = 0; i < nocc; i++) {
        ii = i * nocc + i;
        fii = local.eps_occ[i];
        for (a = 0; a < local.pairdom_nrlen[ii]; a++) {
            value = local.eps_vir[ii][a] - fii;
            for (m = 0; m < nroot; m++) {
                if ((std::fabs(value) < std::fabs(stack[m].value))) {
                    stack_insert(stack, value, i, a, m, nroot);
                    break;
                }
            }
        }
    }

    T1bar = init_array(nso);
    T1tilde = init_array(nso);

    outfile->Printf("\n\tTransitions for local guesses:\n");
    outfile->Printf("\t------------------------------\n");
    for (m = 0; m < nroot; m++) {
        outfile->Printf("\t%3d %3d %14.10f\n", stack[m].i, stack[m].a, stack[m].value);

        memset((void *)T1bar, 0, nso * sizeof(double));
        memset((void *)T1tilde, 0, nso * sizeof(double));

        i = stack[m].i;
        ii = i * nocc + i;

        /* Unit guess vector */
        T1bar[stack[m].a] = 1.0;

        /* Transform this to the canonical MO basis */
        sprintf(lbl, "%s %d", "CME", m);
        global_dpd_->file2_init(&CME, PSIF_EOM_CME, 0, 0, 1, lbl);
        global_dpd_->file2_mat_init(&CME);

        C_DGEMV('n', local.pairdom_len[ii], local.pairdom_nrlen[ii], 1.0, &(local.W[ii][0][0]), local.pairdom_nrlen[ii],
                &(T1bar[0]), 1, 0.0, &(T1tilde[0]), 1);
        C_DGEMV('n', nvir, local.pairdom_len[ii], 1.0, &(local.V[ii][0][0]), local.pairdom_len[ii], &(T1tilde[0]), 1,
                0.0, &(CME.matrix[0][i][0]), 1);

        /* normalize this guess in the MO basis */
        norm = 0.0;
        for (a = 0; a < nvir; a++) {
            norm += CME.matrix[0][i][a] * CME.matrix[0][i][a];
        }
        norm = sqrt(2.0 * norm);
        outfile->Printf("Norm of guess vector %d = %20.14f\n", m, norm);
        for (a = 0; a < nvir; a++) {
            CME.matrix[0][i][a] *= 1.0 / norm;
        }

        global_dpd_->file2_mat_wrt(&CME);
        global_dpd_->file2_mat_close(&CME);

        global_dpd_->file2_close(&CME);
    }

    outfile->Printf("\n");

    free(T1bar);
    free(T1tilde);

    free(stack);

    eom_params.cs_per_irrep[0] = nroot;

    /* Free Local Memory */
    for (i = 0; i < nocc * nocc; i++) {
        free_block(local.W[i]);
        free_block(local.V[i]);
        free(local.eps_vir[i]);
    }
    free(local.W);
    free(local.V);
    free(local.eps_vir);

    free(local.eps_occ);
    free(local.pairdom_len);
    free(local.pairdom_nrlen);
}

void stack_insert(struct onestack *stack, double value, int i, int a, int level, int stacklen) {
    int l;
    struct onestack temp;

    temp = stack[level];

    stack[level].value = value;
    stack[level].i = i;
    stack[level].a = a;

    value = temp.value;
    i = temp.i;
    a = temp.a;

    for (l = level; l < stacklen - 1; l++) {
        temp = stack[l + 1];

        stack[l + 1].value = value;
        stack[l + 1].i = i;
        stack[l + 1].a = a;

        value = temp.value;
        i = temp.i;
        a = temp.a;
    }
}

}  // namespace cceom
}  // namespace psi
