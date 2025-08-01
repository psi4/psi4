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
#include "psi4/libpsio/psio.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"
#include "MOInfo.h"
#include "Params.h"
#include "psi4/cc/ccwave.h"

namespace psi {
namespace ccenergy {

/*
** DIIS: Direct inversion in the iterative subspace routine to
** accelerate convergence of the CCSD amplitude equations.
**
** Substantially improved efficiency of this routine:
** (1) Keeping at most two error vectors in core at once.
** (2) Limiting direct product (overlap) calculation to unique pairs.
** (3) Using LAPACK's linear equation solver DGESV instead of flin.
**
** These improvements have been applied only to RHF cases so far.
**
** -TDC  12/22/01
*/

void CCEnergyWavefunction::diis_RHF(int iter) {
    int nvector = 8; /* Number of error vectors to keep */
    dpdfile2 T1, T1a, T1b;
    dpdbuf4 T2, T2a, T2b;
    psio_address start, end;
    double **B, *C, **vector;
    double product, maximum;

    auto nirreps = moinfo_.nirreps;

    /* Compute the length of a single error vector */
    global_dpd_->file2_init(&T1, PSIF_CC_TAMPS, 0, 0, 1, "tIA");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    auto vector_length = 0;
    for (int h = 0; h < nirreps; h++) {
        vector_length += T1.params->rowtot[h] * T1.params->coltot[h];
        vector_length += T2.params->rowtot[h] * T2.params->coltot[h];
    }
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&T2);

    /* Set the diis cycle value */
    auto diis_cycle = (iter - 1) % nvector;

    /* Build the current error vector and dump it to disk */
    auto error = global_dpd_->dpd_block_matrix(1, vector_length);

    auto word = 0;
    global_dpd_->file2_init(&T1a, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    global_dpd_->file2_mat_init(&T1a);
    global_dpd_->file2_mat_rd(&T1a);
    global_dpd_->file2_init(&T1b, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_mat_init(&T1b);
    global_dpd_->file2_mat_rd(&T1b);
    for (int h = 0; h < nirreps; h++)
        for (int row = 0; row < T1a.params->rowtot[h]; row++)
            for (int col = 0; col < T1a.params->coltot[h]; col++)
                error[0][word++] = T1a.matrix[h][row][col] - T1b.matrix[h][row][col];
    global_dpd_->file2_mat_close(&T1a);
    global_dpd_->file2_close(&T1a);
    global_dpd_->file2_mat_close(&T1b);
    global_dpd_->file2_close(&T1b);

    auto t1_word = word;
    global_dpd_->buf4_init(&T2a, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    for (int h = 0; h < nirreps; h++) {
        global_dpd_->buf4_mat_irrep_init(&T2a, h);
        global_dpd_->buf4_mat_irrep_rd(&T2a, h);
        for (int row = 0; row < T2a.params->rowtot[h]; row++)
            for (int col = 0; col < T2a.params->coltot[h]; col++) error[0][word++] = T2a.matrix[h][row][col];
        global_dpd_->buf4_mat_irrep_close(&T2a, h);
    }
    global_dpd_->buf4_close(&T2a);

    word = t1_word;
    global_dpd_->buf4_init(&T2b, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    for (int h = 0; h < nirreps; h++) {
        global_dpd_->buf4_mat_irrep_init(&T2b, h);
        global_dpd_->buf4_mat_irrep_rd(&T2b, h);
        for (int row = 0; row < T2b.params->rowtot[h]; row++)
            for (int col = 0; col < T2b.params->coltot[h]; col++) error[0][word++] -= T2b.matrix[h][row][col];
        global_dpd_->buf4_mat_irrep_close(&T2b, h);
    }
    global_dpd_->buf4_close(&T2b);

    start = psio_get_address(PSIO_ZERO, sizeof(double) * diis_cycle * vector_length);
    psio_write(PSIF_CC_DIIS_ERR, "DIIS Error Vectors", (char *)error[0], vector_length * sizeof(double), start, &end);

    /* Store the current amplitude vector on disk */
    word = 0;

    global_dpd_->file2_init(&T1a, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    global_dpd_->file2_mat_init(&T1a);
    global_dpd_->file2_mat_rd(&T1a);
    for (int h = 0; h < nirreps; h++)
        for (int row = 0; row < T1a.params->rowtot[h]; row++)
            for (int col = 0; col < T1a.params->coltot[h]; col++) error[0][word++] = T1a.matrix[h][row][col];
    global_dpd_->file2_mat_close(&T1a);
    global_dpd_->file2_close(&T1a);

    global_dpd_->buf4_init(&T2a, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    for (int h = 0; h < nirreps; h++) {
        global_dpd_->buf4_mat_irrep_init(&T2a, h);
        global_dpd_->buf4_mat_irrep_rd(&T2a, h);
        for (int row = 0; row < T2a.params->rowtot[h]; row++)
            for (int col = 0; col < T2a.params->coltot[h]; col++) error[0][word++] = T2a.matrix[h][row][col];
        global_dpd_->buf4_mat_irrep_close(&T2a, h);
    }
    global_dpd_->buf4_close(&T2a);

    start = psio_get_address(PSIO_ZERO, sizeof(double) * diis_cycle * vector_length);
    psio_write(PSIF_CC_DIIS_AMP, "DIIS Amplitude Vectors", (char *)error[0], vector_length * sizeof(double), start,
               &end);

    /* If we haven't run through enough iterations, set the correct dimensions
       for the extrapolation */
    if (!(iter >= (nvector))) {
        if (iter < 2) { /* Leave if we can't extrapolate at all */
            global_dpd_->free_dpd_block(error, 1, vector_length);
            return;
        }
        nvector = iter;
    }

    /* Build B matrix of error vector products */
    vector = global_dpd_->dpd_block_matrix(2, vector_length);
    B = block_matrix(nvector + 1, nvector + 1);
    for (int p = 0; p < nvector; p++) {
        start = psio_get_address(PSIO_ZERO, sizeof(double) * p * vector_length);

        psio_read(PSIF_CC_DIIS_ERR, "DIIS Error Vectors", (char *)vector[0], vector_length * sizeof(double), start,
                  &end);

        product = C_DDOT(vector_length, vector[0], 1, vector[0], 1);
        // dot_arr(vector[0], vector[0], vector_length, &product);

        B[p][p] = product;

        for (int q = 0; q < p; q++) {
            start = psio_get_address(PSIO_ZERO, sizeof(double) * q * vector_length);

            psio_read(PSIF_CC_DIIS_ERR, "DIIS Error Vectors", (char *)vector[1], vector_length * sizeof(double), start,
                      &end);

            // dot_arr(vector[1], vector[0], vector_length, &product);
            product = C_DDOT(vector_length, vector[1], 1, vector[0], 1);

            B[p][q] = B[q][p] = product;
        }
    }
    global_dpd_->free_dpd_block(vector, 2, vector_length);

    for (int p = 0; p < nvector; p++) {
        B[p][nvector] = -1;
        B[nvector][p] = -1;
    }

    B[nvector][nvector] = 0;

    /* Find the maximum value in B and scale all its elements */
    maximum = std::fabs(B[0][0]);
    for (int p = 0; p < nvector; p++)
        for (int q = 0; q < nvector; q++)
            if (std::fabs(B[p][q]) > maximum) maximum = std::fabs(B[p][q]);

    for (int p = 0; p < nvector; p++)
        for (int q = 0; q < nvector; q++) B[p][q] /= maximum;

    /* Build the constant vector */
    C = init_array(nvector + 1);
    C[nvector] = -1;

    /* Solve the linear equations */
    diis_invert_B(B, C, nvector + 1, 1.0E-12);

    /* Build a new amplitude vector from the old ones */
    vector = global_dpd_->dpd_block_matrix(1, vector_length);
    for (int p = 0; p < vector_length; p++) error[0][p] = 0.0;
    for (int p = 0; p < nvector; p++) {
        start = psio_get_address(PSIO_ZERO, sizeof(double) * p * vector_length);

        psio_read(PSIF_CC_DIIS_AMP, "DIIS Amplitude Vectors", (char *)vector[0], vector_length * sizeof(double), start,
                  &end);

        for (int q = 0; q < vector_length; q++) error[0][q] += C[p] * vector[0][q];
    }
    global_dpd_->free_dpd_block(vector, 1, vector_length);

    /* Now place these elements into the DPD amplitude arrays */
    word = 0;
    global_dpd_->file2_init(&T1a, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    global_dpd_->file2_mat_init(&T1a);
    for (int h = 0; h < nirreps; h++)
        for (int row = 0; row < T1a.params->rowtot[h]; row++)
            for (int col = 0; col < T1a.params->coltot[h]; col++) T1a.matrix[h][row][col] = error[0][word++];
    global_dpd_->file2_mat_wrt(&T1a);
    global_dpd_->file2_mat_close(&T1a);
    global_dpd_->file2_close(&T1a);

    global_dpd_->buf4_init(&T2a, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    for (int h = 0; h < nirreps; h++) {
        global_dpd_->buf4_mat_irrep_init(&T2a, h);
        for (int row = 0; row < T2a.params->rowtot[h]; row++)
            for (int col = 0; col < T2a.params->coltot[h]; col++) T2a.matrix[h][row][col] = error[0][word++];
        global_dpd_->buf4_mat_irrep_wrt(&T2a, h);
        global_dpd_->buf4_mat_irrep_close(&T2a, h);
    }
    global_dpd_->buf4_close(&T2a);

    /* Release memory and return */
    /*    free_matrix(vector, nvector); */
    free_block(B);
    free(C);
    global_dpd_->free_dpd_block(error, 1, vector_length);
}
}  // namespace ccenergy
}  // namespace psi
