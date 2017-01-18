/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

/*! \file
    \ingroup ccresponse
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/psifiles.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

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
** Modified for ccresponse, TDC, 5/03
*/

void diis(int iter, const char *pert, int irrep, double omega)
{
  int nvector=8;  /* Number of error vectors to keep */
  int h, nirreps;
  int row, col, word, p, q;
  int diis_cycle;
  int vector_length=0;
  int errcod, *ipiv;
  dpdfile2 T1, T1a, T1b;
  dpdbuf4 T2, T2a, T2b, T2c;
  psio_address start, end, next;
  double **error;
  double **B, *C, **vector;
  double product, determinant, maximum;
  char lbl[32];

  nirreps = moinfo.nirreps;

  if(params.ref == 0) { /** RHF **/
    /* Compute the length of a single error vector */
    global_dpd_->file2_init(&T1, PSIF_CC_MISC, irrep, 0, 1, "XXX");
    global_dpd_->buf4_init(&T2, PSIF_CC_MISC, irrep, 0, 5, 0, 5, 0, "XXX");
    for(h=0; h < nirreps; h++) {
      vector_length += T1.params->rowtot[h] * T1.params->coltot[h^irrep];
      vector_length += T2.params->rowtot[h] * T2.params->coltot[h^irrep];
    }
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_close(&T2);

    /* Set the diis cycle value */
    diis_cycle = (iter-1) % nvector;

    /* Build the current error vector and dump it to disk */
    error = global_dpd_->dpd_block_matrix(1,vector_length);

    word=0;
    sprintf(lbl, "New X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&T1a, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->file2_mat_init(&T1a);
    global_dpd_->file2_mat_rd(&T1a);
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&T1b, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->file2_mat_init(&T1b);
    global_dpd_->file2_mat_rd(&T1b);
    for(h=0; h < nirreps; h++)
      for(row=0; row < T1a.params->rowtot[h]; row++)
        for(col=0; col < T1a.params->coltot[h^irrep]; col++)
          error[0][word++] = T1a.matrix[h][row][col] - T1b.matrix[h][row][col];
    global_dpd_->file2_mat_close(&T1a);
    global_dpd_->file2_close(&T1a);
    global_dpd_->file2_mat_close(&T1b);
    global_dpd_->file2_close(&T1b);

    sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&T2a, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&T2b, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&T2a, h);
      global_dpd_->buf4_mat_irrep_rd(&T2a, h);
      global_dpd_->buf4_mat_irrep_init(&T2b, h);
      global_dpd_->buf4_mat_irrep_rd(&T2b, h);
      for(row=0; row < T2a.params->rowtot[h]; row++)
        for(col=0; col < T2a.params->coltot[h^irrep]; col++)
          error[0][word++] = T2a.matrix[h][row][col] - T2b.matrix[h][row][col];
      global_dpd_->buf4_mat_irrep_close(&T2a, h);
      global_dpd_->buf4_mat_irrep_close(&T2b, h);
    }
    global_dpd_->buf4_close(&T2a);
    global_dpd_->buf4_close(&T2b);

    start = psio_get_address(PSIO_ZERO, diis_cycle*vector_length*sizeof(double));
    sprintf(lbl, "DIIS %s Error Vectors", pert);
    psio_write(PSIF_CC_DIIS_ERR, lbl , (char *) error[0],
               vector_length*sizeof(double), start, &end);

    /* Store the current amplitude vector on disk */
    word=0;

    sprintf(lbl, "New X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&T1a, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->file2_mat_init(&T1a);
    global_dpd_->file2_mat_rd(&T1a);
    for(h=0; h < nirreps; h++)
      for(row=0; row < T1a.params->rowtot[h]; row++)
        for(col=0; col < T1a.params->coltot[h^irrep]; col++)
          error[0][word++] = T1a.matrix[h][row][col];
    global_dpd_->file2_mat_close(&T1a);
    global_dpd_->file2_close(&T1a);

    sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&T2a, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&T2a, h);
      global_dpd_->buf4_mat_irrep_rd(&T2a, h);
      for(row=0; row < T2a.params->rowtot[h]; row++)
        for(col=0; col < T2a.params->coltot[h^irrep]; col++)
          error[0][word++] = T2a.matrix[h][row][col];
      global_dpd_->buf4_mat_irrep_close(&T2a, h);
    }
    global_dpd_->buf4_close(&T2a);

    start = psio_get_address(PSIO_ZERO, diis_cycle*vector_length*sizeof(double));
    sprintf(lbl, "DIIS %s Amplitude Vectors", pert);
    psio_write(PSIF_CC_DIIS_AMP, lbl , (char *) error[0],
               vector_length*sizeof(double), start, &end);

    /* If we haven't run through enough iterations, set the correct dimensions
       for the extrapolation */
    if(!(iter >= (nvector))) {
      if(iter < 2) { /* Leave if we can't extrapolate at all */
        global_dpd_->free_dpd_block(error, 1, vector_length);
        return;
      }
      nvector = iter;
    }

    /* Build B matrix of error vector products */
    vector = global_dpd_->dpd_block_matrix(2, vector_length);
    B = block_matrix(nvector+1,nvector+1);
    for(p=0; p < nvector; p++) {

      start = psio_get_address(PSIO_ZERO, p*vector_length*sizeof(double));

      sprintf(lbl, "DIIS %s Error Vectors", pert);
      psio_read(PSIF_CC_DIIS_ERR, lbl, (char *) vector[0],
                vector_length*sizeof(double), start, &end);

      dot_arr(vector[0], vector[0], vector_length, &product);

      B[p][p] = product;

      for(q=0; q < p; q++) {

        start = psio_get_address(PSIO_ZERO, q*vector_length*sizeof(double));

        sprintf(lbl, "DIIS %s Error Vectors", pert);
        psio_read(PSIF_CC_DIIS_ERR, lbl, (char *) vector[1],
                  vector_length*sizeof(double), start, &end);

        dot_arr(vector[1], vector[0], vector_length, &product);

        B[p][q] = B[q][p] = product;
      }
    }
    global_dpd_->free_dpd_block(vector, 2, vector_length);

    for(p=0; p < nvector; p++) {
      B[p][nvector] = -1;
      B[nvector][p] = -1;
    }

    B[nvector][nvector] = 0;

    /* Find the maximum value in B and scale all its elements */
    maximum = fabs(B[0][0]);
    for(p=0; p < nvector; p++)
      for(q=0; q < nvector; q++)
        if(fabs(B[p][q]) > maximum) maximum = fabs(B[p][q]);

    for(p=0; p < nvector; p++)
      for(q=0; q < nvector; q++)
        B[p][q] /= maximum;

    /* Build the constant vector */
    C = init_array(nvector+1);
    C[nvector] = -1;

    /* Solve the linear equations */
    ipiv = init_int_array(nvector+1);

    errcod = C_DGESV(nvector+1, 1, &(B[0][0]), nvector+1, &(ipiv[0]), &(C[0]), nvector+1);
    if(errcod) {
      throw PsiException("Error in DGESV return in diis",__FILE__,__LINE__);
    }

    /* Build a new amplitude vector from the old ones */
    vector = global_dpd_->dpd_block_matrix(1, vector_length);
    for(p=0; p < vector_length; p++) error[0][p] = 0.0;
    for(p=0; p < nvector; p++) {

      start = psio_get_address(PSIO_ZERO, p*vector_length*sizeof(double));

      sprintf(lbl, "DIIS %s Amplitude Vectors", pert);
      psio_read(PSIF_CC_DIIS_AMP, lbl, (char *) vector[0],
                vector_length*sizeof(double), start, &end);

      for(q=0; q < vector_length; q++)
        error[0][q] += C[p] * vector[0][q];

    }
    global_dpd_->free_dpd_block(vector, 1, vector_length);

    /* Now place these elements into the DPD amplitude arrays */
    word=0;
    sprintf(lbl, "New X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&T1a, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->file2_mat_init(&T1a);
    for(h=0; h < nirreps; h++)
      for(row=0; row < T1a.params->rowtot[h]; row++)
        for(col=0; col < T1a.params->coltot[h^irrep]; col++)
          T1a.matrix[h][row][col] = error[0][word++];
    global_dpd_->file2_mat_wrt(&T1a);
    global_dpd_->file2_mat_close(&T1a);
    global_dpd_->file2_close(&T1a);

    sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&T2a, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    for(h=0; h < nirreps; h++) {
      global_dpd_->buf4_mat_irrep_init(&T2a, h);
      for(row=0; row < T2a.params->rowtot[h]; row++)
        for(col=0; col < T2a.params->coltot[h^irrep]; col++)
          T2a.matrix[h][row][col] = error[0][word++];
      global_dpd_->buf4_mat_irrep_wrt(&T2a, h);
      global_dpd_->buf4_mat_irrep_close(&T2a, h);
    }
    global_dpd_->buf4_close(&T2a);

    /* Release memory and return */
    /*    free_matrix(vector, nvector); */
    free_block(B);
    free(C);
    free(ipiv);
    global_dpd_->free_dpd_block(error, 1, vector_length);
  }

  return;
}

}} // namespace psi::ccresponse
