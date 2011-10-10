/*! \file
    \ingroup ccresponse
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <exception.h>
#include <psifiles.h>
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
    dpd_file2_init(&T1, CC_MISC, irrep, 0, 1, "XXX");
    dpd_buf4_init(&T2, CC_MISC, irrep, 0, 5, 0, 5, 0, "XXX");
    for(h=0; h < nirreps; h++) {
      vector_length += T1.params->rowtot[h] * T1.params->coltot[h^irrep];
      vector_length += T2.params->rowtot[h] * T2.params->coltot[h^irrep];
    }
    dpd_file2_close(&T1);
    dpd_buf4_close(&T2);

    /* Set the diis cycle value */
    diis_cycle = (iter-1) % nvector;

    /* Build the current error vector and dump it to disk */
    error = dpd_block_matrix(1,vector_length);

    word=0;
    sprintf(lbl, "New X_%s_IA (%5.3f)", pert, omega);
    dpd_file2_init(&T1a, CC_OEI, irrep, 0, 1, lbl);
    dpd_file2_mat_init(&T1a);
    dpd_file2_mat_rd(&T1a);
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    dpd_file2_init(&T1b, CC_OEI, irrep, 0, 1, lbl);
    dpd_file2_mat_init(&T1b);
    dpd_file2_mat_rd(&T1b);
    for(h=0; h < nirreps; h++)
      for(row=0; row < T1a.params->rowtot[h]; row++)
        for(col=0; col < T1a.params->coltot[h^irrep]; col++)
          error[0][word++] = T1a.matrix[h][row][col] - T1b.matrix[h][row][col];
    dpd_file2_mat_close(&T1a);
    dpd_file2_close(&T1a);
    dpd_file2_mat_close(&T1b);
    dpd_file2_close(&T1b);

    sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
    dpd_buf4_init(&T2a, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    dpd_buf4_init(&T2b, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&T2a, h);
      dpd_buf4_mat_irrep_rd(&T2a, h);
      dpd_buf4_mat_irrep_init(&T2b, h);
      dpd_buf4_mat_irrep_rd(&T2b, h);
      for(row=0; row < T2a.params->rowtot[h]; row++)
        for(col=0; col < T2a.params->coltot[h^irrep]; col++)
          error[0][word++] = T2a.matrix[h][row][col] - T2b.matrix[h][row][col];
      dpd_buf4_mat_irrep_close(&T2a, h);
      dpd_buf4_mat_irrep_close(&T2b, h);
    }
    dpd_buf4_close(&T2a);
    dpd_buf4_close(&T2b);

    start = psio_get_address(PSIO_ZERO, diis_cycle*vector_length*sizeof(double));
    sprintf(lbl, "DIIS %s Error Vectors", pert);
    psio_write(CC_DIIS_ERR, lbl , (char *) error[0],
               vector_length*sizeof(double), start, &end);

    /* Store the current amplitude vector on disk */
    word=0;

    sprintf(lbl, "New X_%s_IA (%5.3f)", pert, omega);
    dpd_file2_init(&T1a, CC_OEI, irrep, 0, 1, lbl);
    dpd_file2_mat_init(&T1a);
    dpd_file2_mat_rd(&T1a);
    for(h=0; h < nirreps; h++)
      for(row=0; row < T1a.params->rowtot[h]; row++)
        for(col=0; col < T1a.params->coltot[h^irrep]; col++)
          error[0][word++] = T1a.matrix[h][row][col];
    dpd_file2_mat_close(&T1a);
    dpd_file2_close(&T1a);

    sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
    dpd_buf4_init(&T2a, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&T2a, h);
      dpd_buf4_mat_irrep_rd(&T2a, h);
      for(row=0; row < T2a.params->rowtot[h]; row++)
        for(col=0; col < T2a.params->coltot[h^irrep]; col++)
          error[0][word++] = T2a.matrix[h][row][col];
      dpd_buf4_mat_irrep_close(&T2a, h);
    }
    dpd_buf4_close(&T2a);

    start = psio_get_address(PSIO_ZERO, diis_cycle*vector_length*sizeof(double));
    sprintf(lbl, "DIIS %s Amplitude Vectors", pert);
    psio_write(CC_DIIS_AMP, lbl , (char *) error[0],
               vector_length*sizeof(double), start, &end);

    /* If we haven't run through enough iterations, set the correct dimensions
       for the extrapolation */
    if(!(iter >= (nvector))) {
      if(iter < 2) { /* Leave if we can't extrapolate at all */
        dpd_free_block(error, 1, vector_length);
        return;
      }
      nvector = iter;
    }

    /* Build B matrix of error vector products */
    vector = dpd_block_matrix(2, vector_length);
    B = block_matrix(nvector+1,nvector+1);
    for(p=0; p < nvector; p++) {

      start = psio_get_address(PSIO_ZERO, p*vector_length*sizeof(double));

      sprintf(lbl, "DIIS %s Error Vectors", pert);
      psio_read(CC_DIIS_ERR, lbl, (char *) vector[0],
                vector_length*sizeof(double), start, &end);

      dot_arr(vector[0], vector[0], vector_length, &product);

      B[p][p] = product;

      for(q=0; q < p; q++) {

        start = psio_get_address(PSIO_ZERO, q*vector_length*sizeof(double));

        sprintf(lbl, "DIIS %s Error Vectors", pert);
        psio_read(CC_DIIS_ERR, lbl, (char *) vector[1],
                  vector_length*sizeof(double), start, &end);

        dot_arr(vector[1], vector[0], vector_length, &product);

        B[p][q] = B[q][p] = product;
      }
    }
    dpd_free_block(vector, 2, vector_length);

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
    vector = dpd_block_matrix(1, vector_length);
    for(p=0; p < vector_length; p++) error[0][p] = 0.0;
    for(p=0; p < nvector; p++) {

      start = psio_get_address(PSIO_ZERO, p*vector_length*sizeof(double));

      sprintf(lbl, "DIIS %s Amplitude Vectors", pert);
      psio_read(CC_DIIS_AMP, lbl, (char *) vector[0],
                vector_length*sizeof(double), start, &end);

      for(q=0; q < vector_length; q++)
        error[0][q] += C[p] * vector[0][q];

    }
    dpd_free_block(vector, 1, vector_length);

    /* Now place these elements into the DPD amplitude arrays */
    word=0;
    sprintf(lbl, "New X_%s_IA (%5.3f)", pert, omega);
    dpd_file2_init(&T1a, CC_OEI, irrep, 0, 1, lbl);
    dpd_file2_mat_init(&T1a);
    for(h=0; h < nirreps; h++)
      for(row=0; row < T1a.params->rowtot[h]; row++)
        for(col=0; col < T1a.params->coltot[h^irrep]; col++)
          T1a.matrix[h][row][col] = error[0][word++];
    dpd_file2_mat_wrt(&T1a);
    dpd_file2_mat_close(&T1a);
    dpd_file2_close(&T1a);

    sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
    dpd_buf4_init(&T2a, CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&T2a, h);
      for(row=0; row < T2a.params->rowtot[h]; row++)
        for(col=0; col < T2a.params->coltot[h^irrep]; col++)
          T2a.matrix[h][row][col] = error[0][word++];
      dpd_buf4_mat_irrep_wrt(&T2a, h);
      dpd_buf4_mat_irrep_close(&T2a, h);
    }
    dpd_buf4_close(&T2a);

    /* Release memory and return */
    /*    free_matrix(vector, nvector); */
    free_block(B);
    free(C);
    free(ipiv);
    dpd_free_block(error, 1, vector_length);
  }

  return;
}

}} // namespace psi::ccresponse
