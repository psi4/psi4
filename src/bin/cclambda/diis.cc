/*! \file
    \ingroup CCLAMBDA
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <exception.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

/*
** DIIS: Direct inversion in the iterative subspace routine to
** accelerate convergence of the CCSD amplitude equations.
**
** SubstanLially improved efficiency of this routine:
** (1) Keeping at most two error vectors in core at once.
** (2) Limiting direct product (overlap) calculation to unique pairs.
** (3) Using LAPACK's linear equation solver DGESV instead of flin.
**
** These improvements have been applied only to RHF cases so far.
**
** -TDC  12/22/01
*/

void diis(int iter, int L_irr)
{
  int nvector=8;  /* Number of error vectors to keep */
  int h, nirreps;
  int row, col, word, p, q, i;
  int diis_cycle;
  int vector_length=0;
  int errcod, *ipiv;
  dpdfile2 L1, L1a, L1b;
  dpdbuf4 L2, L2a, L2b, L2c;
  psio_address start, end, next;
  double **error;
  double **B, *C, **vector;
  double product, determinant, maximum;

  nirreps = moinfo.nirreps;

  if(params.ref == 0) { /** RHF **/
    /* Compute the length of a single error vector */
    dpd_file2_init(&L1, CC_LAMBDA, L_irr, 0, 1, "LIA");
    dpd_buf4_init(&L2, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    for(h=0; h < nirreps; h++) {
      vector_length += L1.params->rowtot[h] * L1.params->coltot[h^L_irr];
      vector_length += L2.params->rowtot[h] * L2.params->coltot[h^L_irr];
    }
    dpd_file2_close(&L1);
    dpd_buf4_close(&L2);

    /* Set the diis cycle value */
    diis_cycle = (iter-1) % nvector;

    /* Build the current error vector and dump it to disk */
    error = dpd_block_matrix(1,vector_length);

    word=0;
    dpd_file2_init(&L1a, CC_LAMBDA, L_irr, 0, 1, "New LIA");
    /*dpd_file2_print(&L1a,outfile);*/
    dpd_file2_mat_init(&L1a);
    dpd_file2_mat_rd(&L1a);
    dpd_file2_init(&L1b, CC_LAMBDA, L_irr, 0, 1, "LIA");
    /*dpd_file2_print(&L1b,outfile);*/
    dpd_file2_mat_init(&L1b);
    dpd_file2_mat_rd(&L1b);
    for(h=0; h < nirreps; h++)
      for(row=0; row < L1a.params->rowtot[h]; row++)
        for(col=0; col < L1a.params->coltot[h^L_irr]; col++) {
          error[0][word++] = L1a.matrix[h][row][col] - L1b.matrix[h][row][col];
    }

    dpd_file2_mat_close(&L1a);
    dpd_file2_close(&L1a);
    dpd_file2_mat_close(&L1b);
    dpd_file2_close(&L1b);

    dpd_buf4_init(&L2a, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    /*dpd_buf4_print(&L2a,outfile,1);*/
    dpd_buf4_init(&L2b, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    /*dpd_buf4_print(&L2b,outfile,1);*/
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&L2a, h);
      dpd_buf4_mat_irrep_rd(&L2a, h);
      dpd_buf4_mat_irrep_init(&L2b, h);
      dpd_buf4_mat_irrep_rd(&L2b, h);
      for(row=0; row < L2a.params->rowtot[h]; row++)
        for(col=0; col < L2a.params->coltot[h^L_irr]; col++) {
          error[0][word++] = L2a.matrix[h][row][col] - L2b.matrix[h][row][col];
/* fprintf(outfile,"%15.10lf\n", error[0][word-1]); */
    }
      dpd_buf4_mat_irrep_close(&L2a, h);
      dpd_buf4_mat_irrep_close(&L2b, h);
    }
    dpd_buf4_close(&L2a);
    dpd_buf4_close(&L2b);

    start = psio_get_address(PSIO_ZERO, diis_cycle*vector_length*sizeof(double));
    psio_write(CC_DIIS_ERR, "DIIS Error Vectors" , (char *) error[0],
               vector_length*sizeof(double), start, &end);


    /* Store the current amplitude vector on disk */
    word=0;

    dpd_file2_init(&L1a, CC_LAMBDA, L_irr, 0, 1, "New LIA");
    dpd_file2_mat_init(&L1a);
    dpd_file2_mat_rd(&L1a);
    for(h=0; h < nirreps; h++)
      for(row=0; row < L1a.params->rowtot[h]; row++)
        for(col=0; col < L1a.params->coltot[h^L_irr]; col++)
          error[0][word++] = L1a.matrix[h][row][col];
    dpd_file2_mat_close(&L1a);
    dpd_file2_close(&L1a);

    dpd_buf4_init(&L2a, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&L2a, h);
      dpd_buf4_mat_irrep_rd(&L2a, h);
      for(row=0; row < L2a.params->rowtot[h]; row++)
        for(col=0; col < L2a.params->coltot[h^L_irr]; col++)
          error[0][word++] = L2a.matrix[h][row][col];
      dpd_buf4_mat_irrep_close(&L2a, h);
    }
    dpd_buf4_close(&L2a);

    start = psio_get_address(PSIO_ZERO, diis_cycle*vector_length*sizeof(double));
    psio_write(CC_DIIS_AMP, "DIIS Amplitude Vectors" , (char *) error[0],
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

      psio_read(CC_DIIS_ERR, "DIIS Error Vectors", (char *) vector[0],
                vector_length*sizeof(double), start, &end);

      /*
      for(i=0; i < vector_length; i++)
        fprintf(outfile,"E[%d][%d] = %20.15lf\n",p,i,vector[0][i]);
      */

      dot_arr(vector[0], vector[0], vector_length, &product);

      B[p][p] = product;

      for(q=0; q < p; q++) {

        start = psio_get_address(PSIO_ZERO, q*vector_length*sizeof(double));

        psio_read(CC_DIIS_ERR, "DIIS Error Vectors", (char *) vector[1],
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

    /*
    fprintf(outfile,"\nDIIS B:\n");
    print_mat(B,nvector,nvector,outfile);
    */

    /* Build the constant vector */
    C = init_array(nvector+1);
    C[nvector] = -1;

    /* Solve the linear equations */
    ipiv = init_int_array(nvector+1);

    errcod = C_DGESV(nvector+1, 1, &(B[0][0]), nvector+1, &(ipiv[0]), &(C[0]), nvector+1);
    if(errcod) {
      fprintf(outfile, "\nError in DGESV return in diis.\n");
      throw PsiException("cclambda: error", __FILE__, __LINE__);
    }

    /* Build a new amplitude vector from the old ones */
    vector = dpd_block_matrix(1, vector_length);
    for(p=0; p < vector_length; p++) error[0][p] = 0.0;
    for(p=0; p < nvector; p++) {
      /*fprintf(outfile,"C[%d] = %20.15lf\n",p,C[p]);*/

      start = psio_get_address(PSIO_ZERO, p*vector_length*sizeof(double));

      psio_read(CC_DIIS_AMP, "DIIS Amplitude Vectors", (char *) vector[0],
                vector_length*sizeof(double), start, &end);

      for(q=0; q < vector_length; q++)
        error[0][q] += C[p] * vector[0][q];

    }
    dpd_free_block(vector, 1, vector_length);

    /* Now place these elements into the DPD amplitude arrays */
    word=0;
    dpd_file2_init(&L1a, CC_LAMBDA, L_irr, 0, 1, "New LIA");
    dpd_file2_mat_init(&L1a);
    for(h=0; h < nirreps; h++)
      for(row=0; row < L1a.params->rowtot[h]; row++)
        for(col=0; col < L1a.params->coltot[h^L_irr]; col++)
          L1a.matrix[h][row][col] = error[0][word++];
    dpd_file2_mat_wrt(&L1a);
    dpd_file2_mat_close(&L1a);
    dpd_file2_copy(&L1a, CC_LAMBDA, "New Lia"); /* to be removed after spin-adaptation */
    /*dpd_file2_print(&L1a,outfile);*/
    dpd_file2_close(&L1a);

    dpd_buf4_init(&L2a, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&L2a, h);
      for(row=0; row < L2a.params->rowtot[h]; row++)
        for(col=0; col < L2a.params->coltot[h^L_irr]; col++)
          L2a.matrix[h][row][col] = error[0][word++];
      dpd_buf4_mat_irrep_wrt(&L2a, h);
      dpd_buf4_mat_irrep_close(&L2a, h);
    }
    dpd_buf4_close(&L2a);

    /* to be removed after spin-adaptation */
    dpd_buf4_init(&L2a, CC_LAMBDA, L_irr, 2, 7, 0, 5, 1, "New LIjAb");
    dpd_buf4_copy(&L2a, CC_LAMBDA, "New LIJAB");
    dpd_buf4_copy(&L2a, CC_LAMBDA, "New Lijab");
    dpd_buf4_close(&L2a);

    /* Release memory and return */
    /*    free_matrix(vector, nvector); */
    free_block(B);
    free(C);
    free(ipiv);
    dpd_free_block(error, 1, vector_length);
  }
  else if(params.ref == 1) { /** ROHF **/

    /* Compute the length of a single error vector */
    /* RAK changed file nums here from CC_TMP0 */
    dpd_file2_init(&L1, CC_LAMBDA, L_irr, 0, 1, "LIA");
    dpd_buf4_init(&L2a, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&L2b, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    for(h=0; h < nirreps; h++) {
      vector_length += 2 * L1.params->rowtot[h] * L1.params->coltot[h^L_irr];
      vector_length += 2 * L2a.params->rowtot[h] * L2a.params->coltot[h^L_irr];
      vector_length += L2b.params->rowtot[h] * L2b.params->coltot[h^L_irr];
    }
    dpd_file2_close(&L1);
    dpd_buf4_close(&L2a);
    dpd_buf4_close(&L2b);

    /* Set the diis cycle value */
    diis_cycle = (iter-1) % nvector;

    /* Build the current error vector and dump it to disk */
    error = dpd_block_matrix(1,vector_length);
    word=0;
    dpd_file2_init(&L1a, CC_LAMBDA, L_irr, 0, 1, "New LIA");
    dpd_file2_mat_init(&L1a);
    dpd_file2_mat_rd(&L1a);
    dpd_file2_init(&L1b, CC_LAMBDA, L_irr, 0, 1, "LIA");
    dpd_file2_mat_init(&L1b);
    dpd_file2_mat_rd(&L1b);
    for(h=0; h < nirreps; h++)
      for(row=0; row < L1a.params->rowtot[h]; row++)
        for(col=0; col < L1a.params->coltot[h^L_irr]; col++)
          error[0][word++] = L1a.matrix[h][row][col] - L1b.matrix[h][row][col];
    dpd_file2_mat_close(&L1a);
    dpd_file2_close(&L1a);
    dpd_file2_mat_close(&L1b);
    dpd_file2_close(&L1b);

    dpd_file2_init(&L1a, CC_LAMBDA, L_irr, 0, 1, "New Lia");
    dpd_file2_mat_init(&L1a);
    dpd_file2_mat_rd(&L1a);
    dpd_file2_init(&L1b, CC_LAMBDA, L_irr, 0, 1, "Lia");
    dpd_file2_mat_init(&L1b);
    dpd_file2_mat_rd(&L1b);
    for(h=0; h < nirreps; h++)
      for(row=0; row < L1a.params->rowtot[h]; row++)
        for(col=0; col < L1a.params->coltot[h^L_irr]; col++)
          error[0][word++] = L1a.matrix[h][row][col] - L1b.matrix[h][row][col];
    dpd_file2_mat_close(&L1a);
    dpd_file2_close(&L1a);
    dpd_file2_mat_close(&L1b);
    dpd_file2_close(&L1b);

    dpd_buf4_init(&L2a, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_init(&L2b, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&L2a, h);
      dpd_buf4_mat_irrep_rd(&L2a, h);
      dpd_buf4_mat_irrep_init(&L2b, h);
      dpd_buf4_mat_irrep_rd(&L2b, h);
      for(row=0; row < L2a.params->rowtot[h]; row++)
        for(col=0; col < L2a.params->coltot[h^L_irr]; col++)
          error[0][word++] = L2a.matrix[h][row][col] - L2b.matrix[h][row][col];
      dpd_buf4_mat_irrep_close(&L2a, h);
      dpd_buf4_mat_irrep_close(&L2b, h);
    }
    dpd_buf4_close(&L2a);
    dpd_buf4_close(&L2b);

    dpd_buf4_init(&L2a, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
    dpd_buf4_init(&L2b, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "Lijab");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&L2a, h);
      dpd_buf4_mat_irrep_rd(&L2a, h);
      dpd_buf4_mat_irrep_init(&L2b, h);
      dpd_buf4_mat_irrep_rd(&L2b, h);
      for(row=0; row < L2a.params->rowtot[h]; row++)
        for(col=0; col < L2a.params->coltot[h^L_irr]; col++)
          error[0][word++] = L2a.matrix[h][row][col] - L2b.matrix[h][row][col];
      dpd_buf4_mat_irrep_close(&L2a, h);
      dpd_buf4_mat_irrep_close(&L2b, h);
    }
    dpd_buf4_close(&L2a);
    dpd_buf4_close(&L2b);

    dpd_buf4_init(&L2a, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    dpd_buf4_init(&L2b, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&L2a, h);
      dpd_buf4_mat_irrep_rd(&L2a, h);
      dpd_buf4_mat_irrep_init(&L2b, h);
      dpd_buf4_mat_irrep_rd(&L2b, h);
      for(row=0; row < L2a.params->rowtot[h]; row++)
        for(col=0; col < L2a.params->coltot[h^L_irr]; col++)
          error[0][word++] = L2a.matrix[h][row][col] - L2b.matrix[h][row][col];
      dpd_buf4_mat_irrep_close(&L2a, h);
      dpd_buf4_mat_irrep_close(&L2b, h);
    }
    dpd_buf4_close(&L2a);
    dpd_buf4_close(&L2b);

    start = psio_get_address(PSIO_ZERO, diis_cycle*vector_length*sizeof(double));
    psio_write(CC_DIIS_ERR, "DIIS Error Vectors" , (char *) error[0],
               vector_length*sizeof(double), start, &end);

    /* Store the current amplitude vector on disk */
    word=0;
    dpd_file2_init(&L1a, CC_LAMBDA, L_irr, 0, 1, "New LIA");
    dpd_file2_mat_init(&L1a);
    dpd_file2_mat_rd(&L1a);
    for(h=0; h < nirreps; h++)
      for(row=0; row < L1a.params->rowtot[h]; row++)
        for(col=0; col < L1a.params->coltot[h^L_irr]; col++)
          error[0][word++] = L1a.matrix[h][row][col];
    dpd_file2_mat_close(&L1a);
    dpd_file2_close(&L1a);

    dpd_file2_init(&L1a, CC_LAMBDA, L_irr, 0, 1, "New Lia");
    dpd_file2_mat_init(&L1a);
    dpd_file2_mat_rd(&L1a);
    for(h=0; h < nirreps; h++)
      for(row=0; row < L1a.params->rowtot[h]; row++)
        for(col=0; col < L1a.params->coltot[h^L_irr]; col++)
          error[0][word++] = L1a.matrix[h][row][col];
    dpd_file2_mat_close(&L1a);
    dpd_file2_close(&L1a);

    dpd_buf4_init(&L2a, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&L2a, h);
      dpd_buf4_mat_irrep_rd(&L2a, h);
      for(row=0; row < L2a.params->rowtot[h]; row++)
        for(col=0; col < L2a.params->coltot[h^L_irr]; col++)
          error[0][word++] = L2a.matrix[h][row][col];
      dpd_buf4_mat_irrep_close(&L2a, h);
    }
    dpd_buf4_close(&L2a);

    dpd_buf4_init(&L2a, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&L2a, h);
      dpd_buf4_mat_irrep_rd(&L2a, h);
      for(row=0; row < L2a.params->rowtot[h]; row++)
        for(col=0; col < L2a.params->coltot[h^L_irr]; col++)
          error[0][word++] = L2a.matrix[h][row][col];
      dpd_buf4_mat_irrep_close(&L2a, h);
    }
    dpd_buf4_close(&L2a);

    dpd_buf4_init(&L2a, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&L2a, h);
      dpd_buf4_mat_irrep_rd(&L2a, h);
      for(row=0; row < L2a.params->rowtot[h]; row++)
        for(col=0; col < L2a.params->coltot[h^L_irr]; col++)
          error[0][word++] = L2a.matrix[h][row][col];
      dpd_buf4_mat_irrep_close(&L2a, h);
    }
    dpd_buf4_close(&L2a);

    start = psio_get_address(PSIO_ZERO, diis_cycle*vector_length*sizeof(double));
    psio_write(CC_DIIS_AMP, "DIIS Amplitude Vectors" , (char *) error[0],
               vector_length*sizeof(double), start, &end);

    /* If we haven't run through enough iterations, set the correct dimensions
       for the extrapolation */
    if(!(iter >= (nvector))) {
      if(iter < 2) { /* Leave if we can't extrapolate at all */
        free(error);
        return;
      }
      nvector = iter;
    }

    /* Now grab the full set of error vectors from the file */
    vector = init_matrix(nvector, vector_length);
    next = PSIO_ZERO;
    for(p=0; p < nvector; p++)
      psio_read(CC_DIIS_ERR, "DIIS Error Vectors", (char *) vector[p],
                vector_length*sizeof(double), next, &next);

    /* Build B matrix of error vector products */
    B = init_matrix(nvector+1,nvector+1);

    for(p=0; p < nvector; p++)
      for(q=0; q < nvector; q++) {
        dot_arr(vector[p], vector[q], vector_length, &product);
        B[p][q] = product;
      }

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
    flin(B, C, nvector+1, 1, &determinant);

    /* Grab the old amplitude vectors */
    next = PSIO_ZERO;
    for(p=0; p < nvector; p++)
      psio_read(CC_DIIS_AMP, "DIIS Amplitude Vectors", (char *) vector[p],
                vector_length*sizeof(double), next, &next);

    /* Build the new amplitude vector from the old ones */
    for(q=0; q < vector_length; q++) {
      error[0][q] = 0.0;
      for(p=0; p < nvector; p++)
        error[0][q] += C[p] * vector[p][q];
    }

    /* Now place these elements into the DPD amplitude arrays */
    word=0;
    dpd_file2_init(&L1a, CC_LAMBDA, L_irr, 0, 1, "New LIA");
    dpd_file2_mat_init(&L1a);
    for(h=0; h < nirreps; h++)
      for(row=0; row < L1a.params->rowtot[h]; row++)
        for(col=0; col < L1a.params->coltot[h^L_irr]; col++)
          L1a.matrix[h][row][col] = error[0][word++];
    dpd_file2_mat_wrt(&L1a);
    dpd_file2_mat_close(&L1a);
    dpd_file2_close(&L1a);

    dpd_file2_init(&L1a, CC_LAMBDA, L_irr, 0, 1, "New Lia");
    dpd_file2_mat_init(&L1a);
    for(h=0; h < nirreps; h++)
      for(row=0; row < L1a.params->rowtot[h]; row++)
        for(col=0; col < L1a.params->coltot[h^L_irr]; col++)
          L1a.matrix[h][row][col] = error[0][word++];
    dpd_file2_mat_wrt(&L1a);
    dpd_file2_mat_close(&L1a);
    dpd_file2_close(&L1a);

    dpd_buf4_init(&L2a, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&L2a, h);
      for(row=0; row < L2a.params->rowtot[h]; row++)
        for(col=0; col < L2a.params->coltot[h^L_irr]; col++)
          L2a.matrix[h][row][col] = error[0][word++];
      dpd_buf4_mat_irrep_wrt(&L2a, h);
      dpd_buf4_mat_irrep_close(&L2a, h);
    }
    dpd_buf4_close(&L2a);

    dpd_buf4_init(&L2a, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New Lijab");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&L2a, h);
      for(row=0; row < L2a.params->rowtot[h]; row++)
        for(col=0; col < L2a.params->coltot[h^L_irr]; col++)
          L2a.matrix[h][row][col] = error[0][word++];
      dpd_buf4_mat_irrep_wrt(&L2a, h);
      dpd_buf4_mat_irrep_close(&L2a, h);
    }
    dpd_buf4_close(&L2a);

    dpd_buf4_init(&L2a, CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&L2a, h);
      for(row=0; row < L2a.params->rowtot[h]; row++)
        for(col=0; col < L2a.params->coltot[h^L_irr]; col++)
          L2a.matrix[h][row][col] = error[0][word++];
      dpd_buf4_mat_irrep_wrt(&L2a, h);
      dpd_buf4_mat_irrep_close(&L2a, h);
    }
    dpd_buf4_close(&L2a);

    /* Release memory and return */
    free_matrix(vector, nvector);
    free_matrix(B, nvector+1);
    free(C);
    dpd_free_block(error, 1, vector_length);
  } /** ROHF **/
  else if(params.ref == 2) { /** UHF **/

    /* Compute the length of a single error vector */
    dpd_file2_init(&L1a, CC_TMP0, L_irr, 0, 1, "LIA");
    dpd_file2_init(&L1b, CC_TMP0, L_irr, 2, 3, "Lia");
    dpd_buf4_init(&L2a, CC_TMP0, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    dpd_buf4_init(&L2b, CC_TMP0, L_irr, 12, 17, 12, 17, 0, "Lijab");
    dpd_buf4_init(&L2c, CC_TMP0, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    for(h=0; h < nirreps; h++) {
      vector_length += L1a.params->rowtot[h] * L1a.params->coltot[h^L_irr];
      vector_length += L1b.params->rowtot[h] * L1b.params->coltot[h^L_irr];
      vector_length += L2a.params->rowtot[h] * L2a.params->coltot[h^L_irr];
      vector_length += L2b.params->rowtot[h] * L2b.params->coltot[h^L_irr];
      vector_length += L2c.params->rowtot[h] * L2c.params->coltot[h^L_irr];
    }
    dpd_file2_close(&L1a);
    dpd_file2_close(&L1b);
    dpd_buf4_close(&L2a);
    dpd_buf4_close(&L2b);
    dpd_buf4_close(&L2c);

    /* Set the diis cycle value */
    diis_cycle = (iter-1) % nvector;

    /* Build the current error vector and dump it to disk */
    error = dpd_block_matrix(1,vector_length);
    word=0;
    dpd_file2_init(&L1a, CC_LAMBDA, L_irr, 0, 1, "New LIA");
    dpd_file2_mat_init(&L1a);
    dpd_file2_mat_rd(&L1a);
    dpd_file2_init(&L1b, CC_LAMBDA, L_irr, 0, 1, "LIA");
    dpd_file2_mat_init(&L1b);
    dpd_file2_mat_rd(&L1b);
    for(h=0; h < nirreps; h++)
      for(row=0; row < L1a.params->rowtot[h]; row++)
        for(col=0; col < L1a.params->coltot[h^L_irr]; col++)
          error[0][word++] = L1a.matrix[h][row][col] - L1b.matrix[h][row][col];
    dpd_file2_mat_close(&L1a);
    dpd_file2_close(&L1a);
    dpd_file2_mat_close(&L1b);
    dpd_file2_close(&L1b);

    dpd_file2_init(&L1a, CC_LAMBDA, L_irr, 2, 3, "New Lia");
    dpd_file2_mat_init(&L1a);
    dpd_file2_mat_rd(&L1a);
    dpd_file2_init(&L1b, CC_LAMBDA, L_irr, 2, 3, "Lia");
    dpd_file2_mat_init(&L1b);
    dpd_file2_mat_rd(&L1b);
    for(h=0; h < nirreps; h++)
      for(row=0; row < L1a.params->rowtot[h]; row++)
        for(col=0; col < L1a.params->coltot[h^L_irr]; col++)
          error[0][word++] = L1a.matrix[h][row][col] - L1b.matrix[h][row][col];
    dpd_file2_mat_close(&L1a);
    dpd_file2_close(&L1a);
    dpd_file2_mat_close(&L1b);
    dpd_file2_close(&L1b);

    dpd_buf4_init(&L2a, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    dpd_buf4_init(&L2b, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "LIJAB");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&L2a, h);
      dpd_buf4_mat_irrep_rd(&L2a, h);
      dpd_buf4_mat_irrep_init(&L2b, h);
      dpd_buf4_mat_irrep_rd(&L2b, h);
      for(row=0; row < L2a.params->rowtot[h]; row++)
        for(col=0; col < L2a.params->coltot[h^L_irr]; col++)
          error[0][word++] = L2a.matrix[h][row][col] - L2b.matrix[h][row][col];
      dpd_buf4_mat_irrep_close(&L2a, h);
      dpd_buf4_mat_irrep_close(&L2b, h);
    }
    dpd_buf4_close(&L2a);
    dpd_buf4_close(&L2b);

    dpd_buf4_init(&L2a, CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "New Lijab");
    dpd_buf4_init(&L2b, CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "Lijab");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&L2a, h);
      dpd_buf4_mat_irrep_rd(&L2a, h);
      dpd_buf4_mat_irrep_init(&L2b, h);
      dpd_buf4_mat_irrep_rd(&L2b, h);
      for(row=0; row < L2a.params->rowtot[h]; row++)
        for(col=0; col < L2a.params->coltot[h^L_irr]; col++)
          error[0][word++] = L2a.matrix[h][row][col] - L2b.matrix[h][row][col];
      dpd_buf4_mat_irrep_close(&L2a, h);
      dpd_buf4_mat_irrep_close(&L2b, h);
    }
    dpd_buf4_close(&L2a);
    dpd_buf4_close(&L2b);

    dpd_buf4_init(&L2a, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    dpd_buf4_init(&L2b, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&L2a, h);
      dpd_buf4_mat_irrep_rd(&L2a, h);
      dpd_buf4_mat_irrep_init(&L2b, h);
      dpd_buf4_mat_irrep_rd(&L2b, h);
      for(row=0; row < L2a.params->rowtot[h]; row++)
        for(col=0; col < L2a.params->coltot[h^L_irr]; col++)
          error[0][word++] = L2a.matrix[h][row][col] - L2b.matrix[h][row][col];
      dpd_buf4_mat_irrep_close(&L2a, h);
      dpd_buf4_mat_irrep_close(&L2b, h);
    }
    dpd_buf4_close(&L2a);
    dpd_buf4_close(&L2b);

    start = psio_get_address(PSIO_ZERO, diis_cycle*vector_length*sizeof(double));
    psio_write(CC_DIIS_ERR, "DIIS Error[0] Vectors" , (char *) error[0],
               vector_length*sizeof(double), start, &end);

    /* Store the current amplitude vector on disk */
    word=0;
    dpd_file2_init(&L1a, CC_LAMBDA, L_irr, 0, 1, "New LIA");
    dpd_file2_mat_init(&L1a);
    dpd_file2_mat_rd(&L1a);
    for(h=0; h < nirreps; h++)
      for(row=0; row < L1a.params->rowtot[h]; row++)
        for(col=0; col < L1a.params->coltot[h^L_irr]; col++)
          error[0][word++] = L1a.matrix[h][row][col];
    dpd_file2_mat_close(&L1a);
    dpd_file2_close(&L1a);

    dpd_file2_init(&L1a, CC_LAMBDA, L_irr, 2, 3, "New Lia");
    dpd_file2_mat_init(&L1a);
    dpd_file2_mat_rd(&L1a);
    for(h=0; h < nirreps; h++)
      for(row=0; row < L1a.params->rowtot[h]; row++)
        for(col=0; col < L1a.params->coltot[h^L_irr]; col++)
          error[0][word++] = L1a.matrix[h][row][col];
    dpd_file2_mat_close(&L1a);
    dpd_file2_close(&L1a);

    dpd_buf4_init(&L2a, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&L2a, h);
      dpd_buf4_mat_irrep_rd(&L2a, h);
      for(row=0; row < L2a.params->rowtot[h]; row++)
        for(col=0; col < L2a.params->coltot[h^L_irr]; col++)
          error[0][word++] = L2a.matrix[h][row][col];
      dpd_buf4_mat_irrep_close(&L2a, h);
    }
    dpd_buf4_close(&L2a);

    dpd_buf4_init(&L2a, CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "New Lijab");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&L2a, h);
      dpd_buf4_mat_irrep_rd(&L2a, h);
      for(row=0; row < L2a.params->rowtot[h]; row++)
        for(col=0; col < L2a.params->coltot[h^L_irr]; col++)
          error[0][word++] = L2a.matrix[h][row][col];
      dpd_buf4_mat_irrep_close(&L2a, h);
    }
    dpd_buf4_close(&L2a);

    dpd_buf4_init(&L2a, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&L2a, h);
      dpd_buf4_mat_irrep_rd(&L2a, h);
      for(row=0; row < L2a.params->rowtot[h]; row++)
        for(col=0; col < L2a.params->coltot[h^L_irr]; col++)
          error[0][word++] = L2a.matrix[h][row][col];
      dpd_buf4_mat_irrep_close(&L2a, h);
    }
    dpd_buf4_close(&L2a);

    start = psio_get_address(PSIO_ZERO, diis_cycle*vector_length*sizeof(double));
    psio_write(CC_DIIS_AMP, "DIIS Amplitude Vectors" , (char *) error[0],
               vector_length*sizeof(double), start, &end);

    /* If we haven't run through enough iterations, set the correct dimensions
       for the extrapolation */
    if(!(iter >= (nvector))) {
      if(iter < 2) { /* Leave if we can't extrapolate at all */
        free(error[0]);
        return;
      }
      nvector = iter;
    }

    /* Now grab the full set of error[0] vectors from the file */
    vector = init_matrix(nvector, vector_length);
    next = PSIO_ZERO;
    for(p=0; p < nvector; p++)
      psio_read(CC_DIIS_ERR, "DIIS Error[0] Vectors", (char *) vector[p],
                vector_length*sizeof(double), next, &next);

    /* Build B matrix of error[0] vector products */
    B = init_matrix(nvector+1,nvector+1);

    for(p=0; p < nvector; p++)
      for(q=0; q < nvector; q++) {
        dot_arr(vector[p], vector[q], vector_length, &product);
        B[p][q] = product;
      }

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
    flin(B, C, nvector+1, 1, &determinant);

    /* Grab the old amplitude vectors */
    next = PSIO_ZERO;
    for(p=0; p < nvector; p++)
      psio_read(CC_DIIS_AMP, "DIIS Amplitude Vectors", (char *) vector[p],
                vector_length*sizeof(double), next, &next);

    /* Build the new amplitude vector from the old ones */
    for(q=0; q < vector_length; q++) {
      error[0][q] = 0.0;
      for(p=0; p < nvector; p++)
        error[0][q] += C[p] * vector[p][q];
    }

    /* Now place these elements into the DPD amplitude arrays */
    word=0;
    dpd_file2_init(&L1a, CC_LAMBDA, L_irr, 0, 1, "New LIA");
    dpd_file2_mat_init(&L1a);
    for(h=0; h < nirreps; h++)
      for(row=0; row < L1a.params->rowtot[h]; row++)
        for(col=0; col < L1a.params->coltot[h^L_irr]; col++)
          L1a.matrix[h][row][col] = error[0][word++];
    dpd_file2_mat_wrt(&L1a);
    dpd_file2_mat_close(&L1a);
    dpd_file2_close(&L1a);

    dpd_file2_init(&L1a, CC_LAMBDA, L_irr, 2, 3, "New Lia");
    dpd_file2_mat_init(&L1a);
    for(h=0; h < nirreps; h++)
      for(row=0; row < L1a.params->rowtot[h]; row++)
        for(col=0; col < L1a.params->coltot[h^L_irr]; col++)
          L1a.matrix[h][row][col] = error[0][word++];
    dpd_file2_mat_wrt(&L1a);
    dpd_file2_mat_close(&L1a);
    dpd_file2_close(&L1a);

    dpd_buf4_init(&L2a, CC_LAMBDA, L_irr, 2, 7, 2, 7, 0, "New LIJAB");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&L2a, h);
      for(row=0; row < L2a.params->rowtot[h]; row++)
        for(col=0; col < L2a.params->coltot[h^L_irr]; col++)
          L2a.matrix[h][row][col] = error[0][word++];
      dpd_buf4_mat_irrep_wrt(&L2a, h);
      dpd_buf4_mat_irrep_close(&L2a, h);
    }
    dpd_buf4_close(&L2a);

    dpd_buf4_init(&L2a, CC_LAMBDA, L_irr, 12, 17, 12, 17, 0, "New Lijab");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&L2a, h);
      for(row=0; row < L2a.params->rowtot[h]; row++)
        for(col=0; col < L2a.params->coltot[h^L_irr]; col++)
          L2a.matrix[h][row][col] = error[0][word++];
      dpd_buf4_mat_irrep_wrt(&L2a, h);
      dpd_buf4_mat_irrep_close(&L2a, h);
    }
    dpd_buf4_close(&L2a);

    dpd_buf4_init(&L2a, CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    for(h=0; h < nirreps; h++) {
      dpd_buf4_mat_irrep_init(&L2a, h);
      for(row=0; row < L2a.params->rowtot[h]; row++)
        for(col=0; col < L2a.params->coltot[h^L_irr]; col++)
          L2a.matrix[h][row][col] = error[0][word++];
      dpd_buf4_mat_irrep_wrt(&L2a, h);
      dpd_buf4_mat_irrep_close(&L2a, h);
    }
    dpd_buf4_close(&L2a);

    /* Release memory and return */
    free_matrix(vector, nvector);
    free_matrix(B, nvector+1);
    free(C);
    dpd_free_block(error, 1, vector_length);
  } /** UHF **/

  return;
}

}} // namespace psi::cclambda
