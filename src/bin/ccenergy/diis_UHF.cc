/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

/*
** DIIS: Direct inversion in the iterative subspace routine to
** accelerate convergence of the CCSD amplitude equations.
**
** Substantially improved efficiency of this routine:
** (1) Keeping at most two error vectors in core at once.
** (2) Limiting direct product (overlap) calculation to unique pairs.
** (3) Using LAPACK's linear equation solver DGESV instead of flin.
**
** -TDC  12/22/01
**
** Updated for UHF refs, TDC, 6/03
*/

void diis_UHF(int iter)
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

  nirreps = moinfo.nirreps;

  /* Compute the length of a single error vector */
  dpd_file2_init(&T1a, CC_TMP0, 0, 0, 1, "tIA");
  dpd_file2_init(&T1b, CC_TMP0, 0, 2, 3, "tIA");
  dpd_buf4_init(&T2a, CC_TMP0, 0, 2, 7, 2, 7, 0, "tIJAB");
  dpd_buf4_init(&T2b, CC_TMP0, 0, 12, 17, 12, 17, 0, "tIJAB");
  dpd_buf4_init(&T2c, CC_TMP0, 0, 22, 28, 22, 28, 0, "tIjAb");
  for(h=0; h < nirreps; h++) {
    vector_length += T1a.params->rowtot[h] * T1a.params->coltot[h];
    vector_length += T1b.params->rowtot[h] * T1b.params->coltot[h];
    vector_length += T2a.params->rowtot[h] * T2a.params->coltot[h];
    vector_length += T2b.params->rowtot[h] * T2b.params->coltot[h];
    vector_length += T2c.params->rowtot[h] * T2c.params->coltot[h];
  }
  dpd_file2_close(&T1a);
  dpd_file2_close(&T1b);
  dpd_buf4_close(&T2a);
  dpd_buf4_close(&T2b);
  dpd_buf4_close(&T2c);

  /* Set the diis cycle value */
  diis_cycle = (iter-1) % nvector;

  /* Build the current error vector and dump it to disk */
  error = dpd_block_matrix(1,vector_length);

  word=0;
  dpd_file2_init(&T1a, CC_OEI, 0, 0, 1, "New tIA");
  dpd_file2_mat_init(&T1a);
  dpd_file2_mat_rd(&T1a);
  dpd_file2_init(&T1b, CC_OEI, 0, 0, 1, "tIA");
  dpd_file2_mat_init(&T1b);
  dpd_file2_mat_rd(&T1b);
  for(h=0; h < nirreps; h++)
    for(row=0; row < T1a.params->rowtot[h]; row++)
      for(col=0; col < T1a.params->coltot[h]; col++)
	error[0][word++] = T1a.matrix[h][row][col] - T1b.matrix[h][row][col];
  dpd_file2_mat_close(&T1a);
  dpd_file2_close(&T1a);
  dpd_file2_mat_close(&T1b);
  dpd_file2_close(&T1b);

  dpd_file2_init(&T1a, CC_OEI, 0, 2, 3, "New tia");
  dpd_file2_mat_init(&T1a);
  dpd_file2_mat_rd(&T1a);
  dpd_file2_init(&T1b, CC_OEI, 0, 2, 3, "tia");
  dpd_file2_mat_init(&T1b);
  dpd_file2_mat_rd(&T1b);
  for(h=0; h < nirreps; h++)
    for(row=0; row < T1a.params->rowtot[h]; row++)
      for(col=0; col < T1a.params->coltot[h]; col++)
	error[0][word++] = T1a.matrix[h][row][col] - T1b.matrix[h][row][col];
  dpd_file2_mat_close(&T1a);
  dpd_file2_close(&T1a);
  dpd_file2_mat_close(&T1b);
  dpd_file2_close(&T1b);
  
  dpd_buf4_init(&T2a, CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
  dpd_buf4_init(&T2b, CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2a, h);
    dpd_buf4_mat_irrep_rd(&T2a, h);
    dpd_buf4_mat_irrep_init(&T2b, h);
    dpd_buf4_mat_irrep_rd(&T2b, h);
    for(row=0; row < T2a.params->rowtot[h]; row++)
      for(col=0; col < T2a.params->coltot[h]; col++)
	error[0][word++] = T2a.matrix[h][row][col] - T2b.matrix[h][row][col];
    dpd_buf4_mat_irrep_close(&T2a, h);
    dpd_buf4_mat_irrep_close(&T2b, h);
  }
  dpd_buf4_close(&T2a);
  dpd_buf4_close(&T2b);

  dpd_buf4_init(&T2a, CC_TAMPS, 0, 12, 17, 12, 17, 0, "New tijab");
  dpd_buf4_init(&T2b, CC_TAMPS, 0, 12, 17, 12, 17, 0, "tijab");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2a, h);
    dpd_buf4_mat_irrep_rd(&T2a, h);
    dpd_buf4_mat_irrep_init(&T2b, h);
    dpd_buf4_mat_irrep_rd(&T2b, h);
    for(row=0; row < T2a.params->rowtot[h]; row++)
      for(col=0; col < T2a.params->coltot[h]; col++)
	error[0][word++] = T2a.matrix[h][row][col] - T2b.matrix[h][row][col];
    dpd_buf4_mat_irrep_close(&T2a, h);
    dpd_buf4_mat_irrep_close(&T2b, h);
  }
  dpd_buf4_close(&T2a);
  dpd_buf4_close(&T2b);

  dpd_buf4_init(&T2a, CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
  dpd_buf4_init(&T2b, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2a, h);
    dpd_buf4_mat_irrep_rd(&T2a, h);
    dpd_buf4_mat_irrep_init(&T2b, h);
    dpd_buf4_mat_irrep_rd(&T2b, h);
    for(row=0; row < T2a.params->rowtot[h]; row++)
      for(col=0; col < T2a.params->coltot[h]; col++)
	error[0][word++] = T2a.matrix[h][row][col] - T2b.matrix[h][row][col];
    dpd_buf4_mat_irrep_close(&T2a, h);
    dpd_buf4_mat_irrep_close(&T2b, h);
  }
  dpd_buf4_close(&T2a);
  dpd_buf4_close(&T2b);

  start = psio_get_address(PSIO_ZERO, diis_cycle*vector_length*sizeof(double));
  psio_write(CC_DIIS_ERR, "DIIS Error Vectors" , (char *) error[0], 
	     vector_length*sizeof(double), start, &end);

  /* Store the current amplitude vector on disk */
  word=0;

  dpd_file2_init(&T1a, CC_OEI, 0, 0, 1, "New tIA");
  dpd_file2_mat_init(&T1a);
  dpd_file2_mat_rd(&T1a);
  for(h=0; h < nirreps; h++)
    for(row=0; row < T1a.params->rowtot[h]; row++)
      for(col=0; col < T1a.params->coltot[h]; col++)
	error[0][word++] = T1a.matrix[h][row][col];
  dpd_file2_mat_close(&T1a);
  dpd_file2_close(&T1a);

  dpd_file2_init(&T1a, CC_OEI, 0, 2, 3, "New tia");
  dpd_file2_mat_init(&T1a);
  dpd_file2_mat_rd(&T1a);
  for(h=0; h < nirreps; h++)
    for(row=0; row < T1a.params->rowtot[h]; row++)
      for(col=0; col < T1a.params->coltot[h]; col++)
	error[0][word++] = T1a.matrix[h][row][col];
  dpd_file2_mat_close(&T1a);
  dpd_file2_close(&T1a);
  
  dpd_buf4_init(&T2a, CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2a, h);
    dpd_buf4_mat_irrep_rd(&T2a, h);
    for(row=0; row < T2a.params->rowtot[h]; row++)
      for(col=0; col < T2a.params->coltot[h]; col++)
	error[0][word++] = T2a.matrix[h][row][col];
    dpd_buf4_mat_irrep_close(&T2a, h);
  }
  dpd_buf4_close(&T2a);

  dpd_buf4_init(&T2a, CC_TAMPS, 0, 12, 17, 12, 17, 0, "New tijab");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2a, h);
    dpd_buf4_mat_irrep_rd(&T2a, h);
    for(row=0; row < T2a.params->rowtot[h]; row++)
      for(col=0; col < T2a.params->coltot[h]; col++)
	error[0][word++] = T2a.matrix[h][row][col];
    dpd_buf4_mat_irrep_close(&T2a, h);
  }
  dpd_buf4_close(&T2a);

  dpd_buf4_init(&T2a, CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2a, h);
    dpd_buf4_mat_irrep_rd(&T2a, h);
    for(row=0; row < T2a.params->rowtot[h]; row++)
      for(col=0; col < T2a.params->coltot[h]; col++)
	error[0][word++] = T2a.matrix[h][row][col];
    dpd_buf4_mat_irrep_close(&T2a, h);
  }
  dpd_buf4_close(&T2a);

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

  /* Now grab the full set of error[0] vectors from the file */
  vector = dpd_block_matrix(2, vector_length);
  B = block_matrix(nvector+1,nvector+1);
  for(p=0; p < nvector; p++) {

    start = psio_get_address(PSIO_ZERO, p*vector_length*sizeof(double));

    psio_read(CC_DIIS_ERR, "DIIS Error Vectors", (char *) vector[0],
	      vector_length*sizeof(double), start, &end);

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

  /* Build the constant vector */
  C = init_array(nvector+1);
  C[nvector] = -1;

  /* Solve the linear equations */
  ipiv = init_int_array(nvector+1);

  errcod = C_DGESV(nvector+1, 1, &(B[0][0]), nvector+1, &(ipiv[0]), &(C[0]), nvector+1);
  if(errcod) {
    fprintf(outfile, "\nError in DGESV return in diis.\n");
    throw PsiException("Error in DGESV return in diis.", __FILE__, __LINE__);
  }

  /* Build the new amplitude vector from the old ones */
  vector = dpd_block_matrix(1, vector_length);
  for(p=0; p < vector_length; p++) error[0][p] = 0.0;
  for(p=0; p < nvector; p++) {

    start = psio_get_address(PSIO_ZERO, p*vector_length*sizeof(double));

    psio_read(CC_DIIS_AMP, "DIIS Amplitude Vectors", (char *) vector[0], 
	      vector_length*sizeof(double), start, &end);

    for(q=0; q < vector_length; q++) 
      error[0][q] += C[p] * vector[0][q];

  }
  dpd_free_block(vector, 1, vector_length);

  /* Now place these elements into the DPD amplitude arrays */
  word=0;
  dpd_file2_init(&T1a, CC_OEI, 0, 0, 1, "New tIA");
  dpd_file2_mat_init(&T1a);
  for(h=0; h < nirreps; h++)
    for(row=0; row < T1a.params->rowtot[h]; row++)
      for(col=0; col < T1a.params->coltot[h]; col++)
	T1a.matrix[h][row][col] = error[0][word++];
  dpd_file2_mat_wrt(&T1a);
  dpd_file2_mat_close(&T1a);
  dpd_file2_close(&T1a);

  dpd_file2_init(&T1a, CC_OEI, 0, 2, 3, "New tia");
  dpd_file2_mat_init(&T1a);
  for(h=0; h < nirreps; h++)
    for(row=0; row < T1a.params->rowtot[h]; row++)
      for(col=0; col < T1a.params->coltot[h]; col++)
	T1a.matrix[h][row][col] = error[0][word++];
  dpd_file2_mat_wrt(&T1a);
  dpd_file2_mat_close(&T1a);
  dpd_file2_close(&T1a);
  
  dpd_buf4_init(&T2a, CC_TAMPS, 0, 2, 7, 2, 7, 0, "New tIJAB");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2a, h);
    for(row=0; row < T2a.params->rowtot[h]; row++)
      for(col=0; col < T2a.params->coltot[h]; col++)
	T2a.matrix[h][row][col] = error[0][word++];
    dpd_buf4_mat_irrep_wrt(&T2a, h);
    dpd_buf4_mat_irrep_close(&T2a, h);
  }
  dpd_buf4_close(&T2a);

  dpd_buf4_init(&T2a, CC_TAMPS, 0, 12, 17, 12, 17, 0, "New tijab");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2a, h);
    for(row=0; row < T2a.params->rowtot[h]; row++)
      for(col=0; col < T2a.params->coltot[h]; col++)
	T2a.matrix[h][row][col] = error[0][word++];
    dpd_buf4_mat_irrep_wrt(&T2a, h);
    dpd_buf4_mat_irrep_close(&T2a, h);
  }
  dpd_buf4_close(&T2a);

  dpd_buf4_init(&T2a, CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
  for(h=0; h < nirreps; h++) {
    dpd_buf4_mat_irrep_init(&T2a, h);
    for(row=0; row < T2a.params->rowtot[h]; row++)
      for(col=0; col < T2a.params->coltot[h]; col++)
	T2a.matrix[h][row][col] = error[0][word++];
    dpd_buf4_mat_irrep_wrt(&T2a, h);
    dpd_buf4_mat_irrep_close(&T2a, h);
  }
  dpd_buf4_close(&T2a);

  /* Release memory and return */
  free_block(B);
  free(C);
  free(ipiv);
  dpd_free_block(error, 1, vector_length);

  return;
}
}} // namespace psi::ccenergy
