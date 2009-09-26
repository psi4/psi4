/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cmath>
#include "dpd.h"
#include <libqt/qt.h>

namespace psi {
	
/* dpd_buf4_axpy(): Evaluates the standard operation a * X + Y -> Y for dpd
** four-index buffers.
**
** Arguments:
**   dpdbuf4 *BufX: A pointer to the leftmost dpdbuf4.
**   dpdbuf4 *BufY: A pointer to the rightmost (and target)
**                        dpdbuf4.
**   double alpha: The scalar prefactor in the multiplication.
*/

int dpd_buf4_axpy(dpdbuf4 *BufX, dpdbuf4 *BufY, double alpha)
{
  int h, nirreps, my_irrep;
  int row, col, incore, n, nbuckets;
  long int length;
  long int memoryd, rows_per_bucket, rows_left;
  double *X, *Y;

  nirreps = BufX->params->nirreps;
  my_irrep = BufX->file.my_irrep;

#ifdef DPD_TIMER
  timer_on("buf4_axpy");
#endif

  for(h=0; h < nirreps; h++) {

    memoryd = (dpd_memfree()-BufX->file.params->coltot[h^my_irrep])/2; /* use half the memory for each buf4 */
    if(BufX->params->rowtot[h] && BufX->params->coltot[h^my_irrep]) {

      rows_per_bucket = memoryd/BufX->params->coltot[h^my_irrep];

      /* enough memory for the whole matrix? */
      if(rows_per_bucket > BufX->params->rowtot[h]) 
	rows_per_bucket = BufX->params->rowtot[h]; 

      if(!rows_per_bucket) dpd_error("buf4_axpy: Not enough memory for one row!", stderr);

      nbuckets = (int) ceil(((double) BufX->params->rowtot[h])/((double) rows_per_bucket));

      rows_left = BufX->params->rowtot[h] % rows_per_bucket;

      incore = 1;
      if(nbuckets > 1) {
	incore = 0;
#if DPD_DEBUG
        fprintf(stderr, "buf4_axpy: memory information.\n");
	fprintf(stderr, "buf4_axpy: rowtot[%d] = %d\n", h, BufX->params->rowtot[h]);
	fprintf(stderr, "buf4_axpy: nbuckets = %d\n", nbuckets);
	fprintf(stderr, "buf4_axpy: rows_per_bucket = %d\n", rows_per_bucket);
	fprintf(stderr, "buf4_axpy: rows_left = %d\n", rows_left);
	fprintf(stderr, "buf4_axpy: out-of-core algorithm used\n");
#endif
      }
    }
    else incore = 1;

    if(incore) {
      dpd_buf4_mat_irrep_init(BufX, h);
      dpd_buf4_mat_irrep_rd(BufX, h);

      dpd_buf4_mat_irrep_init(BufY, h);
      dpd_buf4_mat_irrep_rd(BufY, h);

      length = ((long) BufX->params->rowtot[h]) * ((long) BufX->params->coltot[h^my_irrep]);
      if(length) {
	X = &(BufX->matrix[h][0][0]);
	Y = &(BufY->matrix[h][0][0]);
	C_DAXPY(length, alpha, X, 1, Y, 1);
      }

      dpd_buf4_mat_irrep_wrt(BufY, h);

      dpd_buf4_mat_irrep_close(BufX, h);
      dpd_buf4_mat_irrep_close(BufY, h);
    }
    else {

      dpd_buf4_mat_irrep_init_block(BufX, h, rows_per_bucket);
      dpd_buf4_mat_irrep_init_block(BufY, h, rows_per_bucket);

      length = ((long) rows_per_bucket) * ((long) BufX->params->coltot[h^my_irrep]);
      X = &(BufX->matrix[h][0][0]);
      Y = &(BufY->matrix[h][0][0]);

      for(n=0; n < (rows_left ? nbuckets-1 : nbuckets); n++) {

	dpd_buf4_mat_irrep_rd_block(BufX, h, n*rows_per_bucket, rows_per_bucket);
	dpd_buf4_mat_irrep_rd_block(BufY, h, n*rows_per_bucket, rows_per_bucket);

	C_DAXPY(length, alpha, X, 1, Y, 1);

	dpd_buf4_mat_irrep_wrt_block(BufY, h, n*rows_per_bucket, rows_per_bucket);
      }

      if(rows_left) {

	length = ((long) rows_left) * ((long) BufX->params->coltot[h^my_irrep]);

	dpd_buf4_mat_irrep_rd_block(BufX, h, n*rows_per_bucket, rows_left);
	dpd_buf4_mat_irrep_rd_block(BufY, h, n*rows_per_bucket, rows_left);

	C_DAXPY(length, alpha, X, 1, Y, 1);

	dpd_buf4_mat_irrep_wrt_block(BufY, h, n*rows_per_bucket, rows_left);
      }

      dpd_buf4_mat_irrep_close_block(BufX, h, rows_per_bucket);
      dpd_buf4_mat_irrep_close_block(BufY, h, rows_per_bucket);
    }
  }

#ifdef DPD_TIMER
  timer_off("buf4_axpy");
#endif

  return 0;
}

}
