/*! \file 
    \ingroup (DPD)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "dpd.h"

extern "C" {
	
/* dpd_buf4_copy(): Copies an existing four-index dpdbuf4 into another file.
**
** Arguments:
**   dpdbuf4 *InBuf: A pointer to the given dpd buffer.
**   int outfilenum: The PSI unit number for the new buffer.
**   char *label: A string labelling for this buffer.
**
** NB: The buffer and file pq/rs parameters are assumed to be
** identical for the copy, obviously.  Hence, anti flag must be off.
**
** Converted to use buf4 only rather than assumptions about file4.
** TDC, September 1999
*/

int dpd_buf4_copy(dpdbuf4 *InBuf, int outfilenum, char *label)
{
  int h, row, col, my_irrep;
  long int rowtot, coltot;
  int nbuckets, incore, n;
  long int memoryd, rows_per_bucket, rows_left, size;
  dpdbuf4 OutBuf;

  my_irrep = InBuf->file.my_irrep;

  dpd_buf4_init(&OutBuf, outfilenum, InBuf->file.my_irrep, InBuf->params->pqnum,
		InBuf->params->rsnum, InBuf->params->pqnum, 
                InBuf->params->rsnum, 0, label);

  for(h=0; h < InBuf->params->nirreps; h++) {

    memoryd = dpd_memfree()/2; /* use half the memory for each buf4 */

    rowtot = InBuf->params->rowtot[h];
    coltot = InBuf->params->coltot[h^my_irrep];

    if(rowtot && coltot) {

      rows_per_bucket = memoryd/coltot;
      /* enough memory for the whole matrix? */
      if(rows_per_bucket > rowtot)
	rows_per_bucket = rowtot;

      if(!rows_per_bucket) dpd_error("buf4_scmcopy: Not enough memory for one row!", stderr);

      nbuckets = (int) ceil(((double) rowtot)/((double) rows_per_bucket));

      rows_left = rowtot % rows_per_bucket;

      incore = 1;
      if(nbuckets > 1) {
	incore = 0;
#if DPD_DEBUG
        fprintf(stderr, "buf4_copy: memory information.\n");
        fprintf(stderr, "buf4_copy: rowtot[%d] = %d.\n", h, InBuf->params->rowtot[h]);
	fprintf(stderr, "buf4_copy: nbuckets = %d\n", nbuckets);
	fprintf(stderr, "buf4_copy: rows_per_bucket = %d\n", rows_per_bucket);
	fprintf(stderr, "buf4_copy: rows_left = %d\n", rows_left);
	fprintf(stderr, "buf4_copy: out-of-core algorithm used\n");
#endif
      }

      if(incore) {


	dpd_buf4_mat_irrep_init(InBuf, h);
	dpd_buf4_mat_irrep_rd(InBuf, h);

	dpd_buf4_mat_irrep_init(&OutBuf, h);

	if(rowtot && coltot) 
	  memcpy((void *) &(OutBuf.matrix[h][0][0]),
		 (const void *) &(InBuf->matrix[h][0][0]),
		 sizeof(double)*rowtot*coltot);

	dpd_buf4_mat_irrep_wrt(&OutBuf, h);

	dpd_buf4_mat_irrep_close(&OutBuf, h);
	dpd_buf4_mat_irrep_close(InBuf, h);
      }
      else {

	dpd_buf4_mat_irrep_init_block(InBuf, h, rows_per_bucket);
	dpd_buf4_mat_irrep_init_block(&OutBuf, h, rows_per_bucket);

	coltot = InBuf->params->coltot[h^my_irrep];
	size = ((long) rows_per_bucket)*((long) coltot);

	for(n=0; n < (rows_left ? nbuckets-1 : nbuckets); n++) {

	  dpd_buf4_mat_irrep_rd_block(InBuf, h, n*rows_per_bucket, rows_per_bucket);

	  memcpy((void *) &(OutBuf.matrix[h][0][0]), (const void *) &(InBuf->matrix[h][0][0]), 
		 ((long) sizeof(double))*size);

	  dpd_buf4_mat_irrep_wrt_block(&OutBuf, h, n*rows_per_bucket, rows_per_bucket);
	}
	if(rows_left) {

	  size = ((long) rows_left) * ((long) coltot);

	  dpd_buf4_mat_irrep_rd_block(InBuf, h, n*rows_per_bucket, rows_left);

	  memcpy((void *) &(OutBuf.matrix[h][0][0]), (const void *) &(InBuf->matrix[h][0][0]), 
		 ((long) sizeof(double))*size);

	  dpd_buf4_mat_irrep_wrt_block(&OutBuf, h, n*rows_per_bucket, rows_left);
	}

	dpd_buf4_mat_irrep_close_block(InBuf, h, rows_per_bucket);
	dpd_buf4_mat_irrep_close_block(&OutBuf, h, rows_per_bucket);

      }
    }

  }

  dpd_buf4_close(&OutBuf);

  return 0;
}

} /* extern "C" */
