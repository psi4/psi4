/*! \file 
    \ingroup (DPD)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <libpsio/psio.h>
#include "dpd.h"

extern "C" {

int dpd_file4_mat_irrep_rd_block(dpdfile4 *File, int irrep, int start_pq,
				int num_pq)
{
  int rowtot, coltot, my_irrep;
  int seek_block;
  psio_address irrep_ptr, next_address;
  long int size;

  my_irrep = File->my_irrep;
  if(File->incore) return 0;  /* We already have this data in core */

  irrep_ptr = File->lfiles[irrep];
  rowtot = num_pq;
  coltot = File->params->coltot[irrep^my_irrep];

  size = ((long) rowtot) * ((long) coltot);

  /* Advance file pointer to current row --- careful about overflows! */
  if(coltot) {
    seek_block = DPD_BIGNUM/(coltot * sizeof(double)); /* no. of rows for which we can compute the address */
    if(seek_block < 1) {
      fprintf(stderr, "\nLIBDPD Error: each row of %s is too long to compute an address.\n",File->label);
      dpd_error("dpd_file4_mat_irrep_rd_block", stderr);
    }
    for(; start_pq > seek_block; start_pq -= seek_block)
      irrep_ptr = psio_get_address(irrep_ptr, seek_block*coltot*sizeof(double));
    irrep_ptr = psio_get_address(irrep_ptr, start_pq*coltot*sizeof(double));
  }

  if(rowtot && coltot)
     psio_read(File->filenum, File->label, (char *) File->matrix[irrep][0],
	       size * ((long) sizeof(double)), irrep_ptr, &next_address);

  return 0;

}

} /* extern "C" */
