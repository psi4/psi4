/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include "dpd.h"

namespace psi {

int dpd_file4_mat_irrep_row_rd(dpdfile4 *File, int irrep, int row)
{
  int coltot, my_irrep, seek_block;
  psio_address row_ptr, next_address;

  if(File->incore) return 0;  /* We already have this data in core */

#ifdef DPD_TIMER
  timer_on("f4_rowrd");
#endif

  my_irrep = File->my_irrep;

  row_ptr = File->lfiles[irrep];
  coltot = File->params->coltot[irrep^my_irrep];

  /* Advance file pointer to current row --- careful about overflows! */
  if(coltot) {
    seek_block = DPD_BIGNUM/(coltot * sizeof(double)); /* no. of rows for which we can compute the address */
    if(seek_block < 1) {
      fprintf(stderr, "\nLIBDPD Error: each row of %s is too long to compute an address.\n",File->label);
      dpd_error("dpd_file4_mat_irrep_row_rd", stderr);
    }
    for(; row > seek_block; row -= seek_block)
      row_ptr = psio_get_address(row_ptr, seek_block*coltot*sizeof(double));
    row_ptr = psio_get_address(row_ptr, row*coltot*sizeof(double));
  }

  if(coltot) 
      psio_read(File->filenum, File->label, (char *) File->matrix[irrep][0],
		coltot*sizeof(double), row_ptr, &next_address);

#ifdef DPD_TIMER
  timer_off("f4_rowrd");
#endif

  return 0;

}

}
