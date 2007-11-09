/*! \file 
    \ingroup (DPD)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <stdlib.h>
#include <libciomr/libciomr.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

extern "C" {

int dpd_trans4_mat_irrep_close(dpdtrans4 *Trans, int irrep)
{
  int h, nirreps, rowtot, coltot, all_buf_irrep;
  long int size;

  all_buf_irrep = Trans->buf.file.my_irrep;
  nirreps = Trans->buf.params->nirreps;
  rowtot = Trans->buf.params->coltot[irrep^all_buf_irrep];
  coltot = Trans->buf.params->rowtot[irrep];
  size = ((long) rowtot) * ((long) coltot);

  /* Free the shift structure for this irrep if used */
  if(Trans->shift.shift_type) {
      for(h=0; h < nirreps; h++)
	  if(Trans->shift.rowtot[irrep][h])
	      free(Trans->shift.matrix[irrep][h]);
      free(Trans->shift.matrix[irrep]);
      Trans->shift.shift_type = 0;
    }

  if(size)
      dpd_free_block(Trans->matrix[irrep], rowtot, coltot);
  
  return 0;
}

} /* extern "C" */
