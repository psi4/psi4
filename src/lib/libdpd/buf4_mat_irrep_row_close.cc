/*! \file 
    \ingroup (DPD)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <libciomr/libciomr.h>
#include "dpd.h"

extern "C" {
	
int dpd_buf4_mat_irrep_row_close(dpdbuf4 *Buf, int irrep)
{
  int all_buf_irrep;
  all_buf_irrep = Buf->file.my_irrep;

  if(Buf->params->coltot[irrep^all_buf_irrep])
    dpd_free_block(Buf->matrix[irrep], 1, Buf->params->coltot[irrep^all_buf_irrep]);

  return 0;
}

} /* extern "C" */