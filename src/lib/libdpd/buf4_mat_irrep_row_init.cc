/*! \file 
    \ingroup (DPD)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <libqt/qt.h>
#include <libciomr/libciomr.h>
#include "dpd.h"

extern "C" {

int dpd_buf4_mat_irrep_row_init(dpdbuf4 *Buf, int irrep)
{
  int all_buf_irrep;
  all_buf_irrep = Buf->file.my_irrep;
#ifdef DPD_TIMER
  timer_on("b4_rowinit");
#endif
  Buf->matrix[irrep] = dpd_block_matrix(1, Buf->params->coltot[irrep^all_buf_irrep]);
#ifdef DPD_TIMER
  timer_off("b4_rowinit");
#endif

  return 0;
}

} /* extern "C" */
