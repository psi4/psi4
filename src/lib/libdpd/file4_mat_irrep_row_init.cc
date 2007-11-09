/*! \file 
    \ingroup (DPD)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "dpd.h"

extern "C" {

int dpd_file4_mat_irrep_row_init(dpdfile4 *File, int irrep)
{
  int my_irrep;

  if(File->incore) return 0;  /* We already have the whole matrix in core */

#ifdef DPD_TIMER
  timer_on("f4_rowinit");
#endif

  my_irrep = File->my_irrep;
  
  File->matrix[irrep] = dpd_block_matrix(1, File->params->coltot[irrep^my_irrep]);

#ifdef DPD_TIMER
  timer_off("f4_rowinit");
#endif

  return 0;
}

} /* extern "C" */
