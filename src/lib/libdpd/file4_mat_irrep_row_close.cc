/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libciomr/libciomr.h>
#include "dpd.h"

namespace psi {

int dpd_file4_mat_irrep_row_close(dpdfile4 *File, int irrep)
{
  int my_irrep;

  if(File->incore) return 0;  /* We're keeping the data in core */

  my_irrep = File->my_irrep;
  
  if(File->params->coltot[irrep])
    dpd_free_block(File->matrix[irrep],1,File->params->coltot[irrep^my_irrep]);

  return 0;
}

} // namespace psi
