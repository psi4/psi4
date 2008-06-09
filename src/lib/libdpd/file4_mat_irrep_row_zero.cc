/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libciomr/libciomr.h>
#include "dpd.h"

namespace psi {

int dpd_file4_mat_irrep_row_zero(dpdfile4 *File, int irrep, int row)
{
  int coltot, my_irrep;

  if(File->incore) return 0;  /* Don't do this if the file is in core */

  my_irrep = File->my_irrep;
  
  coltot = File->params->coltot[irrep^my_irrep];

  if(coltot)
      zero_arr(File->matrix[irrep][0], coltot);

  return 0;

}

} // namespace psi
