/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libciomr/libciomr.h>
#include "dpd.h"

namespace psi {

int dpd_buf4_mat_irrep_row_zero(dpdbuf4 *Buf, int irrep, int row)
{
  int coltot, all_buf_irrep;
  
  all_buf_irrep = Buf->file.my_irrep;
  coltot = Buf->params->coltot[irrep^all_buf_irrep];

  if(coltot)
      zero_arr(Buf->matrix[irrep][0], coltot);

  return 0;

}

}
