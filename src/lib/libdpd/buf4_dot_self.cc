/*! \file 
    \ingroup (DPD)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include "dpd.h"

extern "C" {
	
/* dpd_buf4_dot_self(): Evaluates the sum of the squares of the elements of a
** given dpdbuf4.
**
** Arguments:
**   dpdbuf4 *BufX: A pointer to the dpdbuf4.
*/

double dpd_buf4_dot_self(dpdbuf4 *BufX)
{
  int h, nirreps, my_irrep;
  int row, col;
  double alpha=0.0;

  nirreps = BufX->params->nirreps;
  my_irrep = BufX->file.my_irrep;

  for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(BufX, h);
      dpd_buf4_mat_irrep_rd(BufX, h);

      for(row=0; row < BufX->params->rowtot[h]; row++)
          for(col=0; col < BufX->params->coltot[h^my_irrep]; col++)
              alpha += BufX->matrix[h][row][col] * BufX->matrix[h][row][col];

      dpd_buf4_mat_irrep_close(BufX, h);
    }

  return alpha;
}

} /* extern "C" */