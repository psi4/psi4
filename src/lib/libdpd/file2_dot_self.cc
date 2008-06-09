/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include "dpd.h"

namespace psi {

/* dpd_file2_dot_self(): Evaluates the sum of the squares of the elements of a
** given dpdfile2.
*/

double dpd_file2_dot_self(dpdfile2 *BufX)
{
  int h, nirreps, my_irrep;
  int row, col;
  double alpha=0.0;

  nirreps = BufX->params->nirreps;
  my_irrep = BufX->my_irrep;

  dpd_file2_mat_init(BufX);
  dpd_file2_mat_rd(BufX);

  for(h=0; h < nirreps; h++) {

      for(row=0; row < BufX->params->rowtot[h]; row++)
          for(col=0; col < BufX->params->coltot[h^my_irrep]; col++)
              alpha += BufX->matrix[h][row][col] * BufX->matrix[h][row][col];

    }

  dpd_file2_mat_close(BufX);

  return alpha;

}

} // namespace psi
