/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libqt/qt.h>
#include "dpd.h"

namespace psi {

/* file2_axpy(): Evaluates the standard operation a * X + Y -> Y for
 ** dpdfile2's.
 **
 ** Arguments:
 **   dpdfile2 *FileA: A pointer to the leftmost dpdfile2.
 **   dpdfile2 *FileB: A pointer to the rightmost (and target) dpdfile2.
 **   double alpha: The scalar prefactor in the multiplication.
 **   int transA: A boolean indicating that we should use the transpose of
 **               FileA
 */

int dpd_file2_axpy(dpdfile2 *FileA, dpdfile2 *FileB, double alpha, 
    int transA)
{
  int h, nirreps, my_irrep;
  int row, col;

  nirreps = FileA->params->nirreps;
  my_irrep = FileA->my_irrep;

  dpd_file2_mat_init(FileA);
  dpd_file2_mat_init(FileB);
  dpd_file2_mat_rd(FileA);
  dpd_file2_mat_rd(FileB);

  for(h=0; h < nirreps; h++) {

    if(!transA) {

      for(row=0; row < FileA->params->rowtot[h]; row++)
        for(col=0; col < FileA->params->coltot[h^my_irrep]; col++)
          FileB->matrix[h][row][col] += alpha*FileA->matrix[h][row][col];

    }
    else {
      for(row=0; row < FileB->params->rowtot[h]; row++)
        for(col=0; col < FileB->params->coltot[h^my_irrep]; col++)
          FileB->matrix[h][row][col] += alpha*FileA->matrix[h^my_irrep][col][row];
    }
  }

  dpd_file2_mat_wrt(FileB);
  dpd_file2_mat_close(FileA);
  dpd_file2_mat_close(FileB);

  return 0;
}


} // namespace psi
