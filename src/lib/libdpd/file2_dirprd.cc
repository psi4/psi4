/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libqt/qt.h>
#include "dpd.h"

namespace psi {

/* file2_dirprd(): Computes the direct product between two two-index dpd
** files.
**
** Arguments:
**   dpdfile2 *FileA, *FileB: Pointers to the two-index dpd files.
**  The result is written to FileB.
*/

int dpd_file2_dirprd(dpdfile2 *FileA, dpdfile2 *FileB)
{
  int h, nirreps, my_irrep;

  nirreps = FileA->params->nirreps;
  my_irrep = FileA->my_irrep;

  dpd_file2_mat_init(FileA);
  dpd_file2_mat_init(FileB);
  dpd_file2_mat_rd(FileA);
  dpd_file2_mat_rd(FileB);

  for(h=0; h < nirreps; h++) {

      dirprd_block(FileA->matrix[h], FileB->matrix[h],
		   FileA->params->rowtot[h], FileA->params->coltot[h^my_irrep]);

    }

  dpd_file2_mat_wrt(FileB);
  dpd_file2_mat_close(FileA);
  dpd_file2_mat_close(FileB);

  return 0;
}
      

}
