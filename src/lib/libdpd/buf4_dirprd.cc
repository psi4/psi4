/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libqt/qt.h>
#include "dpd.h"

namespace psi {
	
/* dpd_buf4_dirprd(): Computes the direct product between two dpd four-index
** buffers.
**
** Arguments:
**   dpdbuf4 *BufA, *BufB: Pointers to the dpd four-index buffers.
**  The results is written to FileB.
*/

int dpd_buf4_dirprd(dpdbuf4 *BufA, dpdbuf4 *BufB)
{
  int h, nirreps, my_irrep;

  nirreps = BufA->params->nirreps;
  my_irrep = BufA->file.my_irrep;

  for(h=0; h < nirreps; h++) {

      dpd_buf4_mat_irrep_init(BufA, h);
      dpd_buf4_mat_irrep_init(BufB, h);
      dpd_buf4_mat_irrep_rd(BufA, h);
      dpd_buf4_mat_irrep_rd(BufB, h);

      dirprd_block(BufA->matrix[h], BufB->matrix[h],
		   BufA->params->rowtot[h], BufA->params->coltot[h^my_irrep]);

      dpd_buf4_mat_irrep_wrt(BufB, h);
      dpd_buf4_mat_irrep_close(BufA, h);
      dpd_buf4_mat_irrep_close(BufB, h);
    }

  return 0;
}
      
} // namespace psi
