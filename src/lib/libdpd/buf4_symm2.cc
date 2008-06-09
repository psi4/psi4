/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include "dpd.h"

namespace psi {

/* dpd_buf4_symm2(): Symmetrizes two dpdbuf4's by
** taking, I'(pq,rs) = 1/2 [I1(pq,rs) + I2(pq,rs)] (note the
** indices!).  Users should keep in mind that the first! buffer will
** be overwritten when this function is called.  Also note that this
** routine will NOT check to see if the row and column dimensions of
** the input buffers are identical, which is necessary for this to
** work.
**
** Arguments:
**   dpdbuf4 *Buf1: A pointer to the left dpdbuf4 to be symmetrized.
**   dpdbuf4 *Buf2: A pointer to the right dpdbuf4 to be symmetrized.  */

int dpd_buf4_symm2(dpdbuf4 *Buf1, dpdbuf4 *Buf2)
{
  int h, row, col, all_buf_irrep;
  double value;

  all_buf_irrep = Buf1->file.my_irrep;

  for(h=0; h < Buf1->params->nirreps; h++) {
      dpd_buf4_mat_irrep_init(Buf1, h);
      dpd_buf4_mat_irrep_rd(Buf1, h);

      dpd_buf4_mat_irrep_init(Buf2, h);
      dpd_buf4_mat_irrep_rd(Buf2, h);

      for(row=0; row < Buf1->params->rowtot[h]; row++)
          for(col=0; col < Buf1->params->coltot[h^all_buf_irrep]; col++) {
              value = 0.5*(Buf1->matrix[h][row][col]+Buf2->matrix[h][col][row]);
              Buf1->matrix[h][row][col] = value;
            }

      dpd_buf4_mat_irrep_wrt(Buf1, h);
      dpd_buf4_mat_irrep_close(Buf1, h);
      dpd_buf4_mat_irrep_close(Buf2, h);
    }

  return 0;
}


} // namespace psi
