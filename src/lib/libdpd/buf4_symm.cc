/*! \file 
    \ingroup (DPD)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include "dpd.h"

extern "C" {

/* dpd_buf4_symm(): Symmetrizes a dpdbuf4 by taking,
** I'(pq,rs) = 1/2 [I(pq,rs) + I(rs,pq)].  Users should keep in mind
** that the original buffer will be overwritten when this function is
** called.  Also note that this routine will NOT check to see if the
** row and column dimensions of the input buffer are identical, which
** is necessary for this to work.
**
** Arguments:
**   dpdbuf4 *Buf: A pointer to the dpdbuf4 to be symmetrized.
*/

int dpd_buf4_symm(dpdbuf4 *Buf)
{
  int h, row, col, all_buf_irrep;
  double value;

  all_buf_irrep = Buf->file.my_irrep;

  for(h=0; h < Buf->params->nirreps; h++) {
      dpd_buf4_mat_irrep_init(Buf, h);
      dpd_buf4_mat_irrep_rd(Buf, h);

      for(row=0; row < Buf->params->rowtot[h]; row++)
          for(col=0; col < Buf->params->coltot[h^all_buf_irrep]; col++) {
              value = 0.5*(Buf->matrix[h][row][col]+Buf->matrix[h][col][row]);
              Buf->matrix[h][row][col] = Buf->matrix[h][col][row] = value;
            }

      dpd_buf4_mat_irrep_wrt(Buf, h);
      dpd_buf4_mat_irrep_close(Buf, h);
    }

  return 0;
}


} /* extern "C" */
