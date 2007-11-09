/*! \file 
    \ingroup (DPD)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <libqt/qt.h>
#include "dpd.h"

extern "C" {

int dpd_trans4_mat_irrep_rd(dpdtrans4 *Trans, int irrep)
{
  int pq, rs, all_buf_irrep;
  int rows, cols;
  dpdbuf4 *Buf;
  double *A, *B;

  Buf = &(Trans->buf);
  all_buf_irrep = Buf->file.my_irrep;

#ifdef DPD_TIMER
  timer_on("trans4_rw");
#endif

  /* Loop over rows of transpose */
  /*
    for(pq=0; pq < Trans->buf.params->coltot[irrep]; pq++) {
    for(rs=0; rs < Trans->buf.params->rowtot[irrep]; rs++) {
    Trans->matrix[irrep][pq][rs] = Buf->matrix[irrep][rs][pq];
    }
    }
  */
  /* Transpose using BLAS DCOPY */
  rows = Buf->params->rowtot[irrep];
  cols = Buf->params->coltot[irrep^all_buf_irrep];

  if(rows && cols) {
    for(rs=0; rs < cols; rs++) {
      A = &(Buf->matrix[irrep][0][rs]);
      B = &(Trans->matrix[irrep][rs][0]);
      C_DCOPY(rows, A, cols, B, 1);
    }
  }

#ifdef DPD_TIMER
  timer_off("trans4_rw");
#endif

  return 0;
}

} /* extern "C" */
