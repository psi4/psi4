#include "dpd.h"

namespace psi {

double DPD::buf4_trace(dpdbuf4 *Buf)
{
  double trace = 0.0;
  for(int h=0; h < Buf->params->nirreps; h++) {
    buf4_mat_irrep_init(Buf, h);
    buf4_mat_irrep_rd(Buf, h);
    for(int row=0; row < Buf->params->rowtot[h]; row++) 
      trace += Buf->matrix[h][row][row];
    buf4_mat_irrep_close(Buf, h);
  }

  return trace;
}

}
