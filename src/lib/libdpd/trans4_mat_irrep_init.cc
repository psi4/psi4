/*! \file 
    \ingroup (DPD)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <libciomr/libciomr.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

extern "C" {

int dpd_trans4_mat_irrep_init(dpdtrans4 *Trans, int irrep)
{
  int rowtot, coltot, all_buf_irrep;
  long int size;

  all_buf_irrep = Trans->buf.file.my_irrep;

  rowtot = Trans->buf.params->coltot[irrep^all_buf_irrep];
  coltot = Trans->buf.params->rowtot[irrep];
  size = ((long) rowtot) * ((long) coltot);

  if(size) Trans->matrix[irrep] = dpd_block_matrix(rowtot,coltot);
  

  return 0;
}

} /* extern "C" */
