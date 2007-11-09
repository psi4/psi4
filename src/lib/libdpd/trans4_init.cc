/*! \file 
    \ingroup (DPD)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <stdlib.h>
#include <libciomr/libciomr.h>
#include "dpd.h"

extern "C" {

int dpd_trans4_init(dpdtrans4 *Trans, dpdbuf4 *Buf)
{
  int nirreps;

  nirreps = Buf->params->nirreps;

  /* Assign the input dpdbuf */
  Trans->buf = *Buf;

  Trans->matrix = (double ***) malloc(nirreps * sizeof(double **));

  /* Set up shifted matrix info */
  Trans->shift.shift_type = 0;
  Trans->shift.rowtot = init_int_matrix(nirreps, nirreps);
  Trans->shift.coltot = init_int_matrix(nirreps, nirreps);
  Trans->shift.matrix = (double ****) malloc(nirreps * sizeof(double ***));

  return 0;
}

} /* extern "C" */
