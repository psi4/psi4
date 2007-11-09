/*! \file 
    \ingroup (DPD)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <stdlib.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "dpd.h"

extern "C" {

/* buf_block is the block of buffer to shift (also the row-irrep) */

int dpd_trans4_mat_irrep_shift13(dpdtrans4 *Trans, int buf_block)
{
  int h, i, nirreps, all_buf_irrep;
  int *count;
  long int *dataoff;
  int rowtot, coltot;
  double *data;

#ifdef DPD_TIMER
  timer_on("shift");
#endif

  all_buf_irrep = Trans->buf.file.my_irrep;

  if(Trans->shift.shift_type) {
      fprintf(stderr, "\n\tShift is already on! %d\n",
	      Trans->shift.shift_type);
      exit(PSI_RETURN_FAILURE);
    }
  else Trans->shift.shift_type = 13;

  nirreps = Trans->buf.params->nirreps;
  rowtot = Trans->buf.params->coltot[buf_block^all_buf_irrep];
  coltot = Trans->buf.params->rowtot[buf_block];
  if (rowtot == 0 || coltot == 0) data = 0;
  else data = Trans->matrix[buf_block][0];

  /* Calculate row and column dimensions of each new sub-block */
  for(h=0; h < nirreps; h++) {
      Trans->shift.rowtot[buf_block][h] = Trans->buf.params->rpi[h];
      Trans->shift.coltot[buf_block][h] =
        coltot * Trans->buf.params->spi[h^buf_block^all_buf_irrep];
    }

  /* Malloc the pointers to the rows for the shifted access matrix */
  Trans->shift.matrix[buf_block] = (double ***) malloc(nirreps * sizeof(double **));
  for(h=0; h < nirreps; h++)
      Trans->shift.matrix[buf_block][h] =
	   ((!Trans->shift.rowtot[buf_block][h]) ? NULL :
	    (double **) malloc(Trans->shift.rowtot[buf_block][h] * sizeof(double *)));

  /* Calculate the data offset */
  dataoff = init_long_int_array(nirreps);
  dataoff[0] = 0;
  for(h=1; h < nirreps; h++)
      dataoff[h] = dataoff[h-1] +
		   ((long) Trans->shift.rowtot[buf_block][h-1]) *
		   ((long) Trans->shift.coltot[buf_block][h-1]);
		     

  /* The row counter for each sub-block */
  count = init_int_array(nirreps);

  /* Loop over irreps of isolated index */
  for(h=0; h < nirreps; h++) {
      for(i=0; (i < Trans->shift.rowtot[buf_block][h]) &&
	    Trans->shift.coltot[buf_block][h]; i++,count[h]++) {
          Trans->shift.matrix[buf_block][h][count[h]] =
		    &(data[dataoff[h]+((long) Trans->shift.coltot[buf_block][h])*((long) i)]);
        }
    }

  free(count); free(dataoff);

#ifdef DPD_TIMER
  timer_off("shift");
#endif

  return 0;
}

} /* extern "C" */
