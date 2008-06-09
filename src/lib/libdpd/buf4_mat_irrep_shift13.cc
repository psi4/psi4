/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "dpd.h"

namespace psi {

int dpd_buf4_mat_irrep_shift13(dpdbuf4 *Buf, int buf_block)
{
  int h, i, nirreps, all_buf_irrep;
  int *count;
  int *dataoff;
  int rowtot, coltot;
  double *data;

#ifdef DPD_TIMER
  timer_on("shift");
#endif

  all_buf_irrep = Buf->file.my_irrep;

  if(Buf->shift.shift_type) {
      fprintf(stderr, "\n\tShift is already on! %d\n",
	      Buf->shift.shift_type);
      exit(PSI_RETURN_FAILURE);
    }
  else Buf->shift.shift_type = 13;
  
  Buf->shift.shift_type = 13;

  nirreps = Buf->params->nirreps;
  rowtot = Buf->params->rowtot[buf_block];
  coltot = Buf->params->coltot[buf_block^all_buf_irrep];
  if (rowtot == 0 || coltot == 0) data = 0;
  else data = Buf->matrix[buf_block][0];

  /* Calculate row and column dimensions of each new sub-block */
  for(h=0; h < nirreps; h++) {
      Buf->shift.rowtot[buf_block][h] = Buf->params->ppi[h];
      Buf->shift.coltot[buf_block][h] = coltot * Buf->params->qpi[h^buf_block];
    }

  /* Malloc the pointers to the rows for the shifted access matrix */
  Buf->shift.matrix[buf_block] = (double ***) malloc(nirreps * sizeof(double **));
  for(h=0; h < nirreps; h++)
      Buf->shift.matrix[buf_block][h] =
	   ((!Buf->shift.rowtot[buf_block][h]) ? NULL :
	    (double **) malloc(Buf->shift.rowtot[buf_block][h] * sizeof(double *)));

  /* Calculate the data offset */
  dataoff = init_int_array(nirreps);
  dataoff[0] = 0;
  for(h=1; h < nirreps; h++)
      dataoff[h] = dataoff[h-1] +
		   ((long) Buf->shift.rowtot[buf_block][h-1]) * 
	           ((long) Buf->shift.coltot[buf_block][h-1]);
		     

  /* The row counter for each sub-block */
  count = init_int_array(nirreps);

  /* Loop over irreps of isolated index */
  for(h=0; h < Buf->params->nirreps; h++) {
      for(i=0; (i < Buf->shift.rowtot[buf_block][h]) && Buf->shift.coltot[buf_block][h];
	  i++,count[h]++) {
          Buf->shift.matrix[buf_block][h][count[h]] =
		    &(data[dataoff[h]+((long) (Buf->shift.coltot[buf_block][h]))*((long) i)]);
        }
    }

  free(count); free(dataoff);

#ifdef DPD_TIMER
  timer_off("shift");
#endif

  return 0;
}

} // namespace psi
