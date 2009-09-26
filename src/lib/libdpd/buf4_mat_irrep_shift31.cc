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

int dpd_buf4_mat_irrep_shift31(dpdbuf4 *Buf, int buf_block)
{
  int h, pq, Gr, Gs, r, nirreps, all_buf_irrep, h_pqr;
  int rowtot, coltot,cnt;
  int *count;
  int *blocklen, *rowoff;
  double *data;
  long int pqcol;

#ifdef DPD_TIMER
  timer_on("shift");
#endif

  all_buf_irrep = Buf->file.my_irrep;

  if(Buf->shift.shift_type) {
    fprintf(stderr, "\n\tShift is already on! %d\n",
        Buf->shift.shift_type);
    exit(PSI_RETURN_FAILURE);
  }
  else Buf->shift.shift_type = 31;

  nirreps = Buf->params->nirreps;
  rowtot = Buf->params->rowtot[buf_block];
  coltot = Buf->params->coltot[buf_block^all_buf_irrep];
  if (rowtot == 0 || coltot == 0) data = 0;
  else data = Buf->matrix[buf_block][0];

  /* Calculate row and column dimensions of each new sub-block */
  /* loop over h_pqr */
  for(h_pqr=0; h_pqr < nirreps; h_pqr++) {
    Buf->shift.rowtot[buf_block][h_pqr] = rowtot * Buf->params->rpi[h_pqr^buf_block];
    Buf->shift.coltot[buf_block][h_pqr] = Buf->params->spi[h_pqr^all_buf_irrep];
  }

  /* Malloc the pointers to the rows for the shifted access matrix */
  Buf->shift.matrix[buf_block] = (double ***) malloc(nirreps*sizeof(double **));
  for(h_pqr=0; h_pqr < nirreps; h_pqr++) 
    Buf->shift.matrix[buf_block][h_pqr] =
      ((!Buf->shift.rowtot[buf_block][h_pqr]) ? NULL :
       (double **) malloc(Buf->shift.rowtot[buf_block][h_pqr] * sizeof(double *)));

  /* Calculate the row offsets */
  blocklen = init_int_array(nirreps);
  for(h_pqr=0; h_pqr < nirreps; h_pqr++)
      blocklen[h_pqr] = Buf->params->rpi[h_pqr^buf_block] *
                        Buf->params->spi[h_pqr^all_buf_irrep];

  rowoff = init_int_array(nirreps);
  cnt = 0;
  for (h=0;h<nirreps;++h) {  /* loop over Gr */
    h_pqr = buf_block^h;
    rowoff[h_pqr] = cnt;
    cnt += blocklen[h_pqr];
  }

  /* The row counter for each sub-block */
  count = init_int_array(nirreps);

  /* Loop over rows of original DPD matrix */
  for(pq=0; pq < Buf->params->rowtot[buf_block]; pq++) {
    pqcol = ((long) pq) * ((long) coltot);

    /* Loop over irreps of pqr */
    for(h_pqr=0; h_pqr < nirreps; h_pqr++) {
      Gr = h_pqr^buf_block;
      Gs = h_pqr^all_buf_irrep;

      /* Loop over orbitals in Gr */
      for(r=0; (r < Buf->params->rpi[Gr]) && Buf->params->spi[Gs]; r++) {

        /* Re-assign the row pointer */
        Buf->shift.matrix[buf_block][h_pqr][count[h_pqr]] =
          &(data[pqcol + rowoff[h_pqr] + (r * Buf->params->spi[Gs])]);

        count[h_pqr]++;

      }
    }
  }

  free(count); free(rowoff); free(blocklen);

#ifdef DPD_TIMER
  timer_off("shift");
#endif

  return 0;
}

}
