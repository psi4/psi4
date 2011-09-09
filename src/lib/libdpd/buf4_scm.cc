/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include "dpd.h"
#include <libqt/qt.h>
#include <libpsio/psio.h>
#define EXTERN
#include "dpd.gbl"

namespace psi {

/* dpd_buf4_scm(): Multiplies every element of a four-index dpdbuf by a scalar.
**
** Arguments:
**   dpdbuf4 *InBuf: A pointer to the dpdbuf.
**   double alpha: The scalar.
**
** NB: This function is sometimes called automatically by the contractXXX
** functions to zero-out data that does not yet exist on disk.  In such a
** case, since each symmetry block is handled separately, it is possible 
** for only the first (totally symmetric) block to exist on disk while the
** others have yet to be created.  In such cases, the buf4_mat_irrep_rd()
** request will fail in libpsio, because the correct TOC entry exists
** (created when the first symmetry block was written), but the length of 
** the TOC entry will be too short for the next symmetry block to be added.
** So, to avoid this problem, here we manually check in the beginning to
** see if the TOC entry exists and should be read before the multiplication.
**
** TDC
** June 2000
**
*/

int dpd_buf4_scm(dpdbuf4 *InBuf, double alpha)
{
  int pq;
  long int length, core, memoryd, core_total, rowtot, coltot, maxrows;
  int h, nirreps, new_buf4, all_buf_irrep;
  int row, col;
  int incore;
  double *X;

  nirreps = InBuf->params->nirreps;
  all_buf_irrep = InBuf->file.my_irrep;

#ifdef DPD_TIMER
  timer_on("buf4_scm");
#endif

  /* Look first for the TOC entry on disk */
  if(psio_tocscan(InBuf->file.filenum, InBuf->file.label) == NULL)
    new_buf4 = 1;
  else new_buf4 = 0;

  for(h=0; h < nirreps; h++) {

    memoryd = dpd_main.memory;
    incore = 1; /* default */

    /* Compute the core requirements for the straight contraction */
    core_total = 0;
    /** X terms **/
    coltot = InBuf->params->coltot[h^all_buf_irrep];
    if(coltot) {
      maxrows = DPD_BIGNUM/coltot;
      if(maxrows < 1) {
	fprintf(stderr, "\nLIBDPD Error: cannot compute even the number of rows in buf4_scm.\n");
	dpd_error("buf4_scm", stderr);
      }
    }
    else maxrows = DPD_BIGNUM;
    rowtot = InBuf->params->rowtot[h];
    for(; rowtot > maxrows; rowtot -= maxrows) {
      if(core_total > (core_total + maxrows*coltot)) incore = 0;
      else core_total += maxrows * coltot;
    }
    if(core_total > (core_total + rowtot*coltot)) incore = 0;
    core_total += rowtot * coltot;

    if(core_total > memoryd) incore = 0;

    if(incore) {

      dpd_buf4_mat_irrep_init(InBuf, h);

      if(!new_buf4) dpd_buf4_mat_irrep_rd(InBuf, h);

      length = ((long) InBuf->params->rowtot[h]) * ((long) InBuf->params->coltot[h^all_buf_irrep]);
      if(length) {
	X = &(InBuf->matrix[h][0][0]);
	C_DSCAL(length, alpha, X, 1);
      }

      dpd_buf4_mat_irrep_wrt(InBuf, h);
      dpd_buf4_mat_irrep_close(InBuf, h);
    }
    else {

      dpd_buf4_mat_irrep_row_init(InBuf, h);

      for(pq=0; pq < InBuf->params->rowtot[h]; pq++) {
	if(!new_buf4) dpd_buf4_mat_irrep_row_rd(InBuf, h, pq);

	length = InBuf->params->coltot[h^all_buf_irrep];

	if(length) {
	  X = &(InBuf->matrix[h][0][0]);
	  C_DSCAL(length, alpha, X, 1);
	}
	dpd_buf4_mat_irrep_row_wrt(InBuf, h, pq);
      }    
      dpd_buf4_mat_irrep_row_close(InBuf, h);
    }
  }

#ifdef DPD_TIMER
  timer_off("buf4_scm");
#endif

  return 0;
}

}
