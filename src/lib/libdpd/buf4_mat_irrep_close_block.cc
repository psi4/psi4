/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

namespace psi {
	
/* dpd_buf4_mat_irrep_close_block(): Releases memory for a subblock of
** a matrix for a single irrep of a dpd four-index buffer.
**
** Arguments:
**   dpdbuf4 *Buf: A pointer to the input dpdbuf4.
**   int irrep: The irrep number to be freed.
**
** Note that shift information is freed here as well.
**
*/

int dpd_buf4_mat_irrep_close_block(dpdbuf4 *Buf, int irrep, int num_pq)
{
  int h, nirreps, all_buf_irrep;

  nirreps = Buf->params->nirreps;
  all_buf_irrep = Buf->file.my_irrep;

  /* Free the shift structure for this irrep if used */
  if(Buf->shift.shift_type) {
      for(h=0; h < nirreps; h++) 
	  if(Buf->shift.rowtot[irrep][h])
	      free(Buf->shift.matrix[irrep][h]);
      free(Buf->shift.matrix[irrep]);
      Buf->shift.shift_type = 0;
    }

  if(num_pq && Buf->params->coltot[irrep^all_buf_irrep])
    dpd_free_block(Buf->matrix[irrep], num_pq, Buf->params->coltot[irrep^all_buf_irrep]);

  return 0;
}

}

