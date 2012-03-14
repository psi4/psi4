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
	
/* dpd_buf4_mat_irrep_close(): Releases memory for a matrix for a
** single irrep of a dpd four-index buffer.
**
** Arguments:
**   dpdbuf4 *Buf: A pointer to the input dpdbuf.
**   int irrep: The irrep number to be freed.
**
** Note that shift information is freed here as well.
*/

int dpd_buf4_mat_irrep_close(dpdbuf4 *Buf, int irrep)
{
  int h, nirreps, rowtot, coltot, my_irrep;
  long int size;

  my_irrep = Buf->file.my_irrep;
  rowtot = Buf->params->rowtot[irrep];
  coltot = Buf->params->coltot[irrep^my_irrep];

  size = ((long) rowtot) * ((long) coltot);

  nirreps = Buf->params->nirreps;

  /* Free the shift structure for this irrep if used */
  if(Buf->shift.shift_type) {
      for(h=0; h < nirreps; h++) 
	  if(Buf->shift.rowtot[irrep][h])
	      free(Buf->shift.matrix[irrep][h]);
      free(Buf->shift.matrix[irrep]);
      Buf->shift.shift_type = 0;
    }

  if(size) {
      /* If the file member is already in cache and its ordering is the 
         same as the buffer, then we just copied the pointer in
         buf4_mat_irrep_init(); don't free! */
      if(Buf->file.incore && !(Buf->anti) && 
          (Buf->params->pqnum == Buf->file.params->pqnum) &&
          (Buf->params->rsnum == Buf->file.params->rsnum))
          return 1; // The return value of this function is never checked
      else
          dpd_free_block(Buf->matrix[irrep], rowtot, coltot);
    }

  return 0;
}

}

