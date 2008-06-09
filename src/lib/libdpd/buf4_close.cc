/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include "dpd.h"

namespace psi {
	
/* dpd_buf4_close(): Closes a dpd four-index buffer.
**
** Arguments:
**   dpdbuf4 *Buf: A pointer to the dpdbuf4 to be closed.
*/

int dpd_buf4_close(dpdbuf4 *Buf)
{
  int nirreps;

  nirreps = Buf->params->nirreps;
  
  dpd_file4_close(&(Buf->file));

  free(Buf->matrix);

  free_int_matrix(Buf->shift.rowtot);
  free_int_matrix(Buf->shift.coltot);

  free_int_matrix(Buf->row_offset);
  free_int_matrix(Buf->col_offset);

  free(Buf->shift.matrix);

  return 0;
}

} // namespace psi
