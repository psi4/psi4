/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include "dpd.h"

namespace psi {

int dpd_trans4_close(dpdtrans4 *Trans)
{
  int nirreps;

  nirreps = Trans->buf.params->nirreps;

  free(Trans->matrix);
  
  free_int_matrix(Trans->shift.rowtot);
  free_int_matrix(Trans->shift.coltot);
  free(Trans->shift.matrix);

  return 0;

}

} // namespace psi
