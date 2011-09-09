/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <libpsio/psio.h>
#include "dpd.h"

namespace psi {

int dpd_file4_mat_irrep_wrt(dpdfile4 *File, int irrep)
{
  int rowtot, coltot, my_irrep;
  psio_address irrep_ptr, next_address;
  long int size;

  if(File->incore) {
    dpd_file4_cache_dirty(File);  /* Flag this cache entry for writing */
    return 0;  /* We're keeping this data in core */
  }

  my_irrep = File->my_irrep;
  irrep_ptr = File->lfiles[irrep];
  rowtot = File->params->rowtot[irrep];
  coltot = File->params->coltot[irrep^my_irrep];
  size = ((long) rowtot) * ((long) coltot);

  if(rowtot && coltot)
    psio_write(File->filenum, File->label, (char *) File->matrix[irrep][0],
	       size*((long) sizeof(double)), irrep_ptr, &next_address);

  return 0;

}

}
