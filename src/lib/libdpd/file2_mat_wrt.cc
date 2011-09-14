/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <libpsio/psio.h>
#include "dpd.h"

namespace psi {

int dpd_file2_mat_wrt(dpdfile2 *File)
{
  int h, my_irrep, rowtot, coltot;
  psio_address irrep_ptr, next_address;

  my_irrep = File->my_irrep;

  if(File->incore) {
    dpd_file2_cache_dirty(File); /* Flag this cache entry for writing */
    return 0;  /* We're keeping this data in core */
  }

  for(h=0; h < File->params->nirreps; h++) {
    irrep_ptr = File->lfiles[h];
    rowtot = File->params->rowtot[h];
    coltot = File->params->coltot[h^my_irrep];

    if(rowtot && coltot)
      psio_write(File->filenum, File->label, (char *) File->matrix[h][0],
		 rowtot*coltot*sizeof(double), irrep_ptr, &next_address);
  }

  return 0;
}

}
