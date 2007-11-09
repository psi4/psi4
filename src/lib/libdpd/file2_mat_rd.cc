/*! \file 
    \ingroup (DPD)
    \brief Enter brief description of file here 
*/
#include <libpsio/psio.h>
#include "dpd.h"

extern "C" {

int dpd_file2_mat_rd(dpdfile2 *File)
{
  int h, my_irrep, rowtot, coltot;
  psio_address irrep_ptr, next_address;

  my_irrep = File->my_irrep;

  if(File->incore) return 0; /* We already have this data in core */

  /* If data doesn't actually exist on disk, we just leave */
  if(psio_tocscan(File->filenum, File->label) == NULL) return 1;

  for(h=0; h < File->params->nirreps; h++) {
      irrep_ptr = File->lfiles[h];
      rowtot = File->params->rowtot[h];
      coltot = File->params->coltot[h^my_irrep];

      if(rowtot && coltot)
	  psio_read(File->filenum, File->label, (char *) File->matrix[h][0],
		    rowtot*coltot*sizeof(double), irrep_ptr, &next_address);
    }

  return 0;
}

} /* extern "C" */
