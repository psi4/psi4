/*! \file 
    \ingroup (DPD)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <libpsio/psio.h>
#include <libqt/qt.h>
#include "dpd.h"

extern "C" {

int dpd_file4_mat_irrep_rd(dpdfile4 *File, int irrep)
{
  int rowtot, coltot, my_irrep;
  psio_address irrep_ptr, next_address;
  long int size;

  if(File->incore) return 0;  /* We already have this data in core */

  /* If the data doesn't actually exist on disk, we just leave */
  if(psio_tocscan(File->filenum, File->label) == NULL) return 1;

#ifdef DPD_TIMER
  timer_on("file4_rd");
#endif

  my_irrep = File->my_irrep;
  irrep_ptr = File->lfiles[irrep];
  rowtot = File->params->rowtot[irrep];
  coltot = File->params->coltot[irrep^my_irrep];
  size = ((long) rowtot) * ((long) coltot);

  if(rowtot && coltot)
     psio_read(File->filenum, File->label, (char *) File->matrix[irrep][0],
	       size*((long) sizeof(double)), irrep_ptr, &next_address);

#ifdef DPD_TIMER
  timer_off("file4_rd");
#endif

  return 0;

}

} /* extern "C" */
