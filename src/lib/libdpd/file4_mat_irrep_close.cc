/*! \file 
    \ingroup (DPD)
    \brief Enter brief description of file here 
*/
#include <stdio.h>
#include <libciomr/libciomr.h>
#include "dpd.h"
#define EXTERN
#include "dpd.gbl"

extern "C" {

/* dpd_file4_mat_irrep_close(): Releases memory for a matrix for a
** single irrep of a dpd four-index file.
**
** Arguments:
**   dpdfile4 *File: A pointer to the input dpdfile.
**   int irrep: The irrep number to be freed.
*/

int dpd_file4_mat_irrep_close(dpdfile4 *File, int irrep)
{
  int my_irrep, rowtot, coltot;
  long int size;

  my_irrep = File->my_irrep;

  rowtot = File->params->rowtot[irrep];
  coltot = File->params->coltot[irrep^my_irrep];
  size = ((long) rowtot) * ((long) coltot);

  if(File->incore) return 0;  /* We need to keep the memory */

  if(size) dpd_free_block(File->matrix[irrep], rowtot, coltot);

  return 0;
}

} /* extern "C" */
