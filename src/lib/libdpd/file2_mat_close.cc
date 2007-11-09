/*! \file 
    \ingroup (DPD)
    \brief Enter brief description of file here 
*/
#include "dpd.h"

extern "C" {

int dpd_file2_mat_close(dpdfile2 *File)
{
  int h, my_irrep;

  my_irrep = File->my_irrep;

  if(File->incore) return 0;  /* We need to keep the memory */

  for(h=0; h < File->params->nirreps; h++)
    if(File->params->rowtot[h] && File->params->coltot[h^my_irrep])
      dpd_free_block(File->matrix[h], File->params->rowtot[h],
          File->params->coltot[h^my_irrep]);


  return 0;
}

} /* extern "C" */
