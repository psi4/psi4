/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libqt/qt.h>
#include "dpd.h"

namespace psi {

double dpd_file2_dot(dpdfile2 *FileA, dpdfile2 *FileB)
{
  int h, nirreps, my_irrep;
  double dot;

  nirreps = FileA->params->nirreps;
  my_irrep = FileA->my_irrep;

  dot = 0.0;

  dpd_file2_mat_init(FileA);
  dpd_file2_mat_init(FileB);
  dpd_file2_mat_rd(FileA);
  dpd_file2_mat_rd(FileB);

  for(h=0; h < nirreps; h++) {

      dot += dot_block(FileA->matrix[h], FileB->matrix[h],
		       FileA->params->rowtot[h],
		       FileA->params->coltot[h^my_irrep], 1.0);
      
    }

  dpd_file2_mat_close(FileA);
  dpd_file2_mat_close(FileB);
  
  return dot;

}

}
