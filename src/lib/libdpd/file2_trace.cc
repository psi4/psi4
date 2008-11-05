/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include "dpd.h"

namespace psi {

double dpd_file2_trace(dpdfile2 *InFile)
{
  int h, nirreps;
  int row;
  double trace;

  nirreps = InFile->params->nirreps;

  dpd_file2_mat_init(InFile);
  dpd_file2_mat_rd(InFile);

  trace = 0.0;
  for(h=0; h < nirreps; h++)
      for(row=0; row < InFile->params->rowtot[h]; row++)
          trace += InFile->matrix[h][row][row];

  dpd_file2_mat_wrt(InFile);
  dpd_file2_mat_close(InFile);

  return trace;
}


} // namespace psi
