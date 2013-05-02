/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 *@END LICENSE
 */

/*! \file
    \ingroup DPD
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstring>
#include "dpd.h"

namespace psi {

int dpd_file2_copy(dpdfile2 *InFile, int outfilenum, const char *label)
{
  int h, row, col, my_irrep, rowtot, coltot;
  double ***matrix;
  dpdfile2 OutFile;

  my_irrep = InFile->my_irrep;

  dpd_file2_init(&OutFile, outfilenum, InFile->my_irrep, InFile->params->pnum,
		 InFile->params->qnum, label);

  dpd_file2_mat_init(InFile);
  dpd_file2_mat_rd(InFile);
  dpd_file2_mat_init(&OutFile);
  
  for(h=0; h < OutFile.params->nirreps; h++) {
      rowtot = OutFile.params->rowtot[h];
      coltot = OutFile.params->coltot[h^my_irrep];
      if(rowtot && coltot) 
         memcpy((void *) &(OutFile.matrix[h][0][0]),
                (const void *) &(InFile->matrix[h][0][0]),
                sizeof(double)*rowtot*coltot);
    }
  
  dpd_file2_mat_wrt(&OutFile);
  dpd_file2_mat_close(&OutFile);
  dpd_file2_mat_close(InFile);
  dpd_file2_close(&OutFile);

  return 0;
}

}
