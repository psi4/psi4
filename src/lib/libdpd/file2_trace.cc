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
#include "dpd.h"

namespace psi {

double dpd_file2_trace(dpdfile2 *InFile)
{
  int h, nirreps;
  int row, col;
  double trace;

  nirreps = InFile->params->nirreps;

  dpd_file2_mat_init(InFile);
  dpd_file2_mat_rd(InFile);

  trace = 0.0;
  for(h=0; h < nirreps; h++)
      for(row=0; row < InFile->params->rowtot[h]; row++)
          trace += InFile->matrix[h][row][row];

  dpd_file2_mat_close(InFile);

  return trace;
}


}
