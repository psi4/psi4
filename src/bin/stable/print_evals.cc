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
    \ingroup STABLE
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace stable {

void print_evals(double **evals, int *rank)
{
  int h, i;

  fprintf(outfile, "\t  #  ");
  for(h=0; h < moinfo.nirreps; h++)
      fprintf(outfile, "    %3s  ",moinfo.labels[h]);
  fprintf(outfile, "\n");
  fprintf(outfile, "\t---- ");
  for(h=0; h < moinfo.nirreps; h++)
      fprintf(outfile, "---------");
  fprintf(outfile, "\n");

  for(i=0; i < 5; i++) {
    fprintf(outfile, "\t %2d  ", i);
      for(h=0; h < moinfo.nirreps; h++) {
	  if(rank[h] <= i) fprintf(outfile, "         ");
	  else fprintf(outfile, " %7.4f ", evals[h][i]);
	}
      fprintf(outfile, "\n");
    }

  fprintf(outfile, "\n");
}

}} // namespace psi::stable
