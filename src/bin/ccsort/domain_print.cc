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
    \ingroup CCSORT
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include "Local.h"
#define EXTERN
#include "globals.h"
#include "psi4-dec.h"
namespace psi { namespace ccsort {

void domain_print(int nocc, int natom, int *domain_len, int **domain,
		  double *fR)
{

  int i, j, cnt, max;
  double domain_tot, domain_ave;

  max = 0;
  for(i=0; i < nocc; i++) 
    if(domain_len[i] > max) max = domain_len[i];

  outfile->Printf( "   Orbital  Domain");
  for(i=0; i < max-2; i++) outfile->Printf( "   "); /* formatting junk */
  outfile->Printf( "  Completeness\n");
  outfile->Printf( "   -------  ------");
  for(i=0; i < max-2; i++) outfile->Printf( "---"); /* more formatting junk */
  outfile->Printf( "  ------------\n");
  for(i=0; i < nocc; i++) {
    outfile->Printf( "      %2d    ",i);
    for(j=0,cnt=0; j < natom; j++) if(domain[i][j]) { outfile->Printf( " %2d", j); cnt++; }
    if(cnt < max) for(; cnt < max; cnt++) outfile->Printf( "   ");
    outfile->Printf( "     %7.5f\n", fR[i]);
  }
  domain_tot = 0;
  for(i=0; i < nocc; i++)
    domain_tot += domain_len[i];
  domain_ave = domain_tot/nocc;
  outfile->Printf( "\n   The average domain length is %4.2lf\n", domain_ave);
  
}

}} // namespace psi::ccsort
