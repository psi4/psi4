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
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <exception.h>
#include <psifiles.h>

namespace psi {
extern FILE* outfile;
namespace ccdensity {

void idx_error(const char *message, int p, int q, int r, int s, int pq, int rs,
               int pq_sym, int rs_sym, FILE *outfile)
{

  fprintf(outfile, "\n\tDPD Parameter Error in %s\n", message);
  fprintf(outfile,
          "\t-------------------------------------------------\n");
  fprintf(outfile,
          "\t    p      q      r      s  [   pq]  [   rs] pq_symm rs_symm\n");
  fprintf(outfile,"\t%5d  %5d  %5d  %5d  [%5d]  [%5d]   %1d   %1d\n", p,q,r,s,
          pq,rs,pq_sym,rs_sym);
  throw PsiException("ccdensity: error", __FILE__, __LINE__);
}


}} // namespace psi::ccdensity
