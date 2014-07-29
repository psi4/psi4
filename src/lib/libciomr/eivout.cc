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

/*!
  \file
  \brief Print eigenvectors and eigenvalues to output file
  \ingroup CIOMR
*/
 
#include <cstdio>
#include "psi4-dec.h"
namespace psi {

/*!
** eivout: Print out eigenvectors and eigenvalues to the output file
**
** \param a   = eigenvectors
** \param b   = eigenvalues
** \param m   = rows of a
** \param n   = columns of a
** \param out = output file pointer
**
** Returns: none
**
** \ingroup CIOMR
*/
void eivout(double **a, double *b, int m, int n, FILE *out)
{
  int ii,jj,kk,nn;
  int i,j;

  ii=0;jj=0;
L200:
  ii++;
  jj++;
  kk=10*jj;
  nn=n;
  if (nn > kk) nn=kk;
  psi::fprintf (out,"\n");
  for (i=ii; i <= nn; i++) psi::fprintf(out,"       %5d",i);
  psi::fprintf (out,"\n");
  for (i=0; i < m; i++) {
    psi::fprintf (out,"\n%5d",i+1);
    for (j=ii-1; j < nn; j++) {
      psi::fprintf (out,"%12.7f",a[i][j]);
    }
  }
  psi::fprintf (out,"\n");
  psi::fprintf (out,"\n     ");
  for (j=ii-1; j < nn; j++) {
    psi::fprintf(out,"%12.7f",b[j]);
  }
  psi::fprintf (out,"\n");
  if (n <= kk) {
    fflush(out);
    return;
  }
  ii=kk; goto L200;
}

}

