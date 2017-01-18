/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

/*!
  \file
  \brief Print eigenvectors and eigenvalues to output file
  \ingroup CIOMR
*/

#include <cstdio>
#include "psi4/psi4-dec.h"
#include "psi4/libparallel/ParallelPrinter.h"
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
void eivout(double **a, double *b, int m, int n, std::string out)
{

   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
         std::shared_ptr<OutFile>(new OutFile(out)));
   int ii,jj,kk,nn;
  int i,j;

  ii=0;jj=0;
L200:
  ii++;
  jj++;
  kk=10*jj;
  nn=n;
  if (nn > kk) nn=kk;
  printer->Printf("\n");
  for (i=ii; i <= nn; i++) printer->Printf("       %5d",i);
  printer->Printf("\n");
  for (i=0; i < m; i++) {
    printer->Printf("\n%5d",i+1);
    for (j=ii-1; j < nn; j++) {
      printer->Printf("%12.7f",a[i][j]);
    }
  }
  printer->Printf("\n");
  printer->Printf("\n     ");
  for (j=ii-1; j < nn; j++) {
    printer->Printf("%12.7f",b[j]);
  }
  printer->Printf("\n");
  if (n <= kk) {

    return;
  }
  ii=kk; goto L200;
}

}
