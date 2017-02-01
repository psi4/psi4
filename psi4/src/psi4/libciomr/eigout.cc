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
  \brief Print eigenvectors and eigenvalues
  \ingroup CIOMR
*/

#include <cstdio>
#include "psi4/psi4-dec.h"
#include "psi4/libparallel/ParallelPrinter.h"
namespace psi {

/*!
** eigout(): Print out eigenvectors and eigenvalues.
**
** Prints an n x m matrix of eigenvectors.  Under each eigenvector,
** the corresponding elements of two arrays, b and c, will also be printed.
** This is useful for printing, for example, the SCF eigenvectors with
** their associated eigenvalues (orbital energies) and also the population.
**
** \param a    = matrix of eigenvectors (eigenvectors are columns)
** \param b    = first array to print under eigenvectors (e.g., eigenvalues)
** \param c    = second array to print under eigenvectors (e.g., populations)
** \param m    = number of rows in matrix a
** \param n    = number of columns in matrix a (and length of b and c)
** \param out = file pointer for output
**
** Returns: none
**
** \ingroup CIOMR
*/
void eigout(double **a, double *b, double *c, int m, int n, std::string out)
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
      printer->Printf("\n     ");
      for (j=ii-1; j < nn; j++) {
         printer->Printf("%12.7f",c[j]);
         }
      printer->Printf("\n");
      if (n <= kk) {
         return;
         }
      ii=kk; goto L200;
      }

}
