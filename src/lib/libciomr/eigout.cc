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
  \brief Print eigenvectors and eigenvalues
  \ingroup CIOMR
*/

#include <cstdio>

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
void eigout(double **a, double *b, double *c, int m, int n, FILE *out)
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
      fprintf (out,"\n");
      for (i=ii; i <= nn; i++) fprintf(out,"       %5d",i);
      fprintf (out,"\n");
      for (i=0; i < m; i++) {
         fprintf (out,"\n%5d",i+1);
         for (j=ii-1; j < nn; j++) {
            fprintf (out,"%12.7f",a[i][j]);
            }
         }
      fprintf (out,"\n");
      fprintf (out,"\n     ");
      for (j=ii-1; j < nn; j++) {
         fprintf(out,"%12.7f",b[j]);
         }
      fprintf (out,"\n");
      fprintf (out,"\n     ");
      for (j=ii-1; j < nn; j++) {
         fprintf(out,"%12.7f",c[j]);
         }
      fprintf (out,"\n");
      if (n <= kk) {
         fflush(out);
         return;
         }
      ii=kk; goto L200;
      }

}
