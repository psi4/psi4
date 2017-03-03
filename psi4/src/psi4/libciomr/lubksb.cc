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


namespace psi {

void lubksb(double** a,int n,int* indx,double* b)
   {
      int i,ii=0,ip,j;
      int t=0;
      double sum;

      for (i=0; i < n ; i++) {
         ip = indx[i];
         sum = b[ip];
         b[ip]=b[i];
         if(t) {
            for (j=ii; j <= i-1 ; j++) sum -= a[i][j]*b[j];
            }
         else if(sum) {
            ii=i;
            t++;
            }
         b[i]=sum;
         }
      for (i=n-1; i >= 0 ; i--) {
         sum = b[i];
         for (j=i+1; j < n ; j++) sum -= a[i][j]*b[j];
         b[i] = sum/a[i][i];
         }
      }

}
