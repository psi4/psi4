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
    \ingroup TRANSQT2
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>

/*
** transone(): Transform a packed symmetric matrix.
**
** int m: input matrix row dimension
** int n: output matrix row dimension
** double *input: pointer to input integrals (the lower-triange of a symmetric matrix)
** double *output: pointer to output integrals (the lower-triangle of a symmetric matrix)
** double **C: transformation matrix (rectangular)
** int nc: column dimension of C
** int *order: reordering array for transformed indices
**
** Written for new transqt module
** TDC, 7/06
*/

namespace psi {
  namespace transqt2 {

#define INDEX(i,j) ((i>j) ? ((i*(i+1)/2)+j) : ((j*(j+1)/2)+i))

void transone(int m, int n, double *input, double *output, double **C, int nc, 
	      int *order)
{
  int p, q, pq, dim;
  double **TMP0, **TMP1;

  dim = (m > n) ? m : n;
  TMP0 = block_matrix(dim,dim);
  TMP1 = block_matrix(dim,dim);

  for(p=0,pq=0; p < m; p++)
    for(q=0; q <= p; q++,pq++) 
      TMP0[p][q] = TMP0[q][p] = input[pq];

  if(m && n) {
    C_DGEMM('n','n',m,n,m,1.0,TMP0[0],dim,C[0],nc,0.0,TMP1[0],dim);
    C_DGEMM('t','n',n,n,m,1.0,C[0],nc,TMP1[0],dim,0.0,TMP0[0],dim);
  }

  for(p=0; p < n; p++)
    for(q=0; q <= p; q++) {
      pq = INDEX(order[p],order[q]);
      output[pq] = TMP0[p][q];
    }

  free_block(TMP0);
  free_block(TMP1);
}

  } // namespace transqt2
} // namespace psi
