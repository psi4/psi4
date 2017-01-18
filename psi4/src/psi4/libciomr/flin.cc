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

#include "libciomr.h"
#include <cstdlib>

namespace psi {

extern void ludcmp(double **, int, int *, double *);
extern void lubksb(double **, int, int *, double *);

/*!
** \file
** \brief Linear equation solver for A * x = b
** \ingroup CIOMR
*/ 


/*!
** flin(): solves linear equations A * x = b.
**
** \param a   = coefficient matrix
** \param b   = known vectors
** \param in  = dimension of a(in*in)
** \param im  = number of b vectors
** \param det = pointer to hold determinant of matrix a
**
** Returns: none
**
** \ingroup CIOMR
*/
void flin(double **a, double *b, int in, int im, double *det)
{
  int i,j,*indx;

  indx = (int *) init_array(in);

  ludcmp(a,in,indx,det);

  for (i=0; i < in ; i++) *det *= a[i][i];

  for (j=0; j<im; j++)
    lubksb(a,in,indx,b+j*in);

  free(indx);
}

}
