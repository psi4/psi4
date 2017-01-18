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

/*! \file
**  \ingroup DETCI
**  \brief Contains C code for vector operations
**
** C. David Sherrill
** Center for Computational Quantum Chemistry
** University of Georgia
** 
*/

#include <cstdio>

namespace psi { namespace detci {

/*
** xey
** 
** Perform the operation X[] = Y[] for vectors 'x' and 'y'
** of length 'size'
**
*/
void xey(double *x, double *y, int size) {
   int i;

   for (i=0; i<size; i++) {
      x[i] = y[i];
      }
}

/*
** xeay
**
** Perform the operation X[] = a * Y[] for vectors 'x' and 'y' 
**   (of length 'size') and constant 'a'.
**
** David Sherrill, November 1995
*/
void xeay(double *x, double a, double *y, int size)
{
   int i;

   for (i=0; i<size; i++) {
      x[i] = a * y[i];
      }
}



/*
** xpeay
**
** Perform the operation X[] += A * Y[] for vectors 'x' and 'y' 
**   (of length 'size') and constant 'a'.
**
** David Sherrill, November 1995
*/
void xpeay(double *x, double a, double *y, int size)
{
   int i;

   for (i=0; i<size; i++) {
      x[i] += a * y[i];
      }
}

/*
** xeax
**
** Perform the operation X[] = A * X[] for vector 'X' and constant 'A'
**
** David Sherrill, February 1996
**
*/
void xeax(double *x, double a, int size)
{
   int i;

   for (i=0; i<size; i++) {
      x[i] *= a;
      }
}


/*
** xeaxmy
**
** Perform X[] = A * X[] - Y[]
**
** David Sherrill, February 1996
**
*/
void xeaxmy(double *x, double *y, double a, int size)
{
   int i;

   for (i=0; i<size; i++) {
      x[i] = x[i] * a - y[i];
      }
}

/*
** xeaxpby
**
** Perform X[] = A * X[] - B * Y[]
** 
** David Sherrill, March 1996
**
*/
void xeaxpby(double *x, double *y, double a, double b, int size)
{
   int i;

   for (i=0; i<size; i++) {
      x[i] = a * x[i] + b * y[i];
      }
}
/*
** xexy
**
** Perform X[] = X[] * Y[]
**
** Matt Leininger, September 1998
**
*/
void xexy(double *x, double *y, int size)
{
  int i;
  
  for (i=0; i<size; i++) {
     x[i] *= y[i];
     }
}

/*
** xexmy
**
** Perform X[] = X[] - Y[]
**
** Matt Leininger and Nick Petraco, February 1999
**
*/
void xexmy(double *x, double *y, int size)
{
   int i;

   for (i=0; i<size; i++) {
      x[i] -= y[i];
      }
}
 
/*
** xpey
**
** Perform X[] = X[] + Y[]
**
** Matt Leininger February 1999
**
*/
void xpey(double *x, double *y, int size)
{
   int i;

   for (i=0; i<size; i++) {
      x[i] += y[i];
      }
}

}} // namespace psi::detci
