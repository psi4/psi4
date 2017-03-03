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

#include <cstdio>

/*!
  \file
  \brief read in a matrix from an input stream (deprecated)
  \ingroup QT
*/

namespace psi {
	
/*!
** MAT_IN(): Function to read in a matrix.  Simple version for now.
**
** Parameters:
**    \param fp         =  file pointer to input stream
**    \param array      =  matrix to hold data
**    \param width      =  number of columns to read
**    \param max_length =  maximum number of rows to read
**    \param stat       =  pointer to int to hold status flag 
**                         (0=read ok, 1=error)
**
** Returns: 
**    number of rows read
**    Also modifies stat to = error code (0 = ok, 1 = error)
** \ingroup QT
*/

int mat_in(FILE *fp, double **array, int width, int max_length, int *stat) 
{
   int i=0, j, errcod=0 ;
   int nr ;
   double data ;

   while ( (i < max_length) && (!errcod) ) {
      for (j=0; j<width; j++) {
         nr = fscanf(fp, "%lf", &data) ;
         if (feof(fp)) break ;
         if (nr != 1) {
            errcod = 1 ;
            break ;
            }
         else {
            array[i][j] = data ;
            }
         }
      if (feof(fp)) break ;
      if (!errcod) i++ ;
      }

   *stat = errcod ;
   return(i) ;
}

}
