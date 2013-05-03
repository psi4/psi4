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
** \file
** \brief This file includes the long integer versions of several psi routines
** for handling arrays and matrices of doubles 
**
** David Sherrill, 1996
**
** \ingroup CIOMR
*/

#include <psifiles.h>
#include <cstdio>
#include <cstdlib>
#include <strings.h>

namespace psi {

/*!
** init_long_int_array(): Allocates memory for one-D array of long ints of 
** dimension  'size' and returns pointer to 1st element.  Zeroes all elements.
**
** Just modified the init_int_array() routine to do long int's instead.
**
** Returns: pointer to new array
**
** C. David Sherrill
** \ingroup CIOMR
*/
long int * init_long_int_array(int size)
{
  long int *array;

  if ((array = (long int *) malloc(sizeof(long int)*size))==NULL) {
    fprintf(stderr,"init_array:  trouble allocating memory \n");
    fprintf(stderr,"size = %d\n",size);
    exit(PSI_RETURN_FAILURE);
  }
  bzero(array,sizeof(long int)*size);
  return(array);
}

}

