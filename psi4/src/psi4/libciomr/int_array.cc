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
** \file
** \brief This file includes the integer versions of several psi routines
** for handling arrays and matrices of doubles
**
** David Sherrill, 1996
**
** \ingroup CIOMR
*/

#include "psi4/psifiles.h"
#include <cstdio>
#include <cstdlib>
#include <strings.h>
#include "psi4/psi4-dec.h"
#include "psi4/libparallel/ParallelPrinter.h"
namespace psi {

/*!
** init_int_array(): Allocates memory for one-D array of ints of dimension
** 'size' and returns pointer to 1st element.  Zeroes all elements.
**
** Just modified the init_array() routine to do int's instead.
** This will avoid the temptation to allocate 5 integers by
**    p = (int *) init_array(5/2), which is bad.
**
** \param size = length of array to allocate
**
** Returns: pointer to new array
**
** C. David Sherrill
** \ingroup CIOMR
*/
int * init_int_array(int size)
{
  int *array;

  if ((array = (int *) malloc(sizeof(int)*size))==NULL) {
    outfile->Printf("init_array:  trouble allocating memory \n");
    outfile->Printf("size = %d\n",size);
    exit(PSI_RETURN_FAILURE);
  }
  bzero(array,sizeof(int)*size);
  return(array);
}


/*!
** zero_int_array()
** Zeroes out an array of integers 'size' integers long
**
** \param a    = integer array to zero out
** \param size = number of elements in a to zero
**
** Returns: none
**
** \ingroup CIOMR
*/
void zero_int_array(int *a, int size)
{
   bzero(a,sizeof(int)*size);
}


/*!
** init_int_matrix():
** Function initializes (allocates and clears) a matrix of integers with
** dimensions 'rows' by 'cols' and returns a pointer to it (ptr to first
** row ptr). The matrix layout is blocked, i.e. like produced by block_matrix()
**
** \param rows = number of rows
** \param cols = number of columns
**
** Returns: pointer to first row of newly-allocated integer block matrix
**
** \ingroup CIOMR
*/
int **init_int_matrix(int rows, int cols)
{
   int **array=NULL;
   int i;

   if ((array = (int **) malloc(sizeof(int *)*rows))==NULL) {
     outfile->Printf("init_int_matrix: trouble allocating memory \n");
     outfile->Printf("rows = %d\n", rows);
     exit(PSI_RETURN_FAILURE);
   }

   if ((array[0] = (int *) malloc (sizeof(int)*cols*rows))==NULL) {
     outfile->Printf("init_int_matrix: trouble allocating memory \n");
     outfile->Printf("rows = %d, cols = %d", rows, cols);
     exit(PSI_RETURN_FAILURE) ;
   }
   for (i=1; i<rows; i++) {
     array[i] = array[i-1] + cols;
   }
   bzero(array[0], sizeof(int)*cols*rows);

   return array;
}


/*!
** free_int_matrix():
** Free a matrix of integers.  Pass a pointer to the matrix.
**
** \param array = pointer to integer matrix
**
** \ingroup CIOMR
*/
void free_int_matrix(int **array)
{
  free(array[0]);
  free(array);
}


/*!
** zero_int_matrix():
** Zero a matrix of integers.  Pass the matrix, the number of rows,
** and the number of columns.
**
** \param array = pointer to integer matrix
** \param rows  = number of rows in matrix
** \param cols  = number of columns in matrix
**
** Returns: none
**
** \ingroup CIOMR
*/
void zero_int_matrix(int **array, int rows, int cols)
{
  zero_int_array(array[0], rows*cols);
}


/*!
** print_int_mat():
** Print a matrix of integers.  Pass the matrix, the number of rows and
** columns, and the output file pointer.
**
** \param a   = integer matrix to print
** \param m   = number of rows in matrix
** \param n   = number of columns in matrix
** \param out = FILE pointer to output file
**
** Returns: none
**
** \ingroup CIOMR
*/
void print_int_mat(int **a, int m, int n, std::string out)
{
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
         std::shared_ptr<OutFile>(new OutFile(out)));
   int ii,jj,kk,nn,ll;
  int i,j;

  ii=0;jj=0;
L200:
  ii++;
  jj++;
  kk=10*jj;
  nn=n;
  if (nn > kk) nn=kk;
  ll = 2*(nn-ii+1)+1;
  printer->Printf("\n   ");
  for (i=ii; i <= nn; i++) printer->Printf("   %5d",i);
  printer->Printf("\n");
  for (i=0; i < m; i++) {
    printer->Printf("\n%5d",i+1);
    for (j=ii-1; j < nn; j++) {
      printer->Printf("%8d",a[i][j]);
    }
  }
  printer->Printf("\n");
  if (n <= kk) {
    return;
  }
  ii=kk; goto L200;
}

}
