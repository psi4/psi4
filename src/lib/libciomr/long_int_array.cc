/*!
** \file int_array.cc
** \ingroup (CIOMR)
**
** LONG_INT_ARRAY.C
** This file includes the long int versions of several psi routines
** for handling arrays and matrices of doubles 
**
** TDC, 2003
** based on int_array.c by David Sherrill, 1996
**
*/

#include <psifiles.h>
#include "includes.h"

extern "C" {

/*!
** init_long_int_array():
** Allocates memory for one-D array of long ints of dimension 'size'
** and returns pointer to 1st element.  Zeroes all elements.
**
** \param size = length of array to allocate
**
** Returns: pointer to new array
** \ingroup (CIOMR)
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


/*!
** zero_long_int_array()
** Zeroes out an array of long integers 'size' integers long
**
** \ingroup (CIOMR)
*/
void zero_long_int_array(long int *a, int size)
{
   bzero(a,sizeof(long int)*size) ;
}


/*!
** init_long_int_matrix():
** Function initializes (allocates and clears) a matrix of integers with 
** dimensions 'rows' by 'cols' and returns a pointer to it (ptr to first 
** row ptr). The matrix layout is like that produced by block_matrix().
**
** \ingroup (CIOMR)
*/
long int **init_long_int_matrix(int rows, int cols)
{
   long int **array=NULL ;
   int i ;

   if ((array = (long int **) malloc(sizeof(long int *)*rows))==NULL) {
      fprintf(stderr,"init_long_int_matrix: trouble allocating memory \n") ; 
      fprintf(stderr,"rows = %d\n", rows) ;
      exit(PSI_RETURN_FAILURE) ;
      }

   if ((array[0] = (long int *) malloc (sizeof(long int)*cols*rows))==NULL) {
	   fprintf(stderr,"init_long_int_matrix: trouble allocating memory \n") ; 
	   fprintf(stderr,"row = %d, cols = %d", i, cols) ;
	   exit(PSI_RETURN_FAILURE) ;
   }
   for (i=1; i<rows; i++) {
	   	array[i] = array[i-1] + cols;
   }
   bzero(array[0], sizeof(long int)*cols*rows) ;

   return array;
}


/*!
** free_long_int_matrix():
** Free a matrix of long integers.  Pass a pointer to the matrix.
** \ingroup (CIOMR)
*/
void free_long_int_matrix(long int **array)
{
	free(array[0]) ;
	free(array) ;
}


/*!
** zero_long_int_matrix():
** Zero a matrix of long integers.  Pass the matrix, the number of rows,
** and the number of columns.
** \ingroup (CIOMR)
*/
void zero_long_int_matrix(long int **array, int rows, int cols)
{
   zero_long_int_array(array[0], rows*cols);
}


/*!
** print_long_int_mat():
** Print a matrix of long integers.  Pass the matrix, the number of rows and
** columns, and the output file pointer.
** \ingroup (CIOMR)
*/
void print_long_int_mat(long int **a, int m, int n, FILE *out)
{
  int ii,jj,kk,nn,ll;
  int i,j,k;

  ii=0;jj=0;
L200:
  ii++;
  jj++;
  kk=10*jj;
  nn=n;
  if (nn > kk) nn=kk;
  ll = 2*(nn-ii+1)+1;
  fprintf (out,"\n   ");
  for (i=ii; i <= nn; i++) fprintf(out,"   %5d",i);
  fprintf (out,"\n");
  for (i=0; i < m; i++) {
    fprintf (out,"\n%5d",i+1);
    for (j=ii-1; j < nn; j++) {
      fprintf (out,"%16ld",a[i][j]);
    }
  }
  fprintf (out,"\n");
  if (n <= kk) {
    fflush(out);
    return;
  }
  ii=kk; goto L200;
}

} /* extern "C" */
