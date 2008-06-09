/*!
** \file
** \brief This file includes the integer versions of several psi routines
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
    fprintf(stderr,"init_array:  trouble allocating memory \n");
    fprintf(stderr,"size = %d\n",size);
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
     fprintf(stderr,"init_int_matrix: trouble allocating memory \n"); 
     fprintf(stderr,"rows = %d\n", rows);
     exit(PSI_RETURN_FAILURE);
   }

   if ((array[0] = (int *) malloc (sizeof(int)*cols*rows))==NULL) {
     fprintf(stderr,"init_int_matrix: trouble allocating memory \n"); 
     fprintf(stderr,"rows = %d, cols = %d", rows, cols);
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
void print_int_mat(int **a, int m, int n, FILE *out)
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
      fprintf (out,"%8d",a[i][j]);
    }
  }
  fprintf (out,"\n");
  if (n <= kk) {
    fflush(out);
    return;
  }
  ii=kk; goto L200;
}

}

