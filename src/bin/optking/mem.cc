/*! \file
    \ingroup OPTKING
    \brief memory allocation/zero/free functions
*/

#define EXTERN
#include "globals.h"
#undef EXTERN

namespace psi { namespace optking {

// get array of doubles
double * init_array(long int size) {
  double *array = (double *) malloc(size*(long int)sizeof(double));
  if (array == NULL)
    punt("init_array: memory allocation error.");
  zero_array(array, size);
  return(array);
}

// zero array of doubles
void zero_array(double *a, long int size) {
  memset(a, '\0', sizeof(double)*size);
}

// free array of doubles
void free_array(double *a) {
  if (a != NULL) {
    free(a);
    a = NULL;
  }
}

// get array of integers
int * init_int_array(int size) {
  int *array = (int *) malloc(sizeof(int)*size);
  if (array == NULL) punt("init_int_array : trouble allocating memory");
  zero_int_array(array,size);
  return(array);
}

// zero array of integers
void zero_int_array(int *a, int size) {
   memset(a, '\0', sizeof(int)*size);
}

// free array of integers
void free_int_array(int *a) {
  free(a);
}

// get matrix of doubles
double **init_matrix(long int m, long int n) {
  double **A = NULL;
  double *B = NULL;
  long int i;

  if(m<=0 || n<=0) return((double **) NULL);

  A = (double **) malloc(m * (long int)sizeof(double *));
  B = (double *) malloc(m*n * (long int)sizeof(double));

  if ((A == NULL) || (B == NULL)) punt ("init_matrix : allocation error.");

  zero_array(B, m*n);

  for (i = 0; i < m; i++)
    A[i] = &(B[i*n]);

  return(A);
}

// get unit matrix of doubles
double **unit_matrix(long int m) {
  long int i;
  double **A = init_matrix(m, m);
  for (i=0; i<m; ++i)
    A[i][i] = 1.0;
  return A; 
}

// zero matrix of doubles
void zero_matrix(double **a, long int m, long int n) {
  long int i;
  for (i=0; i < m; i++)
    zero_array(a[i], n);
}

// free matrix of doubles
void free_matrix(double **array) {
  if(array == NULL) return;
  free(array[0]);
  free(array);
}

// get matrix of ints
int **init_int_matrix(int rows, int cols) {
   int **array = NULL;
   int i;

   array = (int **) malloc(sizeof(int *)*rows);
   if (array == NULL)
     punt("init_int_matrix : trouble allocating memory");

   array[0] = (int *) malloc (sizeof(int)*cols*rows);
   if (array[0] == NULL) 
     punt("init_int_matrix: trouble allocating memory");

   zero_int_array(array[0], cols*rows);

   for (i=1; i<rows; i++)
     array[i] = array[i-1] + cols;

   return array;
}

// zero matrix of ints
void zero_int_matrix(int **a, int n, int m) {
  int i;
  for (i=0; i < n; i++)
    zero_int_array(a[i],m);
}

// free matrix of ints
void free_int_matrix(int **array) {
  free(array[0]);
  free(array);
}

}} /* namespace psi::optking */

