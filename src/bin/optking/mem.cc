/*! \file
    \ingroup OPTKING
    \brief memory allocation/zero/free functions
*/

#define EXTERN
#include "globals.h"
#undef EXTERN

namespace psi { //namespace optking {

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

// free array of integers
void free_int_array(int *a) {
  free(a);
}

// get unit matrix of doubles
double **unit_matrix(long int m) {
  long int i;
  double **A = block_matrix(m, m);
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


}//} /* namespace psi::optking */

