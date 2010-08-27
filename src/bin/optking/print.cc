/*! \file
    \ingroup OPT10
    \brief memory allocation
*/

#include "print.h"

namespace opt {


void print_matrix(const FILE *fp, double **A, const int nrow, const int ncol) {
  int i,j,col=0;

  const int max_col = 9;

  for (i=0; i<nrow; ++i) {
    for (j=0; j<ncol; ++j) {
      fprintf(const_cast<FILE *>(fp), "%13.8f", A[i][j]);
      ++col;
      if ((col == max_col) && (j != ncol-1)) {
        fprintf(const_cast<FILE *>(fp), "\n");
        col = 0;
      }
    }
    fprintf(const_cast<FILE *>(fp),"\n");
    col = 0;
  }
}

}

