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

/*! \file    print.cc
    \ingroup optking
    \brief memory allocation
*/

#include "print.h"

namespace opt {


void print_matrix(const FILE *fp, double **A, const int nrow, const int ncol) {
  int i,j,col=0;

  //const int max_col = 12;
  const int max_col = 18;

  for (i=0; i<nrow; ++i) {
    for (j=0; j<ncol; ++j) {
      //fprintf(const_cast<FILE *>(fp), "%13.8f", A[i][j]);
      fprintf(const_cast<FILE *>(fp), "%10.6f", A[i][j]);
      ++col;
      if ((col == max_col) && (j != ncol-1)) {
        fprintf(const_cast<FILE *>(fp), "\n");
        col = 0;
      }
    }
    fprintf(const_cast<FILE *>(fp),"\n");
    col = 0;
  }
  fflush(const_cast<FILE *>(fp));
}

void print_array(const FILE *fp, double *A, const int ncol) {
  int j,col=0;

  const int max_col = 9;

  for (j=0; j<ncol; ++j) {
    fprintf(const_cast<FILE *>(fp), "%13.8f", A[j]);
    ++col;
    if ((col == max_col) && (j != ncol-1)) {
      fprintf(const_cast<FILE *>(fp), "\n");
      col = 0;
    }
  }
  fprintf(const_cast<FILE *>(fp),"\n");
  fflush(const_cast<FILE *>(fp));
}

void print_geom_array(const FILE *fp, double *A, const int natom) {
  int i;

  int cnt = -1;
  for (int i=0; i<natom; ++i) {
    fprintf(const_cast<FILE *>(fp), "\t%13.8f", A[++cnt]);
    fprintf(const_cast<FILE *>(fp), "%13.8f", A[++cnt]);
    fprintf(const_cast<FILE *>(fp), "%13.8f", A[++cnt]);
    fprintf(const_cast<FILE *>(fp), "\n");
  }
  fprintf(const_cast<FILE *>(fp),"\n");
}

}

