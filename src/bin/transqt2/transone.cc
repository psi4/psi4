/*! \file
    \ingroup TRANSQT2
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>

/*
** transone(): Transform a packed symmetric matrix.
**
** int m: input matrix row dimension
** int n: output matrix row dimension
** double *input: pointer to input integrals (the lower-triange of a symmetric matrix)
** double *output: pointer to output integrals (the lower-triangle of a symmetric matrix)
** double **C: transformation matrix (rectangular)
** int nc: column dimension of C
** int *order: reordering array for transformed indices
**
** Written for new transqt module
** TDC, 7/06
*/

namespace psi {
  namespace transqt2 {

#define INDEX(i,j) ((i>j) ? ((i*(i+1)/2)+j) : ((j*(j+1)/2)+i))

void transone(int m, int n, double *input, double *output, double **C, int nc, 
	      int *order)
{
  int p, q, pq, dim;
  double **TMP0, **TMP1;

  dim = (m > n) ? m : n;
  TMP0 = block_matrix(dim,dim);
  TMP1 = block_matrix(dim,dim);

  for(p=0,pq=0; p < m; p++)
    for(q=0; q <= p; q++,pq++) 
      TMP0[p][q] = TMP0[q][p] = input[pq];

  if(m && n) {
    C_DGEMM('n','n',m,n,m,1.0,TMP0[0],dim,C[0],nc,0.0,TMP1[0],dim);
    C_DGEMM('t','n',n,n,m,1.0,C[0],nc,TMP1[0],dim,0.0,TMP0[0],dim);
  }

  for(p=0; p < n; p++)
    for(q=0; q <= p; q++) {
      pq = INDEX(order[p],order[q]);
      output[pq] = TMP0[p][q];
    }

  free_block(TMP0);
  free_block(TMP1);
}

  } // namespace transqt2
} // namespace psi
