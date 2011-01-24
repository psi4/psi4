/*! \file
    \ingroup CCEOM
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

void cleanup(void);

/* This routine transposes matrices and calls lapack dgeev() in libqt *
 * to diagonalize a square nonsymmetric matrix.  The eigenvalues      *
 * are returned in random order.                                      */

void dgeev_eom(int L, double **G, double *lambda, double **alpha) {

  int i, j, lwork, info;
  double *evals_i, *work, **left_evects, tval, temp;

  evals_i = init_array(L); 
  left_evects = block_matrix(L,L);

  work = init_array(20*L);
  lwork = 20*L;          

  for (i=0; i<L; ++i)
    for (j=0; j<i; ++j) {
      temp = G[i][j];
      G[i][j] = G[j][i];
      G[j][i] = temp;
    }

  i = C_DGEEV('V', 'V', L, G[0], L, lambda, evals_i, left_evects[0],
    L, alpha[0], L, work, lwork);

  for (i=0; i<L; ++i)
    for (j=0; j<i; ++j) {
      temp = alpha[i][j];
      alpha[i][j] = alpha[j][i];
      alpha[j][i] = temp;
    }

  free(work);

  tval = 0.0;
  for (i=0; i<L; ++i) {
    tval += fabs(evals_i[i]);
  }
  if (tval > (eom_params.complex_tol)) {
    fprintf(outfile,"\nSum of complex part of HBar eigenvalues %20.15f, %10.2e\n",
      tval, eom_params.complex_tol);
    fflush(outfile);
    /*    exit(1); */
  }
  free(evals_i);
  free_block(left_evects);
  return;
}

}} // namespace psi::cceom
