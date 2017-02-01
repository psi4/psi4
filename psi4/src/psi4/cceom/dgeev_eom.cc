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

/*! \file
    \ingroup CCEOM
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
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
    outfile->Printf("\nSum of complex part of HBar eigenvalues %20.15f, %10.2e\n",
      tval, eom_params.complex_tol);

    /*    exit(1); */
  }
  free(evals_i);
  free_block(left_evects);
  return;
}

}} // namespace psi::cceom
