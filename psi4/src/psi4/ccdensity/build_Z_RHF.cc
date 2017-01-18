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
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <stdlib.h>
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/libciomr/libciomr.h"
#include <math.h>
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

/* build_Z_RHF():  Solve the orbital Z-vector equations for RHF refs:
**
**    sum E,M A(AI,EM) D(orb)(E,M) = -X(A,I)
**
** where A(AI,EM) is the orbital Hessian computed in build_A(), X(A,I)
** is the orbital rotation gradient computed in build_X(), and
** D(orb)(E,M) is the final Z-vector we want.
**
*/

void build_Z_RHF(void)
{
  dpdbuf4 A;
  dpdfile2 X1, D;
  double *X;
  int h, nirreps, a, i, count;

  nirreps = moinfo.nirreps;

  /* Grab only irrep 0 of the orbital Hessian */
  global_dpd_->buf4_init(&A, PSIF_CC_MISC, 0, 11, 11, 11, 11, 0, "A(EM,AI)");
  global_dpd_->buf4_mat_irrep_init(&A, 0);
  global_dpd_->buf4_mat_irrep_rd(&A, 0);

  /* Place all the elements of the orbital rotation gradient, X into a
     linear array, Z */
  global_dpd_->file2_init(&X1, PSIF_CC_OEI, 0, 1, 0, "XAI");
  global_dpd_->file2_mat_init(&X1);
  global_dpd_->file2_mat_rd(&X1);
  X = init_array(A.params->rowtot[0]);
  for(h=0,count=0; h < nirreps; h++)
    for(a=0; a < X1.params->rowtot[h]; a++)
      for(i=0; i < X1.params->coltot[h]; i++)
	X[count++] = -X1.matrix[h][a][i];

  global_dpd_->file2_mat_close(&X1);
  global_dpd_->file2_close(&X1);

  /* Trying out Matt's Pople code --- way to go, Matt! */
  pople(A.matrix[0], X, A.params->rowtot[0], 1, 1e-12, "outfile", 0);

  /* Build the orbital component of Dai */
  global_dpd_->file2_init(&D, PSIF_CC_OEI, 0, 1, 0, "D(orb)(A,I)");
  global_dpd_->file2_mat_init(&D);
  for(h=0,count=0; h < nirreps; h++)
    for(a=0; a < D.params->rowtot[h]; a++)
      for(i=0; i < D.params->coltot[h]; i++)
	D.matrix[h][a][i] = X[count++];
  global_dpd_->file2_mat_wrt(&D);
  global_dpd_->file2_mat_close(&D);
  global_dpd_->file2_close(&D);

  free(X);

  global_dpd_->buf4_mat_irrep_close(&A, 0);
  global_dpd_->buf4_close(&A);
}



}} // namespace psi::ccdensity
