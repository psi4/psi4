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
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

/** AEL() computes the approximate excitation level according to
 ** Stanton and Bartlett, JCP, 98, 1993, 7034.
 ** Trace [rho(excited) - rho(ground)] = AEL
 ** where both densities are expressed in the basis that diagonalizes
 ** the ground-state CCSD density.
 ** I was never able to get these results to agree with those of JFS or
 ** the current ACES2 so I'm not going to use this right now.
 ** --RAK */

namespace psi { namespace ccdensity {

void ael(struct RHO_Params *rho_params)
{
  int dim,i,j,k;
  double **rho_g, *evals, **evects, **tmat, **rho_x, ael, **rho_diff, trace;

  dim = moinfo.nmo - moinfo.nfzv;
  rho_g = block_matrix(dim,dim);
  rho_x = block_matrix(dim,dim);
  rho_diff = block_matrix(dim,dim);
  evals = init_array(dim);
  evects = block_matrix(dim,dim);
  tmat = block_matrix(dim,dim);

  /* read and diagonalize the ground-state rho */
  psio_open(PSIF_MO_OPDM, PSIO_OPEN_OLD);
  psio_read_entry(PSIF_MO_OPDM, rho_params[0].opdm_lbl, (char *) &(rho_g[0][0]), sizeof(double)*dim*dim);
  psio_close(PSIF_MO_OPDM, 1);

  sq_rsp(dim, dim, rho_g, evals, 3, evects, 1.0E-14);
  C_DGEMM('t', 'n', dim, dim, dim, 1.0, &(evects[0][0]), dim, &(rho_g[0][0]), dim, 0.0, &(tmat[0][0]), dim);
  C_DGEMM('n', 'n', dim, dim, dim, 1.0,   &(tmat[0][0]), dim, &(evects[0][0]), dim, 0.0, &(rho_g[0][0]), dim);

  for (i=1; i<params.nstates; ++i) {
    /* read in the excited state density */
    psio_open(PSIF_MO_OPDM, PSIO_OPEN_OLD);
    psio_read_entry(PSIF_MO_OPDM, rho_params[i].opdm_lbl, (char *) &(rho_x[0][0]),
		     sizeof(double)*dim*dim);
    psio_close(PSIF_MO_OPDM, 1);

    /* transform the excited-state density */
    C_DGEMM('t', 'n', dim, dim, dim, 1.0, &(evects[0][0]), dim, &(rho_x[0][0]), dim, 0.0, &(tmat[0][0]), dim);
    C_DGEMM('n', 'n', dim, dim, dim, 1.0,   &(tmat[0][0]), dim, &(evects[0][0]), dim, 0.0, &(rho_x[0][0]), dim);

    /* compute the ith AEL */
    ael = 0.0;
    for (j=0; j<dim; ++j) {
      ael += 0.5 * fabs( rho_x[j][k] - rho_g[j][k] );
    }
    outfile->Printf( "\tAEL %d: %10.7lf\n", i, ael);
  }

  free(evals);
  free_block(evects);
  free_block(tmat);
  free_block(rho_g);
  free_block(rho_x);
  free_block(rho_diff);
  return;
}

}} // namespace psi::ccdensity
