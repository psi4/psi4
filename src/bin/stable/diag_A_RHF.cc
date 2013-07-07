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

/*! \file
    \ingroup STABLE
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace stable {

void diag_A_RHF(void)
{
  int h, dim, i;
  double *eps, **v;
  dpdbuf4 A, B;

  global_dpd_->buf4_init(&A, PSIF_MO_HESS, 0, 11, 11, 11, 11, 0, "A(AI,BJ)");
  global_dpd_->buf4_init(&B, PSIF_MO_HESS, 0, 11, 11, 11, 11, 0, "A(AI,BJ) triplet");
  for(h=0; h < moinfo.nirreps; h++) {

    dim = A.params->rowtot[h];
    eps = init_array(dim);
    v = block_matrix(dim, dim);
    moinfo.rank[h] = dim;

    global_dpd_->buf4_mat_irrep_init(&A, h);
    global_dpd_->buf4_mat_irrep_rd(&A, h);
    sq_rsp(dim, dim, A.matrix[h], eps, 1, v, 1e-12);
    global_dpd_->buf4_mat_irrep_close(&A, h);

    for(i=0; i < MIN0(dim, 5); i++)
      moinfo.A_evals[h][i] = eps[i];

    zero_mat(v, dim, dim);
    zero_arr(eps, dim);

    global_dpd_->buf4_mat_irrep_init(&B, h);
    global_dpd_->buf4_mat_irrep_rd(&B, h);
    sq_rsp(dim, dim, B.matrix[h], eps, 1, v, 1e-12);
    global_dpd_->buf4_mat_irrep_close(&B, h);

    for(i=0; i < MIN0(dim, 5); i++)
      moinfo.A_triplet_evals[h][i] = eps[i];

    free_block(v);
    free(eps);

  }
  global_dpd_->buf4_close(&B);
  global_dpd_->buf4_close(&A);
}

}} // namespace psi::stable
