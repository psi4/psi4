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
    \ingroup RESPONSE
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libdpd/dpd.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace response {

void invert_RPA_RHF(double omega)
{
  int h, h1, nirreps, row, col, dim;
  int Ga, Gi, i, a, ai, aa, ii;
  int lwork, *ipiv, info;
  double *work;
  dpdbuf4 A, B;
  double ***C;

  global_dpd_->buf4_init(&A, PSIF_MO_HESS, 0, 11, 11, 11, 11, 0, "A(AI,BJ)");
  global_dpd_->buf4_init(&B, PSIF_MO_HESS, 0, 11, 11, 11, 11, 0, "B(AI,BJ)");

  C = (double ***) malloc(moinfo.nirreps * sizeof(double **));
  moinfo.RPA_dim = init_int_array(moinfo.nirreps);

  for(h=0; h < moinfo.nirreps; h++) {

    dim = moinfo.RPA_dim[h] = A.params->rowtot[h];

    if(dim) {

      global_dpd_->buf4_mat_irrep_init(&A, h);
      global_dpd_->buf4_mat_irrep_rd(&A, h);
      global_dpd_->buf4_mat_irrep_init(&B, h);
      global_dpd_->buf4_mat_irrep_rd(&B, h);

      C[h] = block_matrix(2*dim, 2*dim);

      for(row=0; row < dim; row++) {
        for(col=0; col < dim; col++) {
          C[h][row][col] = 2 * A.matrix[h][row][col];
          C[h][row+dim][col+dim] = 2 *A.matrix[h][row][col];
          C[h][row][col+dim] = -2 * B.matrix[h][row][col];
          C[h][row+dim][col] = -2 * B.matrix[h][row][col];
        }

        C[h][row][row] += 2 * omega;
        C[h][row+dim][row+dim] -= 2 * omega;
      }

      ipiv = init_int_array(2*dim);
      lwork = 20 * 2*dim;
      work = init_array(lwork);
      info = C_DGETRF(2*dim, 2*dim, &(C[h][0][0]), 2*dim, ipiv);
      info = C_DGETRI(2*dim, &(C[h][0][0]), 2*dim, ipiv, work, lwork);
      if(info) {
        fprintf(outfile, "\n\tDGETRI failed. info = %d. Exiting.\n", info);
        exit(PSI_RETURN_FAILURE);
      }

      free(ipiv);
      free(work);
    }

    global_dpd_->buf4_mat_irrep_close(&A, h);
    global_dpd_->buf4_mat_irrep_close(&B, h);
  }

  global_dpd_->buf4_close(&A);
  global_dpd_->buf4_close(&B);

  moinfo.RPA_inv = C;
}

}} // namespace psi::response
