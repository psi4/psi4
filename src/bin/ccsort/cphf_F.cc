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
    \ingroup CCSORT
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <libqt/qt.h>
#include <libdpd/dpd.h>
#include <libciomr/libciomr.h>
#include <psifiles.h>
#include "MOInfo.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccsort {

void cphf_F(const char *cart)
{
  int irrep, row, a, i, asym, isym, num_ai, info, *ipiv;
  double *vector;
  dpdbuf4 A;
  dpdfile2 mu;

  psio_open(PSIF_MO_HESS, 1);
  global_dpd_->buf4_init(&A, PSIF_MO_HESS, 0, 11, 11, 11, 11, 0, "A(AI,BJ)");

  if (!strcmp(cart,"X")) {
    irrep = moinfo.irrep_x;

    /* sort Mu elements into a single vector for lineq solver */
    global_dpd_->file2_init(&mu, PSIF_CC_OEI, irrep, 0, 1, "Mu_X_IA");
    global_dpd_->file2_mat_init(&mu);
    global_dpd_->file2_mat_rd(&mu);
    num_ai = A.params->rowtot[irrep];
    vector = init_array(num_ai);
    for(row=0; row < num_ai; row++) {
      i = A.params->roworb[irrep][row][0];
      a = A.params->roworb[irrep][row][1];
      isym = A.params->psym[i];
      asym = A.params->qsym[a];
      vector[row] = -mu.matrix[asym][a-A.params->qoff[asym]][i-A.params->poff[isym]];
    }
    global_dpd_->file2_mat_close(&mu);
    global_dpd_->file2_close(&mu);

    /* grab current irrep of MO Hessian */
    global_dpd_->buf4_mat_irrep_init(&A, irrep);
    global_dpd_->buf4_mat_irrep_rd(&A, irrep);


    /* solve CPHF equations */
    ipiv = init_int_array(num_ai);
    info = C_DGESV(num_ai, 1, &(A.matrix[irrep][0][0]), num_ai, ipiv, vector, num_ai);
    if(info) {
      psi::fprintf(outfile, "CCSORT: cphf_F: Error in C_DGESV.  Info = %d.  Cart = X. Exiting.\n", info);
      exit(PSI_RETURN_FAILURE);
    }
    free(ipiv);

    global_dpd_->buf4_mat_irrep_close(&A, irrep);

    /* sort CPHF solution to DPD format */
    global_dpd_->file2_init(&mu, PSIF_CC_OEI, irrep, 1, 0, "CPHF Uf_X_AI");
    global_dpd_->file2_mat_init(&mu);
    for(row=0; row < num_ai; row++) {
      a = A.params->roworb[irrep][row][0];
      i = A.params->roworb[irrep][row][1];
      asym = A.params->psym[a];
      isym = A.params->qsym[i];
      mu.matrix[asym][a-A.params->poff[asym]][i-A.params->qoff[isym]] = vector[row];
    }
    global_dpd_->file2_mat_wrt(&mu);
    global_dpd_->file2_close(&mu);
  }

  if (!strcmp(cart,"Y")) {
    irrep = moinfo.irrep_y;

    /* sort Mu elements into a single vector for lineq solver */
    global_dpd_->file2_init(&mu, PSIF_CC_OEI, irrep, 0, 1, "Mu_Y_IA");
    global_dpd_->file2_mat_init(&mu);
    global_dpd_->file2_mat_rd(&mu);
    num_ai = A.params->rowtot[irrep];
    vector = init_array(num_ai);
    for(row=0; row < num_ai; row++) {
      i = A.params->roworb[irrep][row][0];
      a = A.params->roworb[irrep][row][1];
      isym = A.params->psym[i];
      asym = A.params->qsym[a];
      vector[row] = -mu.matrix[asym][a-A.params->qoff[asym]][i-A.params->poff[isym]];
    }
    global_dpd_->file2_mat_close(&mu);
    global_dpd_->file2_close(&mu);

    /* grab current irrep of MO Hessian */
    global_dpd_->buf4_mat_irrep_init(&A, irrep);
    global_dpd_->buf4_mat_irrep_rd(&A, irrep);


    /* solve CPHF equations */
    ipiv = init_int_array(num_ai);
    info = C_DGESV(num_ai, 1, &(A.matrix[irrep][0][0]), num_ai, ipiv, vector, num_ai);
    if(info) {
      psi::fprintf(outfile, "CCSORT: cphf_F: Error in C_DGESV.  Info = %d.  Cart = Y. Exiting.\n", info);
      exit(PSI_RETURN_FAILURE);
    }
    free(ipiv);

    global_dpd_->buf4_mat_irrep_close(&A, irrep);

    /* sort CPHF solution to DPD format */
    global_dpd_->file2_init(&mu, PSIF_CC_OEI, irrep, 1, 0, "CPHF Uf_Y_AI");
    global_dpd_->file2_mat_init(&mu);
    for(row=0; row < num_ai; row++) {
      a = A.params->roworb[irrep][row][0];
      i = A.params->roworb[irrep][row][1];
      asym = A.params->psym[a];
      isym = A.params->qsym[i];
      mu.matrix[asym][a-A.params->poff[asym]][i-A.params->qoff[isym]] = vector[row];
    }
    global_dpd_->file2_mat_wrt(&mu);
    global_dpd_->file2_close(&mu);
  }

  if (!strcmp(cart,"Z")) {
    irrep = moinfo.irrep_z;

    /* sort Mu elements into a single vector for lineq solver */
    global_dpd_->file2_init(&mu, PSIF_CC_OEI, irrep, 0, 1, "Mu_Z_IA");
    global_dpd_->file2_mat_init(&mu);
    global_dpd_->file2_mat_rd(&mu);
    num_ai = A.params->rowtot[irrep];
    vector = init_array(num_ai);
    for(row=0; row < num_ai; row++) {
      i = A.params->roworb[irrep][row][0];
      a = A.params->roworb[irrep][row][1];
      isym = A.params->psym[i];
      asym = A.params->qsym[a];
      vector[row] = -mu.matrix[asym][a-A.params->qoff[asym]][i-A.params->poff[isym]];
    }
    global_dpd_->file2_mat_close(&mu);
    global_dpd_->file2_close(&mu);

    /* grab current irrep of MO Hessian */
    global_dpd_->buf4_mat_irrep_init(&A, irrep);
    global_dpd_->buf4_mat_irrep_rd(&A, irrep);

    /* solve CPHF equations */
    ipiv = init_int_array(num_ai);
    info = C_DGESV(num_ai, 1, &(A.matrix[irrep][0][0]), num_ai, ipiv, vector, num_ai);
    if(info) {
      psi::fprintf(outfile, "CCSORT: cphf_F: Error in C_DGESV.  Info = %d.  Cart = Z. Exiting.\n", info);
      exit(PSI_RETURN_FAILURE);
    }
    free(ipiv);

    global_dpd_->buf4_mat_irrep_close(&A, irrep);

    /* sort CPHF solution to DPD format */
    global_dpd_->file2_init(&mu, PSIF_CC_OEI, irrep, 1, 0, "CPHF Uf_Z_AI");
    global_dpd_->file2_mat_init(&mu);
    for(row=0; row < num_ai; row++) {
      a = A.params->roworb[irrep][row][0];
      i = A.params->roworb[irrep][row][1];
      asym = A.params->psym[a];
      isym = A.params->qsym[i];
      mu.matrix[asym][a-A.params->poff[asym]][i-A.params->qoff[isym]] = vector[row];
    }
    global_dpd_->file2_mat_wrt(&mu);
    global_dpd_->file2_close(&mu);
  }

  global_dpd_->buf4_close(&A);

  if (!strcmp(cart,"Z"))
    psio_close(PSIF_MO_HESS, 0);
  else
    psio_close(PSIF_MO_HESS, 1);
}

}} // namespace psi::ccsort
