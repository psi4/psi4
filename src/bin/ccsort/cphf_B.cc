/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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

void cphf_B(const char *cart)
{
  int irrep, row, a, i, asym, isym, num_ai, info, *ipiv;
  double *vector, *vector2, polar;
  dpdbuf4 A;
  dpdfile2 L;

  psio_open(PSIF_MO_HESS, 1);
  global_dpd_->buf4_init(&A, PSIF_MO_HESS, 0, 11, 11, 11, 11, 0, "B(AI,BJ)");

  if (!strcmp(cart, "X")) {
    irrep = moinfo.irrep_x;

    /* sort L elements into a single vector for lineq solver */
    global_dpd_->file2_init(&L, PSIF_CC_OEI, irrep, 0, 1, "L_X_IA");
    global_dpd_->file2_mat_init(&L);
    global_dpd_->file2_mat_rd(&L);
    num_ai = A.params->rowtot[irrep];
    vector = init_array(num_ai);
    for(row=0; row < num_ai; row++) {
      a = A.params->roworb[irrep][row][0];
      i = A.params->roworb[irrep][row][1];
      asym = A.params->psym[a];
      isym = A.params->qsym[i];
      vector[row] = L.matrix[isym][i-A.params->qoff[isym]][a-A.params->poff[asym]];
    }
    global_dpd_->file2_mat_close(&L);
    global_dpd_->file2_close(&L);

    /* grab current irrep of MO Hessian */
    global_dpd_->buf4_mat_irrep_init(&A, irrep);
    global_dpd_->buf4_mat_irrep_rd(&A, irrep);


    /* solve CPHF equations */
    ipiv = init_int_array(num_ai);
    info = C_DGESV(num_ai, 1, &(A.matrix[irrep][0][0]), num_ai, ipiv, vector, num_ai);
    if(info) {
      outfile->Printf( "CCSORT: cphf_B: Error in C_DGESV.  Info = %d.  Exiting.\n", info);
      exit(PSI_RETURN_FAILURE);
    }
    free(ipiv);

    global_dpd_->buf4_mat_irrep_close(&A, irrep);

    /* sort CPHF solution to DPD format */
    global_dpd_->file2_init(&L, PSIF_CC_OEI, irrep, 1, 0, "CPHF Ub_X_AI");
    global_dpd_->file2_mat_init(&L);
    for(row=0; row < num_ai; row++) {
      a = A.params->roworb[irrep][row][0];
      i = A.params->roworb[irrep][row][1];
      asym = A.params->psym[a];
      isym = A.params->qsym[i];
      L.matrix[asym][a-A.params->poff[asym]][i-A.params->qoff[isym]] = vector[row];
    }
    global_dpd_->file2_mat_wrt(&L);
    global_dpd_->file2_close(&L);
  }

  if (!strcmp(cart, "Y")) {
    irrep = moinfo.irrep_y;

    /* sort L elements into a single vector for lineq solver */
    global_dpd_->file2_init(&L, PSIF_CC_OEI, irrep, 0, 1, "L_Y_IA");
    global_dpd_->file2_mat_init(&L);
    global_dpd_->file2_mat_rd(&L);
    num_ai = A.params->rowtot[irrep];
    vector = init_array(num_ai);
    for(row=0; row < num_ai; row++) {
      a = A.params->roworb[irrep][row][0];
      i = A.params->roworb[irrep][row][1];
      asym = A.params->psym[a];
      isym = A.params->qsym[i];
      vector[row] = L.matrix[isym][i-A.params->qoff[isym]][a-A.params->poff[asym]];
    }
    global_dpd_->file2_mat_close(&L);
    global_dpd_->file2_close(&L);

    /* grab current irrep of MO Hessian */
    global_dpd_->buf4_mat_irrep_init(&A, irrep);
    global_dpd_->buf4_mat_irrep_rd(&A, irrep);


    /* solve CPHF equations */
    ipiv = init_int_array(num_ai);
    info = C_DGESV(num_ai, 1, &(A.matrix[irrep][0][0]), num_ai, ipiv, vector, num_ai);
    if(info) {
      outfile->Printf( "CCSORT: cphf_B: Error in C_DGESV.  Info = %d.  Exiting.\n", info);
      exit(PSI_RETURN_FAILURE);
    }
    free(ipiv);

    global_dpd_->buf4_mat_irrep_close(&A, irrep);

    /* sort CPHF solution to DPD format */
    global_dpd_->file2_init(&L, PSIF_CC_OEI, irrep, 1, 0, "CPHF Ub_Y_AI");
    global_dpd_->file2_mat_init(&L);
    for(row=0; row < num_ai; row++) {
      a = A.params->roworb[irrep][row][0];
      i = A.params->roworb[irrep][row][1];
      asym = A.params->psym[a];
      isym = A.params->qsym[i];
      L.matrix[asym][a-A.params->poff[asym]][i-A.params->qoff[isym]] = vector[row];
    }
    global_dpd_->file2_mat_wrt(&L);
    global_dpd_->file2_close(&L);
  }

  if (!strcmp(cart, "Z")) {
    irrep = moinfo.irrep_z;

    /* sort L elements into a single vector for lineq solver */
    global_dpd_->file2_init(&L, PSIF_CC_OEI, irrep, 0, 1, "L_Z_IA");
    global_dpd_->file2_mat_init(&L);
    global_dpd_->file2_mat_rd(&L);
    num_ai = A.params->rowtot[irrep];
    vector = init_array(num_ai);
    for(row=0; row < num_ai; row++) {
      a = A.params->roworb[irrep][row][0];
      i = A.params->roworb[irrep][row][1];
      asym = A.params->psym[a];
      isym = A.params->qsym[i];
      vector[row] = L.matrix[isym][i-A.params->qoff[isym]][a-A.params->poff[asym]];
    }
    global_dpd_->file2_mat_close(&L);
    global_dpd_->file2_close(&L);

    /*
      vector2 = init_array(num_ai);
      for(row=0; row < num_ai; row++) vector2[row] = vector[row];
    */

    /* grab current irrep of MO Hessian */
    global_dpd_->buf4_mat_irrep_init(&A, irrep);
    global_dpd_->buf4_mat_irrep_rd(&A, irrep);

    /* solve CPHF equations */
    ipiv = init_int_array(num_ai);
    info = C_DGESV(num_ai, 1, &(A.matrix[irrep][0][0]), num_ai, ipiv, vector, num_ai);
    if(info) {
      outfile->Printf( "CCSORT: cphf_B: Error in C_DGESV.  Info = %d.  Exiting.\n", info);
      exit(PSI_RETURN_FAILURE);
    }
    free(ipiv);

    global_dpd_->buf4_mat_irrep_close(&A, irrep);

    /*
      polar = 0.0;
      for(row=0; row < num_ai; row++) polar += vector2[row] * vector[row];
      polar *= 4.0;
      outfile->Printf( "\talpha_zz = %20.12f\n", polar);
      free(vector2);
    */

    /* sort CPHF solution to DPD format */
    global_dpd_->file2_init(&L, PSIF_CC_OEI, irrep, 1, 0, "CPHF Ub_Z_AI");
    global_dpd_->file2_mat_init(&L);
    for(row=0; row < num_ai; row++) {
      a = A.params->roworb[irrep][row][0];
      i = A.params->roworb[irrep][row][1];
      asym = A.params->psym[a];
      isym = A.params->qsym[i];
      L.matrix[asym][a-A.params->poff[asym]][i-A.params->qoff[isym]] = vector[row];
    }
    global_dpd_->file2_mat_wrt(&L);
    global_dpd_->file2_close(&L);
  }

  global_dpd_->buf4_close(&A);
  if (!strcmp(cart,"Z"))
    psio_close(PSIF_MO_HESS, 0);
  else
    psio_close(PSIF_MO_HESS, 1);
}

}} // namespace psi::ccsort