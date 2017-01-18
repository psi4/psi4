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
#include <cmath>
#include "psi4/libpsio/psio.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

/* This function orthogonalizes the r residual vector with the set of
numCs C vectors and adds the new vector to the C list if its norm is greater
than params.residual_tol */

extern double norm_C(dpdfile2 *CME, dpdfile2 *Cme,
  dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf);

extern void scm_C(dpdfile2 *CME, dpdfile2 *Cme, dpdbuf4 *CMNEF,
  dpdbuf4 *Cmnef, dpdbuf4 *CMnEf, double a);

/* use for ROHF and UHF */
void schmidt_add(dpdfile2 *RIA, dpdfile2 *Ria,
  dpdbuf4 *RIJAB, dpdbuf4 *Rijab, dpdbuf4 *RIjAb, int *numCs, int irrep)
{
  double dotval;
  double norm;
  int i, I;
  dpdfile2 Cme, CME, Cme2, CME2;
  dpdbuf4 CMNEF, Cmnef, CMnEf, CMNEF2, Cmnef2, CMnEf2;
  dpdbuf4 CMnEf_buf;
  char CME_lbl[32], Cme_lbl[32], CMNEF_lbl[32], Cmnef_lbl[32], CMnEf_lbl[32];

  for (i=0; i<*numCs; i++) {
    sprintf(CME_lbl, "%s %d", "CME", i);
    sprintf(Cme_lbl, "%s %d", "Cme", i);
    sprintf(CMNEF_lbl, "%s %d", "CMNEF", i);
    sprintf(Cmnef_lbl, "%s %d", "Cmnef", i);
    sprintf(CMnEf_lbl, "%s %d", "CMnEf", i);

    global_dpd_->file2_init(&CME, PSIF_EOM_CME, irrep, 0, 1, CME_lbl);
    global_dpd_->buf4_init(&CMNEF, PSIF_EOM_CMNEF, irrep, 2, 7, 2, 7, 0, CMNEF_lbl);
    if (params.eom_ref == 1) {
      global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, irrep, 0, 1, Cme_lbl);
      global_dpd_->buf4_init(&Cmnef, PSIF_EOM_Cmnef, irrep, 2, 7, 2, 7, 0, Cmnef_lbl);
      global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, irrep, 0, 5, 0, 5, 0, CMnEf_lbl);
    }
    else if (params.eom_ref == 2) {
      global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, irrep, 2, 3, Cme_lbl);
      global_dpd_->buf4_init(&Cmnef, PSIF_EOM_Cmnef, irrep, 12, 17, 12, 17, 0, Cmnef_lbl);
      global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, irrep, 22, 28, 22, 28, 0, CMnEf_lbl);
    }

    dotval  = global_dpd_->file2_dot(RIA, &CME);
    dotval += global_dpd_->file2_dot(Ria, &Cme);
  //outfile->Printf( "OE Dotval for vector %d = %20.14f\n", i, dotval);
    dotval += global_dpd_->buf4_dot(RIJAB, &CMNEF);
    dotval += global_dpd_->buf4_dot(Rijab, &Cmnef);
    dotval += global_dpd_->buf4_dot(RIjAb, &CMnEf);

  //outfile->Printf( "Dotval for vector %d = %20.14f\n", i, dotval);

    global_dpd_->file2_axpy(&CME, RIA, -1.0*dotval, 0);
    global_dpd_->file2_axpy(&Cme, Ria, -1.0*dotval, 0);
    global_dpd_->buf4_axpy(&CMNEF, RIJAB, -1.0*dotval);
    global_dpd_->buf4_axpy(&Cmnef, Rijab, -1.0*dotval);
    global_dpd_->buf4_axpy(&CMnEf, RIjAb, -1.0*dotval);

    global_dpd_->file2_close(&CME);
    global_dpd_->file2_close(&Cme);
    global_dpd_->buf4_close(&CMNEF);
    global_dpd_->buf4_close(&Cmnef);
    global_dpd_->buf4_close(&CMnEf);
  }

  norm = norm_C(RIA, Ria, RIJAB, Rijab, RIjAb);
  //outfile->Printf( "Norm of residual (TDC) = %20.14f\n", norm);

  if (norm < eom_params.schmidt_add_residual_tol) {
    return;
  }
  else {
    scm_C(RIA, Ria, RIJAB, Rijab, RIjAb, 1.0/norm);
    sprintf(CME_lbl, "%s %d", "CME", *numCs);
    sprintf(Cme_lbl, "%s %d", "Cme", *numCs);
    sprintf(CMNEF_lbl, "%s %d", "CMNEF", *numCs);
    sprintf(Cmnef_lbl, "%s %d", "Cmnef", *numCs);
    sprintf(CMnEf_lbl, "%s %d", "CMnEf", *numCs);

    global_dpd_->file2_copy(RIA, PSIF_EOM_CME, CME_lbl);
    global_dpd_->file2_copy(Ria, PSIF_EOM_Cme, Cme_lbl);
    global_dpd_->buf4_copy(RIJAB, PSIF_EOM_CMNEF, CMNEF_lbl);
    global_dpd_->buf4_copy(Rijab, PSIF_EOM_Cmnef, Cmnef_lbl);
    global_dpd_->buf4_copy(RIjAb, PSIF_EOM_CMnEf, CMnEf_lbl);

    ++(*numCs);
  }
  return;
}

void schmidt_add_RHF(dpdfile2 *RIA, dpdbuf4 *RIjAb, int *numCs, int irrep)
{
  double dotval, norm, R0, C0;
  int i, I;
  dpdfile2 CME;
  dpdbuf4 CMnEf, CAB1, CAB2;
  dpdfile2 R1;
  dpdbuf4 R2a, R2b;
  char CME_lbl[32], Cme_lbl[32], CMNEF_lbl[32], Cmnef_lbl[32], CMnEf_lbl[32], C0_lbl[32];

  if (params.full_matrix) psio_read_entry(PSIF_EOM_R, "R0", (char *) &R0, sizeof(double));

  for (i=0; i<*numCs; i++) {
    /* Spin-adapt the residual */
    global_dpd_->buf4_copy(RIjAb, PSIF_EOM_TMP, "RIjAb");
    global_dpd_->buf4_sort(RIjAb, PSIF_EOM_TMP, pqsr, 0, 5, "RIjbA");

    global_dpd_->buf4_init(&R2a, PSIF_EOM_TMP, irrep, 0, 5, 0, 5, 0, "RIjAb");
    global_dpd_->buf4_init(&R2b, PSIF_EOM_TMP, irrep, 0, 5, 0, 5, 0, "RIjbA");
    global_dpd_->buf4_scm(&R2a, 2.0);
    global_dpd_->buf4_axpy(&R2b, &R2a, -1.0);
    global_dpd_->buf4_close(&R2b);

    sprintf(CME_lbl, "%s %d", "CME", i);
    sprintf(CMnEf_lbl, "%s %d", "CMnEf", i);
    global_dpd_->file2_init(&CME, PSIF_EOM_CME, irrep, 0, 1, CME_lbl);
    global_dpd_->buf4_init(&CMnEf, PSIF_EOM_CMnEf, irrep, 0, 5, 0, 5, 0, CMnEf_lbl);
    dotval  = 2.0 * global_dpd_->file2_dot(RIA, &CME);
 //outfile->Printf( "OE Dotval for vector %d = %20.14f\n", i, dotval);
    dotval += global_dpd_->buf4_dot(&R2a, &CMnEf);
    global_dpd_->buf4_close(&R2a);
		if (params.full_matrix) {
      sprintf(C0_lbl, "%s %d", "C0", i);
			psio_read_entry(PSIF_EOM_CME, C0_lbl, (char *) &C0, sizeof(double));
			dotval += C0 * R0;
		}

 //outfile->Printf( "Dotval for vector %d = %20.14f\n", i, dotval);
		R0 = R0 - 1.0 * dotval * C0;
    global_dpd_->file2_axpy(&CME, RIA, -1.0*dotval, 0);
    global_dpd_->buf4_axpy(&CMnEf, RIjAb, -1.0*dotval);
    global_dpd_->file2_close(&CME);
    global_dpd_->buf4_close(&CMnEf);
  }

  global_dpd_->buf4_sort(RIjAb, PSIF_EOM_TMP, pqsr, 0, 5, "RIjbA");
  global_dpd_->buf4_init(&R2b, PSIF_EOM_TMP, irrep, 0, 5, 0, 5, 0, "RIjbA");

  /* norm = norm_C_rhf(RIA, RIjAb, &R2b); */
  norm  = 2.0 * global_dpd_->file2_dot_self(RIA);
  norm += 2.0 * global_dpd_->buf4_dot_self(RIjAb);
  norm -= global_dpd_->buf4_dot(RIjAb, &R2b);
	if (params.full_matrix)
	  norm += R0 * R0;
  norm = sqrt(norm);

  global_dpd_->buf4_close(&R2b);

  //outfile->Printf( "Norm of residual (TDC) = %20.14f\n", norm);

  if (norm < eom_params.schmidt_add_residual_tol) {
    return;
  }
  else {

    if (params.full_matrix) R0 *= 1.0/norm;
    global_dpd_->file2_scm(RIA, 1.0/norm);
    global_dpd_->buf4_scm(RIjAb, 1.0/norm);

#ifdef EOM_DEBUG
    dpd_buf4_sort(RIjAb, EOM_TMP, pqsr, 0, 5, "RIjbA");
    dpd_buf4_init(&R2b, EOM_TMP, irrep, 0, 5, 0, 5, 0, "RIjbA");
    norm  = 2.0 * dpd_file2_dot_self(RIA);
    norm += 2.0 * dpd_buf4_dot_self(RIjAb);
    norm -= dpd_buf4_dot(RIjAb, &R2b);
		if (params.full_matrix) norm += R0 * R0;
    norm = sqrt(norm);
    outfile->Printf("Norm of final new C in schmidt_add(): %20.15lf\n", norm);
    dpd_buf4_close(&R2b);
#endif

    sprintf(CME_lbl, "%s %d", "CME", *numCs);
    sprintf(CMnEf_lbl, "%s %d", "CMnEf", *numCs);

    global_dpd_->file2_copy(RIA, PSIF_EOM_CME, CME_lbl);
    global_dpd_->buf4_copy(RIjAb, PSIF_EOM_CMnEf, CMnEf_lbl);

    /* Generate AA and BB C2 vectors from AB vector */
    /* C(IJ,AB) = C(ij,ab) = C(Ij,Ab) - C(Ij,bA) */
    global_dpd_->buf4_copy(RIjAb, PSIF_EOM_TMP, "CMnEf");
    global_dpd_->buf4_sort(RIjAb, PSIF_EOM_TMP, pqsr, 0, 5, "CMnfE");

    global_dpd_->buf4_init(&CAB1, PSIF_EOM_TMP, irrep, 0, 5, 0, 5, 0, "CMnEf");
    global_dpd_->buf4_init(&CAB2, PSIF_EOM_TMP, irrep, 0, 5, 0, 5, 0, "CMnfE");
    global_dpd_->buf4_axpy(&CAB2, &CAB1, -1.0);
    global_dpd_->buf4_close(&CAB2);
    global_dpd_->buf4_close(&CAB1);

		if (params.full_matrix) {
      sprintf(C0_lbl, "%s %d", "C0", *numCs);
		  psio_write_entry(PSIF_EOM_CME, C0_lbl, (char *) &R0, sizeof(double));
		}
    ++(*numCs);
  }
  return;
}

}} // namespace psi::cceom
