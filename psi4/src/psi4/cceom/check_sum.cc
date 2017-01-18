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
#include <cstring>
#include "psi4/libpsio/psio.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

extern double norm_C(dpdfile2 *CME, dpdfile2 *Cme,
    dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf);
extern double norm_C_full(double C0, dpdfile2 *CME, dpdfile2 *Cme,
    dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf);
extern double norm_C_rhf(dpdfile2 *CME, dpdbuf4 *CMnEf, dpdbuf4 *CMnfE);
extern double norm_C_rhf_full(double C0, dpdfile2 *CME, dpdbuf4 *CMnEf, dpdbuf4 *CMnfE);

void check_sum(const char *term_lbl, int index, int irrep) {
  int save_params_ref;
  dpdfile2 Sia, SIA;
  dpdbuf4 SIJAB, Sijab, SIjAb, SIjbA;
  static double old_norm=0;
  double norm,dotval,S0;
  char lbl[80];

  if (!strcmp(term_lbl,"reset"))  {
    outfile->Printf("resetting norm\n");
  	old_norm = 0;
  	return;
  }

  /* save_params_ref = params.ref;
     params.ref = 0; */

  if (params.eom_ref == 0) {
    sprintf(lbl, "%s %d", "SIA", index);
    global_dpd_->file2_init(&SIA, PSIF_EOM_SIA, irrep, 0, 1, lbl);
    sprintf(lbl, "%s %d", "SIjAb", index);
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_sort(&SIjAb, PSIF_EOM_SIjAb, pqsr, 0, 5, "SIjbA");
    global_dpd_->buf4_init(&SIjbA, PSIF_EOM_SIjAb, irrep, 0, 5, 0, 5, 0, "SIjbA");

    if (!params.full_matrix) {
      norm = norm_C_rhf(&SIA, &SIjAb, &SIjbA);
		}
    else {
      sprintf(lbl, "%s %d", "S0", index);
      psio_read_entry(PSIF_EOM_SIA, lbl, (char *) &S0, sizeof(double));
      norm = norm_C_rhf_full(S0, &SIA, &SIjAb, &SIjbA);
		}

    global_dpd_->file2_close(&SIA);
    global_dpd_->buf4_close(&SIjAb);
    global_dpd_->buf4_close(&SIjbA);
  }
  else if (params.eom_ref == 1) {
    sprintf(lbl, "%s %d", "SIA", index);
    global_dpd_->file2_init(&SIA, PSIF_EOM_SIA, irrep, 0, 1, lbl);
    sprintf(lbl, "%s %d", "Sia", index);
    global_dpd_->file2_init(&Sia, PSIF_EOM_Sia, irrep, 0, 1, lbl);
    sprintf(lbl, "%s %d", "SIJAB", index);
    global_dpd_->buf4_init(&SIJAB, PSIF_EOM_SIJAB, irrep, 2, 7, 2, 7, 0, lbl);
    sprintf(lbl, "%s %d", "Sijab", index);
    global_dpd_->buf4_init(&Sijab, PSIF_EOM_Sijab, irrep, 2, 7, 2, 7, 0, lbl);
    sprintf(lbl, "%s %d", "SIjAb", index);
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, irrep, 0, 5, 0, 5, 0, lbl);

    norm = norm_C(&SIA, &Sia, &SIJAB, &Sijab, &SIjAb);

    global_dpd_->file2_close(&SIA);
    global_dpd_->file2_close(&Sia);
    global_dpd_->buf4_close(&SIJAB);
    global_dpd_->buf4_close(&Sijab);
    global_dpd_->buf4_close(&SIjAb);
  }
  else if (params.eom_ref == 2) {
    sprintf(lbl, "%s %d", "SIA", index);
    global_dpd_->file2_init(&SIA, PSIF_EOM_SIA, irrep, 0, 1, lbl);
    sprintf(lbl, "%s %d", "Sia", index);
    global_dpd_->file2_init(&Sia, PSIF_EOM_Sia, irrep, 2, 3, lbl);
    sprintf(lbl, "%s %d", "SIJAB", index);
    global_dpd_->buf4_init(&SIJAB, PSIF_EOM_SIJAB, irrep, 2, 7, 2, 7, 0, lbl);
    sprintf(lbl, "%s %d", "Sijab", index);
    global_dpd_->buf4_init(&Sijab, PSIF_EOM_Sijab, irrep, 12, 17, 12, 17, 0, lbl);
    sprintf(lbl, "%s %d", "SIjAb", index);
    global_dpd_->buf4_init(&SIjAb, PSIF_EOM_SIjAb, irrep, 22, 28, 22, 28, 0, lbl);

    norm = norm_C(&SIA, &Sia, &SIJAB, &Sijab, &SIjAb);

    global_dpd_->file2_close(&SIA);
    global_dpd_->file2_close(&Sia);
    global_dpd_->buf4_close(&SIJAB);
    global_dpd_->buf4_close(&Sijab);
    global_dpd_->buf4_close(&SIjAb);
  }

  outfile->Printf("%7s, D(norm sigma)=%15.10lf\n", term_lbl, norm - old_norm);

  old_norm = norm;

  /* params.ref = save_params_ref; */

  return;
}

}} // namespace psi::cceom
