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
    \ingroup CCLAMBDA
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstring>
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

double norm_C(dpdfile2 *CME, dpdfile2 *Cme,
    dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf);

double norm_C_rhf(dpdfile2 *CME, dpdbuf4 *CMnEf, dpdbuf4 *CMnfE);

void check_sum(char *term_lbl, int irrep) {
  dpdfile2 Lia, LIA;
  dpdbuf4 LIJAB, Lijab, LIjAb, LIjbA;
  static double old_norm=0;
  double norm,dotval;
  char lbl[80];

  if (!strcmp(term_lbl,"reset"))  {
    outfile->Printf("resetting norm\n");
  	old_norm = 0;
  	return;
  }

  if (params.ref <= 1) {
    global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, irrep, 0, 1, "New LIA");
    global_dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, irrep, 0, 1, "New Lia");
    global_dpd_->buf4_init(&LIJAB, PSIF_CC_LAMBDA, irrep, 2, 7, 2, 7, 0, "New LIJAB");
    global_dpd_->buf4_init(&Lijab, PSIF_CC_LAMBDA, irrep, 2, 7, 2, 7, 0, "New Lijab");
    global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, irrep, 0, 5, 0, 5, 0, "New LIjAb");

    norm = norm_C(&LIA, &Lia, &LIJAB, &Lijab, &LIjAb);

    global_dpd_->file2_close(&LIA);
    global_dpd_->file2_close(&Lia);
    global_dpd_->buf4_close(&LIJAB);
    global_dpd_->buf4_close(&Lijab);
    global_dpd_->buf4_close(&LIjAb);
  }
  else if (params.ref == 2) {
    global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, irrep, 0, 1, "New LIA");
    global_dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, irrep, 2, 3, "New Lia");
    global_dpd_->buf4_init(&LIJAB, PSIF_CC_LAMBDA, irrep, 2, 7, 2, 7, 0, "New LIJAB");
    global_dpd_->buf4_init(&Lijab, PSIF_CC_LAMBDA, irrep, 12, 17, 12, 17, 0, "New Lijab");
    global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, irrep, 22, 28, 22, 28, 0, "New LIjAb");

    norm = norm_C(&LIA, &Lia, &LIJAB, &Lijab, &LIjAb);

    global_dpd_->file2_close(&LIA);
    global_dpd_->file2_close(&Lia);
    global_dpd_->buf4_close(&LIJAB);
    global_dpd_->buf4_close(&Lijab);
    global_dpd_->buf4_close(&LIjAb);
  }

  outfile->Printf("%7s, D(norm L)=%15.10lf\n", term_lbl, norm - old_norm);

  old_norm = norm;
  return;
}

double norm_C(dpdfile2 *CME, dpdfile2 *Cme,
        dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf)
{
  double norm = 0.0;
  norm += global_dpd_->file2_dot_self(CME);
  norm += global_dpd_->file2_dot_self(Cme);
  norm += global_dpd_->buf4_dot_self(CMNEF);
  norm += global_dpd_->buf4_dot_self(Cmnef);
  norm += global_dpd_->buf4_dot_self(CMnEf);
  norm = sqrt(norm);
  return norm;
}

double norm_C_rhf(dpdfile2 *CME, dpdbuf4 *CMnEf, dpdbuf4 *CMnfE) {
  double norm = 0.0;
  norm = 2.0 * global_dpd_->file2_dot_self(CME);
  norm += 2.0 * global_dpd_->buf4_dot_self(CMnEf);
  norm -= global_dpd_->buf4_dot(CMnEf, CMnfE);
  norm = sqrt(norm);
  return norm;
}


}} // namespace psi::cclambda
