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
/*
**  X_XI_CHECK: check sum for xi
*/

#include <cstdio>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

extern double norm_C(dpdfile2 *CME, dpdfile2 *Cme,
        dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf);

extern double norm_C_rhf(dpdfile2 *CME, dpdbuf4 *CMnEf, dpdbuf4 *CMnfE);

extern void c_cleanSS(dpdfile2 *CME, dpdfile2 *Cme);

void c_clean(dpdfile2 *CME, dpdfile2 *Cme, dpdbuf4 *CMNEF, dpdbuf4 *Cmnef, dpdbuf4 *CMnEf);

void x_xi_check(char *term_lbl) {
  dpdfile2 Xia, XIA;
  dpdbuf4 XIJAB, Xijab, XIjAb, XIjbA;
  static double old_norm=0;
  double norm,dotval;
  char lbl[80];
  int irrep;
  irrep = params.G_irr;

  /*
  if (!strcmp(term_lbl,"reset"))  {
    outfile->Printf("resetting norm\n");
    old_norm = 0;
    return;
  }
  */

  if (params.ref == 0) {
    global_dpd_->file2_init(&XIA, PSIF_EOM_XI, irrep, 0, 1, "XIA");
    global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, irrep, 0, 5, 0, 5, 0, "XIjAb");
    global_dpd_->buf4_sort(&XIjAb, PSIF_EOM_XI, pqsr, 0, 5, "XIjbA");
    global_dpd_->buf4_init(&XIjbA, PSIF_EOM_XI, irrep, 0, 5, 0, 5, 0, "XIjbA");

    norm = norm_C_rhf(&XIA, &XIjAb, &XIjbA);

    global_dpd_->file2_close(&XIA);
    global_dpd_->buf4_close(&XIjAb);
    global_dpd_->buf4_close(&XIjbA);
  }
  else if (params.ref == 1) {
    global_dpd_->file2_init(&XIA, PSIF_EOM_XI, irrep, 0, 1, "XIA");
    global_dpd_->file2_init(&Xia, PSIF_EOM_XI, irrep, 0, 1, "Xia");
    global_dpd_->buf4_init(&XIJAB, PSIF_EOM_XI, irrep, 2, 7, 2, 7, 0, "XIJAB");
    global_dpd_->buf4_init(&Xijab, PSIF_EOM_XI, irrep, 2, 7, 2, 7, 0, "Xijab");
    global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, irrep, 0, 5, 0, 5, 0, "XIjAb");

    c_clean(&XIA, &Xia, &XIJAB, &Xijab, &XIjAb);
    norm = norm_C(&XIA, &Xia, &XIJAB, &Xijab, &XIjAb);

    global_dpd_->file2_close(&XIA);
    global_dpd_->file2_close(&Xia);
    global_dpd_->buf4_close(&XIJAB);
    global_dpd_->buf4_close(&Xijab);
    global_dpd_->buf4_close(&XIjAb);
  }
  else if (params.ref == 2) {
    global_dpd_->file2_init(&XIA, PSIF_EOM_XI, irrep, 0, 1, "XIA");
    global_dpd_->file2_init(&Xia, PSIF_EOM_XI, irrep, 2, 3, "Xia");
    global_dpd_->buf4_init(&XIJAB, PSIF_EOM_XI, irrep, 2, 7, 2, 7, 0, "XIJAB");
    global_dpd_->buf4_init(&Xijab, PSIF_EOM_XI, irrep, 12, 17, 12, 17, 0, "Xijab");
    global_dpd_->buf4_init(&XIjAb, PSIF_EOM_XI, irrep, 22, 28, 22, 28, 0, "XIjAb");

    norm = norm_C(&XIA, &Xia, &XIJAB, &Xijab, &XIjAb);

    global_dpd_->file2_close(&XIA);
    global_dpd_->file2_close(&Xia);
    global_dpd_->buf4_close(&XIJAB);
    global_dpd_->buf4_close(&Xijab);
    global_dpd_->buf4_close(&XIjAb);
  }

  outfile->Printf("%7s, D(norm sigma)=%15.10lf\n", term_lbl, norm - old_norm);

  old_norm = norm;
  return;
}

}} // namespace psi::ccdensity
