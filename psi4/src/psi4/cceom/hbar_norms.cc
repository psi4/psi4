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
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cceom {

void hbar_norms() {
  double tval;
  dpdfile2 FAE, Fae, FMI, Fmi, FME, Fme;
  dpdbuf4 WMBIJ, Wmbij, WMbIj, WmBiJ, W;

  outfile->Printf("\n");

  if ((params.eom_ref == 0) || (params.eom_ref == 1)) {
  global_dpd_->file2_init(&FAE, PSIF_CC_OEI, H_IRR, 1, 1, "FAE");
  global_dpd_->file2_init(&Fae, PSIF_CC_OEI, H_IRR, 1, 1, "Fae");
  tval = global_dpd_->file2_dot_self(&FAE);
  tval += global_dpd_->file2_dot_self(&Fae);
  global_dpd_->file2_close(&Fae);
  global_dpd_->file2_close(&FAE);
  outfile->Printf("Fae   dot Fae   total %15.10lf\n", tval);

  global_dpd_->file2_init(&FMI, PSIF_CC_OEI, H_IRR, 0, 0, "FMI");
  global_dpd_->file2_init(&Fmi, PSIF_CC_OEI, H_IRR, 0, 0, "Fmi");
  tval = global_dpd_->file2_dot_self(&FMI);
  tval += global_dpd_->file2_dot_self(&Fmi); 
  global_dpd_->file2_close(&Fmi);
  global_dpd_->file2_close(&FMI);
  outfile->Printf("Fmi   dot Fmi   total %15.10lf\n", tval);

  global_dpd_->file2_init(&FME, PSIF_CC_OEI, H_IRR, 0, 1, "FME");
  global_dpd_->file2_init(&Fme, PSIF_CC_OEI, H_IRR, 0, 1, "Fme");
  tval = global_dpd_->file2_dot_self(&FME);
  tval += global_dpd_->file2_dot_self(&Fme);
  global_dpd_->file2_close(&Fme);
  global_dpd_->file2_close(&FME);
  outfile->Printf("Fme   dot Fme   total %15.10lf\n", tval);


  global_dpd_->buf4_init(&WMBIJ, PSIF_CC_HBAR, H_IRR, 10, 2, 10, 2, 0, "WMBIJ");
  tval = 2 * global_dpd_->buf4_dot_self(&WMBIJ);
  global_dpd_->buf4_close(&WMBIJ);
  outfile->Printf("WMBIJ dot WMBIJ total %15.10lf\n", tval);

  global_dpd_->buf4_init(&Wmbij, PSIF_CC_HBAR, H_IRR, 10, 2, 10, 2, 0, "Wmbij");
  tval = 2 * global_dpd_->buf4_dot_self(&Wmbij);
  global_dpd_->buf4_close(&Wmbij);
  outfile->Printf("Wmbij dot Wmbij total %15.10lf\n", tval);

  global_dpd_->buf4_init(&WMbIj, PSIF_CC_HBAR, H_IRR, 10, 0, 10, 0, 0, "WMbIj");
  tval = global_dpd_->buf4_dot_self(&WMbIj);
  global_dpd_->buf4_close(&WMbIj);
  outfile->Printf("WMbIj dot WMbIj total %15.10lf\n", tval);

  global_dpd_->buf4_init(&WmBiJ, PSIF_CC_HBAR, H_IRR, 11, 0, 11, 0, 0, "WmBiJ (Bm,Ji)");
  tval = global_dpd_->buf4_dot_self(&WmBiJ);
  global_dpd_->buf4_close(&WmBiJ);
  outfile->Printf("WmBiJ dot WmBiJ total %15.10lf\n", tval);

	  if (params.full_matrix) {
      global_dpd_->file2_init(&FME, PSIF_CC_OEI, H_IRR, 0, 1, "FAI residual");
      tval = global_dpd_->file2_dot_self(&FME);
			global_dpd_->file2_close(&FME);
      outfile->Printf("FAI residual dot FAI residual %15.10lf\n", tval);
	  }
  }
  
  else if (params.eom_ref == 2) {

    global_dpd_->file2_init(&FAE, PSIF_CC_OEI, H_IRR, 1, 1, "FAE");
    global_dpd_->file2_init(&Fae, PSIF_CC_OEI, H_IRR, 3, 3, "Fae");
    tval = global_dpd_->file2_dot_self(&FAE);
    tval += global_dpd_->file2_dot_self(&Fae);
    global_dpd_->file2_close(&Fae);
    global_dpd_->file2_close(&FAE);
    outfile->Printf("Fae   dot Fae   total %15.10lf\n", tval);

    global_dpd_->file2_init(&FMI, PSIF_CC_OEI, H_IRR, 0, 0, "FMI");
    global_dpd_->file2_init(&Fmi, PSIF_CC_OEI, H_IRR, 2, 2, "Fmi");
    tval = global_dpd_->file2_dot_self(&FMI);
    tval += global_dpd_->file2_dot_self(&Fmi); 
    global_dpd_->file2_close(&Fmi);
    global_dpd_->file2_close(&FMI);
    outfile->Printf("Fmi   dot Fmi   total %15.10lf\n", tval);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 20, 30, 20, 30, 0, "WmBeJ (JB,me)"); /* (me,JB) */
    tval = global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, H_IRR, 30, 20, 30, 20, 0, "WMbEj (jb,ME)"); /* (ME,jb) */
    tval += global_dpd_->buf4_dot_self(&W);
    global_dpd_->buf4_close(&W);
    outfile->Printf("WmBeJ and WMbEj dots %15.10lf\n",tval);


    /*
    dpd_file2_init(&FME, CC_OEI, H_IRR, 0, 1, "FME");
    dpd_file2_init(&Fme, CC_OEI, H_IRR, 0, 1, "Fme");
    dpd_buf4_init(&W, CC_HBAR, H_IRR, 2, 21, 2, 21, 0, "WMNIE");
    tval = 2.0 * dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("WMNIE dot WMNIE total %15.10lf\n", tval);

    dpd_buf4_init(&W, CC_HBAR, H_IRR, 12, 31, 12, 31, 0, "Wmnie");
    tval += 2.0 * dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("Wmnie dot Wmnie total %15.10lf\n", tval);

    dpd_buf4_init(&W, CC_HBAR, H_IRR, 22, 25, 22, 25, 0, "WMnIe");
    tval += dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("WMnIe dot WMnIe total %15.10lf\n", tval);

    dpd_buf4_init(&W, CC_HBAR, H_IRR, 23, 26, 23, 26, 0, "WmNiE");
    tval += dpd_buf4_dot_self(&W);
    dpd_buf4_close(&W);
    outfile->Printf("WmNiE dot WmNiE total %15.10lf\n", tval);
    */
  }
  return;
}

}} // namespace psi::cceom
