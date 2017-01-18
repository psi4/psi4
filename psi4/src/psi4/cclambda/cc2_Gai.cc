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
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

void cc2_Gai_build(int L_irr) {
  dpdbuf4 tIJAB, tijab, tiJaB, tIjAb, tijAB, tIJab, t2;
  dpdfile2 G, GAI, Gai, L1, LIA, Lia;
  dpdbuf4 LIJAB, Lijab, LIjAb, LiJaB;

    if(params.ref == 0) {
      global_dpd_->file2_init(&G, PSIF_CC_TMP0, L_irr, 1, 0, "CC2 GAI");

      global_dpd_->file2_init(&L1, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
      global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "2 tIAjb - tIBja");
      global_dpd_->contract422(&t2, &L1, &G, 0, 1, 1, 0);
      global_dpd_->buf4_close(&t2);
      global_dpd_->file2_close(&L1);

      global_dpd_->file2_close(&G);
    }

    else if(params.ref == 1) { /** ROHF **/

      global_dpd_->file2_init(&G, PSIF_CC_TMP0, L_irr, 1, 0, "GAI");
      global_dpd_->file2_init(&G, PSIF_CC_TMP0, L_irr, 4, 3, "Gai");

      /* T2(MJ,AB) * L2(IJ,AB) --> G(M,I) */
      global_dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
      global_dpd_->buf4_init(&LIJAB, PSIF_CC_LAMBDA, L_irr, 0, 7, 2, 7, 0, "LIJAB");
      global_dpd_->contract442(&tIJAB, &LIJAB, &G, 0, 0, 1.0, 0.0);
      global_dpd_->buf4_close(&tIJAB);
      global_dpd_->buf4_close(&LIJAB);

      /* T2(Mj,Ab) * L2(Ij,Ab) --> G(M,I) */
      global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
      global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
      global_dpd_->contract442(&tIjAb, &LIjAb, &G, 0, 0, 1.0, 1.0);
      global_dpd_->buf4_close(&tIjAb);
      global_dpd_->buf4_close(&LIjAb);

      /* T2(mj,ab) * L2(ij,ab) --> G(m,i) */
      global_dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
      global_dpd_->buf4_init(&Lijab, PSIF_CC_LAMBDA, L_irr, 0, 7, 2, 7, 0, "Lijab");
      global_dpd_->contract442(&tijab, &Lijab, &G, 0, 0, 1.0, 0.0);
      global_dpd_->buf4_close(&tijab);
      global_dpd_->buf4_close(&Lijab);

      /* T2(mJ,aB) * L2(iJ,aB) --> G(m,i) */
      global_dpd_->buf4_init(&tiJaB, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
      global_dpd_->buf4_init(&LiJaB, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LiJaB");
      global_dpd_->contract442(&tiJaB, &LiJaB, &G, 0, 0, 1.0, 1.0);
      global_dpd_->buf4_close(&tiJaB);
      global_dpd_->buf4_close(&LiJaB);

      global_dpd_->file2_close(&G);
      global_dpd_->file2_close(&G);
    }

    else if(params.ref == 2) { /** UHF **/

      global_dpd_->file2_init(&GAI, PSIF_CC_TMP0, L_irr, 1, 0, "CC2 GAI");
      global_dpd_->file2_init(&Gai, PSIF_CC_TMP0, L_irr, 3, 2, "CC2 Gai");

      /** AA **/
      global_dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
      global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
      global_dpd_->contract422(&tIJAB, &LIA, &GAI, 0, 1, 1, 0);
      global_dpd_->file2_close(&LIA);
      global_dpd_->buf4_close(&tIJAB);

      global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
      global_dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, L_irr, 2, 3, "Lia");
      global_dpd_->contract422(&tIjAb, &Lia, &GAI, 0, 1, 1, 1);
      global_dpd_->file2_close(&Lia);
      global_dpd_->buf4_close(&tIjAb);

      /** BB **/
      global_dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 30, 30, 30, 30, 0, "tiajb");
      global_dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, L_irr, 2, 3, "Lia");
      global_dpd_->contract422(&tijab, &Lia, &Gai, 0, 1, 1, 0);
      global_dpd_->file2_close(&Lia);
      global_dpd_->buf4_close(&tijab);

      global_dpd_->buf4_init(&tiJaB, PSIF_CC_TAMPS, 0, 30, 20, 30, 20, 0, "tiaJB");
      global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
      global_dpd_->contract422(&tiJaB, &LIA, &Gai, 0, 1, 1, 1);
      global_dpd_->file2_close(&LIA);
      global_dpd_->buf4_close(&tiJaB);

      global_dpd_->file2_close(&Gai);
      global_dpd_->file2_close(&GAI);
    }

  return;
}



}} // namespace psi::cclambda
