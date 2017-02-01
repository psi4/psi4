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
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "Params.h"
#include "ccwave.h"

namespace psi { namespace ccenergy {

void CCEnergyWavefunction::Wmnij_build(void)
{
  dpdbuf4 A_anti, A;
  dpdbuf4 WMNIJ, Wmnij, WMnIj, W;
  dpdfile2 tIA, tia;
  dpdbuf4 Eijka, Eijka_anti, Eaijk, Eaijk_anti;
  dpdbuf4 D_anti, D, tauIJAB, tauijab, tauIjAb;

  if(params_.ref == 0) { /** RHF **/
    global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
    global_dpd_->buf4_copy(&A, PSIF_CC_HBAR, "WMnIj");
    global_dpd_->buf4_close(&A);
  }
  else if(params_.ref == 1) { /** ROHF **/
    global_dpd_->buf4_init(&A_anti, PSIF_CC_AINTS, 0, 2, 2, 0, 0, 1, "A <ij|kl>");
    global_dpd_->buf4_copy(&A_anti, PSIF_CC_HBAR, "WMNIJ");
    global_dpd_->buf4_copy(&A_anti, PSIF_CC_HBAR, "Wmnij");
    global_dpd_->buf4_close(&A_anti);

    global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
    global_dpd_->buf4_copy(&A, PSIF_CC_HBAR, "WMnIj");
    global_dpd_->buf4_close(&A);
  }
  else if(params_.ref == 2) { /*** UHF ***/
    global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 2, 2, 0, 0, 1, "A <IJ|KL>");
    global_dpd_->buf4_copy(&A, PSIF_CC_HBAR, "WMNIJ");
    global_dpd_->buf4_close(&A);

    global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 12, 12, 10, 10, 1, "A <ij|kl>");
    global_dpd_->buf4_copy(&A, PSIF_CC_HBAR, "Wmnij");
    global_dpd_->buf4_close(&A);

    global_dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 22, 22, 22, 22, 0, "A <Ij|Kl>");
    global_dpd_->buf4_copy(&A, PSIF_CC_HBAR, "WMnIj");
    global_dpd_->buf4_close(&A);
  }

  if(params_.ref == 0) { /** RHF **/
    global_dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");

    global_dpd_->buf4_init(&Eaijk, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    global_dpd_->contract244(&tIA, &Eaijk, &WMnIj, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&Eaijk);

    global_dpd_->buf4_init(&Eijka, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    global_dpd_->contract424(&Eijka, &tIA, &WMnIj, 3, 1, 0, 1, 1);
    global_dpd_->buf4_close(&Eijka);

    global_dpd_->file2_close(&tIA);
    global_dpd_->buf4_close(&WMnIj);
  }
  else if(params_.ref == 1) { /** ROHF **/
    global_dpd_->buf4_init(&WMNIJ, PSIF_CC_HBAR, 0, 2, 0, 2, 2, 0, "WMNIJ");
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, 0, 2, 0, 2, 2, 0, "Wmnij");
    global_dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    global_dpd_->buf4_init(&Eijka_anti, PSIF_CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
    global_dpd_->buf4_init(&Eijka, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    global_dpd_->buf4_init(&Eaijk_anti, PSIF_CC_EINTS, 0, 11, 2, 11, 0, 1, "E <ai|jk>");
    global_dpd_->buf4_init(&Eaijk, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 2, 0, 2, 0, 0, "W (MN,IJ)");
    global_dpd_->contract424(&Eijka_anti, &tIA, &W, 3, 1, 0, 1, 0);
    global_dpd_->contract244(&tIA, &Eaijk_anti, &W, 1, 0, 1, 1, 1);
    global_dpd_->buf4_axpy(&W, &WMNIJ, 1);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 2, 0, 2, 0, 0, "W (MN,IJ)");
    global_dpd_->contract424(&Eijka_anti, &tia, &W, 3, 1, 0, 1, 0);
    global_dpd_->contract244(&tia, &Eaijk_anti, &W, 1, 0, 1, 1, 1);
    global_dpd_->buf4_axpy(&W, &Wmnij, 1);
    global_dpd_->buf4_close(&W);

    global_dpd_->contract424(&Eijka, &tia, &WMnIj, 3, 1, 0, 1, 1);
    global_dpd_->contract244(&tIA, &Eaijk, &WMnIj, 1, 0, 1, 1, 1);

    global_dpd_->buf4_close(&Eijka_anti);
    global_dpd_->buf4_close(&Eijka);
    global_dpd_->buf4_close(&Eaijk_anti);
    global_dpd_->buf4_close(&Eaijk);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

    global_dpd_->buf4_close(&WMNIJ);
    global_dpd_->buf4_close(&Wmnij);
    global_dpd_->buf4_close(&WMnIj);
  }
  else if(params_.ref == 2) { /*** UHF ***/

    global_dpd_->buf4_init(&WMNIJ, PSIF_CC_HBAR, 0, 2, 0, 2, 2, 0, "WMNIJ");
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, 0, 12, 10, 12, 12, 0, "Wmnij");
    global_dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, 0, 22, 22, 22, 22, 0, "WMnIj");

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    global_dpd_->buf4_init(&Eijka, PSIF_CC_EINTS, 0, 2, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
    global_dpd_->buf4_init(&Eaijk, PSIF_CC_EINTS, 0, 21, 2, 21, 0, 1, "E <AI|JK>");
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 2, 0, 2, 0, 0, "W (MN,IJ)");
    global_dpd_->contract424(&Eijka, &tIA, &W, 3, 1, 0, 1, 0);
    global_dpd_->contract244(&tIA, &Eaijk, &W, 1, 0, 1, 1, 1);
    global_dpd_->buf4_axpy(&W, &WMNIJ, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Eijka);
    global_dpd_->buf4_close(&Eaijk);

    global_dpd_->buf4_init(&Eijka, PSIF_CC_EINTS, 0, 12, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
    global_dpd_->buf4_init(&Eaijk, PSIF_CC_EINTS, 0, 31, 12, 31, 10, 1, "E <ai|jk>");
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 12, 10, 12, 10, 0, "W (mn,ij)");
    global_dpd_->contract424(&Eijka, &tia, &W, 3, 1, 0, 1, 0);
    global_dpd_->contract244(&tia, &Eaijk, &W, 1, 0, 1, 1, 1);
    global_dpd_->buf4_axpy(&W, &Wmnij, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Eijka);
    global_dpd_->buf4_close(&Eaijk);

    global_dpd_->buf4_init(&Eijka, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    global_dpd_->buf4_init(&Eaijk, PSIF_CC_EINTS, 0, 26, 22, 26, 22, 0, "E <Ai|Jk>");
    global_dpd_->contract424(&Eijka, &tia, &WMnIj, 3, 1, 0, 1, 1);
    global_dpd_->contract244(&tIA, &Eaijk, &WMnIj, 1, 0, 1, 1, 1);
    global_dpd_->buf4_close(&Eijka);
    global_dpd_->buf4_close(&Eaijk);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

    global_dpd_->buf4_close(&WMNIJ);
    global_dpd_->buf4_close(&Wmnij);
    global_dpd_->buf4_close(&WMnIj);

  }

  if(params_.ref == 0) { /** RHF **/
    global_dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->buf4_init(&tauIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    global_dpd_->contract444(&D, &tauIjAb, &WMnIj, 0, 0, 1, 1);
    global_dpd_->buf4_close(&tauIjAb);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&WMnIj);
  }
  else if(params_.ref == 1) { /** ROHF **/
    global_dpd_->buf4_init(&WMNIJ, PSIF_CC_HBAR, 0, 2, 2, 2, 2, 0, "WMNIJ");
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, 0, 2, 2, 2, 2, 0, "Wmnij");
    global_dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");

    global_dpd_->buf4_init(&D_anti, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");

    global_dpd_->buf4_init(&tauIJAB, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    global_dpd_->buf4_init(&tauijab, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
    global_dpd_->buf4_init(&tauIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");

    global_dpd_->contract444(&D_anti, &tauIJAB, &WMNIJ, 0, 0, 1, 1);
    global_dpd_->contract444(&D_anti, &tauijab, &Wmnij, 0, 0, 1, 1);
    global_dpd_->contract444(&D, &tauIjAb, &WMnIj, 0, 0, 1, 1);

    global_dpd_->buf4_close(&tauIJAB);
    global_dpd_->buf4_close(&tauijab);
    global_dpd_->buf4_close(&tauIjAb);

    global_dpd_->buf4_close(&D_anti);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_close(&WMNIJ);
    global_dpd_->buf4_close(&Wmnij);
    global_dpd_->buf4_close(&WMnIj);
  }
  else if(params_.ref == 2) { /*** UHF ***/
    global_dpd_->buf4_init(&WMNIJ, PSIF_CC_HBAR, 0, 2, 2, 2, 2, 0, "WMNIJ");
    global_dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, 0, 12, 12, 12, 12, 0, "Wmnij");
    global_dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, 0, 22, 22, 22, 22, 0, "WMnIj");

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
    global_dpd_->buf4_init(&tauIJAB, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    global_dpd_->contract444(&D, &tauIJAB, &WMNIJ, 0, 0, 1, 1);
    global_dpd_->buf4_close(&tauIJAB);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
    global_dpd_->buf4_init(&tauijab, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
    global_dpd_->contract444(&D, &tauijab, &Wmnij, 0, 0, 1, 1);
    global_dpd_->buf4_close(&tauijab);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->buf4_init(&tauIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
    global_dpd_->contract444(&D, &tauIjAb, &WMnIj, 0, 0, 1, 1);
    global_dpd_->buf4_close(&tauIjAb);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_close(&WMNIJ);
    global_dpd_->buf4_close(&Wmnij);
    global_dpd_->buf4_close(&WMnIj);
  }

}
}} // namespace psi::ccenergy
