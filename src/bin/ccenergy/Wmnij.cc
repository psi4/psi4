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
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

void Wmnij_build(void)
{
  dpdbuf4 A_anti, A;
  dpdbuf4 WMNIJ, Wmnij, WMnIj, W;
  dpdfile2 tIA, tia;
  dpdbuf4 Eijka, Eijka_anti, Eaijk, Eaijk_anti;
  dpdbuf4 D_anti, D, tauIJAB, tauijab, tauIjAb;

  if(params.ref == 0) { /** RHF **/
    dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
    dpd_->buf4_copy(&A, PSIF_CC_HBAR, "WMnIj");
    dpd_->buf4_close(&A);
  }
  else if(params.ref == 1) { /** ROHF **/
    dpd_->buf4_init(&A_anti, PSIF_CC_AINTS, 0, 2, 2, 0, 0, 1, "A <ij|kl>");
    dpd_->buf4_copy(&A_anti, PSIF_CC_HBAR, "WMNIJ");
    dpd_->buf4_copy(&A_anti, PSIF_CC_HBAR, "Wmnij");
    dpd_->buf4_close(&A_anti);

    dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 0, 0, 0, 0, 0, "A <ij|kl>");
    dpd_->buf4_copy(&A, PSIF_CC_HBAR, "WMnIj");
    dpd_->buf4_close(&A);
  }
  else if(params.ref == 2) { /*** UHF ***/
    dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 2, 2, 0, 0, 1, "A <IJ|KL>");
    dpd_->buf4_copy(&A, PSIF_CC_HBAR, "WMNIJ");
    dpd_->buf4_close(&A);

    dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 12, 12, 10, 10, 1, "A <ij|kl>");
    dpd_->buf4_copy(&A, PSIF_CC_HBAR, "Wmnij");
    dpd_->buf4_close(&A);

    dpd_->buf4_init(&A, PSIF_CC_AINTS, 0, 22, 22, 22, 22, 0, "A <Ij|Kl>");
    dpd_->buf4_copy(&A, PSIF_CC_HBAR, "WMnIj");
    dpd_->buf4_close(&A);
  }

  if(params.ref == 0) { /** RHF **/
    dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");

    dpd_->buf4_init(&Eaijk, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    dpd_->contract244(&tIA, &Eaijk, &WMnIj, 1, 0, 1, 1, 1);
    dpd_->buf4_close(&Eaijk);

    dpd_->buf4_init(&Eijka, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_->contract424(&Eijka, &tIA, &WMnIj, 3, 1, 0, 1, 1);
    dpd_->buf4_close(&Eijka);

    dpd_->file2_close(&tIA);
    dpd_->buf4_close(&WMnIj);
  }
  else if(params.ref == 1) { /** ROHF **/  
    dpd_->buf4_init(&WMNIJ, PSIF_CC_HBAR, 0, 2, 0, 2, 2, 0, "WMNIJ");
    dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, 0, 2, 0, 2, 2, 0, "Wmnij");
    dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");

    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    dpd_->buf4_init(&Eijka_anti, PSIF_CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
    dpd_->buf4_init(&Eijka, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    dpd_->buf4_init(&Eaijk_anti, PSIF_CC_EINTS, 0, 11, 2, 11, 0, 1, "E <ai|jk>");
    dpd_->buf4_init(&Eaijk, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");

    dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 2, 0, 2, 0, 0, "W (MN,IJ)");
    dpd_->contract424(&Eijka_anti, &tIA, &W, 3, 1, 0, 1, 0);
    dpd_->contract244(&tIA, &Eaijk_anti, &W, 1, 0, 1, 1, 1);
    dpd_->buf4_axpy(&W, &WMNIJ, 1);
    dpd_->buf4_close(&W);

    dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 2, 0, 2, 0, 0, "W (MN,IJ)");
    dpd_->contract424(&Eijka_anti, &tia, &W, 3, 1, 0, 1, 0);
    dpd_->contract244(&tia, &Eaijk_anti, &W, 1, 0, 1, 1, 1);
    dpd_->buf4_axpy(&W, &Wmnij, 1);
    dpd_->buf4_close(&W);

    dpd_->contract424(&Eijka, &tia, &WMnIj, 3, 1, 0, 1, 1);
    dpd_->contract244(&tIA, &Eaijk, &WMnIj, 1, 0, 1, 1, 1);

    dpd_->buf4_close(&Eijka_anti);
    dpd_->buf4_close(&Eijka);
    dpd_->buf4_close(&Eaijk_anti);
    dpd_->buf4_close(&Eaijk);

    dpd_->file2_close(&tIA);
    dpd_->file2_close(&tia);

    dpd_->buf4_close(&WMNIJ);
    dpd_->buf4_close(&Wmnij);
    dpd_->buf4_close(&WMnIj);
  }
  else if(params.ref == 2) { /*** UHF ***/

    dpd_->buf4_init(&WMNIJ, PSIF_CC_HBAR, 0, 2, 0, 2, 2, 0, "WMNIJ");
    dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, 0, 12, 10, 12, 12, 0, "Wmnij");
    dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, 0, 22, 22, 22, 22, 0, "WMnIj");

    dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    dpd_->buf4_init(&Eijka, PSIF_CC_EINTS, 0, 2, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
    dpd_->buf4_init(&Eaijk, PSIF_CC_EINTS, 0, 21, 2, 21, 0, 1, "E <AI|JK>");
    dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 2, 0, 2, 0, 0, "W (MN,IJ)");
    dpd_->contract424(&Eijka, &tIA, &W, 3, 1, 0, 1, 0);
    dpd_->contract244(&tIA, &Eaijk, &W, 1, 0, 1, 1, 1);
    dpd_->buf4_axpy(&W, &WMNIJ, 1);
    dpd_->buf4_close(&W);
    dpd_->buf4_close(&Eijka);
    dpd_->buf4_close(&Eaijk);

    dpd_->buf4_init(&Eijka, PSIF_CC_EINTS, 0, 12, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
    dpd_->buf4_init(&Eaijk, PSIF_CC_EINTS, 0, 31, 12, 31, 10, 1, "E <ai|jk>");
    dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 12, 10, 12, 10, 0, "W (mn,ij)");
    dpd_->contract424(&Eijka, &tia, &W, 3, 1, 0, 1, 0);
    dpd_->contract244(&tia, &Eaijk, &W, 1, 0, 1, 1, 1);
    dpd_->buf4_axpy(&W, &Wmnij, 1);
    dpd_->buf4_close(&W);
    dpd_->buf4_close(&Eijka);
    dpd_->buf4_close(&Eaijk);

    dpd_->buf4_init(&Eijka, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    dpd_->buf4_init(&Eaijk, PSIF_CC_EINTS, 0, 26, 22, 26, 22, 0, "E <Ai|Jk>");
    dpd_->contract424(&Eijka, &tia, &WMnIj, 3, 1, 0, 1, 1);
    dpd_->contract244(&tIA, &Eaijk, &WMnIj, 1, 0, 1, 1, 1);
    dpd_->buf4_close(&Eijka);
    dpd_->buf4_close(&Eaijk);

    dpd_->file2_close(&tIA);
    dpd_->file2_close(&tia);

    dpd_->buf4_close(&WMNIJ);
    dpd_->buf4_close(&Wmnij);
    dpd_->buf4_close(&WMnIj);

  }

  if(params.ref == 0) { /** RHF **/
    dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_->buf4_init(&tauIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    dpd_->contract444(&D, &tauIjAb, &WMnIj, 0, 0, 1, 1);
    dpd_->buf4_close(&tauIjAb);
    dpd_->buf4_close(&D);
    dpd_->buf4_close(&WMnIj);
  }
  else if(params.ref == 1) { /** ROHF **/
    dpd_->buf4_init(&WMNIJ, PSIF_CC_HBAR, 0, 2, 2, 2, 2, 0, "WMNIJ");
    dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, 0, 2, 2, 2, 2, 0, "Wmnij");
    dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");

    dpd_->buf4_init(&D_anti, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");

    dpd_->buf4_init(&tauIJAB, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    dpd_->buf4_init(&tauijab, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauijab");
    dpd_->buf4_init(&tauIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");

    dpd_->contract444(&D_anti, &tauIJAB, &WMNIJ, 0, 0, 1, 1);
    dpd_->contract444(&D_anti, &tauijab, &Wmnij, 0, 0, 1, 1);
    dpd_->contract444(&D, &tauIjAb, &WMnIj, 0, 0, 1, 1);

    dpd_->buf4_close(&tauIJAB);
    dpd_->buf4_close(&tauijab);
    dpd_->buf4_close(&tauIjAb);

    dpd_->buf4_close(&D_anti);
    dpd_->buf4_close(&D);

    dpd_->buf4_close(&WMNIJ);
    dpd_->buf4_close(&Wmnij);
    dpd_->buf4_close(&WMnIj);
  }
  else if(params.ref == 2) { /*** UHF ***/
    dpd_->buf4_init(&WMNIJ, PSIF_CC_HBAR, 0, 2, 2, 2, 2, 0, "WMNIJ");
    dpd_->buf4_init(&Wmnij, PSIF_CC_HBAR, 0, 12, 12, 12, 12, 0, "Wmnij");
    dpd_->buf4_init(&WMnIj, PSIF_CC_HBAR, 0, 22, 22, 22, 22, 0, "WMnIj");

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
    dpd_->buf4_init(&tauIJAB, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
    dpd_->contract444(&D, &tauIJAB, &WMNIJ, 0, 0, 1, 1);
    dpd_->buf4_close(&tauIJAB);
    dpd_->buf4_close(&D);

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
    dpd_->buf4_init(&tauijab, PSIF_CC_TAMPS, 0, 12, 17, 12, 17, 0, "tauijab");
    dpd_->contract444(&D, &tauijab, &Wmnij, 0, 0, 1, 1);
    dpd_->buf4_close(&tauijab);
    dpd_->buf4_close(&D);

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_->buf4_init(&tauIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");
    dpd_->contract444(&D, &tauIjAb, &WMnIj, 0, 0, 1, 1);
    dpd_->buf4_close(&tauIjAb);
    dpd_->buf4_close(&D);

    dpd_->buf4_close(&WMNIJ);
    dpd_->buf4_close(&Wmnij);
    dpd_->buf4_close(&WMnIj);
  }

}
}} // namespace psi::ccenergy
