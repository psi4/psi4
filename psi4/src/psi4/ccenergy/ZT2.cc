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
#include "psi4/libdpd/dpd.h"
#include "Params.h"
#include "ccwave.h"

namespace psi { namespace ccenergy {

void CCEnergyWavefunction::ZT2(void)
{
  dpdbuf4 ZIJMA, ZIJAM, Zijma, Zijam, ZIjMa, ZIjAm, Z;
  dpdbuf4 newtIJAB, newtijab, newtIjAb, T2;
  dpdfile2 tIA, tia, T1;
  dpdbuf4 t2, X;

  if(params_.ref == 0) { /** RHF **/
    global_dpd_->buf4_init(&X, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "X(Ab,Ij)");
    global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->buf4_init(&Z, PSIF_CC_MISC, 0, 10, 0, 10, 0, 0, "ZMbIj");
    global_dpd_->contract244(&T1, &Z, &X, 0, 0, 0, -1, 0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->file2_close(&T1);
    global_dpd_->buf4_sort_axpy(&X, PSIF_CC_TAMPS, rspq, 0, 5, "New tIjAb", 1);
    global_dpd_->buf4_sort_axpy(&X, PSIF_CC_TAMPS, srqp, 0, 5, "New tIjAb", 1);
    global_dpd_->buf4_close(&X);
  }
  else if(params_.ref == 1) { /** ROHF **/
    global_dpd_->buf4_init(&ZIJMA, PSIF_CC_MISC, 0, 2, 10, 2, 10, 0, "ZIJMA");
    global_dpd_->buf4_init(&ZIJAM, PSIF_CC_MISC, 0, 2, 11, 2, 11, 0, "ZIJAM");
    global_dpd_->buf4_init(&Zijma, PSIF_CC_MISC, 0, 2, 10, 2, 10, 0, "Zijma");
    global_dpd_->buf4_init(&Zijam, PSIF_CC_MISC, 0, 2, 11, 2, 11, 0, "Zijam");
    global_dpd_->buf4_init(&ZIjMa, PSIF_CC_MISC, 0, 0, 10, 0, 10, 0, "ZIjMa");
    global_dpd_->buf4_init(&ZIjAm, PSIF_CC_MISC, 0, 0, 11, 0, 11, 0, "ZIjAm");

    global_dpd_->buf4_init(&newtIJAB, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tIJAB");
    global_dpd_->buf4_init(&newtijab, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tijab");
    global_dpd_->buf4_init(&newtIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    global_dpd_->contract424(&ZIJAM, &tIA, &t2, 3, 0, 0, 1, 0);
    global_dpd_->contract244(&tIA, &ZIJMA, &t2, 0, 2, 1, -1, 1);
    global_dpd_->buf4_axpy(&t2, &newtIJAB, 1);
    global_dpd_->buf4_close(&t2);

    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    global_dpd_->contract424(&Zijam, &tia, &t2, 3, 0, 0, 1, 0);
    global_dpd_->contract244(&tia, &Zijma, &t2, 0, 2, 1, -1, 1);
    global_dpd_->buf4_axpy(&t2, &newtijab, 1);
    global_dpd_->buf4_close(&t2);

    global_dpd_->contract424(&ZIjAm, &tia, &newtIjAb, 3, 0, 0, -1, 1);
    global_dpd_->contract244(&tIA, &ZIjMa, &newtIjAb, 0, 2, 1, -1, 1);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

    global_dpd_->buf4_close(&newtIJAB);
    global_dpd_->buf4_close(&newtijab);
    global_dpd_->buf4_close(&newtIjAb);

    global_dpd_->buf4_close(&ZIJMA);
    global_dpd_->buf4_close(&ZIJAM);
    global_dpd_->buf4_close(&Zijma);
    global_dpd_->buf4_close(&Zijam);
    global_dpd_->buf4_close(&ZIjMa);
    global_dpd_->buf4_close(&ZIjAm);
  }
  else if(params_.ref == 2) { /*** UHF ***/

    global_dpd_->buf4_init(&ZIJMA, PSIF_CC_MISC, 0, 2, 20, 2, 20, 0, "ZIJMA");
    global_dpd_->buf4_init(&ZIJAM, PSIF_CC_MISC, 0, 2, 21, 2, 21, 0, "ZIJAM");
    global_dpd_->buf4_init(&Zijma, PSIF_CC_MISC, 0, 12, 30, 12, 30, 0, "Zijma");
    global_dpd_->buf4_init(&Zijam, PSIF_CC_MISC, 0, 12, 31, 12, 31, 0, "Zijam");
    global_dpd_->buf4_init(&ZIjMa, PSIF_CC_MISC, 0, 22, 24, 22, 24, 0, "ZIjMa");
    global_dpd_->buf4_init(&ZIjAm, PSIF_CC_MISC, 0, 22, 26, 22, 26, 0, "ZIjAm");

    global_dpd_->buf4_init(&newtIJAB, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tIJAB");
    global_dpd_->buf4_init(&newtijab, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "New tijab");
    global_dpd_->buf4_init(&newtIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    global_dpd_->contract424(&ZIJAM, &tIA, &t2, 3, 0, 0, 1, 0);
    global_dpd_->contract244(&tIA, &ZIJMA, &t2, 0, 2, 1, -1, 1);
    global_dpd_->buf4_axpy(&t2, &newtIJAB, 1);
    global_dpd_->buf4_close(&t2);

    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 12, 15, 12, 15, 0, "T (i>j,ab)");
    global_dpd_->contract424(&Zijam, &tia, &t2, 3, 0, 0, 1, 0);
    global_dpd_->contract244(&tia, &Zijma, &t2, 0, 2, 1, -1, 1);
    global_dpd_->buf4_axpy(&t2, &newtijab, 1);
    global_dpd_->buf4_close(&t2);

    global_dpd_->contract424(&ZIjAm, &tia, &newtIjAb, 3, 0, 0, -1, 1);
    global_dpd_->contract244(&tIA, &ZIjMa, &newtIjAb, 0, 2, 1, -1, 1);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

    global_dpd_->buf4_close(&newtIJAB);
    global_dpd_->buf4_close(&newtijab);
    global_dpd_->buf4_close(&newtIjAb);

    global_dpd_->buf4_close(&ZIJMA);
    global_dpd_->buf4_close(&ZIJAM);
    global_dpd_->buf4_close(&Zijma);
    global_dpd_->buf4_close(&Zijam);
    global_dpd_->buf4_close(&ZIjMa);
    global_dpd_->buf4_close(&ZIjAm);

  }
}
}} // namespace psi::ccenergy
