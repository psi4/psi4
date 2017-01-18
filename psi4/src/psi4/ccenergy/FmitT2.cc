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

void CCEnergyWavefunction::FmitT2(void)
{
  dpdfile2 FMIt, Fmit;
  dpdbuf4 newtIJAB, newtijab, newtIjAb;
  dpdbuf4 tIJAB, tijab, tIjAb;
  dpdbuf4 t2;
  dpdbuf4 Z;

  if(params_.ref == 0) { /** RHF **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
    global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->file2_init(&FMIt, PSIF_CC_OEI, 0, 0, 0, "FMIt");
    global_dpd_->contract244(&FMIt, &tIjAb, &Z, 0, 0, 0, 1, 0);
    global_dpd_->file2_close(&FMIt);
    global_dpd_->buf4_close(&tIjAb);
    global_dpd_->buf4_init(&newtIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    global_dpd_->buf4_axpy(&Z, &newtIjAb, -1);
    global_dpd_->buf4_close(&newtIjAb);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_TAMPS, qpsr, 0, 5, "New tIjAb", -1);
    global_dpd_->buf4_close(&Z);
  }
  else if(params_.ref == 1) { /** ROHF **/

    global_dpd_->buf4_init(&newtIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "New tIJAB");
    global_dpd_->buf4_init(&newtijab, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "New tijab");
    global_dpd_->buf4_init(&newtIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

    global_dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
    global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");

    global_dpd_->file2_init(&FMIt, PSIF_CC_OEI, 0, 0, 0, "FMIt");
    global_dpd_->file2_init(&Fmit, PSIF_CC_OEI, 0, 0, 0, "Fmit");

    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    global_dpd_->contract424(&tIJAB, &FMIt, &t2, 1, 0, 1, -1, 0);
    global_dpd_->contract244(&FMIt, &tIJAB, &t2, 0, 0, 0, -1, 1);
    global_dpd_->buf4_axpy(&t2, &newtIJAB, 1);
    global_dpd_->buf4_close(&t2);

    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    global_dpd_->contract424(&tijab, &Fmit, &t2, 1, 0, 1, -1, 0);
    global_dpd_->contract244(&Fmit, &tijab, &t2, 0, 0, 0, -1, 1);
    global_dpd_->buf4_axpy(&t2, &newtijab, 1);
    global_dpd_->buf4_close(&t2);

    global_dpd_->contract424(&tIjAb, &Fmit, &newtIjAb, 1, 0, 1, -1, 1);
    global_dpd_->contract244(&FMIt, &tIjAb, &newtIjAb, 0, 0, 0, -1, 1);

    global_dpd_->file2_close(&FMIt);
    global_dpd_->file2_close(&Fmit);

    global_dpd_->buf4_close(&tIJAB);
    global_dpd_->buf4_close(&tijab);
    global_dpd_->buf4_close(&tIjAb);

    global_dpd_->buf4_close(&newtIJAB);
    global_dpd_->buf4_close(&newtijab);
    global_dpd_->buf4_close(&newtIjAb);
  }
  else if(params_.ref == 2) { /*** UHF ***/

    global_dpd_->buf4_init(&newtIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "New tIJAB");
    global_dpd_->buf4_init(&newtijab, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "New tijab");
    global_dpd_->buf4_init(&newtIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");

    global_dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
    global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");

    global_dpd_->file2_init(&FMIt, PSIF_CC_OEI, 0, 0, 0, "FMIt");
    global_dpd_->file2_init(&Fmit, PSIF_CC_OEI, 0, 2, 2, "Fmit");

    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    global_dpd_->contract424(&tIJAB, &FMIt, &t2, 1, 0, 1, -1, 0);
    global_dpd_->contract244(&FMIt, &tIJAB, &t2, 0, 0, 0, -1, 1);
    global_dpd_->buf4_axpy(&t2, &newtIJAB, 1);
    global_dpd_->buf4_close(&t2);

    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 10, 17, 10, 17, 0, "T (ij,a>b)");
    global_dpd_->contract424(&tijab, &Fmit, &t2, 1, 0, 1, -1, 0);
    global_dpd_->contract244(&Fmit, &tijab, &t2, 0, 0, 0, -1, 1);
    global_dpd_->buf4_axpy(&t2, &newtijab, 1);
    global_dpd_->buf4_close(&t2);

    global_dpd_->contract424(&tIjAb, &Fmit, &newtIjAb, 1, 0, 1, -1, 1);
    global_dpd_->contract244(&FMIt, &tIjAb, &newtIjAb, 0, 0, 0, -1, 1);

    global_dpd_->file2_close(&FMIt);
    global_dpd_->file2_close(&Fmit);

    global_dpd_->buf4_close(&tIJAB);
    global_dpd_->buf4_close(&tijab);
    global_dpd_->buf4_close(&tIjAb);

    global_dpd_->buf4_close(&newtIJAB);
    global_dpd_->buf4_close(&newtijab);
    global_dpd_->buf4_close(&newtIjAb);

  }
}
}} // namespace psi::ccenergy
