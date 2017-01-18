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
#include "psi4/libpsio/psio.h"
#include "psi4/libqt/qt.h"
#include "Params.h"
#include "ccwave.h"

namespace psi { namespace ccenergy {

void CCEnergyWavefunction::cc2_WabeiT2(void) {

  int rowx, colx, rowz, colz, ab;
  int GX, GZ, Ge, Gi, Gj, hxbuf, hzbuf;
  dpdfile2 t1, tia, tIA;
  dpdbuf4 Z, W, X;
  dpdbuf4 t2, t2a, t2b, tIJAB, tijab, tIjAb;

  if(params_.ref == 0) { /** RHF **/

    global_dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "tIA");

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "CC2 ZAbIj");
    global_dpd_->buf4_init(&X, PSIF_CC2_HET1, 0, 5, 11, 5, 11, 0, "CC2 WAbEi");
    global_dpd_->contract244(&t1, &X, &Z, 1, 2, 1, 1, 0);

    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_TAMPS, rspq, 0, 5, "New tIjAb", 1);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_TAMPS, srqp, 0, 5, "New tIjAb", 1);
    global_dpd_->buf4_close(&Z);

    global_dpd_->file2_close(&t1);
  }

  else if(params_.ref == 1) { /** ROHF **/

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    /*** AA ***/
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 11, 7, 11, 7, 0, "CC2 WABEI (EI,A>B)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    global_dpd_->contract244(&tIA, &W, &t2, 1, 0, 0, 1, 0);
    global_dpd_->buf4_sort(&t2, PSIF_CC_TMP0, qprs, 0, 7, "T (JI,A>B)");
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_init(&t2a, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    global_dpd_->buf4_init(&t2b, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (JI,A>B)");
    global_dpd_->buf4_axpy(&t2b, &t2a, -1);
    global_dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "New tIJAB");
    global_dpd_->buf4_axpy(&t2a, &tIJAB, 1);
    global_dpd_->buf4_close(&tIJAB);
    global_dpd_->buf4_close(&t2b);
    global_dpd_->buf4_close(&t2a);
    global_dpd_->buf4_close(&W);

    /*** BB ***/
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 11, 7, 11, 7, 0, "CC2 Wabei (ei,a>b)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    global_dpd_->contract244(&tia, &W, &t2, 1, 0, 0, 1, 0);
    global_dpd_->buf4_sort(&t2, PSIF_CC_TMP0, qprs, 0, 7, "T (JI,A>B)");
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_init(&t2a, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    global_dpd_->buf4_init(&t2b, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (JI,A>B)");
    global_dpd_->buf4_axpy(&t2b, &t2a, -1);
    global_dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "New tijab");
    global_dpd_->buf4_axpy(&t2a, &tijab, 1);
    global_dpd_->buf4_close(&tijab);
    global_dpd_->buf4_close(&t2b);
    global_dpd_->buf4_close(&t2a);
    global_dpd_->buf4_close(&W);

    /*** AB ***/
    global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 11, 5, 11, 5, 0, "CC2 WAbEi (Ei,Ab)");
    global_dpd_->contract244(&tIA, &W, &tIjAb, 1, 0, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&tIjAb);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "CC2 ZjIbA");
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 11, 5, 11, 5, 0, "CC2 WaBeI (eI,aB)");
    global_dpd_->contract244(&tia, &W, &Z, 1, 0, 0, 1, 0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_TAMPS, qpsr, 0, 5, "New tIjAb", 1);
    global_dpd_->buf4_close(&Z);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

  }

  else if(params_.ref == 2) { /*** UHF ***/

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    /*** AA ***/
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 21, 7, 21, 7, 0, "CC2 WABEI (EI,A>B)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    global_dpd_->contract244(&tIA, &W, &t2, 1, 0, 0, 1, 0);
    global_dpd_->buf4_sort(&t2, PSIF_CC_TMP0, qprs, 0, 7, "T (JI,A>B)");
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_init(&t2a, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
    global_dpd_->buf4_init(&t2b, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (JI,A>B)");
    global_dpd_->buf4_axpy(&t2b, &t2a, -1);
    global_dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "New tIJAB");
    global_dpd_->buf4_axpy(&t2a, &tIJAB, 1);
    global_dpd_->buf4_close(&t2b);
    global_dpd_->buf4_close(&t2a);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&tIJAB);

    /*** BB ***/
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 31, 17, 31, 17, 0, "CC2 Wabei (ei,a>b)");
    global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 10, 17, 10, 17, 0, "T (ij,a>b)");
    global_dpd_->contract244(&tia, &W, &t2, 1, 0, 0, 1, 0);
    global_dpd_->buf4_sort(&t2, PSIF_CC_TMP0, qprs, 10, 17, "T (ji,a>b)");
    global_dpd_->buf4_close(&t2);
    global_dpd_->buf4_init(&t2a, PSIF_CC_TMP0, 0, 10, 17, 10, 17, 0, "T (ij,a>b)");
    global_dpd_->buf4_init(&t2b, PSIF_CC_TMP0, 0, 10, 17, 10, 17, 0, "T (ji,a>b)");
    global_dpd_->buf4_axpy(&t2b, &t2a, -1);
    global_dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "New tijab");
    global_dpd_->buf4_axpy(&t2a, &tijab, 1);
    global_dpd_->buf4_close(&t2b);
    global_dpd_->buf4_close(&t2a);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&tijab);

    /*** AB ***/
    global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 26, 28, 26, 28, 0, "CC2 WAbEi (Ei,Ab)");
    global_dpd_->contract244(&tIA, &W, &tIjAb, 1, 0, 0, 1, 1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&tIjAb);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 23, 29, 23, 29, 0, "CC2 ZjIbA");
    global_dpd_->buf4_init(&W, PSIF_CC2_HET1, 0, 25, 29, 25, 29, 0, "CC2 WaBeI (eI,aB)");
    global_dpd_->contract244(&tia, &W, &Z, 1, 0, 0, 1, 0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_TAMPS, qpsr, 22, 28, "New tIjAb", 1);
    global_dpd_->buf4_close(&Z);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

  }

}
}} // namespace psi::ccenergy
