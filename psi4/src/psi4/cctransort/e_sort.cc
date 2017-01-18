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

#include <cstdio>
#include <cstdlib>
#include "psi4/libdpd/dpd.h"

namespace psi { namespace cctransort {

void e_sort(int reference)
{
  dpdbuf4 E;

  if(reference == 2) {  /** UHF **/
    /*** AA ***/
    /* <ij|ka> */
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 21, 0, 21, 0, 0, "E <AI|JK>");
    global_dpd_->buf4_sort(&E, PSIF_CC_EINTS, srqp, 0, 20, "E <IJ|KA>");
    global_dpd_->buf4_close(&E);

    /* <ij||ka> (i>j,ka) */
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 21, 0, 21, 0, 1, "E <AI|JK>");
    global_dpd_->buf4_sort(&E, PSIF_CC_EINTS, srqp, 2, 20, "E <IJ||KA> (I>J,KA)");
    global_dpd_->buf4_close(&E);

    /* <ij||ka> (i>j,ak) */
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 2, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
    global_dpd_->buf4_sort(&E, PSIF_CC_EINTS, pqsr, 2, 21, "E <IJ||KA> (I>J,AK)");
    global_dpd_->buf4_close(&E);

    /*** BB ***/
    /* <ij|ka> */
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 31, 10, 31, 10, 0, "E <ai|jk>");
    global_dpd_->buf4_sort(&E, PSIF_CC_EINTS, srqp, 10, 30, "E <ij|ka>");
    global_dpd_->buf4_close(&E);

    /* <ij||ka> (i>j,ka) */
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 31, 10, 31, 10, 1, "E <ai|jk>");
    global_dpd_->buf4_sort(&E, PSIF_CC_EINTS, srqp, 12, 30, "E <ij||ka> (i>j,ka)");
    global_dpd_->buf4_close(&E);

    /* <ij||ka> (i>j,ak) */
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 12, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
    global_dpd_->buf4_sort(&E, PSIF_CC_EINTS, pqsr, 12, 31, "E <ij||ka> (i>j,ak)");
    global_dpd_->buf4_close(&E);

    /*** AB ***/
    /* <iJ|kA> */
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 26, 22, 26, 22, 0, "E <Ai|Jk>");
    global_dpd_->buf4_sort(&E, PSIF_CC_EINTS, qrsp, 23, 27, "E <iJ|kA>");
    global_dpd_->buf4_close(&E);

    /* <Ij|Ak> */
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
    global_dpd_->buf4_sort(&E, PSIF_CC_EINTS, qpsr, 22, 26, "E <Ij|Ak>");
    global_dpd_->buf4_close(&E);

    /* <iJ|aK> */
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    global_dpd_->buf4_sort(&E, PSIF_CC_EINTS, qpsr, 23, 25, "E <iJ|aK>");
    global_dpd_->buf4_close(&E);

    /* <Ia|Jk> */
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    global_dpd_->buf4_sort(&E, PSIF_CC_EINTS, rspq, 24, 22, "E <Ia|Jk>");
    global_dpd_->buf4_close(&E);

  }
  else {  /** RHF/ROHF **/
    /* <ij|ka> */
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    global_dpd_->buf4_sort(&E, PSIF_CC_EINTS, srqp, 0, 10, "E <ij|ka>");
    global_dpd_->buf4_close(&E);

    /* <ij||ka> (i>j,ka) */
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 1, "E <ai|jk>");
    global_dpd_->buf4_sort(&E, PSIF_CC_EINTS, srqp, 2, 10, "E <ij||ka> (i>j,ka)");
    global_dpd_->buf4_close(&E);

    /* <ij|ka> (ij,ak) */
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    global_dpd_->buf4_sort(&E, PSIF_CC_EINTS, pqsr, 0, 11, "E <ij|ka> (ij,ak)");
    global_dpd_->buf4_close(&E);

    /* <ij||ka> (i>j,ak) */
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
    global_dpd_->buf4_sort(&E, PSIF_CC_EINTS, pqsr, 2, 11, "E <ij||ka> (i>j,ak)");
    global_dpd_->buf4_close(&E);

    /* <ia|jk> */
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    global_dpd_->buf4_sort(&E, PSIF_CC_EINTS, qpsr, 10, 0, "E <ia|jk>");
    global_dpd_->buf4_close(&E);

    /* <ij|ak> */
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    global_dpd_->buf4_sort(&E, PSIF_CC_EINTS, rspq, 0, 11, "E <ij|ak>");
    global_dpd_->buf4_close(&E);
  }
}

}} // namespace psi::cctransort
