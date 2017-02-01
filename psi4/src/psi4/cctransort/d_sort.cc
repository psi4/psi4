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

#include "psi4/libdpd/dpd.h"

namespace psi { namespace cctransort {

void d_sort(int reference)
{
  dpdbuf4 D;

  if(reference == 2) { /*** UHF ***/
    /*** AA ***/
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 7, 0, 5, 1, "D <IJ|AB>");
    global_dpd_->buf4_copy(&D, PSIF_CC_DINTS, "D <IJ||AB> (I>J,A>B)");
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 5, 0, 5, 1, "D <IJ|AB>");
    global_dpd_->buf4_copy(&D, PSIF_CC_DINTS, "D <IJ||AB> (I>J,AB)");
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 7, 0, 5, 1, "D <IJ|AB>");
    global_dpd_->buf4_copy(&D, PSIF_CC_DINTS, "D <IJ||AB> (IJ,A>B)");
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 1, "D <IJ|AB>");
    global_dpd_->buf4_copy(&D, PSIF_CC_DINTS, "D <IJ||AB>");
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ||AB>");
    global_dpd_->buf4_sort(&D, PSIF_CC_DINTS, prqs, 20, 20, "D <IJ||AB> (IA,JB)");
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
    global_dpd_->buf4_sort(&D, PSIF_CC_DINTS, pqsr, 20, 21, "D <IJ||AB> (IA,BJ)");
    global_dpd_->buf4_close(&D);

    /*** BB ***/
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 12, 17, 10, 15, 1, "D <ij|ab>");
    global_dpd_->buf4_copy(&D, PSIF_CC_DINTS, "D <ij||ab> (i>j,a>b)");
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 12, 15, 10, 15, 1, "D <ij|ab>");
    global_dpd_->buf4_copy(&D, PSIF_CC_DINTS, "D <ij||ab> (i>j,ab)");
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 17, 10, 15, 1, "D <ij|ab>");
    global_dpd_->buf4_copy(&D, PSIF_CC_DINTS, "D <ij||ab> (ij,a>b)");
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 15, 10, 15, 1, "D <ij|ab>");
    global_dpd_->buf4_copy(&D, PSIF_CC_DINTS, "D <ij||ab>");
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij||ab>");
    global_dpd_->buf4_sort(&D, PSIF_CC_DINTS, prqs, 30, 30, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_sort(&D, PSIF_CC_DINTS, pqsr, 30, 31, "D <ij||ab> (ia,bj)");
    global_dpd_->buf4_close(&D);

    /*** AB ***/
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->buf4_sort(&D, PSIF_CC_DINTS, qpsr, 23, 29, "D <iJ|aB>");
    global_dpd_->buf4_sort(&D, PSIF_CC_DINTS, psrq, 24, 26, "D <Ij|Ab> (Ib,Aj)");
    global_dpd_->buf4_sort(&D, PSIF_CC_DINTS, prqs, 20, 30, "D <Ij|Ab> (IA,jb)");
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
    global_dpd_->buf4_sort(&D, PSIF_CC_DINTS, rspq, 30, 20, "D <Ij|Ab> (ia,JB)");
    global_dpd_->buf4_sort(&D, PSIF_CC_DINTS, pqsr, 20, 31, "D <Ij|Ab> (IA,bj)");
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
    global_dpd_->buf4_sort(&D, PSIF_CC_DINTS, pqsr, 30, 21, "D <Ij|Ab> (ia,BJ)");
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    global_dpd_->buf4_sort(&D, PSIF_CC_DINTS, psrq, 27, 25, "D <iJ|aB> (iB,aJ)");
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 24, 26, 24, 26, 0, "D <Ij|Ab> (Ib,Aj)");
    global_dpd_->buf4_sort(&D, PSIF_CC_DINTS, pqsr, 24, 27, "D <Ij|Ab> (Ib,jA)");
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 27, 25, 27, 25, 0, "D <iJ|aB> (iB,aJ)");
    global_dpd_->buf4_sort(&D, PSIF_CC_DINTS, pqsr, 27, 24, "D <iJ|aB> (iB,Ja)");
    global_dpd_->buf4_close(&D);
  }
  else {  /*** RHF/ROHF ***/
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 7, 0, 5, 1, "D <ij|ab>");
    global_dpd_->buf4_copy(&D, PSIF_CC_DINTS, "D <ij||ab> (i>j,a>b)");
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 5, 0, 5, 1, "D <ij|ab>");
    global_dpd_->buf4_copy(&D, PSIF_CC_DINTS, "D <ij||ab> (i>j,ab)");
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 7, 0, 5, 1, "D <ij|ab>");
    global_dpd_->buf4_copy(&D, PSIF_CC_DINTS, "D <ij||ab> (ij,a>b)");
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 1, "D <ij|ab>");
    global_dpd_->buf4_copy(&D, PSIF_CC_DINTS, "D <ij||ab>");
    global_dpd_->buf4_close(&D);

    /* <ij|ab> (ia,jb) */
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->buf4_sort(&D, PSIF_CC_DINTS, prqs, 10, 10, "D <ij|ab> (ia,jb)");
    global_dpd_->buf4_close(&D);

    /* <ij|ab> (ai,jb) */
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->buf4_sort(&D, PSIF_CC_DINTS, qprs, 11, 10, "D <ij|ab> (ai,jb)");
    global_dpd_->buf4_close(&D);

    /* <ij|ab> (aj,ib) */
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->buf4_sort(&D, PSIF_CC_DINTS, rqps, 11, 10, "D <ij|ab> (aj,ib)");
    global_dpd_->buf4_close(&D);

    /* <ij|ab> (bi,ja) */
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->buf4_sort(&D, PSIF_CC_DINTS, spqr, 11, 10, "D <ij|ab> (bi,ja)");
    global_dpd_->buf4_close(&D);

    /* <ij||ab> (ia,jb) */
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij||ab>");
    global_dpd_->buf4_sort(&D, PSIF_CC_DINTS, prqs, 10, 10, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_close(&D);

    /* <ij|ab> (ib,ja) */
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->buf4_sort(&D, PSIF_CC_DINTS, psrq, 10, 10, "D <ij|ab> (ib,ja)");
    global_dpd_->buf4_close(&D);

    /* <ij|ab> (ib,aj) */
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    global_dpd_->buf4_sort(&D, PSIF_CC_DINTS, pqsr, 10, 11, "D <ij|ab> (ib,aj)");
    global_dpd_->buf4_close(&D);

    /* <ij|ab> (ia,bj) */
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ia,jb)");
    global_dpd_->buf4_sort(&D, PSIF_CC_DINTS, pqsr, 10, 11, "D <ij|ab> (ia,bj)");
    global_dpd_->buf4_close(&D);

    /* <ij||ab> (ia,bj) */
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij||ab> (ia,jb)");
    global_dpd_->buf4_sort(&D, PSIF_CC_DINTS, pqsr, 10, 11, "D <ij||ab> (ia,bj)");
    global_dpd_->buf4_close(&D);

  }
}

}} // namespace psi::cctransort
