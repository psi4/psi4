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

void c_sort(int reference)
{
  dpdbuf4 C, D;

  if(reference == 2) { /** UHF **/

    /* <ia||jb> = <ia|jb> - <ia|bj> = <ia|jb> - <ij|ba> */

    /*** AA ***/
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA|JB>");
    global_dpd_->buf4_copy(&C, PSIF_CC_CINTS, "C <IA||JB>");
    global_dpd_->buf4_close(&C);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ|AB>");
    global_dpd_->buf4_sort(&D, PSIF_CC_TMP0, psqr, 20, 20, "D <IJ|AB> (IB,JA)");
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_TMP0, 0, 20, 20, 20, 20, 0, "D <IJ|AB> (IB,JA)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
    global_dpd_->buf4_axpy(&D, &C, -1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&C);

    /* <IA||JB> (IA,BJ) (Wmbej.c) */
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
    global_dpd_->buf4_sort(&C, PSIF_CC_CINTS, pqsr, 20, 21, "C <IA||JB> (IA,BJ)");
    global_dpd_->buf4_close(&C);

    /*** BB ***/
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia|jb>");
    global_dpd_->buf4_copy(&C, PSIF_CC_CINTS, "C <ia||jb>");
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij|ab>");
    global_dpd_->buf4_sort(&D, PSIF_CC_TMP0, psqr, 30, 30, "D <ij|ab> (ib,ja)");
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_TMP0, 0, 30, 30, 30, 30, 0, "D <ij|ab> (ib,ja)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
    global_dpd_->buf4_axpy(&D, &C, -1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&C);

    /* <ia||jb> (ia,bj) (Wmbej.c) */
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
    global_dpd_->buf4_sort(&C, PSIF_CC_CINTS, pqsr, 30, 31, "C <ia||jb> (ia,bj)");
    global_dpd_->buf4_close(&C);

    /*** AB ***/

    /* <Ai|Bj> (iA,Bj) (Wmbej.c) */
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 26, 26, 26, 26, 0, "C <Ai|Bj>");
    global_dpd_->buf4_sort(&C, PSIF_CC_CINTS, qpsr, 27, 27, "C <iA|jB>");
    global_dpd_->buf4_sort(&C, PSIF_CC_CINTS, qprs, 27, 26, "C <Ai|Bj> (iA,Bj)");
    global_dpd_->buf4_close(&C);

    /* <Ia|Jb> (Ia,bJ) (Wmbej.c) */
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
    global_dpd_->buf4_sort(&C, PSIF_CC_CINTS, pqsr, 24, 25, "C <Ia|Jb> (Ia,bJ)");
    global_dpd_->buf4_close(&C);

  }
  else { /** RHF/ROHF **/
    /* <ia||jb> = <ia|jb> - <ia|bj> = <ia|jb> - <ij|ba> */
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->buf4_copy(&C, PSIF_CC_CINTS, "C <ia||jb>");
    global_dpd_->buf4_close(&C);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->buf4_sort(&D, PSIF_CC_TMP0, psqr, 10, 10, "D <ij|ab> (ib,ja)");
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&D, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    global_dpd_->buf4_axpy(&D, &C, -1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&C);

    /* <ia|jb> (bi,ja) */
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->buf4_sort(&C, PSIF_CC_CINTS, sprq, 11, 10, "C <ia|jb> (bi,ja)");
    global_dpd_->buf4_close(&C);

    /* <ia||jb> (bi,ja) */
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    global_dpd_->buf4_sort(&C, PSIF_CC_CINTS, sprq, 11, 10, "C <ia||jb> (bi,ja)");
    global_dpd_->buf4_close(&C);

    /* <ia|jb> (ia,bj) */
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->buf4_sort(&C, PSIF_CC_CINTS, pqsr, 10, 11, "C <ia|jb> (ia,bj)");
    global_dpd_->buf4_close(&C);

    /* <ia||jb> (ia,bj) (Wmbej.c) */
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    global_dpd_->buf4_sort(&C, PSIF_CC_CINTS, pqsr, 10, 11, "C <ia||jb> (ia,bj)");
    global_dpd_->buf4_close(&C);

    /* <ai|bj> (cchbar/Wabei_RHF.c) */
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->buf4_sort(&C, PSIF_CC_CINTS, qpsr, 11, 11, "C <ai|bj>");
    global_dpd_->buf4_close(&C);

  }

}

}} // namespace psi::cctransort
