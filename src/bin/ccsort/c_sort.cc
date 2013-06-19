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
    \ingroup CCSORT
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccsort {

void c_sort(void)
{
  dpdbuf4 C, D;

  if(params.ref == 2) { /** UHF **/

    /* <ia||jb> = <ia|jb> - <ia|bj> = <ia|jb> - <ij|ba> */

    /*** AA ***/
    dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA|JB>");
    dpd_->buf4_copy(&C, PSIF_CC_CINTS, "C <IA||JB>");
    dpd_->buf4_close(&C);
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <IJ|AB>");
    dpd_->buf4_sort(&D, PSIF_CC_TMP0, psqr, 20, 20, "D <IJ|AB> (IB,JA)");
    dpd_->buf4_close(&D);
    dpd_->buf4_init(&D, PSIF_CC_TMP0, 0, 20, 20, 20, 20, 0, "D <IJ|AB> (IB,JA)");
    dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
    dpd_->buf4_axpy(&D, &C, -1);
    dpd_->buf4_close(&D);
    dpd_->buf4_close(&C);

    /* <IA||JB> (IA,BJ) (Wmbej.c) */
    dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
    dpd_->buf4_sort(&C, PSIF_CC_CINTS, pqsr, 20, 21, "C <IA||JB> (IA,BJ)");
    dpd_->buf4_close(&C);

    /*** BB ***/
    dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia|jb>");
    dpd_->buf4_copy(&C, PSIF_CC_CINTS, "C <ia||jb>");
    dpd_->buf4_close(&C);

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 15, 10, 15, 0, "D <ij|ab>");
    dpd_->buf4_sort(&D, PSIF_CC_TMP0, psqr, 30, 30, "D <ij|ab> (ib,ja)");
    dpd_->buf4_close(&D);
    dpd_->buf4_init(&D, PSIF_CC_TMP0, 0, 30, 30, 30, 30, 0, "D <ij|ab> (ib,ja)");
    dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
    dpd_->buf4_axpy(&D, &C, -1);
    dpd_->buf4_close(&D);
    dpd_->buf4_close(&C);

    /* <ia||jb> (ia,bj) (Wmbej.c) */
    dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
    dpd_->buf4_sort(&C, PSIF_CC_CINTS, pqsr, 30, 31, "C <ia||jb> (ia,bj)");
    dpd_->buf4_close(&C);

    /*** AB ***/

    /* <Ai|Bj> (iA,Bj) (Wmbej.c) */
    dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 26, 26, 26, 26, 0, "C <Ai|Bj>");
    dpd_->buf4_sort(&C, PSIF_CC_CINTS, qpsr, 27, 27, "C <iA|jB>");
    dpd_->buf4_sort(&C, PSIF_CC_CINTS, qprs, 27, 26, "C <Ai|Bj> (iA,Bj)");
    dpd_->buf4_close(&C);

    /* <Ia|Jb> (Ia,bJ) (Wmbej.c) */
    dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
    dpd_->buf4_sort(&C, PSIF_CC_CINTS, pqsr, 24, 25, "C <Ia|Jb> (Ia,bJ)");
    dpd_->buf4_close(&C);

  }
  else { /** RHF/ROHF **/
    /* <ia||jb> = <ia|jb> - <ia|bj> = <ia|jb> - <ij|ba> */
    dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    dpd_->buf4_copy(&C, PSIF_CC_CINTS, "C <ia||jb>");
    dpd_->buf4_close(&C);

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_->buf4_sort(&D, PSIF_CC_TMP0, psqr, 10, 10, "D <ij|ab> (ib,ja)");
    dpd_->buf4_close(&D);
    dpd_->buf4_init(&D, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    dpd_->buf4_axpy(&D, &C, -1);
    dpd_->buf4_close(&D);
    dpd_->buf4_close(&C);

    /* <ia|jb> (bi,ja) */
    dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    dpd_->buf4_sort(&C, PSIF_CC_CINTS, sprq, 11, 10, "C <ia|jb> (bi,ja)");
    dpd_->buf4_close(&C);

    /* <ia||jb> (bi,ja) */
    dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    dpd_->buf4_sort(&C, PSIF_CC_CINTS, sprq, 11, 10, "C <ia||jb> (bi,ja)");
    dpd_->buf4_close(&C);

    /* <ia|jb> (ia,bj) */
    dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    dpd_->buf4_sort(&C, PSIF_CC_CINTS, pqsr, 10, 11, "C <ia|jb> (ia,bj)");
    dpd_->buf4_close(&C);

    /* <ia||jb> (ia,bj) (Wmbej.c) */
    dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    dpd_->buf4_sort(&C, PSIF_CC_CINTS, pqsr, 10, 11, "C <ia||jb> (ia,bj)");
    dpd_->buf4_close(&C);

    /* <ai|bj> (cchbar/Wabei_RHF.c) */
    dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    dpd_->buf4_sort(&C, PSIF_CC_CINTS, qpsr, 11, 11, "C <ai|bj>");
    dpd_->buf4_close(&C);

  }

}

}} // namespace psi::ccsort
