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
    \ingroup CCHBAR
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <math.h>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cchbar {

/* Wmnie_build(): Computes all contributions to the Wmnie HBAR matrix
** elements.  The spin-orbital expression for this term is:
**
** Wmnie = <mn||ie> + t_i^f <mn||fe>
**
** [cf. Gauss and Stanton, JCP 103, 3561-3577 (1995)]
**
** TDC, June 2002
**
** The storage a naming convention for each of the four spin cases
** are as follows:
** Spin Case    Storage    Name
** ----------   ---------  -------
** WMNIE        (M>N,EI)   "WMNIE (M>N,EI)"
** Wmnie        (m>n,ei)   "Wmnie (m>n,ei)"
** WMnIe        (Mn,eI)    "WMnIe (Mn,eI)"
** WmNiE        (mN,Ei)    "WmNiE (mN,Ei)"
** ----------   ---------  -------
** WMNIE        (M>N,IE)   "WMNIE"
** Wmnie        (m>n,ie)   "Wmnie"
** WMnIe        (Mn,Ie)    "WMnIe"
** WmNiE        (mN,iE)    "WmNiE"
** -------------------------------
**
** Labels have been changed to the above.  Also, all 8 of the above plus
** plus the following for RHF are now stored by cchbar.
** 2WMnIe - WnMIe  (Mn,Ie)  "2WMnIe - WnMIe"
** WMnIe - 2WnMIe  (Mn,Ie)  "WMnIe - 2WnMIe"
** 2WMnIe - WnMIe  (Mn,eI)  "2WMnIe - WnMIe (Mn,eI)"
** WMnIe - 2WnMIe  (Mn,eI)  "WMnIe - 2WnMIe (Mn,eI)"
** RAK, April 2004
*/

void purge_Wmnie(void);

void Wmnie_build(void) {
  dpdbuf4 W, Wmnie, WMNIE, WMnIe, WmNiE, WMniE, WmNIe;
  dpdbuf4 E, Z;
  dpdbuf4 D, D_a;
  dpdfile2 t1, tIA, tia;

  if(params.ref == 0) { /** RHF **/

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    global_dpd_->buf4_copy(&E, PSIF_CC_HBAR, "WMnIe");
    global_dpd_->buf4_close(&E);

    /* D(Mn,Fe) * T(I,F) --> W(Mn,Ie) */
    global_dpd_->buf4_init(&WMnIe, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract244(&t1, &D, &WMnIe, 1, 2, 1, 1, 1);
    global_dpd_->file2_close(&t1);
    global_dpd_->buf4_close(&D);
    /* W(Mn,Ie) --> W(Mn,eI) */
    global_dpd_->buf4_sort(&WMnIe, PSIF_CC_HBAR, pqsr, 0, 11, "WMnIe (Mn,eI)");
    global_dpd_->buf4_sort(&WMnIe, PSIF_CC_HBAR, qpsr, 0, 11, "WMnIe (nM,eI)");
    global_dpd_->buf4_close(&WMnIe);

    /* make spin-combinations */
    global_dpd_->buf4_init(&WMnIe, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
    global_dpd_->buf4_copy(&WMnIe, PSIF_CC_HBAR, "WMnIe - 2WnMIe");
    global_dpd_->buf4_copy(&WMnIe, PSIF_CC_HBAR, "2WMnIe - WnMIe");
    global_dpd_->buf4_close(&WMnIe);

    global_dpd_->buf4_init(&WMnIe, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
    global_dpd_->buf4_sort_axpy(&WMnIe, PSIF_CC_HBAR, qprs, 0, 10, "WMnIe - 2WnMIe", -2.0);
    global_dpd_->buf4_close(&WMnIe);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "2WMnIe - WnMIe");
    global_dpd_->buf4_scm(&W, 2.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&WMnIe, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
    global_dpd_->buf4_sort_axpy(&WMnIe, PSIF_CC_HBAR, qprs, 0, 10, "2WMnIe - WnMIe", -1.0);
    global_dpd_->buf4_close(&WMnIe);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "2WMnIe - WnMIe");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, pqsr, 0, 11, "2WMnIe - WnMIe (Mn,eI)");
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe - 2WnMIe");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, pqsr, 0, 11, "WMnIe - 2WnMIe (Mn,eI)");
    global_dpd_->buf4_close(&W);

  }
  else if(params.ref == 1) { /** ROHF **/

    /* E(M>N,EI) --> W(M>N,EI) */
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 2, 10, 2, 10, 0, "E <ij||ka> (i>j,ka)");
    global_dpd_->buf4_sort(&E, PSIF_CC_HBAR, pqsr, 2, 11, "WMNIE (M>N,EI)");
    global_dpd_->buf4_sort(&E, PSIF_CC_HBAR, pqsr, 2, 11, "Wmnie (m>n,ei)");
    global_dpd_->buf4_close(&E);

    /* D(M>N,EF) * T(I,F) --> W(M>N,EI) */
    global_dpd_->buf4_init(&WMNIE, PSIF_CC_HBAR, 0, 2, 11, 2, 11, 0, "WMNIE (M>N,EI)");
    global_dpd_->buf4_init(&D_a, PSIF_CC_DINTS, 0, 2, 5, 2, 5,0, "D <ij||ab> (i>j,ab)");
    global_dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract424(&D_a,&t1,&WMNIE, 3, 1, 0, -1, 1);
    global_dpd_->file2_close(&t1);
    global_dpd_->buf4_close(&D_a);
    global_dpd_->buf4_close(&WMNIE);

    /* D(m>n,ef) * T(i,f) --> W(m>n,ei) */
    global_dpd_->buf4_init(&Wmnie, PSIF_CC_HBAR, 0, 2, 11, 2, 11, 0, "Wmnie (m>n,ei)");
    global_dpd_->buf4_init(&D_a, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
    global_dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->contract424(&D_a, &t1, &Wmnie, 3, 1, 0, -1, 1);
    global_dpd_->file2_close(&t1);
    global_dpd_->buf4_close(&D_a);
    global_dpd_->buf4_close(&Wmnie);

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 0, 10, 0, 10, 0, "E <ij|ka>");
    global_dpd_->buf4_copy(&E, PSIF_CC_TMP0, "WMnIe (Mn,Ie)");
    global_dpd_->buf4_copy(&E, PSIF_CC_TMP1, "WmNiE (mN,iE)");
    global_dpd_->buf4_close(&E);

    /* D(Mn,Fe) * T(I,F) --> W(Mn,Ie) */
    global_dpd_->buf4_init(&WMnIe, PSIF_CC_TMP0, 0, 0, 10, 0, 10, 0, "WMnIe (Mn,Ie)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract244(&t1, &D, &WMnIe, 1, 2, 1, 1, 1);
    global_dpd_->file2_close(&t1);
    global_dpd_->buf4_close(&D);
    /* W(Mn,Ie) --> W(Mn,eI) */
    global_dpd_->buf4_sort(&WMnIe, PSIF_CC_HBAR, pqsr, 0, 11, "WMnIe (Mn,eI)");
    global_dpd_->buf4_close(&WMnIe);

    /* D(mN,fE) * T(i,f) --> W(mN.iE) */
    global_dpd_->buf4_init(&WmNiE, PSIF_CC_TMP1, 0, 0, 10, 0, 10, 0, "WmNiE (mN,iE)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "tia");
    global_dpd_->contract244(&t1,&D,&WmNiE, 1, 2, 1, 1, 1);
    global_dpd_->file2_close(&t1);
    global_dpd_->buf4_close(&D);
    /* W(mN,iE) --> W(mN,Ei) */
    global_dpd_->buf4_sort(&WmNiE, PSIF_CC_HBAR, pqsr, 0, 11, "WmNiE (mN,Ei)");
    global_dpd_->buf4_close(&WmNiE);

    purge_Wmnie();

    /* also put "normal" sorted versions in CC_HBAR */
    global_dpd_->buf4_init(&WMNIE, PSIF_CC_HBAR, 0, 2, 11, 2, 11, 0, "WMNIE (M>N,EI)");
    global_dpd_->buf4_sort(&WMNIE, PSIF_CC_HBAR, pqsr, 2, 10, "WMNIE");
    global_dpd_->buf4_close(&WMNIE);
    global_dpd_->buf4_init(&Wmnie, PSIF_CC_HBAR, 0, 2, 11, 2, 11, 0, "Wmnie (m>n,ei)");
    global_dpd_->buf4_sort(&Wmnie, PSIF_CC_HBAR, pqsr, 2, 10, "Wmnie");
    global_dpd_->buf4_close(&Wmnie);
    global_dpd_->buf4_init(&WMnIe, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe (Mn,eI)");
    global_dpd_->buf4_sort(&WMnIe, PSIF_CC_HBAR, pqsr, 0, 10, "WMnIe");
    global_dpd_->buf4_close(&WMnIe);
    global_dpd_->buf4_init(&WmNiE, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "WmNiE (mN,Ei)");
    global_dpd_->buf4_sort(&WmNiE, PSIF_CC_HBAR, pqsr, 0, 10, "WmNiE");
    global_dpd_->buf4_close(&WmNiE);
  }
  else if(params.ref == 2) { /** UHF **/

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    /* <M>N||IE> --> W(M>N,EI) */
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 2, 20, 2, 20, 0, "E <IJ||KA> (I>J,KA)");
    global_dpd_->buf4_sort(&E, PSIF_CC_HBAR, pqsr, 2, 21, "WMNIE (M>N,EI)");
    global_dpd_->buf4_close(&E);

    /* <M>N||EF> T(I,F) --> W(M>N,EI) */
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 2, 21, 2, 21, 0, "WMNIE (M>N,EI)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <IJ||AB> (I>J,AB)");
    global_dpd_->contract424(&D, &tIA, &W, 3, 1, 0, -1, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);

    /* <m>n||ie> --> W(m>n,ei) */
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 12, 30, 12, 30, 0, "E <ij||ka> (i>j,ka)");
    global_dpd_->buf4_sort(&E, PSIF_CC_HBAR, pqsr, 12, 31, "Wmnie (m>n,ei)");
    global_dpd_->buf4_close(&E);

    /* <m>n||ef> T(i,f) --> W(m>n,ei) */
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 12, 31, 12, 31, 0, "Wmnie (m>n,ei)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 12, 15, 12, 15, 0, "D <ij||ab> (i>j,ab)");
    global_dpd_->contract424(&D, &tia, &W, 3, 1, 0, -1, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&W);


    /* <Mn|Ie> --> W(Mn,eI) */
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 24, 22, 24, 0, "E <Ij|Ka>");
    global_dpd_->buf4_sort(&E, PSIF_CC_HBAR, pqsr, 22, 25, "WMnIe (Mn,eI)");
    global_dpd_->buf4_close(&E);

    /* Z(nM,eI) = <nM|eF> T(I,F) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 23, 25, 23, 25, 0, "Z(nM,eI)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
    global_dpd_->contract424(&D, &tIA, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);
    /* Z(nM,eI) --> W(Mn,eI) */
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_HBAR, qprs, 22, 25, "WMnIe (Mn,eI)", 1);
    global_dpd_->buf4_close(&Z);

    /* <mN|iE> --> W(mN,Ei) */
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 23, 27, 23, 27, 0, "E <iJ|kA>");
    global_dpd_->buf4_sort(&E, PSIF_CC_HBAR, pqsr, 23, 26, "WmNiE (mN,Ei)");
    global_dpd_->buf4_close(&E);

    /* Z(Nm,Ei) = <Nm|Ef> T(i,f) */
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP1, 0, 22, 26, 22, 26, 0, "Z(Nm,Ei)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->contract424(&D, &tia, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);
    /* Z(Nm,Ei) --> W(mN,Ei) */
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_HBAR, qprs, 23, 26, "WmNiE (mN,Ei)", 1);
    global_dpd_->buf4_close(&Z);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

    /* also put "normal" sorted versions in CC_HBAR */
    global_dpd_->buf4_init(&WMNIE, PSIF_CC_HBAR, 0, 2, 21, 2, 21, 0, "WMNIE (M>N,EI)");
    global_dpd_->buf4_sort(&WMNIE, PSIF_CC_HBAR, pqsr, 2, 20, "WMNIE");
    global_dpd_->buf4_close(&WMNIE);
    global_dpd_->buf4_init(&Wmnie, PSIF_CC_HBAR, 0, 12, 31, 12, 31, 0, "Wmnie (m>n,ei)");
    global_dpd_->buf4_sort(&Wmnie, PSIF_CC_HBAR, pqsr, 12, 30, "Wmnie");
    global_dpd_->buf4_close(&Wmnie);
    global_dpd_->buf4_init(&WMnIe, PSIF_CC_HBAR, 0, 22, 25, 22, 25, 0, "WMnIe (Mn,eI)");
    global_dpd_->buf4_sort(&WMnIe, PSIF_CC_HBAR, pqsr, 22, 24, "WMnIe");
    global_dpd_->buf4_close(&WMnIe);
    global_dpd_->buf4_init(&WmNiE, PSIF_CC_HBAR, 0, 23, 26, 23, 26, 0, "WmNiE (mN,Ei)");
    global_dpd_->buf4_sort(&WmNiE, PSIF_CC_HBAR, pqsr, 23, 27, "WmNiE");
    global_dpd_->buf4_close(&WmNiE);
  }

  return;
}

}} // namespace psi::cchbar
