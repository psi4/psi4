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
    \ingroup CCLAMBDA
    \brief Enter brief description of file here
*/
#include <cstdio>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cclambda {

/* WijmbL2(): Computes the contributions of the Wmnie HBAR matrix
** elements to the Lambda double de-excitation amplitude equations.
** These contributions are given in spin orbitals as:
**
** L_ij^ab = - P(ab) L_m^a Wijmb
**
** where Wijmb = Wmnie is defined as:
**
** Wmnie = <mn||ie> + t_i^f <mn||fe>
**
** [cf. Gauss and Stanton, JCP 103, 3561-3577 (1995)]
**
** All four spin cases are stored in (mn,ei) ordering (see Wmnie.c in
** the CCHBAR code).  The three L_ij^ab spin cases are computed as:
**
** L(IJ,AB) <-- - L(M,A) W(IJ,BM) + L(M,B) W(IJ,AM) (one unique
** contraction)
** L(ij,ab) <-- - L(m,a) W(ij,bm) + L(m,b) W(ij,am) (one unique
** contraction)
** L(Ij,Ab) <-- - L(M,A) W(Ij,bM) - L(m,b) W(jI,Am)
**
** TDC, July 2002
*/

void WijmbL2(int L_irr)
{
  dpdfile2 LIA, Lia;
  dpdbuf4 L2, newLijab, newLIJAB, newLIjAb;
  dpdbuf4 W, WMNIE, Wmnie, WMnIe, WmNiE;
  dpdbuf4 X1, X2, Z, Z1, Z2;

  /* RHS += -P(ab) Lma * Wijmb */
  if(params.ref == 0) { /** RHF **/

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, L_irr, 0, 5, 0, 5, 0, "Z(Ij,bA)");

    global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe (Mn,eI)");
    global_dpd_->contract424(&W, &LIA, &Z, 3, 0, 0, 1, 0);
    global_dpd_->buf4_close(&W);
    global_dpd_->file2_close(&LIA);

    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_LAMBDA, pqsr, 0, 5, "New LIjAb", -1);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_LAMBDA, qprs, 0, 5, "New LIjAb", -1);
    global_dpd_->buf4_close(&Z);
  }
  else if(params.ref == 1) { /** ROHF **/

    global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
    global_dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, L_irr, 0, 1, "Lia");

    global_dpd_->buf4_init(&WMNIE, PSIF_CC_HBAR, 0, 2, 11, 2, 11, 0, "WMNIE (M>N,EI)");
    global_dpd_->buf4_init(&X1, PSIF_CC_TMP1, L_irr, 2, 5, 2, 5, 0, "X(2,5) 1");
    global_dpd_->contract424(&WMNIE, &LIA, &X1, 3, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&WMNIE);
    global_dpd_->buf4_sort(&X1, PSIF_CC_TMP1, pqsr, 2, 5, "X(2,5) 2");
    global_dpd_->buf4_init(&X2, PSIF_CC_TMP1, L_irr, 2, 5, 2, 5, 0, "X(2,5) 2");
    global_dpd_->buf4_axpy(&X2, &X1, -1.0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_init(&newLIJAB, PSIF_CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "New LIJAB");
    global_dpd_->buf4_axpy(&X1, &newLIJAB, 1.0);
    global_dpd_->buf4_close(&newLIJAB);


    global_dpd_->buf4_init(&Wmnie, PSIF_CC_HBAR, 0, 2, 11, 2, 11, 0, "Wmnie (m>n,ei)");
    global_dpd_->buf4_init(&X1, PSIF_CC_TMP1, L_irr, 2, 5, 2, 5, 0, "X(2,5) 1");
    global_dpd_->contract424(&Wmnie, &Lia, &X1, 3, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&Wmnie);
    global_dpd_->buf4_sort(&X1, PSIF_CC_TMP1, pqsr, 2, 5, "X(2,5) 2");
    global_dpd_->buf4_init(&X2, PSIF_CC_TMP1, L_irr, 2, 5, 2, 5, 0, "X(2,5) 2");
    global_dpd_->buf4_axpy(&X2, &X1, -1.0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_init(&newLijab, PSIF_CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "New Lijab");
    global_dpd_->buf4_axpy(&X1, &newLijab, 1.0);
    global_dpd_->buf4_close(&newLijab);

    global_dpd_->buf4_init(&newLIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");

    global_dpd_->buf4_init(&WMnIe, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
    global_dpd_->contract244(&LIA, &WMnIe, &newLIjAb, 0, 2, 1, -1.0, 1.0);
    global_dpd_->buf4_close(&WMnIe);

    global_dpd_->buf4_init(&WmNiE, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "WmNiE (mN,Ei)");
    global_dpd_->buf4_sort(&WmNiE, PSIF_CC_TMP0, qprs, 0, 11, "WmNiE (Nm,Ei)");
    global_dpd_->buf4_close(&WmNiE);

    /* W(Nm,Ei) * L(i,b) --> L(Nm,Eb) */
    global_dpd_->buf4_init(&WmNiE, PSIF_CC_TMP0, 0, 0, 11, 0, 11, 0, "WmNiE (Nm,Ei)");
    global_dpd_->contract424(&WmNiE, &Lia, &newLIjAb, 3, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&WmNiE);

    global_dpd_->buf4_close(&newLIjAb);

    global_dpd_->file2_close(&Lia);
    global_dpd_->file2_close(&LIA);
  }
  else if(params.ref == 2) { /** UHF **/

    global_dpd_->file2_init(&LIA, PSIF_CC_LAMBDA, L_irr, 0, 1, "LIA");
    global_dpd_->file2_init(&Lia, PSIF_CC_LAMBDA, L_irr, 2, 3, "Lia");

    /** W(IJ,AM) L(M,B) --> Z(IJ,AB) **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, L_irr, 2, 5, 2, 5, 0, "Z'(IJ,AB)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 2, 21, 2, 21, 0, "WMNIE (M>N,EI)");
    global_dpd_->contract424(&W, &LIA, &Z, 3, 0, 0, 1, 0);
    global_dpd_->buf4_close(&W);
    /** Z(IJ,AB) --> Z(IJ,BA) **/
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP2, pqsr, 2, 5, "Z'(IJ,BA)");
    global_dpd_->buf4_close(&Z);
    /** Z(IJ,AB) = Z(IJ,AB) - Z(IJ,BA) **/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP2, L_irr, 2, 5, 2, 5, 0, "Z'(IJ,AB)");
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP2, L_irr, 2, 5, 2, 5, 0, "Z'(IJ,BA)");
    global_dpd_->buf4_axpy(&Z2, &Z1, -1);
    global_dpd_->buf4_close(&Z2);
    /** Z(IJ,AB) --> New L(IJ,AB) **/
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "New LIJAB");
    global_dpd_->buf4_axpy(&Z1, &L2, 1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&Z1);


    /** W(ij,am) L(m,b) --> Z(ij,ab) **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, L_irr, 12, 15, 12, 15, 0, "Z'(ij,ab)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 12, 31, 12, 31, 0, "Wmnie (m>n,ei)");
    global_dpd_->contract424(&W, &Lia, &Z, 3, 0, 0, 1, 0);
    global_dpd_->buf4_close(&W);
    /** Z(ij,ab) --> Z(ij,ba) **/
    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP2, pqsr, 12, 15, "Z'(ij,ba)");
    global_dpd_->buf4_close(&Z);
    /** Z(ij,ab) = Z(ij,ab) - Z(ij,ba) **/
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP2, L_irr, 12, 15, 12, 15, 0, "Z'(ij,ab)");
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP2, L_irr, 12, 15, 12, 15, 0, "Z'(ij,ba)");
    global_dpd_->buf4_axpy(&Z2, &Z1, -1);
    global_dpd_->buf4_close(&Z2);
    /** Z(ij,ab) --> New L(ij,ab) **/
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 12, 15, 12, 17, 0, "New Lijab");
    global_dpd_->buf4_axpy(&Z1, &L2, 1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&Z1);


    /** Z(jI,Ab) = W(jI,Am) L(m,b) **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, L_irr, 23, 28, 23, 28, 0, "Z(jI,Ab)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 23, 26, 23, 26, 0, "WmNiE (mN,Ei)");
    global_dpd_->contract424(&W, &Lia, &Z, 3, 0, 0, -1, 0);
    global_dpd_->buf4_close(&W);
    /** Z(jI,Ab) --> New L(Ij,Ab) **/
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_LAMBDA, qprs, 22, 28, "New LIjAb", 1);
    global_dpd_->buf4_close(&Z);

    /** Z(Ij,bA) = W(Ij,bM) L(M,A) **/
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP2, L_irr, 22, 29, 22, 29, 0, "Z(Ij,bA)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 22, 25, 22, 25, 0, "WMnIe (Mn,eI)");
    global_dpd_->contract424(&W, &LIA, &Z, 3, 0, 0, -1, 0);
    global_dpd_->buf4_close(&W);
    /** Z(Ij,bA) --> New L(Ij,Ab) **/
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_LAMBDA, pqsr, 22, 28, "New LIjAb", 1);
    global_dpd_->buf4_close(&Z);

    global_dpd_->file2_close(&Lia);
    global_dpd_->file2_close(&LIA);
  }
}


}} // namespace psi::cclambda
