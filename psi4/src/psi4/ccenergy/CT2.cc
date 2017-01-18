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
#include "psi4/libqt/qt.h"
#include "Params.h"
#include "ccwave.h"

namespace psi { namespace ccenergy {

/* CT2(): Contributions of C-class integrals to T2.
**
** t(ij,ab) <--- P(ij) P(ab) t(m,a) t(i,e) <mb||je>
**
** This term is evaluated in two N^5 steps:
**  (1)  Y(mb,ji) = t(i,e) <mb||je>  (o^3 v^2)
**  (2)  t(ij,ab) <--- P(ij) P(ab) t(m,a) Y(mb,ji) (4 * o^3 v^2)
**
** Spin cases for UHF or ROHF orbitals:
** ------------------------------------
**                 *** AA ***
** + t(M,A) t(I,E) <MB||JE> - t(M,B) t(I,E) <MA||JE>
** - t(M,A) t(J,E) <MB||IE> + t(M,B) t(J,E) <MA||IE>
**
**                 *** BB ***
** + t(m,a) t(i,e) <mb||je> - t(m,b) t(i,e) <ma||je>
** - t(m,a) t(j,e) <me||ie> + t(m,b) t(j,e) <ma||ie>
**
**                 *** AB ***
** - t(M,A) t(I,E) <Mj|Eb> - t(m,b) t(I,E) <mA|jE>
** - t(M,A) t(j,e) <Mb|Ie> - t(m,b) t(j,e) <mI|eA>
**
** For the AA and BB spin cases, only the first term needs to be evaluated,
** while for the AB case, all four terms are different.
**
** This code was rewritten to eliminate all buf4_sort calls involving
** mixing indices between bra and ket (i.e., a "complex" sort).
** The current version requires six pairs of contractions (one each for AA
** and BB and four for AB) and twelve simple sorts (four each for AA, BB,
** and AB).
** TDC
** May 2000
*/

void CCEnergyWavefunction::CT2(void)
{
  dpdfile2 tIA, tia;
  dpdbuf4 Y, C, D, T2new, T2;

  if(params_.ref == 0) { /** RHF **/

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");

    /*** AB ***/

    /* C(mA|jE) * T(I,E) --> Y(mA,jI) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (mA,jI)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->contract424(&C, &tIA, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&C);

    /* T(m,b) * Y(mA,jI) --> T2(bA,jI) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "X(5,0)");
    global_dpd_->contract244(&tIA, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);

    /* T(bA,jI) --> Tnew(Ij,Ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, srqp, 0, 5, "X(0,5) 1");
    global_dpd_->buf4_close(&T2new);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 1");
    global_dpd_->buf4_init(&T2new, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);


    /* C(Mb|Ie) * T(j,e) --> Y(Mb,Ij) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (Mb,Ij)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->contract424(&C, &tIA, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&C);

    /* T(M,A) * Y(Mb,Ij) --> T2(Ab,Ij) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "X(5,0)");
    global_dpd_->contract244(&tIA, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);

    /* T(Ab,Ij) --> Tnew(Ij,Ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, rspq, 0, 5, "X(0,5) 1");
    global_dpd_->buf4_close(&T2new);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 1");
    global_dpd_->buf4_init(&T2new, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);


    /* D(Mb,jE) * T(I,E) --> Y(Mb,jI) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y(Mb,jI)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    global_dpd_->contract424(&D, &tIA, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);

    /* T(M,A) * Y(Mb,jI) --> T2(Ab,jI) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "X(5,0)");
    global_dpd_->contract244(&tIA, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);

    /* T2(Ab,jI) --> Tnew(Ij,Ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, srpq, 0, 5, "X(0,5) 1");
    global_dpd_->buf4_close(&T2new);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 1");
    global_dpd_->buf4_init(&T2new, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);


    /* D(mA,Ie) * T(j,e) --> Y(mA,Ij) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y(mA,Ij)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    global_dpd_->contract424(&D, &tIA, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);

    /* T(m,b) * Y(mA,Ij) --> T2(bA,Ij) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "X(5,0)");
    global_dpd_->contract244(&tIA, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);

    /* T2(bA,Ij) --> Tnew(Ij,Ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, rsqp, 0, 5, "X(0,5) 1");
    global_dpd_->buf4_close(&T2new);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 1");
    global_dpd_->buf4_init(&T2new, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);

    global_dpd_->file2_close(&tIA);

  }
  else if(params_.ref == 1) { /** ROHF **/

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    /*** AA ***/

    /* C(MB||JE) * T(I,E) --> Y(MB,JI) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (MB,JI)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    global_dpd_->contract424(&C, &tIA, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&C);

    /* T(M,A) * Y(MB,JI) --> T(AB,JI) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "X(5,0)");
    global_dpd_->contract244(&tIA, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);

    /* T(AB,JI) --> T(IJ,AB) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, srpq, 0, 5, "X(0,5) 1");
    global_dpd_->buf4_close(&T2new);

    /* P(IJ) P(AB) T2(IJ,AB) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 1");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qprs, 0, 5, "X(0,5) 2");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, pqsr, 0, 5, "X(0,5) 3");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qpsr, 0, 5, "X(0,5) 4");

    /* T2(IJ,AB) - T2(JI,AB) - T2(IJ,BA) - T2(JI,BA) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 2");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 3");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 4");
    global_dpd_->buf4_axpy(&T2, &T2new, 1);
    global_dpd_->buf4_close(&T2);

    /* T2(IJ,AB) --> T2new (IJ,AB) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "New tIJAB");
    global_dpd_->buf4_axpy(&T2new, &T2, 1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);


    /*** BB ***/

    /* C(mb||je) * T(i,e) --> Y(mb,ji) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (MB,JI)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia||jb>");
    global_dpd_->contract424(&C, &tia, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&C);

    /* T(m,a) * Y(mb,ji) --> T(ab,ji) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "X(5,0)");
    global_dpd_->contract244(&tia, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);

    /* T(ab,ji) --> T(ij,ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, srpq, 0, 5, "X(0,5) 1");
    global_dpd_->buf4_close(&T2new);

    /* P(ij) P(ab) T2(ij,ab) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 1");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qprs, 0, 5, "X(0,5) 2");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, pqsr, 0, 5, "X(0,5) 3");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qpsr, 0, 5, "X(0,5) 4");

    /* T2(ij,ab) - T2(ji,ab) - T2(ij,ba) - T2(ji,ba) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 2");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 3");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 4");
    global_dpd_->buf4_axpy(&T2, &T2new, 1);
    global_dpd_->buf4_close(&T2);

    /* T2(ij,ab) --> T2new (ij,ab) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "New tijab");
    global_dpd_->buf4_axpy(&T2new, &T2, 1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);


    /*** AB ***/

    /* C(mA|jE) * T(I,E) --> Y(mA,jI) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (mA,jI)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->contract424(&C, &tIA, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&C);

    /* T(m,b) * Y(mA,jI) --> T2(bA,jI) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "X(5,0)");
    global_dpd_->contract244(&tia, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);

    /* T(bA,jI) --> Tnew(Ij,Ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, srqp, 0, 5, "X(0,5) 1");
    global_dpd_->buf4_close(&T2new);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 1");
    global_dpd_->buf4_init(&T2new, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);


    /* C(Mb|Ie) * T(j,e) --> Y(Mb,Ij) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y (Mb,Ij)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 10, 10, 10, 10, 0, "C <ia|jb>");
    global_dpd_->contract424(&C, &tia, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&C);

    /* T(M,A) * Y(Mb,Ij) --> T2(Ab,Ij) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "X(5,0)");
    global_dpd_->contract244(&tIA, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);

    /* T(Ab,Ij) --> Tnew(Ij,Ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, rspq, 0, 5, "X(0,5) 1");
    global_dpd_->buf4_close(&T2new);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 1");
    global_dpd_->buf4_init(&T2new, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);


    /* D(Mb,jE) * T(I,E) --> Y(Mb,jI) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y(Mb,jI)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    global_dpd_->contract424(&D, &tIA, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);

    /* T(M,A) * Y(Mb,jI) --> T2(Ab,jI) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "X(5,0)");
    global_dpd_->contract244(&tIA, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);

    /* T2(Ab,jI) --> Tnew(Ij,Ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, srpq, 0, 5, "X(0,5) 1");
    global_dpd_->buf4_close(&T2new);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 1");
    global_dpd_->buf4_init(&T2new, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);


    /* D(mA,Ie) * T(j,e) --> Y(mA,Ij) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 10, 0, 10, 0, 0, "Y(mA,Ij)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D <ij|ab> (ib,ja)");
    global_dpd_->contract424(&D, &tia, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);

    /* T(m,b) * Y(mA,Ij) --> T2(bA,Ij) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "X(5,0)");
    global_dpd_->contract244(&tia, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);

    /* T2(bA,Ij) --> Tnew(Ij,Ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, rsqp, 0, 5, "X(0,5) 1");
    global_dpd_->buf4_close(&T2new);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 1");
    global_dpd_->buf4_init(&T2new, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);

    global_dpd_->file2_close(&tIA); global_dpd_->file2_close(&tia);

  }
  else if(params_.ref == 2) { /** UHF **/

    global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    /*** AA ***/

    /* C(MB||JE) * T(I,E) --> Y(MB,JI) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 20, 0, 20, 0, 0, "Y (MB,JI)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 20, 20, 20, 20, 0, "C <IA||JB>");
    global_dpd_->contract424(&C, &tIA, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&C);

    /* T(M,A) * Y(MB,JI) --> T(AB,JI) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "X(5,0)");
    global_dpd_->contract244(&tIA, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);

    /* T(AB,JI) --> T(IJ,AB) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, srpq, 0, 5, "X(0,5) 1");
    global_dpd_->buf4_close(&T2new);

    /* P(IJ) P(AB) T2(IJ,AB) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 1");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qprs, 0, 5, "X(0,5) 2");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, pqsr, 0, 5, "X(0,5) 3");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qpsr, 0, 5, "X(0,5) 4");

    /* T2(IJ,AB) - T2(JI,AB) - T2(IJ,BA) - T2(JI,BA) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 2");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 3");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(0,5) 4");
    global_dpd_->buf4_axpy(&T2, &T2new, 1);
    global_dpd_->buf4_close(&T2);

    /* T2(IJ,AB) --> T2new (IJ,AB) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "New tIJAB");
    global_dpd_->buf4_axpy(&T2new, &T2, 1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);


    /*** BB ***/

    /* C(mb||je) * T(i,e) --> Y(mb,ji) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 30, 10, 30, 10, 0, "Y (mb,ji)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 30, 30, 30, 30, 0, "C <ia||jb>");
    global_dpd_->contract424(&C, &tia, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&C);

    /* T(m,a) * Y(mb,ji) --> T(ab,ji) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 15, 10, 15, 10, 0, "X(15,10)");
    global_dpd_->contract244(&tia, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);

    /* T(ab,ji) --> T(ij,ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, srpq, 10, 15, "X(10,15) 1");
    global_dpd_->buf4_close(&T2new);

    /* P(ij) P(ab) T2(ij,ab) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 10, 15, 10, 15, 0, "X(10,15) 1");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qprs, 10, 15, "X(10,15) 2");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, pqsr, 10, 15, "X(10,15) 3");
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, qpsr, 10, 15, "X(10,15) 4");

    /* T2(ij,ab) - T2(ji,ab) - T2(ij,ba) - T2(ji,ba) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 10, 15, 10, 15, 0, "X(10,15) 2");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 10, 15, 10, 15, 0, "X(10,15) 3");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 10, 15, 10, 15, 0, "X(10,15) 4");
    global_dpd_->buf4_axpy(&T2, &T2new, 1);
    global_dpd_->buf4_close(&T2);

    /* T2(ij,ab) --> T2new (ij,ab) */
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 10, 15, 12, 17, 0, "New tijab");
    global_dpd_->buf4_axpy(&T2new, &T2, 1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);


    /*** AB ***/

    /* C(mA|jE) * T(I,E) --> Y(mA,jI) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 27, 23, 27, 23, 0, "Y (mA,jI)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 27, 27, 27, 27, 0, "C <iA|jB>");
    global_dpd_->contract424(&C, &tIA, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&C);

    /* T(m,b) * Y(mA,jI) --> T2(bA,jI) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 29, 23, 29, 23, 0, "X(29,23)");
    global_dpd_->contract244(&tia, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);

    /* T(bA,jI) --> Tnew(Ij,Ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, srqp, 22, 28, "X(22,28) 1");
    global_dpd_->buf4_close(&T2new);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 22, 28, 22, 28, 0, "X(22,28) 1");
    global_dpd_->buf4_init(&T2new, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);


    /* C(Mb|Ie) * T(j,e) --> Y(Mb,Ij) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 24, 22, 24, 22, 0, "Y (Mb,Ij)");
    global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
    global_dpd_->contract424(&C, &tia, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&C);

    /* T(M,A) * Y(Mb,Ij) --> T2(Ab,Ij) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 28, 22, 28, 22, 0, "X(28,22)");
    global_dpd_->contract244(&tIA, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);

    /* T(Ab,Ij) --> Tnew(Ij,Ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, rspq, 22, 28, "X(22,28) 1");
    global_dpd_->buf4_close(&T2new);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 22, 28, 22, 28, 0, "X(22,28) 1");
    global_dpd_->buf4_init(&T2new, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);


    /* D(Mb,jE) * T(I,E) --> Y(Mb,jI) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 24, 23, 24, 23, 0, "Y(Mb,jI)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 24, 27, 24, 27, 0, "D <Ij|Ab> (Ib,jA)");
    global_dpd_->contract424(&D, &tIA, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);

    /* T(M,A) * Y(Mb,jI) --> T2(Ab,jI) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 28, 23, 28, 23, 0, "X(28,23)");
    global_dpd_->contract244(&tIA, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);

    /* T2(Ab,jI) --> Tnew(Ij,Ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, srpq, 22, 28, "X(22,28) 1");
    global_dpd_->buf4_close(&T2new);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 22, 28, 22, 28, 0, "X(22,28) 1");
    global_dpd_->buf4_init(&T2new, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);


    /* D(mA,Ie) * T(j,e) --> Y(mA,Ij) */
    global_dpd_->buf4_init(&Y, PSIF_CC_TMP0, 0, 27, 22, 27, 22, 0, "Y(mA,Ij)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 27, 24, 27, 24, 0, "D <iJ|aB> (iB,Ja)");
    global_dpd_->contract424(&D, &tia, &Y, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);

    /* T(m,b) * Y(mA,Ij) --> T2(bA,Ij) */
    global_dpd_->buf4_init(&T2new, PSIF_CC_TMP0, 0, 29, 22, 29, 22, 0, "X(29,22)");
    global_dpd_->contract244(&tia, &Y, &T2new, 0, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y);

    /* T2(bA,Ij) --> Tnew(Ij,Ab) */
    global_dpd_->buf4_sort(&T2new, PSIF_CC_TMP0, rsqp, 22, 28, "X(22,28) 1");
    global_dpd_->buf4_close(&T2new);
    global_dpd_->buf4_init(&T2, PSIF_CC_TMP0, 0, 22, 28, 22, 28, 0, "X(22,28) 1");
    global_dpd_->buf4_init(&T2new, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");
    global_dpd_->buf4_axpy(&T2, &T2new, -1);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&T2new);

    global_dpd_->file2_close(&tIA);
    global_dpd_->file2_close(&tia);

  } /*** UHF ***/

}
}} // namespace psi::ccenergy
