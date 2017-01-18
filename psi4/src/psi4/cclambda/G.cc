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

void G_build(int L_irr) {
  dpdbuf4 LIJAB, Lijab, LiJaB, LIjAb, LijAB, LIJab;
  dpdbuf4 tIJAB, tijab, tiJaB, tIjAb, tijAB, tIJab;
  dpdfile2 GAE, Gae, GMI, Gmi;

  if(params.ref == 0) {
    global_dpd_->file2_init(&GMI, PSIF_CC_LAMBDA, L_irr, 0, 0, "GMI");

    /* T(Mj,Ab) * [ 2 L(Ij,Ab) - L(Ij,Ba) ] --> G(M,I) */
    global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract442(&tIjAb, &LIjAb, &GMI, 0, 0, 1, 0);
    global_dpd_->buf4_close(&tIjAb);
    global_dpd_->buf4_close(&LIjAb);

    global_dpd_->file2_close(&GMI);

    global_dpd_->file2_init(&GAE, PSIF_CC_LAMBDA, L_irr, 1, 1, "GAE");

    /* T(Ij,Eb) * [ 2 L(Ij,Ab) - L(Ij,Ba) ] --> G(A,E) */
    global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract442(&LIjAb, &tIjAb, &GAE, 2, 2, -1, 0);
    global_dpd_->buf4_close(&tIjAb);
    global_dpd_->buf4_close(&LIjAb);

    global_dpd_->file2_close(&GAE);
  }
  else if(params.ref == 1) { /** ROHF **/

    global_dpd_->file2_init(&GMI, PSIF_CC_LAMBDA, L_irr, 0, 0, "GMI");
    global_dpd_->file2_init(&Gmi, PSIF_CC_LAMBDA, L_irr, 0, 0, "Gmi");

    /* T2(MJ,AB) * L2(IJ,AB) --> G(M,I) */
    global_dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&LIJAB, PSIF_CC_LAMBDA, L_irr, 0, 7, 2, 7, 0, "LIJAB");
    global_dpd_->contract442(&tIJAB, &LIJAB, &GMI, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&tIJAB);
    global_dpd_->buf4_close(&LIJAB);

    /* T2(Mj,Ab) * L2(Ij,Ab) --> G(M,I) */
    global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->contract442(&tIjAb, &LIjAb, &GMI, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&tIjAb);
    global_dpd_->buf4_close(&LIjAb);

    /* T2(mj,ab) * L2(ij,ab) --> G(m,i) */
    global_dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tijab");
    global_dpd_->buf4_init(&Lijab, PSIF_CC_LAMBDA, L_irr, 0, 7, 2, 7, 0, "Lijab");
    global_dpd_->contract442(&tijab, &Lijab, &Gmi, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&tijab);
    global_dpd_->buf4_close(&Lijab);

    /* T2(mJ,aB) * L2(iJ,aB) --> G(m,i) */
    global_dpd_->buf4_init(&tiJaB, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    global_dpd_->buf4_init(&LiJaB, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->contract442(&tiJaB, &LiJaB, &Gmi, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&tiJaB);
    global_dpd_->buf4_close(&LiJaB);

    global_dpd_->file2_close(&Gmi);
    global_dpd_->file2_close(&GMI);



    global_dpd_->file2_init(&GAE, PSIF_CC_LAMBDA, L_irr, 1, 1, "GAE");
    global_dpd_->file2_init(&Gae, PSIF_CC_LAMBDA, L_irr, 1, 1, "Gae");

    /* T2(IJ,AB) * L2(IJ,EB) --> G(A,E) */
    global_dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&LIJAB, PSIF_CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "LIJAB");
    global_dpd_->contract442(&LIJAB, &tIJAB, &GAE, 2, 2, -1.0, 0.0);
    global_dpd_->buf4_close(&tIJAB);
    global_dpd_->buf4_close(&LIJAB);

    /* T2(Ij,Ab) * L2(Ij,Eb) --> G(A,E) */
    global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LIjAb");
    global_dpd_->contract442(&LIjAb, &tIjAb, &GAE, 2, 2, -1.0, 1.0);
    global_dpd_->buf4_close(&tIjAb);
    global_dpd_->buf4_close(&LIjAb);

    /* T2(ij,ab) * L2(ij,eb) --> G(a,e) */
    global_dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tijab");
    global_dpd_->buf4_init(&Lijab, PSIF_CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "Lijab");
    global_dpd_->contract442(&Lijab, &tijab, &Gae, 2, 2, -1.0, 0.0);
    global_dpd_->buf4_close(&tijab);
    global_dpd_->buf4_close(&Lijab);

    /* T2(iJ,aB) * L2(iJ,eB) --> G(a,e) */
    global_dpd_->buf4_init(&tiJaB, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tiJaB");
    global_dpd_->buf4_init(&LiJaB, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "LiJaB");
    global_dpd_->contract442(&LiJaB, &tiJaB, &Gae, 2, 2, -1.0, 1.0);
    global_dpd_->buf4_close(&tiJaB);
    global_dpd_->buf4_close(&LiJaB);

    global_dpd_->file2_close(&GAE);
    global_dpd_->file2_close(&Gae);
  }
  else if(params.ref == 2) { /** UHF **/

    global_dpd_->file2_init(&GMI, PSIF_CC_LAMBDA, L_irr, 0, 0, "GMI");
    global_dpd_->file2_init(&Gmi, PSIF_CC_LAMBDA, L_irr, 2, 2, "Gmi");

    /* T2(MJ,AB) * L2(IJ,AB) --> G(M,I) */
    global_dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&LIJAB, PSIF_CC_LAMBDA, L_irr, 0, 7, 2, 7, 0, "LIJAB");
    global_dpd_->contract442(&tIJAB, &LIJAB, &GMI, 0, 0, 1, 0);
    global_dpd_->buf4_close(&tIJAB);
    global_dpd_->buf4_close(&LIJAB);

    /* T2(Mj,Ab) * L2(Ij,Ab) --> G(M,I) */
    global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->contract442(&tIjAb, &LIjAb, &GMI, 0, 0, 1, 1);
    global_dpd_->buf4_close(&tIjAb);
    global_dpd_->buf4_close(&LIjAb);

    /* T2(mj,ab) * L2(ij,ab) --> G(m,i) */
    global_dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "tijab");
    global_dpd_->buf4_init(&Lijab, PSIF_CC_LAMBDA, L_irr, 10, 17, 12, 17, 0, "Lijab");
    global_dpd_->contract442(&tijab, &Lijab, &Gmi, 0, 0, 1, 0);
    global_dpd_->buf4_close(&tijab);
    global_dpd_->buf4_close(&Lijab);

    /* T2(mJ,aB) * L2(iJ,aB) --> G(m,i) */
    global_dpd_->buf4_init(&tiJaB, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    global_dpd_->buf4_init(&LiJaB, PSIF_CC_LAMBDA, L_irr, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->contract442(&tiJaB, &LiJaB, &Gmi, 0, 0, 1, 1);
    global_dpd_->buf4_close(&tiJaB);
    global_dpd_->buf4_close(&LiJaB);

    global_dpd_->file2_close(&Gmi);
    global_dpd_->file2_close(&GMI);



    global_dpd_->file2_init(&GAE, PSIF_CC_LAMBDA, L_irr, 1, 1, "GAE");
    global_dpd_->file2_init(&Gae, PSIF_CC_LAMBDA, L_irr, 3, 3, "Gae");

    /* T2(JI,BA) * L2(JI,BE) --> G(A,E) */
    global_dpd_->buf4_init(&tIJAB, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_init(&LIJAB, PSIF_CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "LIJAB");
    global_dpd_->contract442(&LIJAB, &tIJAB, &GAE, 3, 3, -1, 0);
    global_dpd_->buf4_close(&tIJAB);
    global_dpd_->buf4_close(&LIJAB);

    /* T2(jI,bA) * L2(jI,bE) --> G(A,E) */
    global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
    global_dpd_->buf4_init(&LIjAb, PSIF_CC_LAMBDA, L_irr, 23, 29, 23, 29, 0, "LiJaB");
    global_dpd_->contract442(&LIjAb, &tIjAb, &GAE, 3, 3, -1, 1);
    global_dpd_->buf4_close(&tIjAb);
    global_dpd_->buf4_close(&LIjAb);

    /* T2(ji,ba) * L2(ji,be) --> G(a,e) */
    global_dpd_->buf4_init(&tijab, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tijab");
    global_dpd_->buf4_init(&Lijab, PSIF_CC_LAMBDA, L_irr, 12, 15, 12, 17, 0, "Lijab");
    global_dpd_->contract442(&Lijab, &tijab, &Gae, 3, 3, -1, 0);
    global_dpd_->buf4_close(&tijab);
    global_dpd_->buf4_close(&Lijab);

    /* T2(Ji,Ba) * L2(Ji,Be) --> G(a,e) */
    global_dpd_->buf4_init(&tiJaB, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->buf4_init(&LiJaB, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "LIjAb");
    global_dpd_->contract442(&LiJaB, &tiJaB, &Gae, 3, 3, -1, 1);
    global_dpd_->buf4_close(&tiJaB);
    global_dpd_->buf4_close(&LiJaB);

    global_dpd_->file2_close(&GAE);
    global_dpd_->file2_close(&Gae);

  }

  return;
}



}} // namespace psi::cclambda
