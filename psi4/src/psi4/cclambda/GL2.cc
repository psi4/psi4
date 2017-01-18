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

/* GaeL2(): Computes Gae three-body contributions of HBAR matrix
** elements to the Lambda doubles equations. These are written in
** spin-orbitals as:
**
** L_ij^ab <-- <ij||ae> Gbe - <ij||be> Gae
**
** where Gae = -1/2 t_mn^ef L_mn^af
**
** TDC, July 2002
*/

void GaeL2(int L_irr)
{
  dpdbuf4 L2, newLijab, newLIJAB, newLIjAb, newL2;
  dpdbuf4 D, Z;
  dpdfile2 GAE, Gae, G;
  dpdbuf4 X1, X2;

  /* RHS += P(ab)<ij||ae>Gbe */
  if(params.ref == 0) { /** RHF **/
    global_dpd_->file2_init(&G, PSIF_CC_LAMBDA, L_irr, 1, 1, "GAE");

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, L_irr, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract424(&D, &G, &Z, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_LAMBDA, qpsr, 0, 5, "New LIjAb", 1);
    global_dpd_->buf4_init(&newL2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    global_dpd_->buf4_axpy(&Z, &newL2, 1);
    global_dpd_->buf4_close(&newL2);
    global_dpd_->buf4_close(&Z);

    global_dpd_->file2_close(&G);
  }
  else if(params.ref == 1) { /** ROHF **/

    global_dpd_->file2_init(&GAE, PSIF_CC_LAMBDA, L_irr, 1, 1, "GAE");
    global_dpd_->file2_init(&Gae, PSIF_CC_LAMBDA, L_irr, 1, 1, "Gae");

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <ij||ab> (i>j,ab)");
    global_dpd_->buf4_init(&X1, PSIF_CC_TMP1, L_irr, 2, 5, 2, 5, 0, "X(2,5) 1");
    global_dpd_->contract424(&D, &GAE, &X1, 3, 1, 0, 1.0, 0.0);
    global_dpd_->buf4_init(&X2, PSIF_CC_TMP1, L_irr, 2, 5, 2, 5, 0, "X(2,5) 2");
    global_dpd_->contract244(&GAE, &D, &X2, 1, 2, 1, 1.0, 0.0);
    global_dpd_->buf4_axpy(&X1, &X2, 1.0);
    global_dpd_->buf4_close(&X1);
    global_dpd_->buf4_init(&newLIJAB, PSIF_CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "New LIJAB");
    global_dpd_->buf4_axpy(&X2, &newLIJAB, 1.0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&newLIJAB);


    global_dpd_->buf4_init(&X1, PSIF_CC_TMP1, L_irr, 2, 5, 2, 5, 0, "X(2,5) 1");
    global_dpd_->contract424(&D, &Gae, &X1, 3, 1, 0, 1.0, 0.0);
    global_dpd_->buf4_init(&X2, PSIF_CC_TMP1, L_irr, 2, 5, 2, 5, 0, "X(2,5) 2");
    global_dpd_->contract244(&Gae, &D, &X2, 1, 2, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_axpy(&X1, &X2, 1.0);
    global_dpd_->buf4_close(&X1);
    global_dpd_->buf4_init(&newLijab, PSIF_CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "New Lijab");
    global_dpd_->buf4_axpy(&X2, &newLijab, 1.0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&newLijab);


    global_dpd_->buf4_init(&newLIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract424(&D, &Gae, &newLIjAb, 3, 1, 0, 1.0, 1.0);
    global_dpd_->contract244(&GAE, &D, &newLIjAb, 1, 2, 1, 1.0, 1.0);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_close(&newLIjAb);

    global_dpd_->file2_close(&GAE);
    global_dpd_->file2_close(&Gae);
  }
  else if(params.ref == 2) { /** UHF **/

    global_dpd_->file2_init(&GAE, PSIF_CC_LAMBDA, L_irr, 1, 1, "GAE");
    global_dpd_->file2_init(&Gae, PSIF_CC_LAMBDA, L_irr, 3, 3, "Gae");

    /** X(IJ,AB) = <IJ||AE> G(B,E) **/
    global_dpd_->buf4_init(&X1, PSIF_CC_TMP2, L_irr, 2, 5, 2, 5, 0, "X(IJ,AB)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <IJ||AB> (I>J,AB)");
    global_dpd_->contract424(&D, &GAE, &X1, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);
    /** X(IJ,AB) --> X(IJ,BA) **/
    global_dpd_->buf4_sort(&X1, PSIF_CC_TMP2, pqsr, 2, 5, "X(IJ,BA)");
    /** X(IJ,AB) = X(IJ,AB) - X(IJ,BA) **/
    global_dpd_->buf4_init(&X2, PSIF_CC_TMP2, L_irr, 2, 5, 2, 5, 0, "X(IJ,BA)");
    global_dpd_->buf4_axpy(&X2, &X1, -1);
    global_dpd_->buf4_close(&X2);
    /** X(IJ,AB) --> New L(IJ,AB) **/
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 2, 5, 2, 7, 0, "New LIJAB");
    global_dpd_->buf4_axpy(&X1, &L2, 1);
    global_dpd_->buf4_close(&X1);
    global_dpd_->buf4_close(&L2);

    /** X(ij,ab) = <ij||ae> G(b,e) **/
    global_dpd_->buf4_init(&X1, PSIF_CC_TMP2, L_irr, 12, 15, 12, 15, 0, "X(ij,ab)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 12, 15, 12, 15, 0, "D <ij||ab> (i>j,ab)");
    global_dpd_->contract424(&D, &Gae, &X1, 3, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);
    /** X(ij,ab) --> X(ij,ba) **/
    global_dpd_->buf4_sort(&X1, PSIF_CC_TMP2, pqsr, 12, 15, "X(ij,ba)");
    /** X(ij,ab) = X(ij,ab) - X(ij,ba) **/
    global_dpd_->buf4_init(&X2, PSIF_CC_TMP2, L_irr, 12, 15, 12, 15, 0, "X(ij,ba)");
    global_dpd_->buf4_axpy(&X2, &X1, -1);
    global_dpd_->buf4_close(&X2);
    /** X(ij,ab) --> New L(ij,ab) **/
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 12, 15, 12, 17, 0, "New Lijab");
    global_dpd_->buf4_axpy(&X1, &L2, 1);
    global_dpd_->buf4_close(&X1);
    global_dpd_->buf4_close(&L2);

    /** New L(Ij,Ab) = <Ij|Ae> G(b,e) + <Ij|Eb> G(A,E) **/
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->contract424(&D, &Gae, &L2, 3, 1, 0, 1, 1);
    global_dpd_->contract244(&GAE, &D, &L2, 1, 2, 1, 1, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&L2);

    global_dpd_->file2_close(&GAE);
    global_dpd_->file2_close(&Gae);
  }

}

/* GmiL2(): Computes Gmi three-body contributions of HBAR matrix
** elements to the Lambda doubles equations. These are written in
** spin-orbitals as:
**
** L_ij^ab <-- - <im||ab> Gmj + <jm||ab> Gmi
**
** where Gmi = -1/2 t_mn^ef L_in^ef
**
** TDC, July 2002
*/

void GmiL2(int L_irr)
{

  dpdbuf4 L2, newLijab, newLIJAB, newLIjAb, newL2;
  dpdbuf4 D, Z;
  dpdfile2 GMI, Gmi, G;
  dpdbuf4 X1, X2;

  /* RHS -= P(ij) * <im||ab> * Gmj */
  if(params.ref == 0) { /** RHF **/

    global_dpd_->file2_init(&G, PSIF_CC_LAMBDA, L_irr, 0, 0, "GMI");

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, L_irr, 0, 5, 0, 5, 0, "Z(Ij,Ab)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract244(&G, &D, &Z, 0, 0, 0, -1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_LAMBDA, qpsr, 0, 5, "New LIjAb", 1);
    global_dpd_->buf4_init(&newL2, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");
    global_dpd_->buf4_axpy(&Z, &newL2, 1);
    global_dpd_->buf4_close(&newL2);
    global_dpd_->buf4_close(&Z);

    global_dpd_->file2_close(&G);
  }
  else if(params.ref == 1) { /** ROHF **/

    global_dpd_->file2_init(&GMI, PSIF_CC_LAMBDA, L_irr, 0, 0, "GMI");
    global_dpd_->file2_init(&Gmi, PSIF_CC_LAMBDA, L_irr, 0, 0, "Gmi");

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 7, 0, 7, 0, "D <ij||ab> (ij,a>b)");
    global_dpd_->buf4_init(&X1, PSIF_CC_TMP2, L_irr, 0, 7, 0, 7, 0, "X(0,7) 1");
    global_dpd_->contract424(&D, &GMI, &X1, 1, 0, 1, -1.0, 0.0);
    global_dpd_->buf4_init(&X2, PSIF_CC_TMP2, L_irr, 0, 7, 0, 7, 0, "X(0,7) 2");
    global_dpd_->contract244(&GMI, &D, &X2, 0, 0, 0, -1.0, 0.0);
    global_dpd_->buf4_axpy(&X1, &X2, 1.0);
    global_dpd_->buf4_close(&X1);
    global_dpd_->buf4_init(&newLIJAB, PSIF_CC_LAMBDA, L_irr, 0, 7, 2, 7, 0, "New LIJAB");
    global_dpd_->buf4_axpy(&X2, &newLIJAB, 1.0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&newLIJAB);


    global_dpd_->buf4_init(&X1, PSIF_CC_TMP2, L_irr, 0, 7, 0, 7, 0, "X(0,7) 1");
    global_dpd_->contract424(&D, &Gmi, &X1, 1, 0, 1, -1.0, 0.0);
    global_dpd_->buf4_init(&X2, PSIF_CC_TMP2, L_irr, 0, 7, 0, 7, 0, "X(0,7) 2");
    global_dpd_->contract244(&Gmi, &D, &X2, 0, 0, 0, -1.0, 0.0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_axpy(&X1, &X2, 1.0);
    global_dpd_->buf4_close(&X1);
    global_dpd_->buf4_init(&newLijab, PSIF_CC_LAMBDA, L_irr, 0, 7, 2, 7, 0, "New Lijab");
    global_dpd_->buf4_axpy(&X2, &newLijab, 1.0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&newLijab);

    global_dpd_->buf4_init(&newLIjAb, PSIF_CC_LAMBDA, L_irr, 0, 5, 0, 5, 0, "New LIjAb");

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract424(&D, &Gmi, &newLIjAb, 1, 0, 1, -1.0, 1.0);
    global_dpd_->contract244(&GMI, &D, &newLIjAb, 0, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_close(&newLIjAb);

    global_dpd_->file2_close(&Gmi);
    global_dpd_->file2_close(&GMI);
  }
  else if(params.ref == 2) { /** UHF **/

    global_dpd_->file2_init(&GMI, PSIF_CC_LAMBDA, L_irr, 0, 0, "GMI");
    global_dpd_->file2_init(&Gmi, PSIF_CC_LAMBDA, L_irr, 2, 2, "Gmi");

    /** X(IJ,AB) = - G(M,I) <MJ||AB> **/
    global_dpd_->buf4_init(&X1, PSIF_CC_TMP2, L_irr, 0, 7, 0, 7, 0, "X(IJ,AB) C");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 7, 0, 7, 0, "D <IJ||AB> (IJ,A>B)");
    global_dpd_->contract244(&GMI, &D, &X1, 0, 0, 0, -1, 0);
    global_dpd_->buf4_close(&D);
    /** X(IJ,AB) --> X(JI,AB) **/
    global_dpd_->buf4_sort(&X1, PSIF_CC_TMP2, qprs, 0, 7, "X(JI,AB)");
    /** X(IJ,AB) = X(IJ,AB) - X(JI,AB) **/
    global_dpd_->buf4_init(&X2, PSIF_CC_TMP2, L_irr, 0, 7, 0, 7, 0, "X(JI,AB)");
    global_dpd_->buf4_axpy(&X2, &X1, -1);
    global_dpd_->buf4_close(&X2);
    /** X(IJ,AB) --> New L(IJ,AB) **/
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 0, 7, 2, 7, 0, "New LIJAB");
    global_dpd_->buf4_axpy(&X1, &L2, 1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&X1);

    /** X(ij,ab) = - G(m,i) <mj||ab> **/
    global_dpd_->buf4_init(&X1, PSIF_CC_TMP2, L_irr, 10, 17, 10, 17, 0, "X(ij,ab) C");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 17, 10, 17, 0, "D <ij||ab> (ij,a>b)");
    global_dpd_->contract244(&Gmi, &D, &X1, 0, 0, 0, -1, 0);
    global_dpd_->buf4_close(&D);
    /** X(ij,ab) --> X(ji,ab) **/
    global_dpd_->buf4_sort(&X1, PSIF_CC_TMP2, qprs, 10, 17, "X(ji,ab)");
    /** X(ij,ab) = X(ij,ab) - X(ji,ab) **/
    global_dpd_->buf4_init(&X2, PSIF_CC_TMP2, L_irr, 10, 17, 10, 17, 0, "X(ji,ab)");
    global_dpd_->buf4_axpy(&X2, &X1, -1);
    global_dpd_->buf4_close(&X2);
    /** X(ij,ab) --> New L(ij,ab) **/
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 10, 17, 12, 17, 0, "New Lijab");
    global_dpd_->buf4_axpy(&X1, &L2, 1);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&X1);


    /* New L(Ij,Ab) <-- - <Im|Ab> G(m,j) - G(M,I) <Mj|Ab> **/
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMBDA, L_irr, 22, 28, 22, 28, 0, "New LIjAb");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    global_dpd_->contract424(&D, &Gmi, &L2, 1, 0, 1, -1, 1);
    global_dpd_->contract244(&GMI, &D, &L2, 0, 0, 0, -1, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&L2);

    global_dpd_->file2_close(&Gmi);
    global_dpd_->file2_close(&GMI);
  }

}

}} // namespace psi::cclambda
