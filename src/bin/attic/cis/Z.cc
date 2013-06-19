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
    \ingroup CIS
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstring>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cis {

void Z_build(int irrep, int root, enum Spin spin)
{
  char lbl[32];
  dpdfile2 B, B_A, B_B;
  dpdbuf4 X, X1, X2, Z, F, E;

  if(params.ref == 0) { /** RHF **/

    if(spin == singlet)
      sprintf(lbl, "BIA(%d)[%d] singlet", root, irrep);
    else 
      sprintf(lbl, "BIA(%d)[%d] triplet", root, irrep);

    global_dpd_->file2_init(&B, PSIF_CC_OEI, irrep, 0, 1, lbl);

    sprintf(lbl, "ZIjAb[%d]", irrep);
    global_dpd_->buf4_init(&Z, PSIF_CC_MISC, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_scm(&Z, 0.0);

    /* X(Ij,Ab) <-- b(I,C) <Cj|Ab> */
    sprintf(lbl, "XbAjI[%d]", irrep);
    global_dpd_->buf4_init(&X, PSIF_CC_TMP0, irrep, 5, 0, 5, 0, 0, lbl);
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    global_dpd_->contract424(&F, &B, &X, 1, 1, 0, 1, 0);
    global_dpd_->buf4_close(&F);
    sprintf(lbl, "XjIbA[%d]", irrep);
    global_dpd_->buf4_sort(&X, PSIF_CC_TMP0, rspq, 0, 5, lbl);
    global_dpd_->buf4_close(&X);
    sprintf(lbl, "XjIbA[%d]", irrep);
    global_dpd_->buf4_init(&X, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
    if(spin == singlet) global_dpd_->buf4_axpy(&X, &Z, 1);
    else global_dpd_->buf4_axpy(&X, &Z, -1);
    sprintf(lbl, "XIjAb[%d]", irrep);
    global_dpd_->buf4_sort(&X, PSIF_CC_TMP0, qpsr, 0, 5, lbl);
    global_dpd_->buf4_close(&X);
    sprintf(lbl, "XIjAb[%d]", irrep);
    global_dpd_->buf4_init(&X, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_axpy(&X, &Z, 1);
    global_dpd_->buf4_close(&X);

    /* X(Ij,Ab) <-- -<Ij|Ak> B(k,b) */
    sprintf(lbl, "XIjAb[%d]", irrep);
    global_dpd_->buf4_init(&X, PSIF_CC_TMP1, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    global_dpd_->contract424(&E, &B, &X, 1, 0, 0, -1, 0);
    global_dpd_->buf4_close(&E);
    sprintf(lbl, "XjIbA[%d]", irrep);
    global_dpd_->buf4_sort(&X, PSIF_CC_TMP1, qpsr, 0, 5, lbl);

    global_dpd_->buf4_axpy(&X, &Z, 1);
    global_dpd_->buf4_close(&X);
    sprintf(lbl, "XjIbA[%d]", irrep);
    global_dpd_->buf4_init(&X, PSIF_CC_TMP1, irrep, 0, 5, 0, 5, 0, lbl);
    if(spin == singlet) global_dpd_->buf4_axpy(&X, &Z, 1);
    else global_dpd_->buf4_axpy(&X, &Z, -1);
    global_dpd_->buf4_close(&X);

    /* Spin-adapt Z */
    if(spin == singlet) {
      sprintf(lbl, "ZIjbA[%d]", irrep);
      global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, pqsr, 0, 5, lbl);
      sprintf(lbl, "(ZIjAb - 1/2 ZIjbA)[%d]", irrep);
      global_dpd_->buf4_copy(&Z, PSIF_CC_MISC, lbl);
      global_dpd_->buf4_close(&Z);
      sprintf(lbl, "(ZIjAb - 1/2 ZIjbA)[%d]", irrep);
      global_dpd_->buf4_init(&Z, PSIF_CC_MISC, irrep, 0, 5, 0, 5, 0, lbl);
      sprintf(lbl, "ZIjbA[%d]", irrep);
      global_dpd_->buf4_init(&X, PSIF_CC_TMP1, irrep, 0, 5, 0, 5, 0, lbl);
      global_dpd_->buf4_axpy(&X, &Z, -0.5);
      global_dpd_->buf4_close(&X);
    }
    else {
      sprintf(lbl, "ZIjbA[%d]", irrep);
      global_dpd_->buf4_sort(&Z, PSIF_CC_TMP1, pqsr, 0, 5, lbl);
      sprintf(lbl, "(ZIjAb + 1/2 ZIjbA)[%d]", irrep);
      global_dpd_->buf4_copy(&Z, PSIF_CC_MISC, lbl);
      global_dpd_->buf4_close(&Z);
      sprintf(lbl, "(ZIjAb + 1/2 ZIjbA)[%d]", irrep);
      global_dpd_->buf4_init(&Z, PSIF_CC_MISC, irrep, 0, 5, 0, 5, 0, lbl);
      sprintf(lbl, "ZIjbA[%d]", irrep);
      global_dpd_->buf4_init(&X, PSIF_CC_TMP1, irrep, 0, 5, 0, 5, 0, lbl);
      global_dpd_->buf4_axpy(&X, &Z, 0.5);
      global_dpd_->buf4_close(&X);
    }
    global_dpd_->buf4_close(&Z);

    global_dpd_->file2_close(&B);

  }
  else if(params.ref == 2) { /** UHF **/

    sprintf(lbl, "BIA(%d)[%d]", root, irrep);
    global_dpd_->file2_init(&B_A, PSIF_CC_OEI, irrep, 0, 1, lbl);
    sprintf(lbl, "Bia(%d)[%d]", root, irrep);
    global_dpd_->file2_init(&B_B, PSIF_CC_OEI, irrep, 2, 3, lbl);

    /* X(IJ,AB) <-- b(I,C) <CJ||AB> */
    sprintf(lbl, "XIJAB[%d]", irrep);
    global_dpd_->buf4_init(&X, PSIF_CC_TMP0, irrep, 0, 7, 0, 7, 0, lbl);
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 21, 7, 21, 5, 1, "F <AI|BC>");
    global_dpd_->contract244(&B_A, &F, &X, 1, 0, 0, 1, 0);
    global_dpd_->buf4_close(&F);
    sprintf(lbl, "XJIAB[%d]", irrep);
    global_dpd_->buf4_sort(&X, PSIF_CC_TMP0, qprs, 0, 7, lbl);
    global_dpd_->buf4_close(&X);

    sprintf(lbl, "XIJAB[%d]", irrep);
    global_dpd_->buf4_init(&X1, PSIF_CC_TMP0, irrep, 0, 7, 0, 7, 0, lbl);
    sprintf(lbl, "XJIAB[%d]", irrep);
    global_dpd_->buf4_init(&X2, PSIF_CC_TMP0, irrep, 0, 7, 0, 7, 0, lbl);
    global_dpd_->buf4_axpy(&X2, &X1, -1);
    global_dpd_->buf4_close(&X2);
    sprintf(lbl, "ZIJAB[%d]", irrep);
    global_dpd_->buf4_init(&Z, PSIF_CC_MISC, irrep, 0, 7, 2, 7, 0, lbl);
    global_dpd_->buf4_scm(&Z, 0.0);
    global_dpd_->buf4_axpy(&X1, &Z, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&X1);

    /* X(IJ,AB) <-- -<IJ||AK> B(K,B) */
    sprintf(lbl, "XIJAB[%d]", irrep);
    global_dpd_->buf4_init(&X, PSIF_CC_TMP1, irrep, 2, 5, 2, 5, 0, lbl);
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 2, 21, 2, 21, 0, "E <IJ||KA> (I>J,AK)");
    global_dpd_->contract424(&E, &B_A, &X, 3, 0, 0, 1, 0);
    global_dpd_->buf4_close(&E);
    sprintf(lbl, "XIJBA[%d]", irrep);
    global_dpd_->buf4_sort(&X, PSIF_CC_TMP1, pqsr, 2, 5, lbl);
    global_dpd_->buf4_close(&X);

    sprintf(lbl, "XIJAB[%d]", irrep);
    global_dpd_->buf4_init(&X1, PSIF_CC_TMP1, irrep, 2, 5, 2, 5, 0, lbl);
    sprintf(lbl, "XIJBA[%d]", irrep);
    global_dpd_->buf4_init(&X2, PSIF_CC_TMP1, irrep, 2, 5, 2, 5, 0, lbl);
    global_dpd_->buf4_axpy(&X2, &X1, -1);
    global_dpd_->buf4_close(&X2);
    sprintf(lbl, "ZIJAB[%d]", irrep);
    global_dpd_->buf4_init(&Z, PSIF_CC_MISC, irrep, 2, 5, 2, 7, 0, lbl);
    global_dpd_->buf4_axpy(&X1, &Z, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&X1);

    /* X(ij,ab) <-- b(i,c) <cj||ab> */
    sprintf(lbl, "Xijab[%d]", irrep);
    global_dpd_->buf4_init(&X, PSIF_CC_TMP0, irrep, 10, 17, 10, 17, 0, lbl);
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 31, 17, 31, 15, 1, "F <ai|bc>");
    global_dpd_->contract244(&B_B, &F, &X, 1, 0, 0, 1, 0);
    global_dpd_->buf4_close(&F);
    sprintf(lbl, "Xjiab[%d]", irrep);
    global_dpd_->buf4_sort(&X, PSIF_CC_TMP0, qprs, 10, 17, lbl);
    global_dpd_->buf4_close(&X);

    sprintf(lbl, "Xijab[%d]", irrep);
    global_dpd_->buf4_init(&X1, PSIF_CC_TMP0, irrep, 10, 17, 10, 17, 0, lbl);
    sprintf(lbl, "Xjiab[%d]", irrep);
    global_dpd_->buf4_init(&X2, PSIF_CC_TMP0, irrep, 10, 17, 10, 17, 0, lbl);
    global_dpd_->buf4_axpy(&X2, &X1, -1);
    global_dpd_->buf4_close(&X2);
    sprintf(lbl, "Zijab[%d]", irrep);
    global_dpd_->buf4_init(&Z, PSIF_CC_MISC, irrep, 10, 17, 12, 17, 0, lbl);
    global_dpd_->buf4_scm(&Z, 0.0);
    global_dpd_->buf4_axpy(&X1, &Z, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&X1);

    /* X(ij,ab) <-- -<ij||ak> B(k,b) */
    sprintf(lbl, "Xijab[%d]", irrep);
    global_dpd_->buf4_init(&X, PSIF_CC_TMP1, irrep, 12, 15, 12, 15, 0, lbl);
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 12, 31, 12, 31, 0, "E <ij||ka> (i>j,ak)");
    global_dpd_->contract424(&E, &B_B, &X, 3, 0, 0, 1, 0);
    global_dpd_->buf4_close(&E);
    sprintf(lbl, "Xijba[%d]", irrep);
    global_dpd_->buf4_sort(&X, PSIF_CC_TMP1, pqsr, 12, 15, lbl);
    global_dpd_->buf4_close(&X);

    sprintf(lbl, "Xijab[%d]", irrep);
    global_dpd_->buf4_init(&X1, PSIF_CC_TMP1, irrep, 12, 15, 12, 15, 0, lbl);
    sprintf(lbl, "Xijba[%d]", irrep);
    global_dpd_->buf4_init(&X2, PSIF_CC_TMP1, irrep, 12, 15, 12, 15, 0, lbl);
    global_dpd_->buf4_axpy(&X2, &X1, -1);
    global_dpd_->buf4_close(&X2);
    sprintf(lbl, "Zijab[%d]", irrep);
    global_dpd_->buf4_init(&Z, PSIF_CC_MISC, irrep, 12, 15, 12, 17, 0, lbl);
    global_dpd_->buf4_axpy(&X1, &Z, 1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&X1);

    /* Z(Ij,Ab) <-- b(I,C) <Cj|Ab> */
    sprintf(lbl, "ZIjAb[%d]", irrep);
    global_dpd_->buf4_init(&Z, PSIF_CC_MISC, irrep, 22, 28, 22, 28, 0, lbl);
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
    global_dpd_->contract244(&B_A, &F, &Z, 1, 0, 0, 1, 0);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&Z);

    /* X(jI,bA) <-- - b(j,c) <cI|bA> */
    sprintf(lbl, "XIjAb[%d]", irrep);
    global_dpd_->buf4_init(&X, PSIF_CC_TMP0, irrep, 23, 29, 23, 29, 0, lbl);
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
    global_dpd_->contract244(&B_B, &F, &X, 1, 0, 0, 1, 0);
    global_dpd_->buf4_close(&F);
    /* X(jI,bA) --> Z(Ij,Ab) */
    sprintf(lbl, "ZIjAb[%d]", irrep);
    global_dpd_->buf4_sort_axpy(&X, PSIF_CC_MISC, qpsr, 22, 28, lbl, 1);
    global_dpd_->buf4_close(&X);

    /* Z(Ij,Ab) <-- -<Ij|Ak> b(k,b) */
    sprintf(lbl, "ZIjAb[%d]", irrep);
    global_dpd_->buf4_init(&Z, PSIF_CC_MISC, irrep, 22, 28, 22, 28, 0, lbl);
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 22, 26, 22, 26, 0, "E <Ij|Ak>");
    global_dpd_->contract424(&E, &B_B, &Z, 3, 0, 0, -1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&Z);

    /* Z(Ij,Ab) <-- - b(K,A) <Kb|Ij> */
    sprintf(lbl, "ZIjAb[%d]", irrep);
    global_dpd_->buf4_init(&Z, PSIF_CC_MISC, irrep, 22, 28, 22, 28, 0, lbl);
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 24, 22, 24, 22, 0, "E <Ia|Jk>");
    global_dpd_->contract244(&B_A, &E, &Z, 0, 0, 1, -1, 1);
    global_dpd_->buf4_close(&E);
    global_dpd_->buf4_close(&Z);

    global_dpd_->file2_close(&B_A);
    global_dpd_->file2_close(&B_B);
  }

}

}} // namespace psi::cis
