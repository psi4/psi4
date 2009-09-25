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

    dpd_file2_init(&B, CC_OEI, irrep, 0, 1, lbl);

    sprintf(lbl, "ZIjAb[%d]", irrep);
    dpd_buf4_init(&Z, CC_MISC, irrep, 0, 5, 0, 5, 0, lbl);
    dpd_buf4_scm(&Z, 0.0);

    /* X(Ij,Ab) <-- b(I,C) <Cj|Ab> */
    sprintf(lbl, "XbAjI[%d]", irrep);
    dpd_buf4_init(&X, CC_TMP0, irrep, 5, 0, 5, 0, 0, lbl);
    dpd_buf4_init(&F, CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
    dpd_contract424(&F, &B, &X, 1, 1, 0, 1, 0);
    dpd_buf4_close(&F);
    sprintf(lbl, "XjIbA[%d]", irrep);
    dpd_buf4_sort(&X, CC_TMP0, rspq, 0, 5, lbl);
    dpd_buf4_close(&X);
    sprintf(lbl, "XjIbA[%d]", irrep);
    dpd_buf4_init(&X, CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
    if(spin == singlet) dpd_buf4_axpy(&X, &Z, 1);
    else dpd_buf4_axpy(&X, &Z, -1);
    sprintf(lbl, "XIjAb[%d]", irrep);
    dpd_buf4_sort(&X, CC_TMP0, qpsr, 0, 5, lbl);
    dpd_buf4_close(&X);
    sprintf(lbl, "XIjAb[%d]", irrep);
    dpd_buf4_init(&X, CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
    dpd_buf4_axpy(&X, &Z, 1);
    dpd_buf4_close(&X);

    /* X(Ij,Ab) <-- -<Ij|Ak> B(k,b) */
    sprintf(lbl, "XIjAb[%d]", irrep);
    dpd_buf4_init(&X, CC_TMP1, irrep, 0, 5, 0, 5, 0, lbl);
    dpd_buf4_init(&E, CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    dpd_contract424(&E, &B, &X, 1, 0, 0, -1, 0);
    dpd_buf4_close(&E);
    sprintf(lbl, "XjIbA[%d]", irrep);
    dpd_buf4_sort(&X, CC_TMP1, qpsr, 0, 5, lbl);

    dpd_buf4_axpy(&X, &Z, 1);
    dpd_buf4_close(&X);
    sprintf(lbl, "XjIbA[%d]", irrep);
    dpd_buf4_init(&X, CC_TMP1, irrep, 0, 5, 0, 5, 0, lbl);
    if(spin == singlet) dpd_buf4_axpy(&X, &Z, 1);
    else dpd_buf4_axpy(&X, &Z, -1);
    dpd_buf4_close(&X);

    /* Spin-adapt Z */
    if(spin == singlet) {
      sprintf(lbl, "ZIjbA[%d]", irrep);
      dpd_buf4_sort(&Z, CC_TMP1, pqsr, 0, 5, lbl);
      sprintf(lbl, "(ZIjAb - 1/2 ZIjbA)[%d]", irrep);
      dpd_buf4_copy(&Z, CC_MISC, lbl);
      dpd_buf4_close(&Z);
      sprintf(lbl, "(ZIjAb - 1/2 ZIjbA)[%d]", irrep);
      dpd_buf4_init(&Z, CC_MISC, irrep, 0, 5, 0, 5, 0, lbl);
      sprintf(lbl, "ZIjbA[%d]", irrep);
      dpd_buf4_init(&X, CC_TMP1, irrep, 0, 5, 0, 5, 0, lbl);
      dpd_buf4_axpy(&X, &Z, -0.5);
      dpd_buf4_close(&X);
    }
    else {
      sprintf(lbl, "ZIjbA[%d]", irrep);
      dpd_buf4_sort(&Z, CC_TMP1, pqsr, 0, 5, lbl);
      sprintf(lbl, "(ZIjAb + 1/2 ZIjbA)[%d]", irrep);
      dpd_buf4_copy(&Z, CC_MISC, lbl);
      dpd_buf4_close(&Z);
      sprintf(lbl, "(ZIjAb + 1/2 ZIjbA)[%d]", irrep);
      dpd_buf4_init(&Z, CC_MISC, irrep, 0, 5, 0, 5, 0, lbl);
      sprintf(lbl, "ZIjbA[%d]", irrep);
      dpd_buf4_init(&X, CC_TMP1, irrep, 0, 5, 0, 5, 0, lbl);
      dpd_buf4_axpy(&X, &Z, 0.5);
      dpd_buf4_close(&X);
    }
    dpd_buf4_close(&Z);

    dpd_file2_close(&B);

  }
  else if(params.ref == 2) { /** UHF **/

    sprintf(lbl, "BIA(%d)[%d]", root, irrep);
    dpd_file2_init(&B_A, CC_OEI, irrep, 0, 1, lbl);
    sprintf(lbl, "Bia(%d)[%d]", root, irrep);
    dpd_file2_init(&B_B, CC_OEI, irrep, 2, 3, lbl);

    /* X(IJ,AB) <-- b(I,C) <CJ||AB> */
    sprintf(lbl, "XIJAB[%d]", irrep);
    dpd_buf4_init(&X, CC_TMP0, irrep, 0, 7, 0, 7, 0, lbl);
    dpd_buf4_init(&F, CC_FINTS, 0, 21, 7, 21, 5, 1, "F <AI|BC>");
    dpd_contract244(&B_A, &F, &X, 1, 0, 0, 1, 0);
    dpd_buf4_close(&F);
    sprintf(lbl, "XJIAB[%d]", irrep);
    dpd_buf4_sort(&X, CC_TMP0, qprs, 0, 7, lbl);
    dpd_buf4_close(&X);

    sprintf(lbl, "XIJAB[%d]", irrep);
    dpd_buf4_init(&X1, CC_TMP0, irrep, 0, 7, 0, 7, 0, lbl);
    sprintf(lbl, "XJIAB[%d]", irrep);
    dpd_buf4_init(&X2, CC_TMP0, irrep, 0, 7, 0, 7, 0, lbl);
    dpd_buf4_axpy(&X2, &X1, -1);
    dpd_buf4_close(&X2);
    sprintf(lbl, "ZIJAB[%d]", irrep);
    dpd_buf4_init(&Z, CC_MISC, irrep, 0, 7, 2, 7, 0, lbl);
    dpd_buf4_scm(&Z, 0.0);
    dpd_buf4_axpy(&X1, &Z, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&X1);

    /* X(IJ,AB) <-- -<IJ||AK> B(K,B) */
    sprintf(lbl, "XIJAB[%d]", irrep);
    dpd_buf4_init(&X, CC_TMP1, irrep, 2, 5, 2, 5, 0, lbl);
    dpd_buf4_init(&E, CC_EINTS, 0, 2, 21, 2, 21, 0, "E <IJ||KA> (I>J,AK)");
    dpd_contract424(&E, &B_A, &X, 3, 0, 0, 1, 0);
    dpd_buf4_close(&E);
    sprintf(lbl, "XIJBA[%d]", irrep);
    dpd_buf4_sort(&X, CC_TMP1, pqsr, 2, 5, lbl);
    dpd_buf4_close(&X);

    sprintf(lbl, "XIJAB[%d]", irrep);
    dpd_buf4_init(&X1, CC_TMP1, irrep, 2, 5, 2, 5, 0, lbl);
    sprintf(lbl, "XIJBA[%d]", irrep);
    dpd_buf4_init(&X2, CC_TMP1, irrep, 2, 5, 2, 5, 0, lbl);
    dpd_buf4_axpy(&X2, &X1, -1);
    dpd_buf4_close(&X2);
    sprintf(lbl, "ZIJAB[%d]", irrep);
    dpd_buf4_init(&Z, CC_MISC, irrep, 2, 5, 2, 7, 0, lbl);
    dpd_buf4_axpy(&X1, &Z, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&X1);

    /* X(ij,ab) <-- b(i,c) <cj||ab> */
    sprintf(lbl, "Xijab[%d]", irrep);
    dpd_buf4_init(&X, CC_TMP0, irrep, 10, 17, 10, 17, 0, lbl);
    dpd_buf4_init(&F, CC_FINTS, 0, 31, 17, 31, 15, 1, "F <ai|bc>");
    dpd_contract244(&B_B, &F, &X, 1, 0, 0, 1, 0);
    dpd_buf4_close(&F);
    sprintf(lbl, "Xjiab[%d]", irrep);
    dpd_buf4_sort(&X, CC_TMP0, qprs, 10, 17, lbl);
    dpd_buf4_close(&X);

    sprintf(lbl, "Xijab[%d]", irrep);
    dpd_buf4_init(&X1, CC_TMP0, irrep, 10, 17, 10, 17, 0, lbl);
    sprintf(lbl, "Xjiab[%d]", irrep);
    dpd_buf4_init(&X2, CC_TMP0, irrep, 10, 17, 10, 17, 0, lbl);
    dpd_buf4_axpy(&X2, &X1, -1);
    dpd_buf4_close(&X2);
    sprintf(lbl, "Zijab[%d]", irrep);
    dpd_buf4_init(&Z, CC_MISC, irrep, 10, 17, 12, 17, 0, lbl);
    dpd_buf4_scm(&Z, 0.0);
    dpd_buf4_axpy(&X1, &Z, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&X1);

    /* X(ij,ab) <-- -<ij||ak> B(k,b) */
    sprintf(lbl, "Xijab[%d]", irrep);
    dpd_buf4_init(&X, CC_TMP1, irrep, 12, 15, 12, 15, 0, lbl);
    dpd_buf4_init(&E, CC_EINTS, 0, 12, 31, 12, 31, 0, "E <ij||ka> (i>j,ak)");
    dpd_contract424(&E, &B_B, &X, 3, 0, 0, 1, 0);
    dpd_buf4_close(&E);
    sprintf(lbl, "Xijba[%d]", irrep);
    dpd_buf4_sort(&X, CC_TMP1, pqsr, 12, 15, lbl);
    dpd_buf4_close(&X);

    sprintf(lbl, "Xijab[%d]", irrep);
    dpd_buf4_init(&X1, CC_TMP1, irrep, 12, 15, 12, 15, 0, lbl);
    sprintf(lbl, "Xijba[%d]", irrep);
    dpd_buf4_init(&X2, CC_TMP1, irrep, 12, 15, 12, 15, 0, lbl);
    dpd_buf4_axpy(&X2, &X1, -1);
    dpd_buf4_close(&X2);
    sprintf(lbl, "Zijab[%d]", irrep);
    dpd_buf4_init(&Z, CC_MISC, irrep, 12, 15, 12, 17, 0, lbl);
    dpd_buf4_axpy(&X1, &Z, 1);
    dpd_buf4_close(&Z);
    dpd_buf4_close(&X1);

    /* Z(Ij,Ab) <-- b(I,C) <Cj|Ab> */
    sprintf(lbl, "ZIjAb[%d]", irrep);
    dpd_buf4_init(&Z, CC_MISC, irrep, 22, 28, 22, 28, 0, lbl);
    dpd_buf4_init(&F, CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
    dpd_contract244(&B_A, &F, &Z, 1, 0, 0, 1, 0);
    dpd_buf4_close(&F);
    dpd_buf4_close(&Z);

    /* X(jI,bA) <-- - b(j,c) <cI|bA> */
    sprintf(lbl, "XIjAb[%d]", irrep);
    dpd_buf4_init(&X, CC_TMP0, irrep, 23, 29, 23, 29, 0, lbl);
    dpd_buf4_init(&F, CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
    dpd_contract244(&B_B, &F, &X, 1, 0, 0, 1, 0);
    dpd_buf4_close(&F);
    /* X(jI,bA) --> Z(Ij,Ab) */
    sprintf(lbl, "ZIjAb[%d]", irrep);
    dpd_buf4_sort_axpy(&X, CC_MISC, qpsr, 22, 28, lbl, 1);
    dpd_buf4_close(&X);

    /* Z(Ij,Ab) <-- -<Ij|Ak> b(k,b) */
    sprintf(lbl, "ZIjAb[%d]", irrep);
    dpd_buf4_init(&Z, CC_MISC, irrep, 22, 28, 22, 28, 0, lbl);
    dpd_buf4_init(&E, CC_EINTS, 0, 22, 26, 22, 26, 0, "E <Ij|Ak>");
    dpd_contract424(&E, &B_B, &Z, 3, 0, 0, -1, 1);
    dpd_buf4_close(&E);
    dpd_buf4_close(&Z);

    /* Z(Ij,Ab) <-- - b(K,A) <Kb|Ij> */
    sprintf(lbl, "ZIjAb[%d]", irrep);
    dpd_buf4_init(&Z, CC_MISC, irrep, 22, 28, 22, 28, 0, lbl);
    dpd_buf4_init(&E, CC_EINTS, 0, 24, 22, 24, 22, 0, "E <Ia|Jk>");
    dpd_contract244(&B_A, &E, &Z, 0, 0, 1, -1, 1);
    dpd_buf4_close(&E);
    dpd_buf4_close(&Z);

    dpd_file2_close(&B_A);
    dpd_file2_close(&B_B);
  }

}

}} // namespace psi::cis
