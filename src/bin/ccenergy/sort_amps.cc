/*! \file
    \ingroup CCENERGY
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

void sort_amps(void)
{
  dpdbuf4 t2, t2B;

  if(params.ref == 0) { /** RHF **/
    /* T(iJ,aB) */
    dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_sort(&t2, CC_TAMPS, qpsr, 0, 5, "tiJaB");
    dpd_buf4_close(&t2);

    /* TIjAb (IA,jb) */
    dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_sort(&t2, CC_TAMPS, prqs, 10, 10, "tIAjb");
    dpd_buf4_close(&t2);

    /* TIjAb (ij,JB) */
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_buf4_sort(&t2, CC_TAMPS, rspq, 10, 10, "tiaJB");
    dpd_buf4_close(&t2);

    /* TIjAb (Ib,jA) */
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_buf4_sort(&t2, CC_TAMPS, psrq, 10, 10, "tIbjA");
    dpd_buf4_close(&t2);

    /* TIjAb (jA,Ib) */
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    dpd_buf4_sort(&t2, CC_TAMPS, rspq, 10, 10, "tjAIb");
    dpd_buf4_close(&t2);

    /* 2 tIjAb - tIjBa */
    dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_scmcopy(&t2, CC_TAMPS, "2 tIjAb - tIjBa", 2);
    dpd_buf4_sort_axpy(&t2, CC_TAMPS, pqsr, 0, 5, "2 tIjAb - tIjBa", -1);
    dpd_buf4_close(&t2);

    /* 2 T(IA,jb) - t(IB,ja) */
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_buf4_scmcopy(&t2, CC_TAMPS, "2 tIAjb - tIBja", 2);
    dpd_buf4_close(&t2);
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "2 tIAjb - tIBja");
    dpd_buf4_init(&t2B, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    dpd_buf4_axpy(&t2B, &t2, -1);
    dpd_buf4_close(&t2B);
    dpd_buf4_close(&t2);

  }
  else if(params.ref == 1) { /** ROHF **/
    /* T(iJ,aB) */
    dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_sort(&t2, CC_TAMPS, qpsr, 0, 5, "tiJaB");
    dpd_buf4_close(&t2);

    /* TIJAB (IA,JB) */
    dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    dpd_buf4_sort(&t2, CC_TAMPS, prqs, 10, 10, "tIAJB");
    dpd_buf4_close(&t2);

    /* Tijab (ia,jb) */
    dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tijab");
    dpd_buf4_sort(&t2, CC_TAMPS, prqs, 10, 10, "tiajb");
    dpd_buf4_close(&t2);

    /* TIjAb (IA,jb) */
    dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_sort(&t2, CC_TAMPS, prqs, 10, 10, "tIAjb");
    dpd_buf4_close(&t2);

    /* TIjAb (ij,JB) */
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_buf4_sort(&t2, CC_TAMPS, rspq, 10, 10, "tiaJB");
    dpd_buf4_close(&t2);

    /* TIjAb (Ib,jA) */
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_buf4_sort(&t2, CC_TAMPS, psrq, 10, 10, "tIbjA");
    dpd_buf4_close(&t2);
    /* TIjAb (jA,Ib) */
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    dpd_buf4_sort(&t2, CC_TAMPS, rspq, 10, 10, "tjAIb");
    dpd_buf4_close(&t2);
  }
  else if(params.ref == 2) { /*** UHF ***/

    /* TIJAB (IA,JB) */
    dpd_buf4_init(&t2, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    dpd_buf4_sort(&t2, CC_TAMPS, prqs, 20, 20, "tIAJB");
    dpd_buf4_close(&t2);

    /* Tijab (ia,jb) */
    dpd_buf4_init(&t2, CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    dpd_buf4_sort(&t2, CC_TAMPS, prqs, 30, 30, "tiajb");
    dpd_buf4_close(&t2);

    /* TIjAb (IA,jb) */
    dpd_buf4_init(&t2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_buf4_sort(&t2, CC_TAMPS, prqs, 20, 30, "tIAjb");
    dpd_buf4_close(&t2);

    dpd_buf4_init(&t2, CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    /* TIjAb (jb,IA) */
    dpd_buf4_sort(&t2, CC_TAMPS, rspq, 30, 20, "tiaJB");
    /* TIjAb (Ib,jA) (Wmbej.c) */
    dpd_buf4_sort(&t2, CC_TAMPS, psrq, 24, 27, "tIbjA");
    dpd_buf4_close(&t2);

    /* TiJaB (iB,Ja) (Wmbej.c) */
    dpd_buf4_init(&t2, CC_TAMPS, 0, 24, 27, 24, 27, 0, "tIbjA");
    dpd_buf4_sort(&t2, CC_TAMPS, rspq, 27, 24, "tiBJa");
    dpd_buf4_close(&t2);

    /* T(iJ,aB) */
    dpd_buf4_init(&t2, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_buf4_sort(&t2, CC_TAMPS, qpsr, 23, 29, "tiJaB");
    dpd_buf4_close(&t2);

  }
}
}} // namespace psi::ccenergy
