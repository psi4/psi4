/*! \file
    \ingroup MP2
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#define EXTERN
#include "globals.h"

namespace psi{ namespace mp2{

void sort_amps(void)
{
  dpdbuf4 T;
  dpdbuf4 T2AB1, T2AB2;

  if(params.ref == 0) { /** RHF **/

    dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_scmcopy(&T, CC_TAMPS, "2 tIjAb - tIjBa", 2);
    dpd_buf4_sort_axpy(&T, CC_TAMPS, pqsr, 0, 5, "2 tIjAb - tIjBa", -1);
    dpd_buf4_close(&T);

    dpd_buf4_init(&T, CC_TAMPS, 0, 2, 7, 0, 5, 1, "tIjAb");
    dpd_buf4_copy(&T, CC_TAMPS, "tIJAB");
    dpd_buf4_copy(&T, CC_TAMPS, "tijab");
    dpd_buf4_close(&T);

    /* T(iJ,aB) */
    dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_sort(&T, CC_TAMPS, qpsr, 0, 5, "tiJaB");
    dpd_buf4_close(&T);

    /* TIjAb (IA,jb) */
    dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_sort(&T, CC_TAMPS, prqs, 10, 10, "tIAjb");
    dpd_buf4_close(&T);

    /* TIjAb (ij,JB) */
    dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_buf4_sort(&T, CC_TAMPS, rspq, 10, 10, "tiaJB");
    dpd_buf4_close(&T);

    /* TIjAb (Ib,jA) */
    dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_buf4_sort(&T, CC_TAMPS, psrq, 10, 10, "tIbjA");
    dpd_buf4_close(&T);

    /* TIjAb (jA,Ib) */
    dpd_buf4_init(&T, CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    dpd_buf4_sort(&T, CC_TAMPS, rspq, 10, 10, "tjAIb");
    dpd_buf4_close(&T);

  }
  else if(params.ref == 2) { /*** UHF ***/

    /* TIJAB (IA,JB) */
    dpd_buf4_init(&T, CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    dpd_buf4_sort(&T, CC_TAMPS, prqs, 20, 20, "tIAJB");
    dpd_buf4_close(&T);

    /* Tijab (ia,jb) */
    dpd_buf4_init(&T, CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    dpd_buf4_sort(&T, CC_TAMPS, prqs, 30, 30, "tiajb");
    dpd_buf4_close(&T);

    /* TIjAb (IA,jb) */
    dpd_buf4_init(&T, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_buf4_sort(&T, CC_TAMPS, prqs, 20, 30, "tIAjb");
    dpd_buf4_close(&T);

    dpd_buf4_init(&T, CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    /* TIjAb (jb,IA) */
    dpd_buf4_sort(&T, CC_TAMPS, rspq, 30, 20, "tiaJB");
    /* TIjAb (Ib,jA) (Wmbej.c) */
    dpd_buf4_sort(&T, CC_TAMPS, psrq, 24, 27, "tIbjA");
    dpd_buf4_close(&T);

    /* TiJaB (iB,Ja) (Wmbej.c) */
    dpd_buf4_init(&T, CC_TAMPS, 0, 24, 27, 24, 27, 0, "tIbjA");
    dpd_buf4_sort(&T, CC_TAMPS, rspq, 27, 24, "tiBJa");
    dpd_buf4_close(&T);

    /* T(iJ,aB) */
    dpd_buf4_init(&T, CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_buf4_sort(&T, CC_TAMPS, qpsr, 23, 29, "tiJaB");
    dpd_buf4_close(&T);

  }
}

}} /* End namespaces */
