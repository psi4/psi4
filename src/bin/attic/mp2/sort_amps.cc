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

    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_scmcopy(&T, PSIF_CC_TAMPS, "2 tIjAb - tIjBa", 2);
    global_dpd_->buf4_sort_axpy(&T, PSIF_CC_TAMPS, pqsr, 0, 5, "2 tIjAb - tIjBa", -1);
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 7, 0, 5, 1, "tIjAb");
    global_dpd_->buf4_copy(&T, PSIF_CC_TAMPS, "tIJAB");
    global_dpd_->buf4_copy(&T, PSIF_CC_TAMPS, "tijab");
    global_dpd_->buf4_close(&T);

    /* T(iJ,aB) */
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_sort(&T, PSIF_CC_TAMPS, qpsr, 0, 5, "tiJaB");
    global_dpd_->buf4_close(&T);

    /* TIjAb (IA,jb) */
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_sort(&T, PSIF_CC_TAMPS, prqs, 10, 10, "tIAjb");
    global_dpd_->buf4_close(&T);

    /* TIjAb (ij,JB) */
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->buf4_sort(&T, PSIF_CC_TAMPS, rspq, 10, 10, "tiaJB");
    global_dpd_->buf4_close(&T);

    /* TIjAb (Ib,jA) */
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    global_dpd_->buf4_sort(&T, PSIF_CC_TAMPS, psrq, 10, 10, "tIbjA");
    global_dpd_->buf4_close(&T);

    /* TIjAb (jA,Ib) */
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    global_dpd_->buf4_sort(&T, PSIF_CC_TAMPS, rspq, 10, 10, "tjAIb");
    global_dpd_->buf4_close(&T);

  }
  else if(params.ref == 2) { /*** UHF ***/

    /* TIJAB (IA,JB) */
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    global_dpd_->buf4_sort(&T, PSIF_CC_TAMPS, prqs, 20, 20, "tIAJB");
    global_dpd_->buf4_close(&T);

    /* Tijab (ia,jb) */
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    global_dpd_->buf4_sort(&T, PSIF_CC_TAMPS, prqs, 30, 30, "tiajb");
    global_dpd_->buf4_close(&T);

    /* TIjAb (IA,jb) */
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->buf4_sort(&T, PSIF_CC_TAMPS, prqs, 20, 30, "tIAjb");
    global_dpd_->buf4_close(&T);

    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    /* TIjAb (jb,IA) */
    global_dpd_->buf4_sort(&T, PSIF_CC_TAMPS, rspq, 30, 20, "tiaJB");
    /* TIjAb (Ib,jA) (Wmbej.c) */
    global_dpd_->buf4_sort(&T, PSIF_CC_TAMPS, psrq, 24, 27, "tIbjA");
    global_dpd_->buf4_close(&T);

    /* TiJaB (iB,Ja) (Wmbej.c) */
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 24, 27, 24, 27, 0, "tIbjA");
    global_dpd_->buf4_sort(&T, PSIF_CC_TAMPS, rspq, 27, 24, "tiBJa");
    global_dpd_->buf4_close(&T);

    /* T(iJ,aB) */
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    global_dpd_->buf4_sort(&T, PSIF_CC_TAMPS, qpsr, 23, 29, "tiJaB");
    global_dpd_->buf4_close(&T);

  }
}

}} /* End namespaces */
