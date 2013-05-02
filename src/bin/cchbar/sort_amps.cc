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
    \ingroup CCHBAR
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cchbar {

void sort_amps(void)
{
  dpdbuf4 t2;

  if(params.ref == 0 || params.ref == 1) { /** RHF or ROHF **/
    /* T(iJ,aB) */
    dpd_buf4_init(&t2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_sort(&t2, PSIF_CC_TAMPS, qpsr, 0, 5, "tiJaB");
    dpd_buf4_close(&t2);

    /* TIJAB (IA,JB) */
    dpd_buf4_init(&t2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    dpd_buf4_sort(&t2, PSIF_CC_TAMPS, prqs, 10, 10, "tIAJB");
    dpd_buf4_close(&t2);

    /* Tijab (ia,jb) */
    dpd_buf4_init(&t2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tijab");
    dpd_buf4_sort(&t2, PSIF_CC_TAMPS, prqs, 10, 10, "tiajb");
    dpd_buf4_close(&t2);

    /* TIjAb (IA,jb) */
    dpd_buf4_init(&t2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    dpd_buf4_sort(&t2, PSIF_CC_TAMPS, prqs, 10, 10, "tIAjb");
    dpd_buf4_close(&t2);

    /* TIjAb (ij,JB) */
    dpd_buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_buf4_sort(&t2, PSIF_CC_TAMPS, rspq, 10, 10, "tiaJB");
    dpd_buf4_close(&t2);

    /* TIjAb (Ib,jA) */
    dpd_buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIAjb");
    dpd_buf4_sort(&t2, PSIF_CC_TAMPS, psrq, 10, 10, "tIbjA");
    dpd_buf4_close(&t2);

    /* TIjAb (jA,Ib) */
    dpd_buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 10, 10, 10, 0, "tIbjA");
    dpd_buf4_sort(&t2, PSIF_CC_TAMPS, rspq, 10, 10, "tjAIb");
    dpd_buf4_sort(&t2, PSIF_CC_TAMPS, pqsr, 10, 11, "tIbAj");
    dpd_buf4_close(&t2);

  }
  else if(params.ref == 2) { /** UHF **/

    /* T(iJ,aB) */
    dpd_buf4_init(&t2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_buf4_sort(&t2, PSIF_CC_TAMPS, qpsr, 23, 29, "tiJaB");
    dpd_buf4_close(&t2);

    /* TIJAB (IA,JB) */
    dpd_buf4_init(&t2, PSIF_CC_TAMPS, 0, 0, 5, 2, 7, 0, "tIJAB");
    dpd_buf4_sort(&t2, PSIF_CC_TAMPS, prqs, 20, 20, "tIAJB");
    dpd_buf4_close(&t2);

    /* Tijab (ia,jb) */
    dpd_buf4_init(&t2, PSIF_CC_TAMPS, 0, 10, 15, 12, 17, 0, "tijab");
    dpd_buf4_sort(&t2, PSIF_CC_TAMPS, prqs, 30, 30, "tiajb");
    dpd_buf4_close(&t2);

    /* TIjAb (IA,jb) */
    dpd_buf4_init(&t2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tIjAb");
    dpd_buf4_sort(&t2, PSIF_CC_TAMPS, prqs, 20, 30, "tIAjb");
    dpd_buf4_close(&t2);

    /* TIjAb (ij,JB) */
    dpd_buf4_init(&t2, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    dpd_buf4_sort(&t2, PSIF_CC_TAMPS, rspq, 30, 20, "tiaJB");
    dpd_buf4_close(&t2);

    /* TIjAb (Ib,jA) */
    dpd_buf4_init(&t2, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
    dpd_buf4_sort(&t2, PSIF_CC_TAMPS, psrq, 24, 27, "tIbjA");
    dpd_buf4_close(&t2);

    dpd_buf4_init(&t2, PSIF_CC_TAMPS, 0, 24, 27, 24, 27, 0, "tIbjA");
    /* TIjAb (jA,Ib) */
    dpd_buf4_sort(&t2, PSIF_CC_TAMPS, rspq, 27, 24, "tjAIb");
    /* TIjAb (Ib,Aj) */
    dpd_buf4_sort(&t2, PSIF_CC_TAMPS, pqsr, 24, 26, "tIbAj");
    dpd_buf4_close(&t2);

    /* TiJaB (iB,aJ) */
    dpd_buf4_init(&t2, PSIF_CC_TAMPS, 0, 27, 24, 27, 24, 0, "tjAIb");
    dpd_buf4_sort(&t2, PSIF_CC_TAMPS, pqsr, 27, 25, "tiBaJ");
    dpd_buf4_close(&t2);

  }
}

}} // namespace psi::cchbar
