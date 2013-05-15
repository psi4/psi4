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

void ET2(void)
{
  dpdfile2 tIA, tia;
  dpdbuf4 newtIJAB, newtijab, newtIjAb;
  dpdbuf4 E, t2, t2a, t2b;

  if(params.ref == 0) { /** RHF **/
    dpd_buf4_init(&newtIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

    dpd_file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");

    dpd_buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    dpd_contract424(&E, &tIA, &newtIjAb, 1, 0, 0, -1, 1);
    dpd_buf4_close(&E);

    dpd_buf4_init(&E, PSIF_CC_EINTS, 0, 10, 0, 10, 0, 0, "E <ia|jk>");
    dpd_contract244(&tIA, &E, &newtIjAb, 0, 0, 1, -1, 1);
    dpd_buf4_close(&E);

    dpd_file2_close(&tIA);

    dpd_buf4_close(&newtIjAb);
  }
  else if(params.ref == 1) { /** ROHF **/

    dpd_buf4_init(&newtIJAB, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tIJAB");
    dpd_buf4_init(&newtijab, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tijab");
    dpd_buf4_init(&newtIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

    dpd_file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

    /*** AA ***/
    dpd_buf4_init(&E, PSIF_CC_EINTS, 0, 11, 2, 11, 0, 1, "E <ai|jk>");
    dpd_buf4_init(&t2, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    dpd_contract424(&E, &tIA, &t2, 1, 0, 0, -1, 0);
    dpd_buf4_sort(&t2, PSIF_CC_TMP0, pqsr, 2, 5, "T (I>J,BA)");
    dpd_buf4_close(&t2);
    dpd_buf4_init(&t2a, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    dpd_buf4_init(&t2b, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,BA)");
    dpd_buf4_axpy(&t2b, &t2a, -1);
    dpd_buf4_close(&t2b);
    dpd_buf4_axpy(&t2a, &newtIJAB, 1);
    dpd_buf4_close(&t2a);
    dpd_buf4_close(&E);

    /*** BB ***/
    dpd_buf4_init(&E, PSIF_CC_EINTS, 0, 11, 2, 11, 0, 1, "E <ai|jk>");
    dpd_buf4_init(&t2, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    dpd_contract424(&E, &tia, &t2, 1, 0, 0, -1, 0);
    dpd_buf4_sort(&t2, PSIF_CC_TMP0, pqsr, 2, 5, "T (I>J,BA)");
    dpd_buf4_close(&t2);
    dpd_buf4_init(&t2a, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    dpd_buf4_init(&t2b, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,BA)");
    dpd_buf4_axpy(&t2b, &t2a, -1);
    dpd_buf4_close(&t2b);
    dpd_buf4_axpy(&t2a, &newtijab, 1);
    dpd_buf4_close(&t2a);
    dpd_buf4_close(&E);

    /*** AB ***/

    dpd_buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    dpd_contract424(&E, &tia, &newtIjAb, 1, 0, 0, -1, 1);
    dpd_buf4_close(&E);
    dpd_buf4_init(&E, PSIF_CC_EINTS, 0, 10, 0, 10, 0, 0, "E <ia|jk>");
    dpd_contract244(&tIA, &E, &newtIjAb, 0, 0, 1, -1, 1);
    dpd_buf4_close(&E);

    dpd_file2_close(&tIA); 
    dpd_file2_close(&tia);

    dpd_buf4_close(&newtIJAB);
    dpd_buf4_close(&newtijab);
    dpd_buf4_close(&newtIjAb);

  }
  else if(params.ref == 2) { /*** UHF ***/

    dpd_buf4_init(&newtIJAB, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tIJAB");
    dpd_buf4_init(&newtijab, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "New tijab");
    dpd_buf4_init(&newtIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");

    dpd_file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

    /*** AA ***/
    dpd_buf4_init(&E, PSIF_CC_EINTS, 0, 21, 2, 21, 0, 1, "E <AI|JK>");
    dpd_buf4_init(&t2, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    dpd_contract424(&E, &tIA, &t2, 1, 0, 0, -1, 0);
    dpd_buf4_sort(&t2, PSIF_CC_TMP0, pqsr, 2, 5, "T (I>J,BA)");
    dpd_buf4_close(&t2);
    dpd_buf4_init(&t2a, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,AB)");
    dpd_buf4_init(&t2b, PSIF_CC_TMP0, 0, 2, 5, 2, 5, 0, "T (I>J,BA)");
    dpd_buf4_axpy(&t2b, &t2a, -1);
    dpd_buf4_close(&t2b);
    dpd_buf4_axpy(&t2a, &newtIJAB, 1);
    dpd_buf4_close(&t2a);
    dpd_buf4_close(&E);

    /*** BB ***/
    dpd_buf4_init(&E, PSIF_CC_EINTS, 0, 31, 12, 31, 10, 1, "E <ai|jk>");
    dpd_buf4_init(&t2, PSIF_CC_TMP0, 0, 12, 15, 12, 15, 0, "T (i>j,ab)");
    dpd_contract424(&E, &tia, &t2, 1, 0, 0, -1, 0);
    dpd_buf4_sort(&t2, PSIF_CC_TMP0, pqsr, 12, 15, "T (i>j,ba)");
    dpd_buf4_close(&t2);
    dpd_buf4_init(&t2a, PSIF_CC_TMP0, 0, 12, 15, 12, 15, 0, "T (i>j,ab)");
    dpd_buf4_init(&t2b, PSIF_CC_TMP0, 0, 12, 15, 12, 15, 0, "T (i>j,ba)");
    dpd_buf4_axpy(&t2b, &t2a, -1);
    dpd_buf4_close(&t2b);
    dpd_buf4_axpy(&t2a, &newtijab, 1);
    dpd_buf4_close(&t2a);
    dpd_buf4_close(&E);

    /*** AB ***/

    dpd_buf4_init(&E, PSIF_CC_EINTS, 0, 22, 26, 22, 26, 0, "E <Ij|Ak>");
    dpd_contract424(&E, &tia, &newtIjAb, 3, 0, 0, -1, 1);
    dpd_buf4_close(&E);
    dpd_buf4_init(&E, PSIF_CC_EINTS, 0, 24, 22, 24, 22, 0, "E <Ia|Jk>");
    dpd_contract244(&tIA, &E, &newtIjAb, 0, 0, 1, -1, 1);
    dpd_buf4_close(&E);

    dpd_file2_close(&tIA); 
    dpd_file2_close(&tia);

    dpd_buf4_close(&newtIJAB);
    dpd_buf4_close(&newtijab);
    dpd_buf4_close(&newtIjAb);

  }
}
}} // namespace psi::ccenergy
