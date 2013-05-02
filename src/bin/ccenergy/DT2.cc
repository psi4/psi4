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
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

void DT2(void)
{
  dpdbuf4 D;

  if(params.ref == 0) { /*** RHF ***/
    dpd_buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_buf4_copy(&D, PSIF_CC_TAMPS, "New tIjAb");
    dpd_buf4_close(&D);
  }
  else if(params.ref == 1) { /*** ROHF ***/

    dpd_buf4_init(&D, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <ij||ab> (i>j,a>b)");
    dpd_buf4_copy(&D, PSIF_CC_TAMPS, "New tIJAB");
    dpd_buf4_copy(&D, PSIF_CC_TAMPS, "New tijab");
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_buf4_copy(&D, PSIF_CC_TAMPS, "New tIjAb");
    dpd_buf4_close(&D);
  }
  else if(params.ref == 2) { /*** UHF ***/

    dpd_buf4_init(&D, PSIF_CC_DINTS, 0, 2, 7, 2, 7, 0, "D <IJ||AB> (I>J,A>B)");
    dpd_buf4_copy(&D, PSIF_CC_TAMPS, "New tIJAB");
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, PSIF_CC_DINTS, 0, 12, 17, 12, 17, 0, "D <ij||ab> (i>j,a>b)");
    dpd_buf4_copy(&D, PSIF_CC_TAMPS, "New tijab");
    dpd_buf4_close(&D);

    dpd_buf4_init(&D, PSIF_CC_DINTS, 0, 22, 28, 22, 28, 0, "D <Ij|Ab>");
    dpd_buf4_copy(&D, PSIF_CC_TAMPS, "New tIjAb");
    dpd_buf4_close(&D);

  }
}
}} // namespace psi::ccenergy
