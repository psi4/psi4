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
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccenergy {

void local_filter_T1(dpdfile2 *T1);
void dijabT2(void);

/* apply denominators to t1 and t2 */

void denom(void)
{
  dpdfile2 newtIA, dIA, tIA, newtia, dia, tia;

  if (params.ref == 0) {
    dpd_file2_init(&newtIA, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    dpd_file2_copy(&newtIA, PSIF_CC_OEI, "New tIA Increment");
    dpd_file2_close(&newtIA);

    dpd_file2_init(&newtIA, PSIF_CC_OEI, 0, 0, 1, "New tIA Increment");
    if(params.local && local.filter_singles) {
      local_filter_T1(&newtIA);
    }
    else {
      dpd_file2_init(&dIA, PSIF_CC_OEI, 0, 0, 1, "dIA");
      dpd_file2_dirprd(&dIA, &newtIA);
      dpd_file2_close(&dIA);
    }
    dpd_file2_close(&newtIA);

    /* Add the new increment to the old tIA to get the New tIA */
    dpd_file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
    dpd_file2_copy(&tIA, PSIF_CC_OEI, "New tIA");
    dpd_file2_close(&tIA);
    dpd_file2_init(&newtIA, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    dpd_file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "New tIA Increment");
    dpd_file2_axpy(&tIA, &newtIA, 1, 0);
    dpd_file2_close(&tIA);
    dpd_file2_close(&newtIA);
  }
  else if (params.ref == 1) {
    dpd_file2_init(&newtIA, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    dpd_file2_init(&dIA, PSIF_CC_OEI, 0, 0, 1, "dIA");
    dpd_file2_dirprd(&dIA, &newtIA);
    dpd_file2_close(&dIA);
    dpd_file2_close(&newtIA); 

    dpd_file2_init(&newtia, PSIF_CC_OEI, 0, 0, 1, "New tia");
    dpd_file2_init(&dia, PSIF_CC_OEI, 0, 0, 1, "dia");
    dpd_file2_dirprd(&dia, &newtia);
    dpd_file2_close(&dia);
    dpd_file2_close(&newtia);
  }
  else if (params.ref == 2) {
    dpd_file2_init(&newtIA, PSIF_CC_OEI, 0, 0, 1, "New tIA");
    dpd_file2_init(&dIA, PSIF_CC_OEI, 0, 0, 1, "dIA");
    dpd_file2_dirprd(&dIA, &newtIA);
    dpd_file2_close(&dIA);
    dpd_file2_close(&newtIA);

    dpd_file2_init(&newtia, PSIF_CC_OEI, 0, 2, 3, "New tia");
    dpd_file2_init(&dia, PSIF_CC_OEI, 0, 2, 3, "dia");
    dpd_file2_dirprd(&dia, &newtia);
    dpd_file2_close(&dia);
    dpd_file2_close(&newtia);
  }

  dijabT2(); 

  return;
}
}} // namespace psi::ccenergy
