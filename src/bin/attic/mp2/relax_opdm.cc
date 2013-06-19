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
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

namespace psi{ namespace mp2{

void relax_opdm(void)
{
  dpdfile2 D1, D2;

  if(params.ref == 0) {   
    dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 0, 1, "DAI");
    dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 1, 0, "D(orb)(A,I)");
    dpd_->file2_axpy(&D2, &D1, 1.0, 1);
    dpd_->file2_close(&D2);
    dpd_->file2_close(&D1);

    dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 0, 1, "DIA");
    dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 1, 0, "D(orb)(A,I)");
    dpd_->file2_axpy(&D2, &D1, 1.0, 1);
    dpd_->file2_close(&D2);
    dpd_->file2_close(&D1);

    dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 0, 1, "Dai");
    dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 1, 0, "D(orb)(a,i)");
    dpd_->file2_axpy(&D2, &D1, 1.0, 1);
    dpd_->file2_close(&D2);
    dpd_->file2_close(&D1);

    dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 0, 1, "Dia");
    dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 1, 0, "D(orb)(a,i)");
    dpd_->file2_axpy(&D2, &D1, 1.0, 1);
    dpd_->file2_close(&D2);
    dpd_->file2_close(&D1);
  }
  else if(params.ref == 2) {

    dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 0, 1, "DAI");
    dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 1, 0, "D(orb)(A,I)");
    dpd_->file2_axpy(&D2, &D1, 1.0, 1);
    dpd_->file2_close(&D2);
    dpd_->file2_close(&D1);

    dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 0, 1, "DIA");
    dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 1, 0, "D(orb)(A,I)");
    dpd_->file2_axpy(&D2, &D1, 1.0, 1);
    dpd_->file2_close(&D2);
    dpd_->file2_close(&D1);

    dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 2, 3, "Dai");
    dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 3, 2, "D(orb)(a,i)");
    dpd_->file2_axpy(&D2, &D1, 1.0, 1);
    dpd_->file2_close(&D2);
    dpd_->file2_close(&D1);

    dpd_->file2_init(&D1, PSIF_CC_OEI, 0, 2, 3, "Dia");
    dpd_->file2_init(&D2, PSIF_CC_OEI, 0, 3, 2, "D(orb)(a,i)");
    dpd_->file2_axpy(&D2, &D1, 1.0, 1);
    dpd_->file2_close(&D2);
    dpd_->file2_close(&D1);
  }
}

}} /* End namespaces */
