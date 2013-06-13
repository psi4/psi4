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
    \ingroup CCSORT
    \brief Enter brief description of file here 
*/
#include <cstdio>
#include <cstdlib>
#include <libciomr/libciomr.h>
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccsort {

void d_spinad(void)
{
  dpdbuf4 D, D1;

  if(params.ref == 0) { /*** RHF ***/
    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    dpd_->buf4_scmcopy(&D, PSIF_CC_DINTS, "D 2<ij|ab> - <ij|ba>", 2);
    dpd_->buf4_sort_axpy(&D, PSIF_CC_DINTS, pqsr, 0, 5, "D 2<ij|ab> - <ij|ba>", -1);
    dpd_->buf4_close(&D);

    dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    dpd_->buf4_sort(&D, PSIF_CC_DINTS, prqs, 10, 10, "D 2<ij|ab> - <ij|ba> (ia,jb)");
    dpd_->buf4_close(&D);

  }
}

}} // namespace psi::ccsort
