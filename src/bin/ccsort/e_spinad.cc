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

void e_spinad(void)
{
  dpdbuf4 E, E1;

  if(params.ref == 0) { /*** RHF ***/
    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E <ai|jk>");
    global_dpd_->buf4_scmcopy(&E, PSIF_CC_EINTS, "E 2<ai|jk> - <ai|kj>", 2);
    global_dpd_->buf4_sort(&E, PSIF_CC_TMP0, pqsr, 11, 0, "E <ai|kj>");
    global_dpd_->buf4_close(&E);

    global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 11, 0, 11, 0, 0, "E 2<ai|jk> - <ai|kj>");
    global_dpd_->buf4_init(&E1, PSIF_CC_TMP0, 0, 11, 0, 11, 0, 0, "E <ai|kj>");
    global_dpd_->buf4_axpy(&E1, &E, -1);
    global_dpd_->buf4_close(&E1);
    global_dpd_->buf4_close(&E);

  }
}

}} // namespace psi::ccsort
