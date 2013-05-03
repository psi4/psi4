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
    \ingroup ccresponse
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccresponse {

void sort_lamps(void)
{
  dpdbuf4 L;

  /* RAK fixing this for new cclambda, assuming A1 ground lambda? */
  dpd_buf4_init(&L, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "LIjAb 0 -1");
  dpd_buf4_scmcopy(&L, PSIF_CC_LAMPS, "2 LIjAb - LIjBa", 2);
  dpd_buf4_sort_axpy(&L, PSIF_CC_LAMPS, pqsr, 0, 5, "2 LIjAb - LIjBa", -1);
  dpd_buf4_close(&L);
}  

}} // namespace psi::ccresponse
