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
#include <cmath>
#include <libdpd/dpd.h>
#define EXTERN
#include "globals.h"

namespace psi{ namespace mp2{

void rhf_twopdm(void);
void uhf_twopdm(void);
void rhf_sf_twopdm(void);
void uhf_sf_twopdm(void);

void twopdm(void)
{
  if(params.ref == 0) rhf_sf_twopdm();
  else if(params.ref == 2) uhf_sf_twopdm();
}

void rhf_twopdm(void)
{
  dpdbuf4 T;

  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  global_dpd_->buf4_copy(&T, PSIF_CC_GAMMA, "GIjAb");
  global_dpd_->buf4_close(&T);
}

void uhf_twopdm(void)
{

}

void rhf_sf_twopdm(void)
{
  dpdbuf4 T;

  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tIJAB");
  global_dpd_->buf4_copy(&T, PSIF_CC_GAMMA, "GIJAB");
  global_dpd_->buf4_close(&T);

  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tijab");
  global_dpd_->buf4_copy(&T, PSIF_CC_GAMMA, "Gijab");
  global_dpd_->buf4_close(&T);

  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
  global_dpd_->buf4_copy(&T, PSIF_CC_GAMMA, "GIjAb");
  global_dpd_->buf4_close(&T);
}

void uhf_sf_twopdm(void)
{

}

}} /* End namespaces */
