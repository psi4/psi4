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
    \ingroup RESPONSE
    \brief Enter brief description of file here 
*/
#include <libdpd/dpd.h>
#include <psifiles.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace response {

/* build_B_RHF(): Builds the RHF B matrix for RPA calculations.
** In spin orbitals, the B matrix is:
**
** B(ai,bj) = <ij||ab>
**
** RHF references and singlet eigenstates:
**  B(AI,BJ) = 2 <IJ|AB> - <IJ|BA>
**
** RHF references and triplet eigenstates:
**  B(AI,BJ) = <IJ|AB>
**
** TDC, March 2003
*/

void build_B_RHF(void)
{
  dpdbuf4 D;

  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
  global_dpd_->buf4_sort(&D, PSIF_MO_HESS, rpsq, 11, 11, "B(AI,BJ)");
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>"); 
  global_dpd_->buf4_sort(&D, PSIF_MO_HESS, rpsq, 11, 11, "B(AI,BJ) triplet");
  global_dpd_->buf4_close(&D);
}

}} // namespace psi::response
