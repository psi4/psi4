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
    \ingroup DETCAS
    \brief Enter brief description of file here 
*/
#include <cstdlib>
#include <cstdio>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include "globaldefs.h"
#include "structs.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace detci {

/*
** cleanup()
**
** This function frees any allocated global variables
**
*/
void mcscf_cleanup(void)
{
  int i;
  
  free(MCSCF_CalcInfo.docc);
  free(MCSCF_CalcInfo.socc);
  free(MCSCF_CalcInfo.frozen_docc);
  free(MCSCF_CalcInfo.frozen_uocc);
  free(MCSCF_CalcInfo.rstr_docc);
  free(MCSCF_CalcInfo.rstr_uocc);
  // free(MCSCF_CalcInfo.orbsym);
  // free(CalcInfo.reorder);
  // free(CalcInfo.order);
  free(MCSCF_CalcInfo.ci2relpitz);
  // free(MCSCF_CalcInfo.first);
  // free(MCSCF_CalcInfo.last);
  // free(MCSCF_CalcInfo.fstact);
  // free(MCSCF_CalcInfo.lstact);
  // free(MCSCF_CalcInfo.active);
  // free_int_matrix(CalcInfo.ras_opi);
  free_int_matrix(MCSCF_CalcInfo.fzc_orbs);
  free_int_matrix(MCSCF_CalcInfo.fzv_orbs);
  for (i=0; i<MAX_RAS_SPACES; i++) 
    free_int_matrix(MCSCF_CalcInfo.ras_orbs[i]);
  free(MCSCF_CalcInfo.ras_orbs);

  // for (i=0; i<CalcInfo.nirreps; i++) 
  //   free(CalcInfo.labels[i]);

  for (i=0; i<CalcInfo.nirreps; i++) {
    if (CalcInfo.orbs_per_irr[i]) 
      free_block(MCSCF_CalcInfo.mo_coeffs[i]);
  }
  free(MCSCF_CalcInfo.mo_coeffs);

  free(MCSCF_CalcInfo.onel_ints);
  free(MCSCF_CalcInfo.onel_ints_bare);
  free(MCSCF_CalcInfo.twoel_ints);
  free_block(MCSCF_CalcInfo.opdm);
  free(MCSCF_CalcInfo.tpdm);
  free_block(MCSCF_CalcInfo.lag);
  free(MCSCF_CalcInfo.F_act);
  free(MCSCF_CalcInfo.mo_grad);
  if (MCSCF_CalcInfo.mo_hess_diag != NULL) free(MCSCF_CalcInfo.mo_hess_diag);
  if (MCSCF_CalcInfo.mo_hess != NULL) free_block(MCSCF_CalcInfo.mo_hess);
  free(MCSCF_CalcInfo.theta_cur);
  free(MCSCF_CalcInfo.theta_step);
  // free(CalcInfo.orbs_per_irr);
}

}} // end namespace psi::detci

