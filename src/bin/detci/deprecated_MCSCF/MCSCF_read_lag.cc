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
/*
** READ_LAG.C
**
** Read the lagrangian
**
** C. David Sherrill
** University of California, Berkeley
**
** April 1998
** Updated to new libpsio libraries 8/03
*/

#include <cstdlib>
#include <cstdio>
#include <libiwl/iwl.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include "globaldefs.h"
#include "structs.h"
#define EXTERN
#include "globals.h"
#include "psi4-dec.h"


namespace psi { namespace detci {

void read_lagrangian(void)
{
  int nmo;

  nmo = MCSCF_CalcInfo.nmo;
  
  MCSCF_CalcInfo.lag = block_matrix(nmo, nmo);

  psio_open(MCSCF_Parameters.lag_file, PSIO_OPEN_OLD);  
  psio_read_entry(MCSCF_Parameters.lag_file, "MO-basis Lagrangian", 
    (char *) MCSCF_CalcInfo.lag[0], nmo*nmo*sizeof(double));

  if (MCSCF_Parameters.print_lvl > 3) {
    outfile->Printf("Lagrangian matrix\n");
    print_mat(MCSCF_CalcInfo.lag, nmo, nmo, "outfile");
  }

  psio_close(MCSCF_Parameters.lag_file, MCSCF_Parameters.lag_erase ? 0 : 1);

} 

}} // end namespace psi::detci

