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
** REF_ORBS.C
** 
** This file contains routines pertaining to the set of reference orbitals,
** C_0, from which the orbital rotation angles are defined.
**
** C. David Sherrill
** University of California, Berkeley
** May 1998
*/

#include <cstdlib>
#include <cstdio>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libpsio/psio.hpp>
#include <libpsio/psio.h>
#include <psifiles.h>
#include "MCSCF_globaldefs.h"
#include "MCSCF_globals.h"
#include "psi4-dec.h"


namespace psi { namespace detcas {

/*
** read_ref_orbs()
**
** This function reads in the ``reference'' orbitals, C_0, from a special
**  file.  If this file does not exist, then it is created, and the
**  current orbitals in file30 are placed there, and the array of orbital
**  rotation angles is reset to zero.
**
** Returns: 1 if read is successful, otherwise 0
*/
int read_ref_orbs(void)
{
  int h, ir_orbs;
  char orb_key[80];

  if (psio_tocentry_exists(PSIF_DETCAS, "Orbs Irrep  0")){ 
    psio_open(PSIF_DETCAS, PSIO_OPEN_OLD);
    for (h=0; h<CalcInfo.nirreps; h++) {
      ir_orbs = CalcInfo.orbs_per_irr[h];
      sprintf(orb_key, "Orbs Irrep %2d", h);
      psio_read_entry(PSIF_DETCAS, orb_key, (char *) CalcInfo.mo_coeffs[h][0],
                      ir_orbs*ir_orbs*sizeof(double));
    }
    psio_close(PSIF_DETCAS, 1);
    return(1);
  }
  else {
    if (Params.print_lvl) 
      outfile->Printf("No previous orbitals ... using new reference orbitals\n");
    return(0);
  }
}


/*
** write_ref_orbs()
**
** This function initializes the set of reference orbitals.  The current
**  orbitals in file 30 are read and written to the reference orbital
**  disk file.
**
** Returns: 1 if successful, otherwise 0
*/
int write_ref_orbs(void)
{
  int h, ir_orbs;
  char orb_key[80];

  psio_open(PSIF_DETCAS, PSIO_OPEN_OLD);
  for (h=0; h<CalcInfo.nirreps; h++) {
    ir_orbs = CalcInfo.orbs_per_irr[h];
    sprintf(orb_key, "Orbs Irrep %2d", h);
    psio_write_entry(PSIF_DETCAS, orb_key, (char *) CalcInfo.mo_coeffs[h][0],
                    ir_orbs*ir_orbs*sizeof(double));
  }
  psio_close(PSIF_DETCAS, 1);
  return(1);
}

}} // end namespace psi::detcas

