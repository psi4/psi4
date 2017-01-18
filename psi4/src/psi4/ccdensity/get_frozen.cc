/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
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
 * @END LICENSE
 */

/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace ccdensity {

/*
** get_frozen(): Routine to get symmetry and orbital reordering arrays
** which include frozen orbitals.
**
** TDC, March 2000.
*/

void get_frozen(void)
{
  int i, nfzc;

  frozen.nfzc = moinfo.nfzc;
  frozen.nfzv = moinfo.nfzv;

  nfzc = moinfo.nfzc;

  /* Get full orbital list occpi and virtpi arrays */
  frozen.occpi = init_int_array(moinfo.nirreps);
  frozen.virtpi = init_int_array(moinfo.nirreps);
  psio_read_entry(PSIF_CC_INFO, "All Occ Orbs Per Irrep",
		  (char *) frozen.occpi, sizeof(int)*moinfo.nirreps);
  psio_read_entry(PSIF_CC_INFO, "All Virt Orbs Per Irrep",
		  (char *) frozen.virtpi, sizeof(int)*moinfo.nirreps);

  /* Get CC->QT and QT->CC full occupied and virtual reordering arrays */
  frozen.qt_occ = init_int_array(moinfo.nmo);
  frozen.qt_vir = init_int_array(moinfo.nmo);
  frozen.allcc_occ = init_int_array(moinfo.nmo);
  frozen.allcc_vir = init_int_array(moinfo.nmo);
  psio_read_entry(PSIF_CC_INFO, "CC->QT All Occ Order",
		   (char *) frozen.qt_occ, sizeof(int)*moinfo.nmo);
  psio_read_entry(PSIF_CC_INFO, "CC->QT All Virt Order",
		   (char *) frozen.qt_vir, sizeof(int)*moinfo.nmo);
  psio_read_entry(PSIF_CC_INFO, "QT->CC All Occ Order",
		   (char *) frozen.allcc_occ, sizeof(int)*moinfo.nmo);
  psio_read_entry(PSIF_CC_INFO, "QT->CC All Virt Order",
		   (char *) frozen.allcc_vir, sizeof(int)*moinfo.nmo);

  /* Build active-only cc_occ and cc_vir, inserting -1s for frozen orbitals */
  /* This is used in resort_gamma */
  frozen.cc_occ = init_int_array(moinfo.nmo);
  frozen.cc_vir = init_int_array(moinfo.nmo);
  for(i=0; i < moinfo.nmo; i++) {
      frozen.cc_occ[i] = -1;
      frozen.cc_vir[i] = -1;
    }
  for(i=0; i < moinfo.nactive; i++) {
      frozen.cc_occ[i+nfzc] = moinfo.cc_occ[i];
      frozen.cc_vir[i+nfzc] = moinfo.cc_vir[i];
    }

  /* Get full orbital list symmetry arrays */
  frozen.occ_sym = init_int_array(moinfo.nmo);
  frozen.vir_sym = init_int_array(moinfo.nmo);
  psio_read_entry(PSIF_CC_INFO, "All Occ Orb Symmetry",
		  (char *) frozen.occ_sym, sizeof(int)*moinfo.nmo);
  psio_read_entry(PSIF_CC_INFO, "All Virt Orb Symmetry",
		  (char *) frozen.vir_sym, sizeof(int)*moinfo.nmo);

  /* Get full orbtial list offset arrays for occupied and virtual */
  frozen.occ_off = init_int_array(moinfo.nirreps);
  frozen.vir_off = init_int_array(moinfo.nirreps);
  psio_read_entry(PSIF_CC_INFO, "All Occ Orb Offsets",
		   (char *) frozen.occ_off, sizeof(int)*moinfo.nirreps);
  psio_read_entry(PSIF_CC_INFO, "All Virt Orb Offsets",
		   (char *) frozen.vir_off, sizeof(int)*moinfo.nirreps);

  /* Get boolean arrays for orbital classification routines */
  frozen.occ = init_int_array(moinfo.nmo);
  frozen.vir = init_int_array(moinfo.nmo);
  frozen.socc = init_int_array(moinfo.nmo);
  psio_read_entry(PSIF_CC_INFO, "All Occ Orbital Boolean",
		  (char *) frozen.occ, sizeof(int)*moinfo.nmo);
  psio_read_entry(PSIF_CC_INFO, "All Virt Orbital Boolean",
		  (char *) frozen.vir, sizeof(int)*moinfo.nmo);
  psio_read_entry(PSIF_CC_INFO, "All Socc Orbital Boolean",
		  (char *) frozen.socc, sizeof(int)*moinfo.nmo);
}

}} // namespace psi::ccdensity
