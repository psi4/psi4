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

/*!
  \file
  \brief Check if wavefunction is coupled-cluster type
  \ingroup QT
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include "psi4/psifiles.h"

namespace psi {

/*!
** cc_wfn(): Checks if the given wavefunction string is a coupled-cluster
** type and returns 1 if yes and 0 if no.
**
** Note: "coupled-cluster type" means it is handled by PSI like the
** coupled-cluster codes, not necessarily that it is literally a
** coupled-cluster wavefunction
**
** \param *wfn = wavefunction string
**
** Returns: 1 if the WFN is a CC method, 0 otherwise
**
** \ingroup QT
*/
int cc_wfn(const char *wfn)
{
  if ( !strcmp(wfn, "CCSD")     || !strcmp(wfn, "CCSD_T") ||
       !strcmp(wfn, "BCCD")     || !strcmp(wfn, "BCCD_T") ||
       !strcmp(wfn, "CC2")      || !strcmp(wfn, "CC3")    ||
       !strcmp(wfn, "EOM_CCSD") || !strcmp(wfn, "LEOM_CCSD") ||
       !strcmp(wfn, "EOM_CC2")  || !strcmp(wfn, "EOM_CC3") ||
       !strcmp(wfn, "CIS")      || !strcmp(wfn, "CCSD_AT")) {
    return 1;
  }
  else {
    return 0;
  }
}

/*!
** cc_wfn(): Checks if the given wavefunction string is a coupled-cluster
** type and returns 1 if yes and 0 if no.
**
** Note: "coupled-cluster type" means it is handled by PSI like the
** coupled-cluster codes, not necessarily that it is literally a
** coupled-cluster wavefunction
**
** \param wfn = wavefunction string
**
** Returns: 1 if the WFN is a CC method, 0 otherwise
**
** \ingroup QT
*/
int cc_wfn(std::string wfn)
{
  return cc_wfn(wfn.c_str());
}

}
