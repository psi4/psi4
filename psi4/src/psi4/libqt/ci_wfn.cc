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
  \brief Check if wavefunction is CI-type
  \ingroup QT
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>
#include "psi4/psifiles.h"

namespace psi {

/*!
** ci_wfn(): Examine the wavefunction type and return 1 if a CI/MCSCF-type,
** otherwise 0
**
** \param *wfn = wavefunction string
**
** Returns: 1 if a CI/MCSCF-type wavefunction, otherwise 0
**
** \ingroup QT
*/
int ci_wfn(char *wfn)
{

  if (strcmp(wfn, "CI")==0     || strcmp(wfn, "DETCAS")==0 ||
      strcmp(wfn, "CASSCF")==0 || strcmp(wfn, "RASSCF")==0 ||
      strcmp(wfn, "DETCI")==0 || strcmp(wfn, "MCSCF")==0 ||
      strcmp(wfn, "OOCCD")==0 || strcmp(wfn,"ZAPTN")==0)
  {
    return(1);
  }
  else  {
    return 0;
  }
}

/*!
** ci_wfn(): Examine the wavefunction type and return 1 if a CI/MCSCF-type,
** otherwise 0
**
** \param wfn = wavefunction string
**
** Returns: 1 if a CI/MCSCF-type wavefunction, otherwise 0
**
** \ingroup QT
*/
int ci_wfn(std::string wfn)
{

  if ((wfn == "CI")    || (wfn == "DETCAS") ||
      (wfn =="CASSCF") || (wfn =="RASSCF")  ||
      (wfn =="DETCI")  || (wfn =="MCSCF")   ||
      (wfn =="OOCCD")  || (wfn == "ZAPTN"))
  {
    return(1);
  }
  else  {
    return 0;
  }
}

}
