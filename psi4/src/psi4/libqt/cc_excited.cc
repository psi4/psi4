/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This file is part of Psi4.
 *
 * Psi4 is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * Psi4 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License along
 * with Psi4; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*!
  \file
  \brief Determine if a wavefunction is a CC-excited wavefunction type
  \ingroup QT
*/

#include "psi4/psifiles.h"
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/exception.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <string>

namespace psi {

/*!
** cc_excited(): This function takes a WFN string and returns 1 if the WFN
** is an excited-state method and 0 if the WFN is a ground-state method.
**
** \param *wfn = wavefunction string
**
** Returns: 1 if an excited state method, else 0
** \ingroup QT
*/
int cc_excited(const char *wfn) {
    if (!strcmp(wfn, "CCSD") || !strcmp(wfn, "CCSD_T") || !strcmp(wfn, "BCCD") || !strcmp(wfn, "BCCD_T") ||
        !strcmp(wfn, "CC2") || !strcmp(wfn, "CC3") || !strcmp(wfn, "CCSD_MVD") || !strcmp(wfn, "CCSD_AT")) {
        return 0;
    } else if (!strcmp(wfn, "EOM_CCSD") || !strcmp(wfn, "EOM_CC2") || !strcmp(wfn, "EOM_CC3")) {
        return 1;
    } else {
        std::string str = "Invalid value of input keyword WFN: ";
        str += wfn;
        throw PsiException(str, __FILE__, __LINE__);
    }
}

/*!
** cc_excited(): This function takes a WFN string and returns 1 if the WFN
** is an excited-state method and 0 if the WFN is a ground-state method.
**
** \param wfn = wavefunction string
**
** Returns: 1 if an excited state method, else 0
** \ingroup QT
*/
int cc_excited(std::string wfn) { return cc_excited(wfn.c_str()); }

}  // namespace psi
