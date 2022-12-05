/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#include "psi4/libmints/wavefunction.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccdensity {

/* SORTI(): Place all the components of the Lagrangian into a large
** matrix, I (moinfo.I), which we also symmetrize by computing Ipq =
** 1/2 (Ipq + Iqp).  This matrix is later written to disk in dump()
** for subsequent backtransformation.  Note that some of the
** components of the Lagrangian computed into the IIJ, Iij, IIA, and
** Iia matrices remain non-symmetric (e.g., IIJ neq IJI).  I re-used
** my sortone.c code here, so don't let some of the variable names
** confuse you. */

void sortI_RHF(Wavefunction&);
void sortI_ROHF(Wavefunction&);
void sortI_UHF(Wavefunction&);

void sortI(Wavefunction& wfn) {
    if (params.ref == 0)
        return sortI_RHF(wfn);
    else if (params.ref == 1)
        return sortI_ROHF(wfn);
    else if (params.ref == 2)
        return sortI_UHF(wfn);
}

}  // namespace ccdensity
}  // namespace psi
