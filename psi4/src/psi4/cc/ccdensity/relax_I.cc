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
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccdensity {

/* relax_I(): Add the orbital-response contributions from the
** one-electron density matrix to the I(I,J) and I(I,A) blocks of the
** Lagrangian.  These terms arise from the first-order CPHF
** equations.  I *think* the following code is general enough to deal
** with both RHF and ROHF cases.
** */

void relax_I_RHF();
void relax_I_ROHF();
void relax_I_UHF();

void relax_I() {
    if (params.ref == 0)
        relax_I_RHF();
    else if (params.ref == 1)
        relax_I_ROHF();
    else if (params.ref == 2)
        relax_I_UHF();
}

}  // namespace ccdensity
}  // namespace psi
