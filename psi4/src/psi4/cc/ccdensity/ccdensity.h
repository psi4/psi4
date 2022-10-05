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

#ifndef CCDENSITY_CCDENSITY_H
#define CCDENSITY_CCDENSITY_H

/*! \file
    \ingroup CCDENSITY
    \brief What will eventually be an honest-to-goodness header file for this module.
           Refactoring needed, so pardon the mess, and help out if you can.
*/

#include "psi4/cc/ccwave.h"
#include "psi4/libmints/typedefs.h"
#include "psi4/libmints/wavefunction.h"

namespace psi {
namespace ccdensity {
// Save a scalar describing a ground->excited transition.
void scalar_saver_ground(ccenergy::CCEnergyWavefunction& wfn, struct TD_Params *S, const std::string suffix, double val);
// Save a scalar describing a excited->excited transition.
void scalar_saver_excited(ccenergy::CCEnergyWavefunction& wfn, struct TD_Params *S, struct TD_Params *U, const std::string suffix, double val);
// Save a density.
void density_saver(ccenergy::CCEnergyWavefunction& wfn, struct RHO_Params *S, const std::string suffix, SharedMatrix val);

}
}

#endif
