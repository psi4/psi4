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

#include "psi4/libmints/potentialint.h"

namespace psi{

PCMPotentialInt::PCMPotentialInt(std::vector<SphericalTransform>& trans,
std::shared_ptr<BasisSet> bs1, std::shared_ptr<BasisSet> bs2, int deriv):
    PotentialInt(trans, bs1, bs1)
{
    UNUSED(bs2);
    UNUSED(deriv);
    // We don't want to transform the integrals from Cartesian (6d, 10f, ...) to Pure (5d, 7f, ...)
    // for each external charge.  It'll be better to backtransform the density / Fock matrices to
    // the Cartesian basis, if necessary.
    force_cartesian_ = true;
}

} //Namespace
