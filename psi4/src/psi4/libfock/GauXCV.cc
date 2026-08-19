/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

#include "GauXCV.h"

#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/psi4-dec.h"

namespace psi {

GauXCRV::GauXCRV(std::shared_ptr<SuperFunctional> functional, std::shared_ptr<BasisSet> primary, Options& options)
    : RV(functional, primary, options) {}

GauXCRV::~GauXCRV() {}

void GauXCRV::compute_V(std::vector<SharedMatrix> ret) {
    throw PSIEXCEPTION("DFT_V_ALGORITHM=GAUXC is not implemented yet.");
}

void GauXCRV::print_header() const { RV::print_header(); }

GauXCUV::GauXCUV(std::shared_ptr<SuperFunctional> functional, std::shared_ptr<BasisSet> primary, Options& options)
    : UV(functional, primary, options) {}

GauXCUV::~GauXCUV() {}

void GauXCUV::compute_V(std::vector<SharedMatrix> ret) {
    throw PSIEXCEPTION("DFT_V_ALGORITHM=GAUXC is not implemented yet.");
}

void GauXCUV::print_header() const { UV::print_header(); }

}  // namespace psi
