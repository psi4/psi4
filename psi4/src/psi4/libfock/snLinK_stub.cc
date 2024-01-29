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

#error We shouldn't be here!

#include "jk.h"
#include "SplitJK.h"
#include "psi4/libpsi4util/PsiOutStream.h"

using namespace psi;

// "stub" implementation for snLinK when GauXC isn't installed. 
// Just throws on calling the constructor.
namespace psi {

snLinK::snLinK(std::shared_ptr<BasisSet> primary, Options& options) : SplitJK(primary, options) {
    throw PSIEXCEPTION("snLinK was requested in SCF_TYPE, but GauXC is not installed! Please recompile Psi4 with ENABLE_GAUXC=True.");
}

snLinK::~snLinK() {}

void snLinK::print_header() const {
    return;
}

void snLinK::build_G_component(std::vector<std::shared_ptr<Matrix>>& D, std::vector<std::shared_ptr<Matrix>>& K,
    std::vector<std::shared_ptr<TwoBodyAOInt> >& eri_computers) {
    return;
}

}  // namespace psi
