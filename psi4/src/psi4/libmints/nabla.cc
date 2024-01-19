/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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
#include "psi4/libmints/nabla.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/molecule.h"

#include <stdexcept>
#include "psi4/libciomr/libciomr.h"

#include <libint2/engine.h>

using namespace psi;
;

// to compute the Nabla derivatives
NablaInt::NablaInt(std::vector<SphericalTransform>& spherical_transforms, std::shared_ptr<BasisSet> bs1,
                   std::shared_ptr<BasisSet> bs2, int nderiv)
    : OneBodyAOInt(spherical_transforms, bs1, bs2, nderiv) {

    if (nderiv > 0) {
        throw std::runtime_error("NablaInt: does not support derivatives.");
    }

    int max_am = std::max(basis1()->max_am(), basis2()->max_am());
    int max_nprim = std::max(basis1()->max_nprimitive(), basis2()->max_nprimitive());

    // We're returning Ax,Ay,Az and ignoring Bx,By,Bz
    set_chunks(3);

    engine0_ = std::make_unique<libint2::Engine>(libint2::Operator::overlap, max_nprim, max_am, 1);

    buffer_ = nullptr;
    buffers_.resize(nchunk_);
}

NablaInt::~NablaInt() { }
