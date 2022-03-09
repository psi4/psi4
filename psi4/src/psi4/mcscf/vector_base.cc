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

#include "vector_base.h"

#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/psi4-dec.h"

#include <cstring>
#include <iostream>

extern FILE* outfile;

namespace psi {
namespace mcscf {

extern MemoryManager* memory_manager;

VectorBase::VectorBase(int elements) : elements_(elements), vector_(nullptr) { allocate1(double, vector_, elements_); }

VectorBase::~VectorBase() { release1(vector_); }

void VectorBase::print() {
    outfile->Printf("\n  ");
    for (size_t i = 0; i < elements_; ++i) {
        outfile->Printf("%10.6f", vector_[i]);
    }
}

void VectorBase::copy(VectorBase& source) {
    for (size_t i = 0; i < elements_; ++i) vector_[i] = source.vector_[i];
}

}  // namespace mcscf
}  // namespace psi
