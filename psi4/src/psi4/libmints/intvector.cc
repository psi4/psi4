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

#include <cstdlib>
#include <cstring>
#include <numeric>
#include "psi4/libqt/qt.h"
#include "vector.h"
#include "psi4/libpsi4util/PsiOutStream.h"
using namespace psi;

void IntVector::copy(const IntVector &other) {
    dimpi_ = other.dimpi_;
    v_ = other.v_;
    assign_pointer_offsets();
}

void IntVector::set(int *vec) {
    int h, i, ij;

    ij = 0;
    for (h = 0; h < nirrep(); ++h) {
        for (i = 0; i < dimpi_[h]; ++i) {
            vector_[h][i] = vec[ij++];
        }
    }
}

void IntVector::print(std::string out, const char *extra) const {
    int h;
    std::shared_ptr<psi::PsiOutStream> printer = (out == "outfile" ? outfile : std::make_shared<PsiOutStream>(out));
    if (extra == nullptr) {
        printer->Printf("\n # %s #\n", name_.c_str());
    } else {
        printer->Printf("\n # %s %s #\n", name_.c_str(), extra);
    }
    for (h = 0; h < nirrep(); ++h) {
        printer->Printf(" Irrep: %d\n", h + 1);
        for (int i = 0; i < dimpi_[h]; ++i) printer->Printf("   %4d: %10d\n", i + 1, vector_[h][i]);
        printer->Printf("\n");
    }
}

IntVector IntVector::iota(const Dimension &dim) {
    IntVector vec(dim);
    for (int h = 0; h < dim.n(); h++) {
        std::iota(vec.pointer(h), vec.pointer(h) + dim[h], 0);
    }
    return vec;
}

