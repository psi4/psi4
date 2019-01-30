/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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
#include "psi4/libqt/qt.h"
#include "matrix.h"
#include "vector.h"
#include "psi4/libpsi4util/PsiOutStream.h"
using namespace psi;

IntVector::IntVector() {
    vector_ = nullptr;
    dimpi_ = nullptr;
    nirrep_ = 0;
    name_ = "";
}

IntVector::IntVector(const IntVector &c) {
    vector_ = nullptr;
    nirrep_ = c.nirrep_;
    dimpi_ = new int[nirrep_];
    for (int h = 0; h < nirrep_; ++h) dimpi_[h] = c.dimpi_[h];
    alloc();
    copy_from(c.vector_);
    name_ = c.name_;
}

IntVector::IntVector(int nirreps, int *dimpi) {
    vector_ = nullptr;
    nirrep_ = nirreps;
    dimpi_ = new int[nirrep_];
    for (int h = 0; h < nirrep_; ++h) dimpi_[h] = dimpi[h];
    alloc();
}
IntVector::IntVector(int dim) {
    vector_ = nullptr;
    nirrep_ = 1;
    dimpi_ = new int[nirrep_];
    dimpi_[0] = dim;
    alloc();
}
IntVector::IntVector(const std::string &name, int nirreps, int *dimpi) {
    vector_ = nullptr;
    nirrep_ = nirreps;
    dimpi_ = new int[nirrep_];
    for (int h = 0; h < nirrep_; ++h) dimpi_[h] = dimpi[h];
    alloc();
    name_ = name;
}
IntVector::IntVector(const std::string &name, int dim) {
    vector_ = nullptr;
    nirrep_ = 1;
    dimpi_ = new int[nirrep_];
    dimpi_[0] = dim;
    alloc();
    name_ = name;
}

IntVector::~IntVector() {
    release();
    if (dimpi_) delete[] dimpi_;
}

void IntVector::init(int nirreps, int *dimpi) {
    if (dimpi_) delete[] dimpi_;
    nirrep_ = nirreps;
    dimpi_ = new int[nirrep_];
    for (int h = 0; h < nirrep_; ++h) dimpi_[h] = dimpi[h];
    alloc();
}

void IntVector::alloc() {
    if (vector_) release();

    vector_ = (int **)malloc(sizeof(int *) * nirrep_);
    for (int h = 0; h < nirrep_; ++h) {
        if (dimpi_[h]) {
            vector_[h] = new int[dimpi_[h]];
            memset(vector_[h], 0, dimpi_[h] * sizeof(int));
        }
    }
}

void IntVector::release() {
    if (!vector_) return;

    for (int h = 0; h < nirrep_; ++h) {
        if (dimpi_[h]) delete[](vector_[h]);
    }
    free(vector_);
    vector_ = nullptr;
}

void IntVector::copy_from(int **c) {
    size_t size;
    for (int h = 0; h < nirrep_; ++h) {
        size = dimpi_[h] * sizeof(int);
        if (size) memcpy(&(vector_[h][0]), &(c[h][0]), size);
    }
}

void IntVector::copy(const IntVector *rhs) {
    if (nirrep_ != rhs->nirrep_) {
        release();
        if (dimpi_) delete[] dimpi_;
        nirrep_ = rhs->nirrep_;
        dimpi_ = new int[nirrep_];
        for (int h = 0; h < nirrep_; ++h) dimpi_[h] = rhs->dimpi_[h];
        alloc();
    }
    copy_from(rhs->vector_);
}
void IntVector::copy(const IntVector &rhs) {
    if (nirrep_ != rhs.nirrep_) {
        release();
        if (dimpi_) delete[] dimpi_;
        nirrep_ = rhs.nirrep_;
        dimpi_ = new int[nirrep_];
        for (int h = 0; h < nirrep_; ++h) dimpi_[h] = rhs.dimpi_[h];
        alloc();
    }
    copy_from(rhs.vector_);
}

void IntVector::set(int *vec) {
    int h, i, ij;

    ij = 0;
    for (h = 0; h < nirrep_; ++h) {
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
    for (h = 0; h < nirrep_; ++h) {
        printer->Printf(" Irrep: %d\n", h + 1);
        for (int i = 0; i < dimpi_[h]; ++i) printer->Printf("   %4d: %10d\n", i + 1, vector_[h][i]);
        printer->Printf("\n");
    }
}

int *IntVector::to_block_vector() {
    size_t size = 0;
    for (int h = 0; h < nirrep_; ++h) size += dimpi_[h];

    auto *temp = new int[size];
    size_t offset = 0;
    for (int h = 0; h < nirrep_; ++h) {
        for (int i = 0; i < dimpi_[h]; ++i) {
            temp[i + offset] = vector_[h][i];
        }
        offset += dimpi_[h];
    }

    return temp;
}
