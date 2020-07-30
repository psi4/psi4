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

#include "dimension.h"

#include <string>

#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/exception.h"

namespace psi {
Dimension::Dimension(Dimension::value_type n, const std::string& name) : name_(name), blocks_(n, 0) {}

Dimension::Dimension(const storage& v) : blocks_(v) {}

Dimension::Dimension(const Dimension& other) : name_(other.name_), blocks_(other.blocks_) {}

void Dimension::init(Dimension::value_type n, const std::string& name) {
    name_ = name;
    blocks_.assign(n, 0);
}

void Dimension::print() const {
    outfile->Printf("  %s (n = %zu): ", name_.c_str(), n());
    for (auto ni : blocks_) {
        outfile->Printf("%d  ", ni);
    }
    outfile->Printf("\n");
}

Dimension& Dimension::operator=(const Dimension& other) {
    name_ = other.name_;
    blocks_ = other.blocks_;
    return *this;
}

Dimension& Dimension::operator=(const Dimension::value_type* other) {
    for (size_t i = 0, maxi = n(); i < maxi; ++i) blocks_[i] = other[i];

    return *this;
}

Dimension& Dimension::operator+=(const Dimension& b) {
    if (n() == b.n()) {
        for (auto i = 0; i < n(); ++i) blocks_[i] += b.blocks_[i];
    } else {
        std::string msg = "Dimension operator+=: adding operators of different size (" + std::to_string(n()) + " and " +
                          std::to_string(b.n()) + ")";
        throw PSIEXCEPTION(msg);
    }

    return *this;
}

Dimension& Dimension::operator-=(const Dimension& b) {
    if (n() == b.n()) {
        for (auto i = 0; i < n(); ++i) blocks_[i] -= b.blocks_[i];
    } else {
        std::string msg = "Dimension operator-=: subtracting operators of different size (" + std::to_string(n()) +
                          " and " + std::to_string(b.n()) + ")";
        throw PSIEXCEPTION(msg);
    }

    return *this;
}

PSI_API bool operator==(const Dimension& a, const Dimension& b) { return (a.blocks_ == b.blocks_); }

PSI_API bool operator!=(const Dimension& a, const Dimension& b) { return !operator==(a, b); }

PSI_API Dimension operator+(const Dimension& a, const Dimension& b) {
    Dimension result = a;
    if (a.n() == b.n()) {
        for (auto i = 0; i < a.n(); ++i) result[i] += b[i];
    } else {
        std::string msg = "Dimension operator+: adding operators of different size (" + std::to_string(a.n()) +
                          " and " + std::to_string(b.n()) + ")";
        throw PSIEXCEPTION(msg);
    }

    return result;
}

PSI_API Dimension operator-(const Dimension& a, const Dimension& b) {
    Dimension result = a;
    if (a.n() == b.n()) {
        for (auto i = 0; i < a.n(); ++i) result[i] -= b[i];
    } else {
        std::string msg = "Dimension operator-: subtracting operators of different size (" + std::to_string(a.n()) +
                          " and " + std::to_string(b.n()) + ")";
        throw PSIEXCEPTION(msg);
    }
    return result;
}

Slice::Slice(const Dimension& begin, const Dimension& end) : begin_(begin), end_(end) { validate_slice(); }

Slice::Slice(const Slice& other) : begin_(other.begin()), end_(other.end()) { validate_slice(); }

Slice& Slice::operator+=(const Dimension& increment) {
    begin_ += increment;
    end_ += increment;
    validate_slice();
    return *this;
}

bool Slice::validate_slice() {
    bool valid = true;
    std::string msg;
    if (begin_.n() != end_.n()) {
        valid = false;
        msg = "Invalid Slice: begin and end Dimension objects have different size.";
        begin_.print();
        end_.print();
        throw PSIEXCEPTION(msg);
    }

    // Check that begin[h] >= 0 and end[h] >= begin[h]
    for (auto h = 0; h < begin_.n(); h++) {
        if (begin_[h] < 0) {
            valid = false;
            msg = "Invalid Slice: element " + std::to_string(h) + " of begin Dimension object is less than zero (" +
                  std::to_string(begin_[h]) + ")";
            break;
        }
        if (end_[h] < begin_[h]) {
            valid = false;
            msg = "Invalid Slice: element " + std::to_string(h) +
                  " of (end - begin) Dimension object is less than zero (" + std::to_string(end_[h] - begin_[h]) + ")";
            break;
        }
    }
    if (!valid) {
        begin_.print();
        end_.print();
        throw PSIEXCEPTION(msg);
    }
    return valid;
}
}  // namespace psi
