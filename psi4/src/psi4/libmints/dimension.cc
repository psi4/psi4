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

#include <string.h>
#include "dimension.h"
#include "psi4/psi4-dec.h"
#include "psi4/libparallel/PsiOutStream.h"
#include "psi4/libpsi4util/exception.h"

#include <string.h>

namespace psi {

Dimension::Dimension()
    : name_("(empty)")
{
}

Dimension::Dimension(int n, const std::string &name)
    : name_(name), blocks_(n,0)
{
}

Dimension::Dimension(const std::vector<int> &v)
    : blocks_(v)
{
}

Dimension::Dimension(const Dimension &other)
    : name_(other.name_), blocks_(other.blocks_)
{
}

Dimension::~Dimension()
{
}

void Dimension::init(int n, const std::string& name)
{
    name_ = name;
    blocks_.assign(n, 0);
}

int Dimension::sum() const
{
    int s = 0;
    for (int ni : blocks_){
        s += ni;
    }
    return s;
}

int Dimension::max() const
{
    int s = 0;
    for (int ni : blocks_){
        s = std::max(ni, s);
    }
    return s;
}

void Dimension::zero()
{
    for (int& ni : blocks_){
        ni = 0;
    }
}

void Dimension::print() const
{
    outfile->Printf( "  %s (n = %d): ", name_.c_str(), n());
    for (int ni : blocks_){
        outfile->Printf( "%d  ", ni);
    }
    outfile->Printf( "\n");
}

Dimension& Dimension::operator =(const Dimension& other)
{
    name_ = other.name_;
    blocks_ = other.blocks_;
    return *this;
}

Dimension& Dimension::operator =(const int* other)
{
    for (int i = 0, maxi = n(); i < maxi; ++i)
        blocks_[i] = other[i];

    return *this;
}

Dimension& Dimension::operator+=(const Dimension& b)
{
    if (n() == b.n()){
        for (int i = 0, maxi = n(); i < maxi; ++i)
            blocks_[i] += b.blocks_[i];
    }else{
        std::string msg = "Dimension operator+=: adding operators of different size ("
                + std::to_string(n()) + " and " + std::to_string(b.n()) + ")";
        throw PSIEXCEPTION(msg);
    }

    return *this;
}

Dimension& Dimension::operator-=(const Dimension& b)
{
    if (n() == b.n()){
        for (int i = 0, maxi = n(); i < maxi; ++i)
            blocks_[i] -= b.blocks_[i];
    }else{
        std::string msg = "Dimension operator-=: subtracting operators of different size ("
                + std::to_string(n()) + " and " + std::to_string(b.n()) + ")";
        throw PSIEXCEPTION(msg);
    }

    return *this;
}

bool operator==(const Dimension& a, const Dimension& b) {
    return (a.blocks_ == b.blocks_);
}

bool operator!=(const Dimension& a, const Dimension& b) {
    return !operator==(a, b);
}

Dimension operator+(const Dimension& a, const Dimension& b) {
    Dimension result = a;
    if (a.n() == b.n()){
        for (int i = 0, maxi = a.n(); i < maxi; ++i)
            result[i] += b[i];
    }else{
        std::string msg = "Dimension operator+: adding operators of different size ("
                + std::to_string(a.n()) + " and " + std::to_string(b.n()) + ")";
        throw PSIEXCEPTION(msg);
    }

    return result;
}

Dimension operator-(const Dimension& a, const Dimension& b) {
    Dimension result = a;
    if (a.n() == b.n()){
        for (int i = 0, maxi = a.n(); i < maxi; ++i)
            result[i] -= b[i];
    }else{
        std::string msg = "Dimension operator-: subtracting operators of different size ("
                + std::to_string(a.n()) + " and " + std::to_string(b.n()) + ")";
        throw PSIEXCEPTION(msg);
    }
    return result;
}

Slice::Slice(const Dimension& begin,const Dimension& end)
    : begin_(begin), end_(end)
{
    validate_slice();
}

Slice::Slice(const Slice& other)
     : begin_(other.begin()), end_(other.end())
{
    validate_slice();
}

Slice& Slice::operator+=(const Dimension& increment)
{
    begin_ += increment;
    end_ += increment;
    validate_slice();
    return *this;
}

bool Slice::validate_slice(){
    bool valid = true;
    if (begin_.n() != end_.n()){
        valid = false;
        std::string msg = "Invalid Slice: begin and end Dimension objects have different size.";
        throw PSIEXCEPTION(msg);
    }
    // Check that begin[h] >= 0 and end[h] >= begin[h]
    for (int h = 0, max_h = begin_.n(); h < max_h; h++){
        if (begin_[h] < 0){
            valid = false;
            std::string msg = "Invalid Slice: element "
                    + std::to_string(h)
                    + " of begin Dimension object is less than zero ("
                    + std::to_string(begin_[h]) + ")";
            throw PSIEXCEPTION(msg);
        }
        if (end_[h] < begin_[h]){
            valid = false;
            std::string msg = "Invalid Slice: element "
                    + std::to_string(h)
                    + " of (end - begin) Dimension object is less than zero ("
                    + std::to_string(end_[h] - begin_[h]) + ")";
            throw PSIEXCEPTION(msg);
        }
    }
    if (not valid){
        begin_.print();
        end_.print();
    }
    return valid;
}

}
