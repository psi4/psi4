/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

#include <string.h>
#include "dimension.h"
#include "psi4-dec.h"
#include <string.h>

namespace psi {

Dimension::Dimension()
    : name_("(empty)"), n_(0), blocks_(0)
{

}

Dimension::Dimension(int n, const std::string &name)
    : name_(name), n_(n)
{
    blocks_ = new int[n_];
    ::memset(blocks_, 0, sizeof(int)*n_);
}

Dimension::Dimension(const std::vector<int> &v)
    : n_(static_cast<int>(v.size()))
{
    blocks_ = new int[v.size()];
    std::copy( v.begin(), v.end(), blocks_);
}

Dimension::Dimension(const Dimension &other)
    : name_(other.name_), n_(other.n_)
{
    blocks_ = new int[other.n_];
    for (int i=0; i<n_; ++i)
        blocks_[i] = other.blocks_[i];
}

Dimension::~Dimension()
{
    delete[] blocks_;
}

void Dimension::init(const std::string& name, int n)
{
    name_ = name;
    n_ = n;
    delete[] blocks_;
    blocks_ = new int[n_];
    ::memset(blocks_, 0, sizeof(int)*n_);
}

int Dimension::sum() const
{
    int s = 0;
    for (int i = 0; i < n_; i++) {
        s += blocks_[i];
    }
    return s;
}

int Dimension::max() const
{
    int s = 0;
    for (int i = 0; i < n_; i++) {
        s = (s > blocks_[i] ? s : blocks_[i]);
    }
    return s;
}

void Dimension::print() const
{
    psi::fprintf(outfile, "  %s (n = %d): ", name_.c_str(), n_);
    for (int i=0; i<n(); ++i) {
        psi::fprintf(outfile, "%d  ", blocks_[i]);
    }
    psi::fprintf(outfile, "\n");
}

Dimension& Dimension::operator =(const Dimension& other)
{
    name_ = other.name_;

    if (n_ < other.n_) {
        delete blocks_;
        blocks_ = new int[other.n_];
    }
    n_ = other.n_;
    for (int i=0; i<n_; ++i)
        blocks_[i] = other.blocks_[i];

    return *this;
}

Dimension& Dimension::operator =(const int* other)
{
    for (int i=0; i<n_; ++i)
        blocks_[i] = other[i];

    return *this;
}

Dimension& Dimension::operator+=(const Dimension& b)
{
    for (int i=0; i<n_; ++i)
        blocks_[i] += b.blocks_[i];

    return *this;
}

Dimension& Dimension::operator-=(const Dimension& b)
{
    for (int i=0; i<n_; ++i)
        blocks_[i] -= b.blocks_[i];

    return *this;
}

bool operator==(const Dimension& a, const Dimension& b) {
    if (a.n() != b.n())
        return false;
    for (int i=0; i<a.n(); ++i)
        if (a[i] != b[i])
            return false;
    return true;
}

bool operator!=(const Dimension& a, const Dimension& b) {
    return !operator==(a, b);
}

Dimension operator+(const Dimension& a, const Dimension& b) {
    Dimension result = a;
    for (int i=0; i<a.n(); ++i)
        result[i] += b[i];
    return result;
}

Dimension operator-(const Dimension& a, const Dimension& b) {
    Dimension result = a;
    for (int i=0; i<a.n(); ++i)
        result[i] -= b[i];
    return result;
}

}
