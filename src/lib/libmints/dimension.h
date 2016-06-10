/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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

#ifndef _psi_src_lib_libmints_dimension_h_
#define _psi_src_lib_libmints_dimension_h_

#include <string>
#include <cstdio>
#include <vector>

namespace psi {



class Dimension
{
    std::string name_;
    int n_;
    int *blocks_;

public:
    Dimension();
    Dimension(const Dimension& other);
    Dimension(int n, const std::string& name = "");
    Dimension(const std::vector<int>& other);
    ~Dimension();

    /// Assignment operator
    Dimension& operator=(const Dimension& other);

    /// Assignment operator, this one can be very dangerous
    Dimension& operator=(const int* other);

    Dimension& operator+=(const Dimension& b);
    Dimension& operator-=(const Dimension& b);

    /// Re-initializes the object.
    void init(const std::string&, int n);

    /// Return the rank
    const int& n() const { return n_; }

    /// Return the name of the dimension.
    const std::string& name() const { return name_; }

    /// Set the name of the dimension.
    void set_name(const std::string& nme) { name_ = nme; }

    /// Blocks access
    int& operator[](int i) { return blocks_[i]; }
    const int& operator[](int i) const { return blocks_[i]; }

    /// Casting operator to int*
    operator int*() const { return blocks_; }
    /// Casting operator to const int*
    operator const int*() const { return blocks_; }

    /// Return the sum of constituent dimensions
    int sum() const;
    int max() const;

    void print() const;

    // Only used for python
    const int& get(int i) const { return blocks_[i]; }
    void set(int i, int val) { blocks_[i] = val; }
};

bool operator==(const Dimension& a, const Dimension& b);
bool operator!=(const Dimension& a, const Dimension& b);
Dimension operator+(const Dimension& a, const Dimension& b);
Dimension operator-(const Dimension& a, const Dimension& b);

}

#endif // _psi_src_lib_libmints_dimension_h_