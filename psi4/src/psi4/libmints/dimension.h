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

#ifndef _psi_src_lib_libmints_dimension_h_
#define _psi_src_lib_libmints_dimension_h_

#include <string>
#include <cstdio>
#include <vector>

#include "psi4/pragma.h"

namespace psi {

class PSI_API Dimension {
   private:
    std::string name_;
    std::vector<int> blocks_;

   public:
    Dimension();
    Dimension(int n, const std::string& name = "");
    Dimension(const std::vector<int>& other);

    /// Assignment operator, this one can be very dangerous
    Dimension& operator=(const int* other);

    Dimension& operator+=(const Dimension& b);
    Dimension& operator-=(const Dimension& b);

    /// Re-initializes the object
    void init(int n, const std::string& name = "");

    /// Return the rank
    int n() const { return static_cast<int>(blocks_.size()); }

    /// Return the name of the dimension
    const std::string& name() const { return name_; }

    /// Set the name of the dimension
    void set_name(const std::string& name) { name_ = name; }

    /// Blocks access
    int& operator[](int i) { return blocks_[i]; }
    const int& operator[](int i) const { return blocks_[i]; }
    const std::vector<int>& blocks() const { return blocks_; }

    /// Casting operator to int*
    operator int*() { return blocks_.data(); }
    /// Casting operator to const int*
    operator const int*() const { return blocks_.data(); }

    /// Return the sum of constituent dimensions
    int sum() const;
    int max() const;

    /// Zero all the elements
    void zero();

    /// Fill all elements in blocks_ with given value
    void fill(int v);

    void print() const;

    // Only used for python
    const int& get(int i) const { return blocks_[i]; }
    void set(int i, int val) { blocks_[i] = val; }

    PSI_API friend bool operator==(const Dimension& a, const Dimension& b);
    PSI_API friend bool operator!=(const Dimension& a, const Dimension& b);
    PSI_API friend Dimension operator+(const Dimension& a, const Dimension& b);
    PSI_API friend Dimension operator-(const Dimension& a, const Dimension& b);
};

/*! \ingroup MINTS
 *  \class Slice
 *  \brief Slicing for Matrices and Vectors objects.
 *
 *  Slices are pairs of Dimension objects used to manipulate parts of vectors
 *  and matrices.
 *
 *  Slices can be constructed from Dimension objects:
 *
 *      Dimension begin;
 *      Dimension end;
 *      Slice slice(begin,end);
 *
 *  or can be implicitly constructed from pairs of Dimension objects. E.g.:
 *
 *      Dimension begin;
 *      Dimension end;
 *      SharedVector v;
 *      v->get_block({begin,end}); // same as v->get_block(slice);
 */
class PSI_API Slice {
    Dimension begin_;
    Dimension end_;

   public:
    /// Creator
    /// Rules: begin must satisfly begin[h] >= 0 for all h
    ///        end must satisfly end[h] >= begin[h] for all h
    Slice(const Dimension& begin, const Dimension& end);
    /// Copy constructor
    Slice(const Slice& other);

    /// Get the first element of this slice
    const Dimension& begin() const { return begin_; }
    /// Get the past-the-end element of this slice
    const Dimension& end() const { return end_; }
    /// Increment the beginning and end of this slice
    Slice& operator+=(const Dimension& increment);

   private:
    /// Check if this Slice is acceptable
    bool validate_slice();
};
}  // namespace psi

#endif  // _psi_src_lib_libmints_dimension_h_
