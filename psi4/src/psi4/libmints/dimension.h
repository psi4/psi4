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

#ifndef _psi_src_lib_libmints_dimension_h_
#define _psi_src_lib_libmints_dimension_h_

#include <string>
#include <cstdio>
#include <vector>

#include "psi4/pragma.h"

namespace psi {

/// @brief Dimension object
/// \ingroup MINTS
class PSI_API Dimension {
   private:
    std::string name_;
    std::vector<int> blocks_;

   public:
    /// @brief Constructs an empty Dimension object, with the name initialized to "(empty)"
    Dimension();

    /// @brief Constructs a Dimension object with a specified number of block numbers and optionally a name. If no name
    /// is given it will be a zero-length string.
    /// @param n : number of blocks
    /// @param name : (optional) name associated with this Dimension object
    Dimension(size_t n, const std::string& name = "");

    /// @brief Constructs a Dimension object from an std::vector<int> object, leaving the name a zero-length string.
    /// @param other : object to copy the block numbers from
    Dimension(const std::vector<int>& other);

    /// @brief Assignment operator, this one can be very dangerous
    PSI_DEPRECATED(
        "The assignment operator for psi::Dimension is being deprecated. Unless someone speaks up, 1.10 may be the "
        "last release to have it.")
    Dimension& operator=(const int* other);

    Dimension& operator+=(const Dimension& b);
    Dimension& operator-=(const Dimension& b);

    /// @brief Re-initializes the object. If no name is given it will be a zero-length string.
    /// @param n : number of blocks
    /// @param name : (optional) name associated with this Dimension object
    void init(size_t n, const std::string& name = "");

    /// @brief Return the rank (number of block numbers)
    size_t n() const { return blocks_.size(); }

    /// @brief Return the name of the Dimension object
    const std::string& name() const { return name_; }

    /// @brief Set the name of the Dimension object
    void set_name(const std::string& name) { name_ = name; }

    /// @brief Access a block number at a partcular index. Not bounds-checked.
    int& operator[](size_t i) { return blocks_[i]; }

    /// @brief Access a block number at a partcular index. Not bounds-checked.
    const int& operator[](size_t i) const { return blocks_[i]; }

    /// @brief Get a const reference to the std::vector storing the block numbers inside the Dimension object.
    const std::vector<int>& blocks() const { return blocks_; }

    /// @brief Casting operator to int*
    PSI_DEPRECATED(
        "Cast-to-pointer operators for psi::Dimension are being deprecated. Unless someone speaks up, 1.10 may be the "
        "last release to have them.")
    operator int*() { return blocks_.data(); }

    /// @brief Casting operator to const int*
    PSI_DEPRECATED(
        "Cast-to-pointer operators for psi::Dimension are being deprecated. Unless someone speaks up, 1.10 may be the "
        "last release to have them.")
    operator const int*() const { return blocks_.data(); }

    /// @brief Return the sum of constituent dimensions (block numbers)
    int sum() const;

    /// @brief Return the maximum of constituent dimensions (block numbers)
    int max() const;

    /// @brief Zero all the dimensions (block numbers)
    void zero();

    /// @brief Fill all elements in blocks_ with given value
    void fill(int v);

    void print() const;

    // Only used for python
    const int& get(size_t i) const { return blocks_[i]; }
    void set(size_t i, int val) { blocks_[i] = val; }

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
