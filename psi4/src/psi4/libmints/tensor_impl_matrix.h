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

#pragma once

#include <memory>
#include <string>

#include "dimension.h"
#include "tensor_impl.h"

namespace psi {
template <typename T, size_t Rank>
class Tensor;

template <typename T, size_t Rank>
using SharedTensor = std::shared_ptr<Tensor<T, Rank>>;

template <typename T>
using Matrix_ = Tensor<T, 2>;

template <typename T>
using SharedMatrix_ = SharedTensor<T, 2>;

namespace detail {
template <typename T>
struct RankDependentImpl<Matrix_<T>> {
    T get(size_t h, size_t i, size_t j) const { return static_cast<const Matrix_<T>*>(this)->store_[h](i, j); }
    T get(size_t i, size_t j) const { return this->get(0, i, j); }

    void set(size_t h, size_t i, size_t j, T val) { static_cast<Matrix_<T>*>(this)->store_[h](i, j) = val; }
    void set(size_t i, size_t j, T val) { this->set(0, i, j, val); }

    /*! C++ string representation of type */
    std::string cxxClassName() const noexcept { return "Matrix<" + detail::Type2String<T>::full() + ">"; }
    /*! Python string representation of type */
    static std::string pyClassName() noexcept { return "Matrix_" + detail::Type2String<T>::suffix(); }

    /// Returns the dimension array for the rows of a rank-2 tensor aka a matrix
    const Dimension& rowspi() const { return static_cast<const Matrix_<T>*>(this)->axes_dimpi_[0]; }
    /*! Returns the number of rows in given irrep
     *  \param[in] h
     */
    size_t rows(size_t h = 0) const { return static_cast<const Matrix_<T>*>(this)->axes_dimpi_[0][h]; }

    /// Returns the dimension array for the rows of a rank-2 tensor aka a matrix
    const Dimension& colspi() const { return static_cast<const Matrix_<T>*>(this)->axes_dimpi_[1]; }
    /*! Returns the number of columns in given irrep
     *  \param[in] h
     */
    size_t cols(size_t h = 0) const { return static_cast<const Matrix_<T>*>(this)->axes_dimpi_[1][h]; }
};
}  // namespace detail
}  // namespace psi
