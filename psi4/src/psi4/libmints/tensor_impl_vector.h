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

#include <xtensor/xadapt.hpp>

#include "dimension.h"
#include "tensor_impl.h"
#include "vector.h"

namespace psi {
template <typename T, size_t Rank>
class Tensor;

template <typename T>
using Vector_ = Tensor<T, 1>;

template <typename T>
using SharedVector_ = std::shared_ptr<Vector_<T>>;

namespace detail {
template <typename T>
struct RankDependentImpl<Vector_<T>> {
    T get(size_t h, size_t i) const { return static_cast<const Vector_<T>*>(this)->store_[h](i); }
    T get(size_t i) const { return this->get(0, i); }

    void set(size_t h, size_t i, T val) { static_cast<Vector_<T>*>(this)->store_[h](i) = val; }
    void set(size_t i, T val) { this->set(0, i, val); }

    /// Returns the dimension array for a rank-1 tensor aka a vector
    const Dimension& dimpi() const { return static_cast<const Vector_<T>*>(this)->axes_dimpi_[0]; }

    /*! C++ string representation of type */
    std::string cxxClassName() const noexcept { return "Vector<" + detail::Type2String<T>::full() + ">"; }
    /*! Python string representation of type */
    static std::string pyClassName() noexcept { return "Vector_" + detail::Type2String<T>::suffix(); }
};
}  // namespace detail

/*! Conversion from Vector to Vector_<double>
 *  \param[in] v
 */
template <typename T = double>
auto transmute(const Vector& v) -> Vector_<T> {
    auto vec = Vector_<T>(v.name(), v.dimpi(), 0);
    // Adapt memory owned by vector class
    for (int h = 0; h < v.nirrep(); ++h) {
        std::vector<std::size_t> shape = {v.dim(h)};
        vec[h] = xt::adapt(v.pointer(h), v.dim(h), xt::no_ownership(), shape);
    }
    return vec;
}

/*! Conversion from SharedVector to SharedVector_<double>
 *  \param[in] v
 */
template <typename T = double>
auto transmute(const std::shared_ptr<Vector>& v) -> SharedVector_<T> {
    auto vec = std::make_shared<Vector_<T>>(v->name(), v->dimpi(), 0);
    // Adapt memory owned by Vector class
    for (int h = 0; h < v->nirrep(); ++h) {
        std::vector<std::size_t> shape = {v->dim(h)};
        vec->block(h) = xt::adapt(v->pointer(h), v->dim(h), xt::no_ownership(), shape);
    }
    return vec;
}

/*! Conversion from Vector_<double> to Vector
 *  \param[in] v
 */
template <typename T = double>
auto transmute(const Vector_<T>& v) -> Vector {
    auto vec = Vector(v.label(), v.dimpi());
    // Element-by-element copy
    for (int h = 0; h < v.nirrep(); ++h) {
        for (auto i = 0; i < v.dim(h); ++i) {
            vec.set(h, i, v.get(h, i));
        }
    }
    return vec;
}

/*! Conversion from SharedVector_<double> to SharedVector
 *  \param[in] v
 */
template <typename T = double>
auto transmute(const SharedVector_<T>& v) -> SharedVector {
    auto vec = std::make_shared<Vector>(v->label(), v->dimpi());
    // Element-by-element copy
    for (int h = 0; h < v->nirrep(); ++h) {
        for (auto i = 0; i < v->dim(h); ++i) {
            vec->set(h, i, v->get(h, i));
        }
    }
    return vec;
}
}  // namespace psi
