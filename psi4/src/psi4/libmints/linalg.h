/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
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

#include <iostream>
#include <array>
#include <string>
#include <type_traits>
#include <vector>

#include <xtensor/xtensor.hpp>
//#include <xtensor-blas/xlinalg.hpp>

#include "dimension.h"

namespace psi {
namespace detail {
/*! Type trait for tensorisable types
 *
 * These are used for SFINAE in the Tensor class.
 * Currently any arithmetic type can be used, which excludes complex numbers.
 */
template <typename T>
struct is_tensorisable : std::integral_constant<bool, std::is_arithmetic<T>::value> {};

// Note for posterity this can be declared inline in C++17
template <typename T>
constexpr bool is_tensorisable_v = is_tensorisable<T>::value;

using rank1 = std::integral_constant<size_t, 1>;

template <size_t Rank>
struct is_rank1 : std::integral_constant<bool, std::is_same<std::integral_constant<size_t, Rank>, rank1>::value> {};

// Note for posterity this can be declared inline in C++17
template <size_t Rank>
constexpr bool is_rank1_v = is_rank1<Rank>::value;

using rank2 = std::integral_constant<size_t, 2>;

template <size_t Rank>
struct is_rank2 : std::integral_constant<bool, std::is_same<std::integral_constant<size_t, Rank>, rank2>::value> {};

// Note for posterity this can be declared inline in C++17
template <size_t Rank>
constexpr bool is_rank2_v = is_rank2<Rank>::value;
}  // namespace detail

/*! Basic linear algebra storage object
 * \tparam T the underlying numerical type
 * \tparam Rank rank of the object, i.e. 1 for a vector, 2 for a matrix, etc.
 */
template <typename T, size_t Rank>
using Storage = xt::xtensor<T, Rank>;

/*! Symmetry-blocked linear algebra storage
 * \tparam T the underlying numerical type
 * \tparam Rank rank of the object, i.e. 1 for a vector, 2 for a matrix, etc.
 * \tparam Blocks number of symmetry blocks, i.e. 1 for C1, 4 for C2v, etc.
 *
 * \note To ensure symmetry blocking, we use a std::array. Each block has the same storage type.
 */
template <typename T, size_t Rank, size_t Blocks>
using BlockStorage = std::array<Storage<T, Rank>, Blocks>;

template <size_t Rank>
using Shape = std::array<size_t, Rank>;

template <typename T, size_t Rank>
class Tensor {
   public:
    /*! Main CTOR
     *  \param[in] name
     *  \param[in] blocks
     *  \param[in] axes_dimpi
     *
     * axes_dimpi.size() is the number of axes in each block of the storage object,
     * hence it must be equal to the Rank of the Tensor. A thing that can be statically checked
     */
    template <typename T_ = T, typename = std::enable_if_t<detail::is_tensorisable_v<T_>>>
    explicit Tensor(const std::string& name, size_t blocks, const std::array<Dimension, Rank>& axes_dimpi)
        : name_(name), axes_dimpi_(axes_dimpi), store_(blocks) {
        for (int h = 0; h < store_.size(); ++h) {
            Shape<Rank> shape{0};
            for (int ax = 0; ax < Rank; ++ax) shape[ax] = axes_dimpi[ax][h];
            store_[h] = xt::zeros<T>(shape);
        }
    }

    // FIXME int vs. size_t in these CTORs
    /*! Unnamed 1-irrep rank-1 CTOR
     *  \param[in] dim
     */
    template <size_t Rank_ = Rank, typename = std::enable_if_t<detail::is_rank1_v<Rank_>>>
    explicit Tensor(int dim) : Tensor("", 1, std::array<Dimension, Rank>{Dimension(std::vector<int>{dim})}) {}

    /*! Named 1-irrep rank-1 CTOR
     *  \param[in] name
     *  \param[in] dim
     */
    template <size_t Rank_ = Rank, typename = std::enable_if_t<detail::is_rank1_v<Rank_>>>
    explicit Tensor(const std::string& name, int dim)
        : Tensor(name, 1, std::array<Dimension, Rank>{Dimension(std::vector<int>{dim})}) {}

    /*! Unnamed 1-irrep rank-2 CTOR
     *  \param[in] rows
     *  \param[in] cols
     */
    template <size_t Rank_ = Rank, typename = std::enable_if_t<detail::is_rank2_v<Rank_>>>
    explicit Tensor(int rows, int cols)
        : Tensor("", 1,
                 std::array<Dimension, Rank>{Dimension(std::vector<int>{rows}), Dimension(std::vector<int>{cols})}) {}

    /*! Named 1-irrep rank-2 CTOR
     *  \param[in] name
     *  \param[in] rows
     *  \param[in] cols
     */
    template <size_t Rank_ = Rank, typename = std::enable_if_t<detail::is_rank2_v<Rank_>>>
    explicit Tensor(const std::string& name, int rows, int cols)
        : Tensor(name, 1,
                 std::array<Dimension, Rank>{Dimension(std::vector<int>{rows}), Dimension(std::vector<int>{cols})}) {}

    size_t dim() const {
        return std::accumulate(std::begin(store_), std::end(store_), 0,
                               [](size_t s, auto& b) { return (s + b.size()); });
    }

    size_t nirrep() const { return store_.size(); }
    std::string name() const { return name_; }
    void set_name(const std::string& name) { name_ = name; }

    /*! Returns Dimension object for given axis
     * \param[in] axis
     */
    const Dimension& axes_dimpi(size_t axis) const { return axes_dimpi_[axis]; }

    /// Returns the dimension array for a rank-1 tensor aka a vector
    template <size_t Rank_ = Rank, typename = std::enable_if_t<detail::is_rank1_v<Rank_>>>
    const Dimension& dimpi() const {
        return axes_dimpi_[0];
    }

    /*
    /// Returns the dimension array for the rows of a rank-2 tensor aka a matrix
    template <typename = typename std::enable_if_t<detail::is_rank2_v<Rank>>>
    const Dimension& rowspi() const {
        return axes_dimpi_[0];
    }

    /// Returns the dimension array for the rows of a rank-2 tensor aka a matrix
    template <typename = typename std::enable_if_t<detail::is_rank2_v<Rank>>>
    const Dimension& colspi() const {
        return axes_dimpi_[1];
    }
    */

   protected:
    std::string name_;
    std::array<Dimension, Rank> axes_dimpi_;
    std::vector<Storage<T, Rank>> store_{};
};

template <typename T, size_t Rank, size_t Blocks>
T norm(const BlockStorage<T, Rank, Blocks>& store);

// template <typename T, size_t Rank, size_t Blocks>
// T dot(const BlockStorage<T, Rank, Blocks>& bs1, const BlockStorage<T, Rank, Blocks>& bs2) {
//    T retval = 0.0;
//    for (size_t h = 0; h < Blocks; ++h) {
//        retval += xt::linalg::vdot(bs1[h], bs2[h]);
//    }
//    return retval;
//}
template <typename T, size_t Rank, size_t Blocks>
void zero(BlockStorage<T, Rank, Blocks>& bs) {
    std::for_each(std::begin(bs), std::end(bs), [](auto& b) { xt::zeros<T>(b); });
}
}  // namespace psi
