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

#include <array>

#include <xtensor/xtensor.hpp>
//#include <xtensor-blas/xlinalg.hpp>

namespace psi {
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
