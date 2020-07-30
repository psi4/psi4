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
#include <vector>

#include "tensor.h"

namespace psi {
template <typename T>
using Shared = std::shared_ptr<T>;

class BasisSet;
class IntegralFactory;
class OneBodyAOInt;
class TwoBodyAOInt;

/*! Methods returning objects with xtensor-based storage (implementation in xthelper.cc) */
template <typename Mints, typename T = double, typename = std::enable_if_t<detail::is_tensorisable_v<T>>>
class xtHelper {
   public:
    auto ao_overlap_() const -> SharedMatrix_<T>;
    auto ao_overlap_(Shared<BasisSet>, Shared<BasisSet>) const -> SharedMatrix_<T>;

    auto ao_kinetic_() const -> SharedMatrix_<T>;
    auto ao_kinetic_(Shared<BasisSet>, Shared<BasisSet>) const -> SharedMatrix_<T>;

    auto ao_potential_() const -> SharedMatrix_<T>;
    auto ao_potential_(Shared<BasisSet>, Shared<BasisSet>) const -> SharedMatrix_<T>;

    /// AO ERI Integrals (Full matrix, not recommended for large systems)
    auto ao_eri_(Shared<IntegralFactory> = nullptr) const -> SharedTensor<T, 4>;
    auto ao_eri_(Shared<BasisSet> bs1, Shared<BasisSet> bs2, Shared<BasisSet> bs3, Shared<BasisSet> bs4) const
        -> SharedTensor<T, 4>;

    /// AO ERI Shell
    auto ao_eri_shell_(int M, int N, int P, int Q) -> SharedTensor<T, 4>;

    /// SO Overlap Integrals
    auto so_overlap_(bool include_perturbations = true) -> SharedMatrix_<T>;
    /// SO Kinetic Integrals
    auto so_kinetic_(bool include_perturbations = true) -> SharedMatrix_<T>;
    /// SO Potential Integrals
    auto so_potential_(bool include_perturbations = true) -> SharedMatrix_<T>;

   private:
    auto one_body_ao_computer(std::vector<Shared<OneBodyAOInt>> ints, SharedMatrix_<T> out, bool symm) const -> void;

    auto ao_helper(const std::string& label, Shared<TwoBodyAOInt> ints) const -> SharedTensor<T, 4>;
    auto ao_shell_getter(const std::string& label, Shared<TwoBodyAOInt> ints, int M, int N, int P, int Q) const
        -> SharedTensor<T, 4>;
};
}  // namespace psi
