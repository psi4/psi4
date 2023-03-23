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

/*!
  \file
  \brief Contains some probability functions
  \ingroup QT
*/

#include "qt.h"
#include "psi4/libpsi4util/exception.h"
namespace psi {

/// @brief Returns n!
/// @param n : number to take factorial of
/// @return n factorial, as a 64-bit int (since n! can get very large)
/// \ingroup QT
uint64_t factorial(const uint64_t n) {
    if (n <= 1) return 1;
    return (n * factorial(n - 1));
}

/// @brief Calculates the number of ways to choose k objects from n objects, or "n choose k"
/// @param n : number of objects in total
/// @param k : number of objects taken at a time
/// @return number of combinations of n objects taken k at a time ("n choose k")
/// \ingroup QT
uint64_t combinations(const uint64_t n, const uint64_t k) {
    if (k > n) throw PSIEXCEPTION("Cannot compute n choose k if k > n!");
    return factorial(n) / (factorial(k) * factorial(n - k));
}
}  // namespace psi
