/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

#ifndef _psi4_libmints_libcint_util_h
#define _psi4_libmints_libcint_util_h

#include <vector>

namespace psi {

class BasisSet;

/// Helpers shared by the libcint one- and two-electron integral backends. Only
/// plain ints and doubles cross this header, so it stays free of <cint.h>.
namespace libcint_util {

/// Number of cartesian functions of angular momentum l.
inline int ncart(int l) { return (l + 1) * (l + 2) / 2; }

/// psi4 / libint2-Gaussian spherical slot for magnetic quantum number m
/// (m == 0 -> 0, m > 0 -> 2m-1, m < 0 -> -2m). All relative signs are +1
/// (machine-verified against libint2 in Gaussian ordering, l = 0..4).
inline int psi4_slot(int m) { return m == 0 ? 0 : (m > 0 ? 2 * m - 1 : -2 * m); }

/// Magnetic quantum number of libcint spherical slot i for an l-shell. libcint
/// orders l >= 2 as m = -l..+l (so m = i - l), but p functions (l == 1) are the
/// special cartesian-like order (px, py, pz) = m(+1, -1, 0). Verified against
/// libint2 (Gaussian ordering) to ~1e-14 for l = 0..4.
inline int libcint_m(int i, int l) {
    if (l == 1) return i == 0 ? +1 : (i == 1 ? -1 : 0);
    return i - l;
}

/// psi4 slot of libcint function i of an l-shell, cartesian (identity: libcint's
/// cartesian order is psi4's CartesianIter order) or spherical.
inline int psi4_index(int i, int l, bool cart) { return cart ? i : psi4_slot(libcint_m(i, l)); }

/// Factor converting a Libint2 contraction coefficient into the libcint coefficient
/// describing the same function, for angular momentum l (see the definition).
double libint2_to_libcint(int l, bool cart);

/// Append one psi4 basis set's atoms (its own molecule's centers) and shells to
/// libcint's atm/bas/env arrays, with exactly the coefficients the Libint2 path
/// uses. coef_ratio[l] is libint2_to_libcint(l, cart). Records per shell whether
/// psi4 wants it spherical in bas_pure. Returns the bas index of the basis's shell 0.
int append_basis(const BasisSet &bs, const std::vector<double> &coef_ratio, std::vector<int> &atm,
                 std::vector<int> &bas, std::vector<double> &env, std::vector<char> &bas_pure);

}  // namespace libcint_util
}  // namespace psi

#endif
