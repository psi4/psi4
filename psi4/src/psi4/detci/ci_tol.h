/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

/*! \file
    \ingroup DETCI
    \brief Enter brief description of file here
*/

#pragma once

namespace psi {
namespace detci {
constexpr auto S_MAX = 0.999995; /* Max overlap for mitrush vects and collapsed vecs */
constexpr auto HD_MIN = 1.0E-4;  /* Minimum diagonal element of preconditioner       */
constexpr auto ZERO = 1e-10;
constexpr auto MPn_ZERO = 1e-14;
constexpr auto SA_NORM_TOL = 1.0E-5; /* norm of schmidt orthogonalized vector */
constexpr auto MPn_NORM_TOL = 1.0E-12;
}  // namespace detci
}  // namespace psi
