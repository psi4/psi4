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
#pragma once

#include <array>
#include <vector>

namespace mdintegrals {
using Point = std::array<double, 3>;
void fill_E_matrix(int maxam1, int maxam2, const Point& P, const Point& A, const Point& B, double a, double b,
                   std::vector<double>& Ex, std::vector<double>& Ey, std::vector<double>& Ez);
std::vector<std::array<int, 4>> generate_am_components_cca(int am);

inline Point point_diff(const Point& A, const Point& B) { return {A[0] - B[0], A[1] - B[1], A[2] - B[2]}; }
inline int address_3d(int i, int j, int k, int dim2, int dim3) { return k + dim3 * (j + dim2 * i); }
}  // namespace mdintegrals