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

// Similar to https://github.com/psi4/psi4/pull/2076/files

#ifndef NUMINT_HELPER_H
#define NUMINT_HELPER_H
#include "psi4/libmints/typedefs.h"
#include "psi4/pragma.h"
#include <vector>
#include <map>
#include <unordered_map>
#include <string>

namespace psi {
class BasisSet;
class DFTGrid;

class NumIntHelper {
   protected:
    int print_;
    int nthread_;

    /// Local grid object
    std::shared_ptr<DFTGrid> numint_grid_;

   public:
    std::shared_ptr<DFTGrid> numint_grid() const { return numint_grid_; }

    NumIntHelper(std::shared_ptr<DFTGrid> numint_grid);

    /// Compute an integral \int -\rho(r) f(r) where f is a vector-valued function.
    /// f is represented for each block of points of the DFT integration grid as
    /// the matrix (n_data, n_points). The output has (n_data)
    SharedVector density_integral(const std::vector<SharedMatrix>& grid_data, const SharedMatrix& D) const;

    /// Same as density_integral, but don't sum over the atoms. Output is a matrix
    /// (n_atoms, n_data).
    SharedMatrix dd_density_integral(const std::vector<SharedMatrix>& grid_data, const SharedMatrix& D) const;

    /// Compute an integral \int \chi_\mu(r) \chi_\nu(r) f(r) where f is a scalar function.
    /// f is represented for each block of points of the DFT integration grid as
    /// the vector (n_points). The output has (n_bas, n_bas)
    SharedMatrix potential_integral(const std::vector<SharedVector>& grid_data) const;
};

}  // namespace psi
#endif
