/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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

//
// Created by Justin Turney on 1/5/16.
//

#include "convert.h"
#include <ambit/tensor.h>
//#include <tensor/core/core.h>

namespace ambit {

namespace helpers {

namespace psi4 {

void PSI_API convert(const psi::Matrix &matrix, ambit::Tensor *target) {
    if (target->rank() != 2) throw std::runtime_error("convert(psi::Matrix, ambit::Tensor): Tensor is not rank 2");

    if (matrix.nirrep() != 1)
        throw std::runtime_error(
            "convert(psi::Matrix, ambit::Tensor): Matrix "
            "appears to have symmetry (nirrep != 1)");

    if (matrix.rowdim() != target->dim(0))
        throw std::runtime_error(
            "convert(psi::Matrix, ambit::Tensor): Matrix "
            "and Tensor do not have the same number of "
            "rows (dim(0))");
    if (matrix.coldim() != target->dim(1))
        throw std::runtime_error(
            "convert(psi::Matrix, ambit::Tensor): Matrix "
            "and Tensor do not have the same number of "
            "columns (dim(1))");

    size_t row = target->dim(0);
    size_t col = target->dim(1);

    Tensor local_tensor = Tensor::build(CoreTensor, "Local Data", {row, col});

    if (row && col) {
        // copy data from SharedMatrix to local_tensor
        std::copy(matrix.pointer()[0], matrix.pointer()[0] + (row * col), local_tensor.data().begin());
    }

    // Splice data into the target tensor
    (*target)() = local_tensor();
}

void PSI_API convert(const psi::Vector &vector, ambit::Tensor *target) {
    if (target->rank() != 1) throw std::runtime_error("convert(psi::Matrix, ambit::Tensor): Tensor is not rank 1");

    if (vector.nirrep() != 1)
        throw std::runtime_error(
            "convert(psi::Matrix, ambit::Tensor): Matrix "
            "appears to have symmetry (nirrep != 1)");

    if (vector.dim() != target->dim(0))
        throw std::runtime_error(
            "convert(psi::Matrix, ambit::Tensor): Matrix "
            "and Tensor do not have the same number of "
            "elements (dim(0))");

    size_t row = target->dim(0);

    Tensor local_tensor = Tensor::build(CoreTensor, "Local Data", {row});

    // copy data from SharedMatrix to local_tensor
    std::copy(vector.pointer(), vector.pointer() + (row), local_tensor.data().begin());

    // Splice data into the target tensor
    (*target)() = local_tensor();
}

}  // namespace psi4

}  // namespace helpers

}  // namespace ambit
