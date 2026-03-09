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

#ifndef AMBIT_CONVERTER_H
#define AMBIT_CONVERTER_H

#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"

namespace ambit {

class Tensor;

namespace helpers {

namespace psi4 {

void convert(const psi::Matrix &matrix, ambit::Tensor *target);

void convert(const psi::Vector &vector, ambit::Tensor *target);

}  // namespace psi4

}  // namespace helpers

}  // namespace ambit

#endif  // AMBIT_CONVERTER_H
