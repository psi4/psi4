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

//
// Created by Justin Turney on 12/17/15.
//

#ifndef AMBIT_INTEGRALS_H
#define AMBIT_INTEGRALS_H

#include "psi4/libmints/onebody.h"
#include "psi4/libmints/twobody.h"

namespace ambit {

class Tensor;

namespace helpers {

namespace psi4 {

void integrals(psi::OneBodyAOInt &integral, ambit::Tensor *target);

void integrals(psi::TwoBodyAOInt &integral, ambit::Tensor *target);

}  // namespace psi4

}  // namespace helpers

}  // namespace ambit

#endif  // AMBIT_INTEGRALS_H
