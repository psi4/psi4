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

#include "psi4/libpsi4util/libpsi4util.h"
#include <cstdio>

#include "blas.h"
#include "matrix.h"
#include "psi4/psi4-dec.h"
namespace psi {

namespace psimrcc {

double* CCOperation::local_work = nullptr;
double* CCOperation::out_of_core_buffer = nullptr;
double CCOperation::zero_timing = 0.0;
double CCOperation::numerical_timing = 0.0;
double CCOperation::contract_timing = 0.0;
double CCOperation::tensor_timing = 0.0;
double CCOperation::dot_timing = 0.0;
double CCOperation::plus_timing = 0.0;
double CCOperation::product_timing = 0.0;
double CCOperation::division_timing = 0.0;
double CCOperation::sort_timing = 0.0;
double CCOperation::PartA_timing = 0.0;
double CCOperation::PartB_timing = 0.0;
double CCOperation::PartC_timing = 0.0;

CCOperation::CCOperation(double in_factor, std::string in_assignment, std::string in_reindexing,
                         std::string in_operation, CCMatrix* in_A_Matrix, CCMatrix* in_B_Matrix, CCMatrix* in_C_Matrix,
                         double* work, double* buffer)
    : factor(in_factor),
      assignment(in_assignment),
      reindexing(in_reindexing),
      operation(in_operation),
      A_Matrix(in_A_Matrix),
      B_Matrix(in_B_Matrix),
      C_Matrix(in_C_Matrix) {
    if (in_B_Matrix != nullptr && in_A_Matrix->wfn() != in_B_Matrix->wfn()) {
        throw PSIEXCEPTION("\nError in CCOperation construction. Matrices A and B must have the same wavefunction.");
    } else if (in_C_Matrix != nullptr && in_A_Matrix->wfn() != in_C_Matrix->wfn()) {
        throw PSIEXCEPTION("\nError in CCOperation construction. Matrices A and C must have the same wavefunction.");
    } else {
        wfn_ = in_A_Matrix->wfn();
    }
    local_work = work;
    out_of_core_buffer = buffer;
}

CCOperation::~CCOperation() {}

void CCOperation::print() {
    if (reindexing.size()) outfile->Printf("\n\tReindexing = %s", reindexing.c_str());
    outfile->Printf("\n\tNumericalFactor = %lf", factor);
    outfile->Printf("\tAssignment = %s", assignment.c_str());
    outfile->Printf("\tOperation = %s", operation.c_str());
    outfile->Printf("\n\tA = %s", A_Matrix->get_label().c_str());
    if (B_Matrix != nullptr) outfile->Printf("\tB = %s", B_Matrix->get_label().c_str());
    if (C_Matrix != nullptr) outfile->Printf("\tC = %s", C_Matrix->get_label().c_str());
}

void CCOperation::print_operation() {
    outfile->Printf("%s", A_Matrix->get_label().c_str());
    outfile->Printf(" %s", assignment.c_str());
    if (reindexing.size()) outfile->Printf(" %s", reindexing.c_str());
    outfile->Printf(" %lf", factor);
    if (B_Matrix != nullptr) outfile->Printf(" %s", B_Matrix->get_label().c_str());
    outfile->Printf(" %s", operation.c_str());
    if (C_Matrix != nullptr) outfile->Printf(" %s", C_Matrix->get_label().c_str());
}

}  // namespace psimrcc
}  // namespace psi
