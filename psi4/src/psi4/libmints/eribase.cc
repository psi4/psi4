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

#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libmints/eri.h"
#include "psi4/libmints/integral.h"
#include "psi4/libmints/basisset.h"
#include "psi4/libmints/fjt.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <libint2/shell.h>
#include <libint2/engine.h>

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>

#define MAX(a, b) ((a) > (b) ? (a) : (b))

// libderiv computes 9 of the 12 total derivatives. It computes 3 of the
// centers we handle the 4th.
#define ERI_1DER_NTYPE (9)
// libderiv computes the second derivatives using the first derivatives.
// The first derivatives are provided when second derivatives are asked
// for.
#define ERI_2DER_NTYPE (ERI_1DER_NTYPE + 45)
;
using namespace psi;

namespace {

unsigned char ntypes[] = {1, ERI_1DER_NTYPE, ERI_2DER_NTYPE};

#ifdef ENABLE_Libint1t
/**
 * @brief Takes care of the changing the results buffer for any reordering done for libderiv.
 * @param permutation How the shells were permuted to satisfy Libint angular momentum requirements.
 * @param libderiv_ The libderiv buffer that contains the integrals
 * @param source_ The destination for the integrals.
 * @param size Total number of integrals computed.
 */
static void handle_reordering1(PermutedOrder permutation, Libderiv_t &libderiv_, double *source_, size_t size) {
    switch (permutation) {
        case ABCD:
            // Ax
            memcpy(source_ + 0 * size, libderiv_.ABCD[0], sizeof(double) * size);
            // Ay
            memcpy(source_ + 1 * size, libderiv_.ABCD[1], sizeof(double) * size);
            // Az
            memcpy(source_ + 2 * size, libderiv_.ABCD[2], sizeof(double) * size);
            // Cx
            memcpy(source_ + 3 * size, libderiv_.ABCD[6], sizeof(double) * size);
            // Cy
            memcpy(source_ + 4 * size, libderiv_.ABCD[7], sizeof(double) * size);
            // Cz
            memcpy(source_ + 5 * size, libderiv_.ABCD[8], sizeof(double) * size);
            // Dx
            memcpy(source_ + 6 * size, libderiv_.ABCD[9], sizeof(double) * size);
            // Dy
            memcpy(source_ + 7 * size, libderiv_.ABCD[10], sizeof(double) * size);
            // Dz
            memcpy(source_ + 8 * size, libderiv_.ABCD[11], sizeof(double) * size);
            break;
        case BACD:
            // Ax
            C_DAXPY(size, -1.0, libderiv_.ABCD[0], 1, source_ + 0 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[6], 1, source_ + 0 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[9], 1, source_ + 0 * size, 1);
            // Ay
            C_DAXPY(size, -1.0, libderiv_.ABCD[1], 1, source_ + 1 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[7], 1, source_ + 1 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[10], 1, source_ + 1 * size, 1);
            // Az
            C_DAXPY(size, -1.0, libderiv_.ABCD[2], 1, source_ + 2 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[8], 1, source_ + 2 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[11], 1, source_ + 2 * size, 1);
            // Cx
            memcpy(source_ + 3 * size, libderiv_.ABCD[6], sizeof(double) * size);
            // Cy
            memcpy(source_ + 4 * size, libderiv_.ABCD[7], sizeof(double) * size);
            // Cz
            memcpy(source_ + 5 * size, libderiv_.ABCD[8], sizeof(double) * size);
            // Dx
            memcpy(source_ + 6 * size, libderiv_.ABCD[9], sizeof(double) * size);
            // Dy
            memcpy(source_ + 7 * size, libderiv_.ABCD[10], sizeof(double) * size);
            // Dz
            memcpy(source_ + 8 * size, libderiv_.ABCD[11], sizeof(double) * size);
            break;
        case ABDC:
            // Ax
            memcpy(source_ + 0 * size, libderiv_.ABCD[0], sizeof(double) * size);
            // Ay
            memcpy(source_ + 1 * size, libderiv_.ABCD[1], sizeof(double) * size);
            // Az
            memcpy(source_ + 2 * size, libderiv_.ABCD[2], sizeof(double) * size);
            // Cx
            memcpy(source_ + 3 * size, libderiv_.ABCD[9], sizeof(double) * size);
            // Cy
            memcpy(source_ + 4 * size, libderiv_.ABCD[10], sizeof(double) * size);
            // Cz
            memcpy(source_ + 5 * size, libderiv_.ABCD[11], sizeof(double) * size);
            // Dx
            memcpy(source_ + 6 * size, libderiv_.ABCD[6], sizeof(double) * size);
            // Dy
            memcpy(source_ + 7 * size, libderiv_.ABCD[7], sizeof(double) * size);
            // Dz
            memcpy(source_ + 8 * size, libderiv_.ABCD[8], sizeof(double) * size);
            break;
        case BADC:
            // Ax
            C_DAXPY(size, -1.0, libderiv_.ABCD[0], 1, source_ + 0 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[6], 1, source_ + 0 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[9], 1, source_ + 0 * size, 1);
            // Ay
            C_DAXPY(size, -1.0, libderiv_.ABCD[1], 1, source_ + 1 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[7], 1, source_ + 1 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[10], 1, source_ + 1 * size, 1);
            // Az
            C_DAXPY(size, -1.0, libderiv_.ABCD[2], 1, source_ + 2 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[8], 1, source_ + 2 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[11], 1, source_ + 2 * size, 1);
            // Cx
            memcpy(source_ + 3 * size, libderiv_.ABCD[9], sizeof(double) * size);
            // Cy
            memcpy(source_ + 4 * size, libderiv_.ABCD[10], sizeof(double) * size);
            // Cz
            memcpy(source_ + 5 * size, libderiv_.ABCD[11], sizeof(double) * size);
            // Dx
            memcpy(source_ + 6 * size, libderiv_.ABCD[6], sizeof(double) * size);
            // Dy
            memcpy(source_ + 7 * size, libderiv_.ABCD[7], sizeof(double) * size);
            // Dz
            memcpy(source_ + 8 * size, libderiv_.ABCD[8], sizeof(double) * size);
            break;
        case CDAB:
            // Ax
            memcpy(source_ + 0 * size, libderiv_.ABCD[6], sizeof(double) * size);
            // Ay
            memcpy(source_ + 1 * size, libderiv_.ABCD[7], sizeof(double) * size);
            // Az
            memcpy(source_ + 2 * size, libderiv_.ABCD[8], sizeof(double) * size);
            // Cx
            memcpy(source_ + 3 * size, libderiv_.ABCD[0], sizeof(double) * size);
            // Cy
            memcpy(source_ + 4 * size, libderiv_.ABCD[1], sizeof(double) * size);
            // Cz
            memcpy(source_ + 5 * size, libderiv_.ABCD[2], sizeof(double) * size);
            // Dx
            C_DAXPY(size, -1.0, libderiv_.ABCD[0], 1, source_ + 6 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[6], 1, source_ + 6 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[9], 1, source_ + 6 * size, 1);
            // Dy
            C_DAXPY(size, -1.0, libderiv_.ABCD[1], 1, source_ + 7 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[7], 1, source_ + 7 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[10], 1, source_ + 7 * size, 1);
            // Dz
            C_DAXPY(size, -1.0, libderiv_.ABCD[2], 1, source_ + 8 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[8], 1, source_ + 8 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[11], 1, source_ + 8 * size, 1);
            break;
        case CDBA:
            // Ax
            memcpy(source_ + 0 * size, libderiv_.ABCD[9], sizeof(double) * size);
            // Ay
            memcpy(source_ + 1 * size, libderiv_.ABCD[10], sizeof(double) * size);
            // Az
            memcpy(source_ + 2 * size, libderiv_.ABCD[11], sizeof(double) * size);
            // Cx
            memcpy(source_ + 3 * size, libderiv_.ABCD[0], sizeof(double) * size);
            // Cy
            memcpy(source_ + 4 * size, libderiv_.ABCD[1], sizeof(double) * size);
            // Cz
            memcpy(source_ + 5 * size, libderiv_.ABCD[2], sizeof(double) * size);
            // Dx
            C_DAXPY(size, -1.0, libderiv_.ABCD[0], 1, source_ + 6 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[6], 1, source_ + 6 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[9], 1, source_ + 6 * size, 1);
            // Dy
            C_DAXPY(size, -1.0, libderiv_.ABCD[1], 1, source_ + 7 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[7], 1, source_ + 7 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[10], 1, source_ + 7 * size, 1);
            // Dz
            C_DAXPY(size, -1.0, libderiv_.ABCD[2], 1, source_ + 8 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[8], 1, source_ + 8 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[11], 1, source_ + 8 * size, 1);
            break;
        case DCAB:
            // Ax
            memcpy(source_ + 0 * size, libderiv_.ABCD[6], sizeof(double) * size);
            // Ay
            memcpy(source_ + 1 * size, libderiv_.ABCD[7], sizeof(double) * size);
            // Az
            memcpy(source_ + 2 * size, libderiv_.ABCD[8], sizeof(double) * size);
            // Cx
            C_DAXPY(size, -1.0, libderiv_.ABCD[0], 1, source_ + 3 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[6], 1, source_ + 3 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[9], 1, source_ + 3 * size, 1);
            // Cy
            C_DAXPY(size, -1.0, libderiv_.ABCD[1], 1, source_ + 4 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[7], 1, source_ + 4 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[10], 1, source_ + 4 * size, 1);
            // Cz
            C_DAXPY(size, -1.0, libderiv_.ABCD[2], 1, source_ + 5 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[8], 1, source_ + 5 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[11], 1, source_ + 5 * size, 1);
            // Dx
            memcpy(source_ + 6 * size, libderiv_.ABCD[0], sizeof(double) * size);
            // Dy
            memcpy(source_ + 7 * size, libderiv_.ABCD[1], sizeof(double) * size);
            // Dz
            memcpy(source_ + 8 * size, libderiv_.ABCD[2], sizeof(double) * size);
            break;
        case DCBA:
            // Ax
            memcpy(source_ + 0 * size, libderiv_.ABCD[9], sizeof(double) * size);
            // Ay
            memcpy(source_ + 1 * size, libderiv_.ABCD[10], sizeof(double) * size);
            // Az
            memcpy(source_ + 2 * size, libderiv_.ABCD[11], sizeof(double) * size);
            // Cx
            C_DAXPY(size, -1.0, libderiv_.ABCD[0], 1, source_ + 3 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[6], 1, source_ + 3 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[9], 1, source_ + 3 * size, 1);
            // Cy
            C_DAXPY(size, -1.0, libderiv_.ABCD[1], 1, source_ + 4 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[7], 1, source_ + 4 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[10], 1, source_ + 4 * size, 1);
            // Cz
            C_DAXPY(size, -1.0, libderiv_.ABCD[2], 1, source_ + 5 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[8], 1, source_ + 5 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[11], 1, source_ + 5 * size, 1);
            // Dx
            memcpy(source_ + 6 * size, libderiv_.ABCD[0], sizeof(double) * size);
            // Dy
            memcpy(source_ + 7 * size, libderiv_.ABCD[1], sizeof(double) * size);
            // Dz
            memcpy(source_ + 8 * size, libderiv_.ABCD[2], sizeof(double) * size);
            break;

        default:
            throw PSIEXCEPTION("Illegal permutation in handle_reordering code");
    }
}

/**
 * @brief Takes care of the changing the results buffer for any reordering done for libderiv.
 * @param permutation How the shells were permuted to satisfy Libint angular momentum requirements.
 * @param libderiv_ The libderiv buffer that contains the integrals
 * @param source_ The destination for the integrals.
 * @param size Total number of integrals computed.
 */
static void handle_reordering12(PermutedOrder permutation, Libderiv_t &libderiv_, double *source_, size_t size) {
    switch (permutation) {
        case ABCD:
            // Ax
            memcpy(source_ + 0 * size, libderiv_.ABCD[0], sizeof(double) * size);
            // Ay
            memcpy(source_ + 1 * size, libderiv_.ABCD[1], sizeof(double) * size);
            // Az
            memcpy(source_ + 2 * size, libderiv_.ABCD[2], sizeof(double) * size);
            // Cx
            memcpy(source_ + 3 * size, libderiv_.ABCD[6], sizeof(double) * size);
            // Cy
            memcpy(source_ + 4 * size, libderiv_.ABCD[7], sizeof(double) * size);
            // Cz
            memcpy(source_ + 5 * size, libderiv_.ABCD[8], sizeof(double) * size);
            // Dx
            memcpy(source_ + 6 * size, libderiv_.ABCD[9], sizeof(double) * size);
            // Dy
            memcpy(source_ + 7 * size, libderiv_.ABCD[10], sizeof(double) * size);
            // Dz
            memcpy(source_ + 8 * size, libderiv_.ABCD[11], sizeof(double) * size);
            // AxAx
            memcpy(source_ + 9 * size, libderiv_.ABCD[12], sizeof(double) * size);
            // AxAy
            memcpy(source_ + 10 * size, libderiv_.ABCD[13], sizeof(double) * size);
            // AxAz
            memcpy(source_ + 11 * size, libderiv_.ABCD[14], sizeof(double) * size);
            // AxCx
            memcpy(source_ + 12 * size, libderiv_.ABCD[18], sizeof(double) * size);
            // AxCy
            memcpy(source_ + 13 * size, libderiv_.ABCD[19], sizeof(double) * size);
            // AxCz
            memcpy(source_ + 14 * size, libderiv_.ABCD[20], sizeof(double) * size);
            // AxDx
            memcpy(source_ + 15 * size, libderiv_.ABCD[21], sizeof(double) * size);
            // AxDy
            memcpy(source_ + 16 * size, libderiv_.ABCD[22], sizeof(double) * size);
            // AxDz
            memcpy(source_ + 17 * size, libderiv_.ABCD[23], sizeof(double) * size);
            // AyAy
            memcpy(source_ + 18 * size, libderiv_.ABCD[25], sizeof(double) * size);
            // AyAz
            memcpy(source_ + 19 * size, libderiv_.ABCD[26], sizeof(double) * size);
            // AyCx
            memcpy(source_ + 20 * size, libderiv_.ABCD[30], sizeof(double) * size);
            // AyCy
            memcpy(source_ + 21 * size, libderiv_.ABCD[31], sizeof(double) * size);
            // AyCz
            memcpy(source_ + 22 * size, libderiv_.ABCD[32], sizeof(double) * size);
            // AyDx
            memcpy(source_ + 23 * size, libderiv_.ABCD[33], sizeof(double) * size);
            // AyDy
            memcpy(source_ + 24 * size, libderiv_.ABCD[34], sizeof(double) * size);
            // AyDz
            memcpy(source_ + 25 * size, libderiv_.ABCD[35], sizeof(double) * size);
            // AzAz
            memcpy(source_ + 26 * size, libderiv_.ABCD[38], sizeof(double) * size);
            // AzCx
            memcpy(source_ + 27 * size, libderiv_.ABCD[42], sizeof(double) * size);
            // AzCy
            memcpy(source_ + 28 * size, libderiv_.ABCD[43], sizeof(double) * size);
            // AzCz
            memcpy(source_ + 29 * size, libderiv_.ABCD[44], sizeof(double) * size);
            // AzDx
            memcpy(source_ + 30 * size, libderiv_.ABCD[45], sizeof(double) * size);
            // AzDy
            memcpy(source_ + 31 * size, libderiv_.ABCD[46], sizeof(double) * size);
            // AzDz
            memcpy(source_ + 32 * size, libderiv_.ABCD[47], sizeof(double) * size);
            // CxCx
            memcpy(source_ + 33 * size, libderiv_.ABCD[90], sizeof(double) * size);
            // CxCy
            memcpy(source_ + 34 * size, libderiv_.ABCD[91], sizeof(double) * size);
            // CxCz
            memcpy(source_ + 35 * size, libderiv_.ABCD[92], sizeof(double) * size);
            // CxDx
            memcpy(source_ + 36 * size, libderiv_.ABCD[93], sizeof(double) * size);
            // CxDy
            memcpy(source_ + 37 * size, libderiv_.ABCD[94], sizeof(double) * size);
            // CxDz
            memcpy(source_ + 38 * size, libderiv_.ABCD[95], sizeof(double) * size);
            // CyCy
            memcpy(source_ + 39 * size, libderiv_.ABCD[103], sizeof(double) * size);
            // CyCz
            memcpy(source_ + 40 * size, libderiv_.ABCD[104], sizeof(double) * size);
            // CyDx
            memcpy(source_ + 41 * size, libderiv_.ABCD[105], sizeof(double) * size);
            // CyDy
            memcpy(source_ + 42 * size, libderiv_.ABCD[106], sizeof(double) * size);
            // CyDz
            memcpy(source_ + 43 * size, libderiv_.ABCD[107], sizeof(double) * size);
            // CzCz
            memcpy(source_ + 44 * size, libderiv_.ABCD[116], sizeof(double) * size);
            // CzDx
            memcpy(source_ + 45 * size, libderiv_.ABCD[117], sizeof(double) * size);
            // CzDy
            memcpy(source_ + 46 * size, libderiv_.ABCD[118], sizeof(double) * size);
            // CzDz
            memcpy(source_ + 47 * size, libderiv_.ABCD[119], sizeof(double) * size);
            // DxDx
            memcpy(source_ + 48 * size, libderiv_.ABCD[129], sizeof(double) * size);
            // DxDy
            memcpy(source_ + 49 * size, libderiv_.ABCD[130], sizeof(double) * size);
            // DxDz
            memcpy(source_ + 50 * size, libderiv_.ABCD[131], sizeof(double) * size);
            // DyDy
            memcpy(source_ + 51 * size, libderiv_.ABCD[142], sizeof(double) * size);
            // DyDz
            memcpy(source_ + 52 * size, libderiv_.ABCD[143], sizeof(double) * size);
            // DzDz
            memcpy(source_ + 53 * size, libderiv_.ABCD[155], sizeof(double) * size);
            break;
        case BACD:
            // Ax
            C_DAXPY(size, -1.0, libderiv_.ABCD[0], 1, source_ + 0 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[6], 1, source_ + 0 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[9], 1, source_ + 0 * size, 1);
            // Ay
            C_DAXPY(size, -1.0, libderiv_.ABCD[1], 1, source_ + 1 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[7], 1, source_ + 1 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[10], 1, source_ + 1 * size, 1);
            // Az
            C_DAXPY(size, -1.0, libderiv_.ABCD[2], 1, source_ + 2 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[8], 1, source_ + 2 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[11], 1, source_ + 2 * size, 1);
            // Cx
            memcpy(source_ + 3 * size, libderiv_.ABCD[6], sizeof(double) * size);
            // Cy
            memcpy(source_ + 4 * size, libderiv_.ABCD[7], sizeof(double) * size);
            // Cz
            memcpy(source_ + 5 * size, libderiv_.ABCD[8], sizeof(double) * size);
            // Dx
            memcpy(source_ + 6 * size, libderiv_.ABCD[9], sizeof(double) * size);
            // Dy
            memcpy(source_ + 7 * size, libderiv_.ABCD[10], sizeof(double) * size);
            // Dz
            memcpy(source_ + 8 * size, libderiv_.ABCD[11], sizeof(double) * size);
            // AxAx
            C_DAXPY(size, 1.0, libderiv_.ABCD[12], 1, source_ + 9 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[18], 1, source_ + 9 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[90], 1, source_ + 9 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[21], 1, source_ + 9 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[93], 1, source_ + 9 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[129], 1, source_ + 9 * size, 1);
            // AxAy
            C_DAXPY(size, 1.0, libderiv_.ABCD[13], 1, source_ + 10 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[19], 1, source_ + 10 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[22], 1, source_ + 10 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[30], 1, source_ + 10 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[91], 1, source_ + 10 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[94], 1, source_ + 10 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[33], 1, source_ + 10 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[105], 1, source_ + 10 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[130], 1, source_ + 10 * size, 1);
            // AxAz
            C_DAXPY(size, 1.0, libderiv_.ABCD[14], 1, source_ + 11 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[20], 1, source_ + 11 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[23], 1, source_ + 11 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[42], 1, source_ + 11 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[92], 1, source_ + 11 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[95], 1, source_ + 11 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[45], 1, source_ + 11 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[117], 1, source_ + 11 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[131], 1, source_ + 11 * size, 1);
            // AxCx
            C_DAXPY(size, -1.0, libderiv_.ABCD[18], 1, source_ + 12 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[90], 1, source_ + 12 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[93], 1, source_ + 12 * size, 1);
            // AxCy
            C_DAXPY(size, -1.0, libderiv_.ABCD[19], 1, source_ + 13 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[91], 1, source_ + 13 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[105], 1, source_ + 13 * size, 1);
            // AxCz
            C_DAXPY(size, -1.0, libderiv_.ABCD[20], 1, source_ + 14 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[92], 1, source_ + 14 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[117], 1, source_ + 14 * size, 1);
            // AxDx
            C_DAXPY(size, -1.0, libderiv_.ABCD[21], 1, source_ + 15 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[93], 1, source_ + 15 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[129], 1, source_ + 15 * size, 1);
            // AxDy
            C_DAXPY(size, -1.0, libderiv_.ABCD[22], 1, source_ + 16 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[94], 1, source_ + 16 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[130], 1, source_ + 16 * size, 1);
            // AxDz
            C_DAXPY(size, -1.0, libderiv_.ABCD[23], 1, source_ + 17 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[95], 1, source_ + 17 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[131], 1, source_ + 17 * size, 1);
            // AyAy
            C_DAXPY(size, 1.0, libderiv_.ABCD[25], 1, source_ + 18 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[31], 1, source_ + 18 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[103], 1, source_ + 18 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[34], 1, source_ + 18 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[106], 1, source_ + 18 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[142], 1, source_ + 18 * size, 1);
            // AyAz
            C_DAXPY(size, 1.0, libderiv_.ABCD[26], 1, source_ + 19 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[32], 1, source_ + 19 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[35], 1, source_ + 19 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[43], 1, source_ + 19 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[104], 1, source_ + 19 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[107], 1, source_ + 19 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[46], 1, source_ + 19 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[118], 1, source_ + 19 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[143], 1, source_ + 19 * size, 1);
            // AyCx
            C_DAXPY(size, -1.0, libderiv_.ABCD[30], 1, source_ + 20 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[91], 1, source_ + 20 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[94], 1, source_ + 20 * size, 1);
            // AyCy
            C_DAXPY(size, -1.0, libderiv_.ABCD[31], 1, source_ + 21 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[103], 1, source_ + 21 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[106], 1, source_ + 21 * size, 1);
            // AyCz
            C_DAXPY(size, -1.0, libderiv_.ABCD[32], 1, source_ + 22 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[104], 1, source_ + 22 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[118], 1, source_ + 22 * size, 1);
            // AyDx
            C_DAXPY(size, -1.0, libderiv_.ABCD[33], 1, source_ + 23 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[105], 1, source_ + 23 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[130], 1, source_ + 23 * size, 1);
            // AyDy
            C_DAXPY(size, -1.0, libderiv_.ABCD[34], 1, source_ + 24 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[106], 1, source_ + 24 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[142], 1, source_ + 24 * size, 1);
            // AyDz
            C_DAXPY(size, -1.0, libderiv_.ABCD[35], 1, source_ + 25 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[107], 1, source_ + 25 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[143], 1, source_ + 25 * size, 1);
            // AzAz
            C_DAXPY(size, 1.0, libderiv_.ABCD[38], 1, source_ + 26 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[44], 1, source_ + 26 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[116], 1, source_ + 26 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[47], 1, source_ + 26 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[119], 1, source_ + 26 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[155], 1, source_ + 26 * size, 1);
            // AzCx
            C_DAXPY(size, -1.0, libderiv_.ABCD[42], 1, source_ + 27 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[92], 1, source_ + 27 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[95], 1, source_ + 27 * size, 1);
            // AzCy
            C_DAXPY(size, -1.0, libderiv_.ABCD[43], 1, source_ + 28 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[104], 1, source_ + 28 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[107], 1, source_ + 28 * size, 1);
            // AzCz
            C_DAXPY(size, -1.0, libderiv_.ABCD[44], 1, source_ + 29 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[116], 1, source_ + 29 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[119], 1, source_ + 29 * size, 1);
            // AzDx
            C_DAXPY(size, -1.0, libderiv_.ABCD[45], 1, source_ + 30 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[117], 1, source_ + 30 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[131], 1, source_ + 30 * size, 1);
            // AzDy
            C_DAXPY(size, -1.0, libderiv_.ABCD[46], 1, source_ + 31 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[118], 1, source_ + 31 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[143], 1, source_ + 31 * size, 1);
            // AzDz
            C_DAXPY(size, -1.0, libderiv_.ABCD[47], 1, source_ + 32 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[119], 1, source_ + 32 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[155], 1, source_ + 32 * size, 1);
            // CxCx
            memcpy(source_ + 33 * size, libderiv_.ABCD[90], sizeof(double) * size);
            // CxCy
            memcpy(source_ + 34 * size, libderiv_.ABCD[91], sizeof(double) * size);
            // CxCz
            memcpy(source_ + 35 * size, libderiv_.ABCD[92], sizeof(double) * size);
            // CxDx
            memcpy(source_ + 36 * size, libderiv_.ABCD[93], sizeof(double) * size);
            // CxDy
            memcpy(source_ + 37 * size, libderiv_.ABCD[94], sizeof(double) * size);
            // CxDz
            memcpy(source_ + 38 * size, libderiv_.ABCD[95], sizeof(double) * size);
            // CyCy
            memcpy(source_ + 39 * size, libderiv_.ABCD[103], sizeof(double) * size);
            // CyCz
            memcpy(source_ + 40 * size, libderiv_.ABCD[104], sizeof(double) * size);
            // CyDx
            memcpy(source_ + 41 * size, libderiv_.ABCD[105], sizeof(double) * size);
            // CyDy
            memcpy(source_ + 42 * size, libderiv_.ABCD[106], sizeof(double) * size);
            // CyDz
            memcpy(source_ + 43 * size, libderiv_.ABCD[107], sizeof(double) * size);
            // CzCz
            memcpy(source_ + 44 * size, libderiv_.ABCD[116], sizeof(double) * size);
            // CzDx
            memcpy(source_ + 45 * size, libderiv_.ABCD[117], sizeof(double) * size);
            // CzDy
            memcpy(source_ + 46 * size, libderiv_.ABCD[118], sizeof(double) * size);
            // CzDz
            memcpy(source_ + 47 * size, libderiv_.ABCD[119], sizeof(double) * size);
            // DxDx
            memcpy(source_ + 48 * size, libderiv_.ABCD[129], sizeof(double) * size);
            // DxDy
            memcpy(source_ + 49 * size, libderiv_.ABCD[130], sizeof(double) * size);
            // DxDz
            memcpy(source_ + 50 * size, libderiv_.ABCD[131], sizeof(double) * size);
            // DyDy
            memcpy(source_ + 51 * size, libderiv_.ABCD[142], sizeof(double) * size);
            // DyDz
            memcpy(source_ + 52 * size, libderiv_.ABCD[143], sizeof(double) * size);
            // DzDz
            memcpy(source_ + 53 * size, libderiv_.ABCD[155], sizeof(double) * size);
            break;
        case ABDC:
            // Ax
            memcpy(source_ + 0 * size, libderiv_.ABCD[0], sizeof(double) * size);
            // Ay
            memcpy(source_ + 1 * size, libderiv_.ABCD[1], sizeof(double) * size);
            // Az
            memcpy(source_ + 2 * size, libderiv_.ABCD[2], sizeof(double) * size);
            // Cx
            memcpy(source_ + 3 * size, libderiv_.ABCD[9], sizeof(double) * size);
            // Cy
            memcpy(source_ + 4 * size, libderiv_.ABCD[10], sizeof(double) * size);
            // Cz
            memcpy(source_ + 5 * size, libderiv_.ABCD[11], sizeof(double) * size);
            // Dx
            memcpy(source_ + 6 * size, libderiv_.ABCD[6], sizeof(double) * size);
            // Dy
            memcpy(source_ + 7 * size, libderiv_.ABCD[7], sizeof(double) * size);
            // Dz
            memcpy(source_ + 8 * size, libderiv_.ABCD[8], sizeof(double) * size);
            // AxAx
            memcpy(source_ + 9 * size, libderiv_.ABCD[12], sizeof(double) * size);
            // AxAy
            memcpy(source_ + 10 * size, libderiv_.ABCD[13], sizeof(double) * size);
            // AxAz
            memcpy(source_ + 11 * size, libderiv_.ABCD[14], sizeof(double) * size);
            // AxCx
            memcpy(source_ + 12 * size, libderiv_.ABCD[21], sizeof(double) * size);
            // AxCy
            memcpy(source_ + 13 * size, libderiv_.ABCD[22], sizeof(double) * size);
            // AxCz
            memcpy(source_ + 14 * size, libderiv_.ABCD[23], sizeof(double) * size);
            // AxDx
            memcpy(source_ + 15 * size, libderiv_.ABCD[18], sizeof(double) * size);
            // AxDy
            memcpy(source_ + 16 * size, libderiv_.ABCD[19], sizeof(double) * size);
            // AxDz
            memcpy(source_ + 17 * size, libderiv_.ABCD[20], sizeof(double) * size);
            // AyAy
            memcpy(source_ + 18 * size, libderiv_.ABCD[25], sizeof(double) * size);
            // AyAz
            memcpy(source_ + 19 * size, libderiv_.ABCD[26], sizeof(double) * size);
            // AyCx
            memcpy(source_ + 20 * size, libderiv_.ABCD[33], sizeof(double) * size);
            // AyCy
            memcpy(source_ + 21 * size, libderiv_.ABCD[34], sizeof(double) * size);
            // AyCz
            memcpy(source_ + 22 * size, libderiv_.ABCD[35], sizeof(double) * size);
            // AyDx
            memcpy(source_ + 23 * size, libderiv_.ABCD[30], sizeof(double) * size);
            // AyDy
            memcpy(source_ + 24 * size, libderiv_.ABCD[31], sizeof(double) * size);
            // AyDz
            memcpy(source_ + 25 * size, libderiv_.ABCD[32], sizeof(double) * size);
            // AzAz
            memcpy(source_ + 26 * size, libderiv_.ABCD[38], sizeof(double) * size);
            // AzCx
            memcpy(source_ + 27 * size, libderiv_.ABCD[45], sizeof(double) * size);
            // AzCy
            memcpy(source_ + 28 * size, libderiv_.ABCD[46], sizeof(double) * size);
            // AzCz
            memcpy(source_ + 29 * size, libderiv_.ABCD[47], sizeof(double) * size);
            // AzDx
            memcpy(source_ + 30 * size, libderiv_.ABCD[42], sizeof(double) * size);
            // AzDy
            memcpy(source_ + 31 * size, libderiv_.ABCD[43], sizeof(double) * size);
            // AzDz
            memcpy(source_ + 32 * size, libderiv_.ABCD[44], sizeof(double) * size);
            // CxCx
            memcpy(source_ + 33 * size, libderiv_.ABCD[129], sizeof(double) * size);
            // CxCy
            memcpy(source_ + 34 * size, libderiv_.ABCD[130], sizeof(double) * size);
            // CxCz
            memcpy(source_ + 35 * size, libderiv_.ABCD[131], sizeof(double) * size);
            // CxDx
            memcpy(source_ + 36 * size, libderiv_.ABCD[93], sizeof(double) * size);
            // CxDy
            memcpy(source_ + 37 * size, libderiv_.ABCD[105], sizeof(double) * size);
            // CxDz
            memcpy(source_ + 38 * size, libderiv_.ABCD[117], sizeof(double) * size);
            // CyCy
            memcpy(source_ + 39 * size, libderiv_.ABCD[142], sizeof(double) * size);
            // CyCz
            memcpy(source_ + 40 * size, libderiv_.ABCD[143], sizeof(double) * size);
            // CyDx
            memcpy(source_ + 41 * size, libderiv_.ABCD[94], sizeof(double) * size);
            // CyDy
            memcpy(source_ + 42 * size, libderiv_.ABCD[106], sizeof(double) * size);
            // CyDz
            memcpy(source_ + 43 * size, libderiv_.ABCD[118], sizeof(double) * size);
            // CzCz
            memcpy(source_ + 44 * size, libderiv_.ABCD[155], sizeof(double) * size);
            // CzDx
            memcpy(source_ + 45 * size, libderiv_.ABCD[95], sizeof(double) * size);
            // CzDy
            memcpy(source_ + 46 * size, libderiv_.ABCD[107], sizeof(double) * size);
            // CzDz
            memcpy(source_ + 47 * size, libderiv_.ABCD[119], sizeof(double) * size);
            // DxDx
            memcpy(source_ + 48 * size, libderiv_.ABCD[90], sizeof(double) * size);
            // DxDy
            memcpy(source_ + 49 * size, libderiv_.ABCD[91], sizeof(double) * size);
            // DxDz
            memcpy(source_ + 50 * size, libderiv_.ABCD[92], sizeof(double) * size);
            // DyDy
            memcpy(source_ + 51 * size, libderiv_.ABCD[103], sizeof(double) * size);
            // DyDz
            memcpy(source_ + 52 * size, libderiv_.ABCD[104], sizeof(double) * size);
            // DzDz
            memcpy(source_ + 53 * size, libderiv_.ABCD[116], sizeof(double) * size);
            break;
        case BADC:
            // Ax
            C_DAXPY(size, -1.0, libderiv_.ABCD[0], 1, source_ + 0 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[6], 1, source_ + 0 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[9], 1, source_ + 0 * size, 1);
            // Ay
            C_DAXPY(size, -1.0, libderiv_.ABCD[1], 1, source_ + 1 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[7], 1, source_ + 1 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[10], 1, source_ + 1 * size, 1);
            // Az
            C_DAXPY(size, -1.0, libderiv_.ABCD[2], 1, source_ + 2 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[8], 1, source_ + 2 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[11], 1, source_ + 2 * size, 1);
            // Cx
            memcpy(source_ + 3 * size, libderiv_.ABCD[9], sizeof(double) * size);
            // Cy
            memcpy(source_ + 4 * size, libderiv_.ABCD[10], sizeof(double) * size);
            // Cz
            memcpy(source_ + 5 * size, libderiv_.ABCD[11], sizeof(double) * size);
            // Dx
            memcpy(source_ + 6 * size, libderiv_.ABCD[6], sizeof(double) * size);
            // Dy
            memcpy(source_ + 7 * size, libderiv_.ABCD[7], sizeof(double) * size);
            // Dz
            memcpy(source_ + 8 * size, libderiv_.ABCD[8], sizeof(double) * size);
            // AxAx
            C_DAXPY(size, 1.0, libderiv_.ABCD[12], 1, source_ + 9 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[18], 1, source_ + 9 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[90], 1, source_ + 9 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[21], 1, source_ + 9 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[93], 1, source_ + 9 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[129], 1, source_ + 9 * size, 1);
            // AxAy
            C_DAXPY(size, 1.0, libderiv_.ABCD[13], 1, source_ + 10 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[19], 1, source_ + 10 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[22], 1, source_ + 10 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[30], 1, source_ + 10 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[91], 1, source_ + 10 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[94], 1, source_ + 10 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[33], 1, source_ + 10 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[105], 1, source_ + 10 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[130], 1, source_ + 10 * size, 1);
            // AxAz
            C_DAXPY(size, 1.0, libderiv_.ABCD[14], 1, source_ + 11 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[20], 1, source_ + 11 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[23], 1, source_ + 11 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[42], 1, source_ + 11 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[92], 1, source_ + 11 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[95], 1, source_ + 11 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[45], 1, source_ + 11 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[117], 1, source_ + 11 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[131], 1, source_ + 11 * size, 1);
            // AxCx
            C_DAXPY(size, -1.0, libderiv_.ABCD[21], 1, source_ + 12 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[93], 1, source_ + 12 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[129], 1, source_ + 12 * size, 1);
            // AxCy
            C_DAXPY(size, -1.0, libderiv_.ABCD[22], 1, source_ + 13 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[94], 1, source_ + 13 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[130], 1, source_ + 13 * size, 1);
            // AxCz
            C_DAXPY(size, -1.0, libderiv_.ABCD[23], 1, source_ + 14 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[95], 1, source_ + 14 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[131], 1, source_ + 14 * size, 1);
            // AxDx
            C_DAXPY(size, -1.0, libderiv_.ABCD[18], 1, source_ + 15 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[90], 1, source_ + 15 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[93], 1, source_ + 15 * size, 1);
            // AxDy
            C_DAXPY(size, -1.0, libderiv_.ABCD[19], 1, source_ + 16 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[91], 1, source_ + 16 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[105], 1, source_ + 16 * size, 1);
            // AxDz
            C_DAXPY(size, -1.0, libderiv_.ABCD[20], 1, source_ + 17 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[92], 1, source_ + 17 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[117], 1, source_ + 17 * size, 1);
            // AyAy
            C_DAXPY(size, 1.0, libderiv_.ABCD[25], 1, source_ + 18 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[31], 1, source_ + 18 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[103], 1, source_ + 18 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[34], 1, source_ + 18 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[106], 1, source_ + 18 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[142], 1, source_ + 18 * size, 1);
            // AyAz
            C_DAXPY(size, 1.0, libderiv_.ABCD[26], 1, source_ + 19 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[32], 1, source_ + 19 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[35], 1, source_ + 19 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[43], 1, source_ + 19 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[104], 1, source_ + 19 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[107], 1, source_ + 19 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[46], 1, source_ + 19 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[118], 1, source_ + 19 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[143], 1, source_ + 19 * size, 1);
            // AyCx
            C_DAXPY(size, -1.0, libderiv_.ABCD[33], 1, source_ + 20 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[105], 1, source_ + 20 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[130], 1, source_ + 20 * size, 1);
            // AyCy
            C_DAXPY(size, -1.0, libderiv_.ABCD[34], 1, source_ + 21 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[106], 1, source_ + 21 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[142], 1, source_ + 21 * size, 1);
            // AyCz
            C_DAXPY(size, -1.0, libderiv_.ABCD[35], 1, source_ + 22 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[107], 1, source_ + 22 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[143], 1, source_ + 22 * size, 1);
            // AyDx
            C_DAXPY(size, -1.0, libderiv_.ABCD[30], 1, source_ + 23 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[91], 1, source_ + 23 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[94], 1, source_ + 23 * size, 1);
            // AyDy
            C_DAXPY(size, -1.0, libderiv_.ABCD[31], 1, source_ + 24 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[103], 1, source_ + 24 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[106], 1, source_ + 24 * size, 1);
            // AyDz
            C_DAXPY(size, -1.0, libderiv_.ABCD[32], 1, source_ + 25 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[104], 1, source_ + 25 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[118], 1, source_ + 25 * size, 1);
            // AzAz
            C_DAXPY(size, 1.0, libderiv_.ABCD[38], 1, source_ + 26 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[44], 1, source_ + 26 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[116], 1, source_ + 26 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[47], 1, source_ + 26 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[119], 1, source_ + 26 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[155], 1, source_ + 26 * size, 1);
            // AzCx
            C_DAXPY(size, -1.0, libderiv_.ABCD[45], 1, source_ + 27 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[117], 1, source_ + 27 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[131], 1, source_ + 27 * size, 1);
            // AzCy
            C_DAXPY(size, -1.0, libderiv_.ABCD[46], 1, source_ + 28 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[118], 1, source_ + 28 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[143], 1, source_ + 28 * size, 1);
            // AzCz
            C_DAXPY(size, -1.0, libderiv_.ABCD[47], 1, source_ + 29 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[119], 1, source_ + 29 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[155], 1, source_ + 29 * size, 1);
            // AzDx
            C_DAXPY(size, -1.0, libderiv_.ABCD[42], 1, source_ + 30 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[92], 1, source_ + 30 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[95], 1, source_ + 30 * size, 1);
            // AzDy
            C_DAXPY(size, -1.0, libderiv_.ABCD[43], 1, source_ + 31 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[104], 1, source_ + 31 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[107], 1, source_ + 31 * size, 1);
            // AzDz
            C_DAXPY(size, -1.0, libderiv_.ABCD[44], 1, source_ + 32 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[116], 1, source_ + 32 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[119], 1, source_ + 32 * size, 1);
            // CxCx
            memcpy(source_ + 33 * size, libderiv_.ABCD[129], sizeof(double) * size);
            // CxCy
            memcpy(source_ + 34 * size, libderiv_.ABCD[130], sizeof(double) * size);
            // CxCz
            memcpy(source_ + 35 * size, libderiv_.ABCD[131], sizeof(double) * size);
            // CxDx
            memcpy(source_ + 36 * size, libderiv_.ABCD[93], sizeof(double) * size);
            // CxDy
            memcpy(source_ + 37 * size, libderiv_.ABCD[105], sizeof(double) * size);
            // CxDz
            memcpy(source_ + 38 * size, libderiv_.ABCD[117], sizeof(double) * size);
            // CyCy
            memcpy(source_ + 39 * size, libderiv_.ABCD[142], sizeof(double) * size);
            // CyCz
            memcpy(source_ + 40 * size, libderiv_.ABCD[143], sizeof(double) * size);
            // CyDx
            memcpy(source_ + 41 * size, libderiv_.ABCD[94], sizeof(double) * size);
            // CyDy
            memcpy(source_ + 42 * size, libderiv_.ABCD[106], sizeof(double) * size);
            // CyDz
            memcpy(source_ + 43 * size, libderiv_.ABCD[118], sizeof(double) * size);
            // CzCz
            memcpy(source_ + 44 * size, libderiv_.ABCD[155], sizeof(double) * size);
            // CzDx
            memcpy(source_ + 45 * size, libderiv_.ABCD[95], sizeof(double) * size);
            // CzDy
            memcpy(source_ + 46 * size, libderiv_.ABCD[107], sizeof(double) * size);
            // CzDz
            memcpy(source_ + 47 * size, libderiv_.ABCD[119], sizeof(double) * size);
            // DxDx
            memcpy(source_ + 48 * size, libderiv_.ABCD[90], sizeof(double) * size);
            // DxDy
            memcpy(source_ + 49 * size, libderiv_.ABCD[91], sizeof(double) * size);
            // DxDz
            memcpy(source_ + 50 * size, libderiv_.ABCD[92], sizeof(double) * size);
            // DyDy
            memcpy(source_ + 51 * size, libderiv_.ABCD[103], sizeof(double) * size);
            // DyDz
            memcpy(source_ + 52 * size, libderiv_.ABCD[104], sizeof(double) * size);
            // DzDz
            memcpy(source_ + 53 * size, libderiv_.ABCD[116], sizeof(double) * size);
            break;
        case CDAB:
            // Ax
            memcpy(source_ + 0 * size, libderiv_.ABCD[6], sizeof(double) * size);
            // Ay
            memcpy(source_ + 1 * size, libderiv_.ABCD[7], sizeof(double) * size);
            // Az
            memcpy(source_ + 2 * size, libderiv_.ABCD[8], sizeof(double) * size);
            // Cx
            memcpy(source_ + 3 * size, libderiv_.ABCD[0], sizeof(double) * size);
            // Cy
            memcpy(source_ + 4 * size, libderiv_.ABCD[1], sizeof(double) * size);
            // Cz
            memcpy(source_ + 5 * size, libderiv_.ABCD[2], sizeof(double) * size);
            // Dx
            C_DAXPY(size, -1.0, libderiv_.ABCD[0], 1, source_ + 6 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[6], 1, source_ + 6 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[9], 1, source_ + 6 * size, 1);
            // Dy
            C_DAXPY(size, -1.0, libderiv_.ABCD[1], 1, source_ + 7 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[7], 1, source_ + 7 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[10], 1, source_ + 7 * size, 1);
            // Dz
            C_DAXPY(size, -1.0, libderiv_.ABCD[2], 1, source_ + 8 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[8], 1, source_ + 8 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[11], 1, source_ + 8 * size, 1);
            // AxAx
            memcpy(source_ + 9 * size, libderiv_.ABCD[90], sizeof(double) * size);
            // AxAy
            memcpy(source_ + 10 * size, libderiv_.ABCD[91], sizeof(double) * size);
            // AxAz
            memcpy(source_ + 11 * size, libderiv_.ABCD[92], sizeof(double) * size);
            // AxCx
            memcpy(source_ + 12 * size, libderiv_.ABCD[18], sizeof(double) * size);
            // AxCy
            memcpy(source_ + 13 * size, libderiv_.ABCD[30], sizeof(double) * size);
            // AxCz
            memcpy(source_ + 14 * size, libderiv_.ABCD[42], sizeof(double) * size);
            // AxDx
            C_DAXPY(size, -1.0, libderiv_.ABCD[18], 1, source_ + 15 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[90], 1, source_ + 15 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[93], 1, source_ + 15 * size, 1);
            // AxDy
            C_DAXPY(size, -1.0, libderiv_.ABCD[30], 1, source_ + 16 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[91], 1, source_ + 16 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[94], 1, source_ + 16 * size, 1);
            // AxDz
            C_DAXPY(size, -1.0, libderiv_.ABCD[42], 1, source_ + 17 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[92], 1, source_ + 17 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[95], 1, source_ + 17 * size, 1);
            // AyAy
            memcpy(source_ + 18 * size, libderiv_.ABCD[103], sizeof(double) * size);
            // AyAz
            memcpy(source_ + 19 * size, libderiv_.ABCD[104], sizeof(double) * size);
            // AyCx
            memcpy(source_ + 20 * size, libderiv_.ABCD[19], sizeof(double) * size);
            // AyCy
            memcpy(source_ + 21 * size, libderiv_.ABCD[31], sizeof(double) * size);
            // AyCz
            memcpy(source_ + 22 * size, libderiv_.ABCD[43], sizeof(double) * size);
            // AyDx
            C_DAXPY(size, -1.0, libderiv_.ABCD[19], 1, source_ + 23 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[91], 1, source_ + 23 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[105], 1, source_ + 23 * size, 1);
            // AyDy
            C_DAXPY(size, -1.0, libderiv_.ABCD[31], 1, source_ + 24 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[103], 1, source_ + 24 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[106], 1, source_ + 24 * size, 1);
            // AyDz
            C_DAXPY(size, -1.0, libderiv_.ABCD[43], 1, source_ + 25 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[104], 1, source_ + 25 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[107], 1, source_ + 25 * size, 1);
            // AzAz
            memcpy(source_ + 26 * size, libderiv_.ABCD[116], sizeof(double) * size);
            // AzCx
            memcpy(source_ + 27 * size, libderiv_.ABCD[20], sizeof(double) * size);
            // AzCy
            memcpy(source_ + 28 * size, libderiv_.ABCD[32], sizeof(double) * size);
            // AzCz
            memcpy(source_ + 29 * size, libderiv_.ABCD[44], sizeof(double) * size);
            // AzDx
            C_DAXPY(size, -1.0, libderiv_.ABCD[20], 1, source_ + 30 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[92], 1, source_ + 30 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[117], 1, source_ + 30 * size, 1);
            // AzDy
            C_DAXPY(size, -1.0, libderiv_.ABCD[32], 1, source_ + 31 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[104], 1, source_ + 31 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[118], 1, source_ + 31 * size, 1);
            // AzDz
            C_DAXPY(size, -1.0, libderiv_.ABCD[44], 1, source_ + 32 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[116], 1, source_ + 32 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[119], 1, source_ + 32 * size, 1);
            // CxCx
            memcpy(source_ + 33 * size, libderiv_.ABCD[12], sizeof(double) * size);
            // CxCy
            memcpy(source_ + 34 * size, libderiv_.ABCD[13], sizeof(double) * size);
            // CxCz
            memcpy(source_ + 35 * size, libderiv_.ABCD[14], sizeof(double) * size);
            // CxDx
            C_DAXPY(size, -1.0, libderiv_.ABCD[12], 1, source_ + 36 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[18], 1, source_ + 36 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[21], 1, source_ + 36 * size, 1);
            // CxDy
            C_DAXPY(size, -1.0, libderiv_.ABCD[13], 1, source_ + 37 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[19], 1, source_ + 37 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[22], 1, source_ + 37 * size, 1);
            // CxDz
            C_DAXPY(size, -1.0, libderiv_.ABCD[14], 1, source_ + 38 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[20], 1, source_ + 38 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[23], 1, source_ + 38 * size, 1);
            // CyCy
            memcpy(source_ + 39 * size, libderiv_.ABCD[25], sizeof(double) * size);
            // CyCz
            memcpy(source_ + 40 * size, libderiv_.ABCD[26], sizeof(double) * size);
            // CyDx
            C_DAXPY(size, -1.0, libderiv_.ABCD[13], 1, source_ + 41 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[30], 1, source_ + 41 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[33], 1, source_ + 41 * size, 1);
            // CyDy
            C_DAXPY(size, -1.0, libderiv_.ABCD[25], 1, source_ + 42 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[31], 1, source_ + 42 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[34], 1, source_ + 42 * size, 1);
            // CyDz
            C_DAXPY(size, -1.0, libderiv_.ABCD[26], 1, source_ + 43 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[32], 1, source_ + 43 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[35], 1, source_ + 43 * size, 1);
            // CzCz
            memcpy(source_ + 44 * size, libderiv_.ABCD[38], sizeof(double) * size);
            // CzDx
            C_DAXPY(size, -1.0, libderiv_.ABCD[14], 1, source_ + 45 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[42], 1, source_ + 45 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[45], 1, source_ + 45 * size, 1);
            // CzDy
            C_DAXPY(size, -1.0, libderiv_.ABCD[26], 1, source_ + 46 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[43], 1, source_ + 46 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[46], 1, source_ + 46 * size, 1);
            // CzDz
            C_DAXPY(size, -1.0, libderiv_.ABCD[38], 1, source_ + 47 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[44], 1, source_ + 47 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[47], 1, source_ + 47 * size, 1);
            // DxDx
            C_DAXPY(size, 1.0, libderiv_.ABCD[12], 1, source_ + 48 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[18], 1, source_ + 48 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[90], 1, source_ + 48 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[21], 1, source_ + 48 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[93], 1, source_ + 48 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[129], 1, source_ + 48 * size, 1);
            // DxDy
            C_DAXPY(size, 1.0, libderiv_.ABCD[13], 1, source_ + 49 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[19], 1, source_ + 49 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[22], 1, source_ + 49 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[30], 1, source_ + 49 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[91], 1, source_ + 49 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[94], 1, source_ + 49 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[33], 1, source_ + 49 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[105], 1, source_ + 49 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[130], 1, source_ + 49 * size, 1);
            // DxDz
            C_DAXPY(size, 1.0, libderiv_.ABCD[14], 1, source_ + 50 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[20], 1, source_ + 50 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[23], 1, source_ + 50 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[42], 1, source_ + 50 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[92], 1, source_ + 50 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[95], 1, source_ + 50 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[45], 1, source_ + 50 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[117], 1, source_ + 50 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[131], 1, source_ + 50 * size, 1);
            // DyDy
            C_DAXPY(size, 1.0, libderiv_.ABCD[25], 1, source_ + 51 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[31], 1, source_ + 51 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[103], 1, source_ + 51 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[34], 1, source_ + 51 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[106], 1, source_ + 51 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[142], 1, source_ + 51 * size, 1);
            // DyDz
            C_DAXPY(size, 1.0, libderiv_.ABCD[26], 1, source_ + 52 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[32], 1, source_ + 52 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[35], 1, source_ + 52 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[43], 1, source_ + 52 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[104], 1, source_ + 52 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[107], 1, source_ + 52 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[46], 1, source_ + 52 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[118], 1, source_ + 52 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[143], 1, source_ + 52 * size, 1);
            // DzDz
            C_DAXPY(size, 1.0, libderiv_.ABCD[38], 1, source_ + 53 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[44], 1, source_ + 53 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[116], 1, source_ + 53 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[47], 1, source_ + 53 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[119], 1, source_ + 53 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[155], 1, source_ + 53 * size, 1);
            break;
        case CDBA:
            // Ax
            memcpy(source_ + 0 * size, libderiv_.ABCD[9], sizeof(double) * size);
            // Ay
            memcpy(source_ + 1 * size, libderiv_.ABCD[10], sizeof(double) * size);
            // Az
            memcpy(source_ + 2 * size, libderiv_.ABCD[11], sizeof(double) * size);
            // Cx
            memcpy(source_ + 3 * size, libderiv_.ABCD[0], sizeof(double) * size);
            // Cy
            memcpy(source_ + 4 * size, libderiv_.ABCD[1], sizeof(double) * size);
            // Cz
            memcpy(source_ + 5 * size, libderiv_.ABCD[2], sizeof(double) * size);
            // Dx
            C_DAXPY(size, -1.0, libderiv_.ABCD[0], 1, source_ + 6 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[6], 1, source_ + 6 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[9], 1, source_ + 6 * size, 1);
            // Dy
            C_DAXPY(size, -1.0, libderiv_.ABCD[1], 1, source_ + 7 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[7], 1, source_ + 7 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[10], 1, source_ + 7 * size, 1);
            // Dz
            C_DAXPY(size, -1.0, libderiv_.ABCD[2], 1, source_ + 8 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[8], 1, source_ + 8 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[11], 1, source_ + 8 * size, 1);
            // AxAx
            memcpy(source_ + 9 * size, libderiv_.ABCD[129], sizeof(double) * size);
            // AxAy
            memcpy(source_ + 10 * size, libderiv_.ABCD[130], sizeof(double) * size);
            // AxAz
            memcpy(source_ + 11 * size, libderiv_.ABCD[131], sizeof(double) * size);
            // AxCx
            memcpy(source_ + 12 * size, libderiv_.ABCD[21], sizeof(double) * size);
            // AxCy
            memcpy(source_ + 13 * size, libderiv_.ABCD[33], sizeof(double) * size);
            // AxCz
            memcpy(source_ + 14 * size, libderiv_.ABCD[45], sizeof(double) * size);
            // AxDx
            C_DAXPY(size, -1.0, libderiv_.ABCD[21], 1, source_ + 15 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[93], 1, source_ + 15 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[129], 1, source_ + 15 * size, 1);
            // AxDy
            C_DAXPY(size, -1.0, libderiv_.ABCD[33], 1, source_ + 16 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[105], 1, source_ + 16 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[130], 1, source_ + 16 * size, 1);
            // AxDz
            C_DAXPY(size, -1.0, libderiv_.ABCD[45], 1, source_ + 17 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[117], 1, source_ + 17 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[131], 1, source_ + 17 * size, 1);
            // AyAy
            memcpy(source_ + 18 * size, libderiv_.ABCD[142], sizeof(double) * size);
            // AyAz
            memcpy(source_ + 19 * size, libderiv_.ABCD[143], sizeof(double) * size);
            // AyCx
            memcpy(source_ + 20 * size, libderiv_.ABCD[22], sizeof(double) * size);
            // AyCy
            memcpy(source_ + 21 * size, libderiv_.ABCD[34], sizeof(double) * size);
            // AyCz
            memcpy(source_ + 22 * size, libderiv_.ABCD[46], sizeof(double) * size);
            // AyDx
            C_DAXPY(size, -1.0, libderiv_.ABCD[22], 1, source_ + 23 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[94], 1, source_ + 23 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[130], 1, source_ + 23 * size, 1);
            // AyDy
            C_DAXPY(size, -1.0, libderiv_.ABCD[34], 1, source_ + 24 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[106], 1, source_ + 24 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[142], 1, source_ + 24 * size, 1);
            // AyDz
            C_DAXPY(size, -1.0, libderiv_.ABCD[46], 1, source_ + 25 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[118], 1, source_ + 25 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[143], 1, source_ + 25 * size, 1);
            // AzAz
            memcpy(source_ + 26 * size, libderiv_.ABCD[155], sizeof(double) * size);
            // AzCx
            memcpy(source_ + 27 * size, libderiv_.ABCD[23], sizeof(double) * size);
            // AzCy
            memcpy(source_ + 28 * size, libderiv_.ABCD[35], sizeof(double) * size);
            // AzCz
            memcpy(source_ + 29 * size, libderiv_.ABCD[47], sizeof(double) * size);
            // AzDx
            C_DAXPY(size, -1.0, libderiv_.ABCD[23], 1, source_ + 30 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[95], 1, source_ + 30 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[131], 1, source_ + 30 * size, 1);
            // AzDy
            C_DAXPY(size, -1.0, libderiv_.ABCD[35], 1, source_ + 31 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[107], 1, source_ + 31 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[143], 1, source_ + 31 * size, 1);
            // AzDz
            C_DAXPY(size, -1.0, libderiv_.ABCD[47], 1, source_ + 32 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[119], 1, source_ + 32 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[155], 1, source_ + 32 * size, 1);
            // CxCx
            memcpy(source_ + 33 * size, libderiv_.ABCD[12], sizeof(double) * size);
            // CxCy
            memcpy(source_ + 34 * size, libderiv_.ABCD[13], sizeof(double) * size);
            // CxCz
            memcpy(source_ + 35 * size, libderiv_.ABCD[14], sizeof(double) * size);
            // CxDx
            C_DAXPY(size, -1.0, libderiv_.ABCD[12], 1, source_ + 36 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[18], 1, source_ + 36 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[21], 1, source_ + 36 * size, 1);
            // CxDy
            C_DAXPY(size, -1.0, libderiv_.ABCD[13], 1, source_ + 37 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[19], 1, source_ + 37 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[22], 1, source_ + 37 * size, 1);
            // CxDz
            C_DAXPY(size, -1.0, libderiv_.ABCD[14], 1, source_ + 38 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[20], 1, source_ + 38 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[23], 1, source_ + 38 * size, 1);
            // CyCy
            memcpy(source_ + 39 * size, libderiv_.ABCD[25], sizeof(double) * size);
            // CyCz
            memcpy(source_ + 40 * size, libderiv_.ABCD[26], sizeof(double) * size);
            // CyDx
            C_DAXPY(size, -1.0, libderiv_.ABCD[13], 1, source_ + 41 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[30], 1, source_ + 41 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[33], 1, source_ + 41 * size, 1);
            // CyDy
            C_DAXPY(size, -1.0, libderiv_.ABCD[25], 1, source_ + 42 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[31], 1, source_ + 42 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[34], 1, source_ + 42 * size, 1);
            // CyDz
            C_DAXPY(size, -1.0, libderiv_.ABCD[26], 1, source_ + 43 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[32], 1, source_ + 43 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[35], 1, source_ + 43 * size, 1);
            // CzCz
            memcpy(source_ + 44 * size, libderiv_.ABCD[38], sizeof(double) * size);
            // CzDx
            C_DAXPY(size, -1.0, libderiv_.ABCD[14], 1, source_ + 45 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[42], 1, source_ + 45 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[45], 1, source_ + 45 * size, 1);
            // CzDy
            C_DAXPY(size, -1.0, libderiv_.ABCD[26], 1, source_ + 46 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[43], 1, source_ + 46 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[46], 1, source_ + 46 * size, 1);
            // CzDz
            C_DAXPY(size, -1.0, libderiv_.ABCD[38], 1, source_ + 47 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[44], 1, source_ + 47 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[47], 1, source_ + 47 * size, 1);
            // DxDx
            C_DAXPY(size, 1.0, libderiv_.ABCD[12], 1, source_ + 48 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[18], 1, source_ + 48 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[90], 1, source_ + 48 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[21], 1, source_ + 48 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[93], 1, source_ + 48 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[129], 1, source_ + 48 * size, 1);
            // DxDy
            C_DAXPY(size, 1.0, libderiv_.ABCD[13], 1, source_ + 49 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[19], 1, source_ + 49 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[22], 1, source_ + 49 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[30], 1, source_ + 49 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[91], 1, source_ + 49 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[94], 1, source_ + 49 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[33], 1, source_ + 49 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[105], 1, source_ + 49 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[130], 1, source_ + 49 * size, 1);
            // DxDz
            C_DAXPY(size, 1.0, libderiv_.ABCD[14], 1, source_ + 50 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[20], 1, source_ + 50 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[23], 1, source_ + 50 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[42], 1, source_ + 50 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[92], 1, source_ + 50 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[95], 1, source_ + 50 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[45], 1, source_ + 50 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[117], 1, source_ + 50 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[131], 1, source_ + 50 * size, 1);
            // DyDy
            C_DAXPY(size, 1.0, libderiv_.ABCD[25], 1, source_ + 51 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[31], 1, source_ + 51 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[103], 1, source_ + 51 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[34], 1, source_ + 51 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[106], 1, source_ + 51 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[142], 1, source_ + 51 * size, 1);
            // DyDz
            C_DAXPY(size, 1.0, libderiv_.ABCD[26], 1, source_ + 52 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[32], 1, source_ + 52 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[35], 1, source_ + 52 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[43], 1, source_ + 52 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[104], 1, source_ + 52 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[107], 1, source_ + 52 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[46], 1, source_ + 52 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[118], 1, source_ + 52 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[143], 1, source_ + 52 * size, 1);
            // DzDz
            C_DAXPY(size, 1.0, libderiv_.ABCD[38], 1, source_ + 53 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[44], 1, source_ + 53 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[116], 1, source_ + 53 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[47], 1, source_ + 53 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[119], 1, source_ + 53 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[155], 1, source_ + 53 * size, 1);
            break;
        case DCAB:
            // Ax
            memcpy(source_ + 0 * size, libderiv_.ABCD[6], sizeof(double) * size);
            // Ay
            memcpy(source_ + 1 * size, libderiv_.ABCD[7], sizeof(double) * size);
            // Az
            memcpy(source_ + 2 * size, libderiv_.ABCD[8], sizeof(double) * size);
            // Cx
            C_DAXPY(size, -1.0, libderiv_.ABCD[0], 1, source_ + 3 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[6], 1, source_ + 3 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[9], 1, source_ + 3 * size, 1);
            // Cy
            C_DAXPY(size, -1.0, libderiv_.ABCD[1], 1, source_ + 4 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[7], 1, source_ + 4 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[10], 1, source_ + 4 * size, 1);
            // Cz
            C_DAXPY(size, -1.0, libderiv_.ABCD[2], 1, source_ + 5 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[8], 1, source_ + 5 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[11], 1, source_ + 5 * size, 1);
            // Dx
            memcpy(source_ + 6 * size, libderiv_.ABCD[0], sizeof(double) * size);
            // Dy
            memcpy(source_ + 7 * size, libderiv_.ABCD[1], sizeof(double) * size);
            // Dz
            memcpy(source_ + 8 * size, libderiv_.ABCD[2], sizeof(double) * size);
            // AxAx
            memcpy(source_ + 9 * size, libderiv_.ABCD[90], sizeof(double) * size);
            // AxAy
            memcpy(source_ + 10 * size, libderiv_.ABCD[91], sizeof(double) * size);
            // AxAz
            memcpy(source_ + 11 * size, libderiv_.ABCD[92], sizeof(double) * size);
            // AxCx
            C_DAXPY(size, -1.0, libderiv_.ABCD[18], 1, source_ + 12 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[90], 1, source_ + 12 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[93], 1, source_ + 12 * size, 1);
            // AxCy
            C_DAXPY(size, -1.0, libderiv_.ABCD[30], 1, source_ + 13 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[91], 1, source_ + 13 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[94], 1, source_ + 13 * size, 1);
            // AxCz
            C_DAXPY(size, -1.0, libderiv_.ABCD[42], 1, source_ + 14 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[92], 1, source_ + 14 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[95], 1, source_ + 14 * size, 1);
            // AxDx
            memcpy(source_ + 15 * size, libderiv_.ABCD[18], sizeof(double) * size);
            // AxDy
            memcpy(source_ + 16 * size, libderiv_.ABCD[30], sizeof(double) * size);
            // AxDz
            memcpy(source_ + 17 * size, libderiv_.ABCD[42], sizeof(double) * size);
            // AyAy
            memcpy(source_ + 18 * size, libderiv_.ABCD[103], sizeof(double) * size);
            // AyAz
            memcpy(source_ + 19 * size, libderiv_.ABCD[104], sizeof(double) * size);
            // AyCx
            C_DAXPY(size, -1.0, libderiv_.ABCD[19], 1, source_ + 20 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[91], 1, source_ + 20 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[105], 1, source_ + 20 * size, 1);
            // AyCy
            C_DAXPY(size, -1.0, libderiv_.ABCD[31], 1, source_ + 21 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[103], 1, source_ + 21 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[106], 1, source_ + 21 * size, 1);
            // AyCz
            C_DAXPY(size, -1.0, libderiv_.ABCD[43], 1, source_ + 22 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[104], 1, source_ + 22 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[107], 1, source_ + 22 * size, 1);
            // AyDx
            memcpy(source_ + 23 * size, libderiv_.ABCD[19], sizeof(double) * size);
            // AyDy
            memcpy(source_ + 24 * size, libderiv_.ABCD[31], sizeof(double) * size);
            // AyDz
            memcpy(source_ + 25 * size, libderiv_.ABCD[43], sizeof(double) * size);
            // AzAz
            memcpy(source_ + 26 * size, libderiv_.ABCD[116], sizeof(double) * size);
            // AzCx
            C_DAXPY(size, -1.0, libderiv_.ABCD[20], 1, source_ + 27 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[92], 1, source_ + 27 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[117], 1, source_ + 27 * size, 1);
            // AzCy
            C_DAXPY(size, -1.0, libderiv_.ABCD[32], 1, source_ + 28 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[104], 1, source_ + 28 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[118], 1, source_ + 28 * size, 1);
            // AzCz
            C_DAXPY(size, -1.0, libderiv_.ABCD[44], 1, source_ + 29 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[116], 1, source_ + 29 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[119], 1, source_ + 29 * size, 1);
            // AzDx
            memcpy(source_ + 30 * size, libderiv_.ABCD[20], sizeof(double) * size);
            // AzDy
            memcpy(source_ + 31 * size, libderiv_.ABCD[32], sizeof(double) * size);
            // AzDz
            memcpy(source_ + 32 * size, libderiv_.ABCD[44], sizeof(double) * size);
            // CxCx
            C_DAXPY(size, 1.0, libderiv_.ABCD[12], 1, source_ + 33 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[18], 1, source_ + 33 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[90], 1, source_ + 33 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[21], 1, source_ + 33 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[93], 1, source_ + 33 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[129], 1, source_ + 33 * size, 1);
            // CxCy
            C_DAXPY(size, 1.0, libderiv_.ABCD[13], 1, source_ + 34 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[19], 1, source_ + 34 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[22], 1, source_ + 34 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[30], 1, source_ + 34 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[91], 1, source_ + 34 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[94], 1, source_ + 34 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[33], 1, source_ + 34 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[105], 1, source_ + 34 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[130], 1, source_ + 34 * size, 1);
            // CxCz
            C_DAXPY(size, 1.0, libderiv_.ABCD[14], 1, source_ + 35 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[20], 1, source_ + 35 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[23], 1, source_ + 35 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[42], 1, source_ + 35 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[92], 1, source_ + 35 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[95], 1, source_ + 35 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[45], 1, source_ + 35 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[117], 1, source_ + 35 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[131], 1, source_ + 35 * size, 1);
            // CxDx
            C_DAXPY(size, -1.0, libderiv_.ABCD[12], 1, source_ + 36 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[18], 1, source_ + 36 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[21], 1, source_ + 36 * size, 1);
            // CxDy
            C_DAXPY(size, -1.0, libderiv_.ABCD[13], 1, source_ + 37 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[30], 1, source_ + 37 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[33], 1, source_ + 37 * size, 1);
            // CxDz
            C_DAXPY(size, -1.0, libderiv_.ABCD[14], 1, source_ + 38 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[42], 1, source_ + 38 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[45], 1, source_ + 38 * size, 1);
            // CyCy
            C_DAXPY(size, 1.0, libderiv_.ABCD[25], 1, source_ + 39 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[31], 1, source_ + 39 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[103], 1, source_ + 39 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[34], 1, source_ + 39 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[106], 1, source_ + 39 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[142], 1, source_ + 39 * size, 1);
            // CyCz
            C_DAXPY(size, 1.0, libderiv_.ABCD[26], 1, source_ + 40 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[32], 1, source_ + 40 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[35], 1, source_ + 40 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[43], 1, source_ + 40 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[104], 1, source_ + 40 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[107], 1, source_ + 40 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[46], 1, source_ + 40 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[118], 1, source_ + 40 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[143], 1, source_ + 40 * size, 1);
            // CyDx
            C_DAXPY(size, -1.0, libderiv_.ABCD[13], 1, source_ + 41 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[19], 1, source_ + 41 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[22], 1, source_ + 41 * size, 1);
            // CyDy
            C_DAXPY(size, -1.0, libderiv_.ABCD[25], 1, source_ + 42 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[31], 1, source_ + 42 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[34], 1, source_ + 42 * size, 1);
            // CyDz
            C_DAXPY(size, -1.0, libderiv_.ABCD[26], 1, source_ + 43 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[43], 1, source_ + 43 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[46], 1, source_ + 43 * size, 1);
            // CzCz
            C_DAXPY(size, 1.0, libderiv_.ABCD[38], 1, source_ + 44 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[44], 1, source_ + 44 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[116], 1, source_ + 44 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[47], 1, source_ + 44 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[119], 1, source_ + 44 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[155], 1, source_ + 44 * size, 1);
            // CzDx
            C_DAXPY(size, -1.0, libderiv_.ABCD[14], 1, source_ + 45 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[20], 1, source_ + 45 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[23], 1, source_ + 45 * size, 1);
            // CzDy
            C_DAXPY(size, -1.0, libderiv_.ABCD[26], 1, source_ + 46 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[32], 1, source_ + 46 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[35], 1, source_ + 46 * size, 1);
            // CzDz
            C_DAXPY(size, -1.0, libderiv_.ABCD[38], 1, source_ + 47 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[44], 1, source_ + 47 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[47], 1, source_ + 47 * size, 1);
            // DxDx
            memcpy(source_ + 48 * size, libderiv_.ABCD[12], sizeof(double) * size);
            // DxDy
            memcpy(source_ + 49 * size, libderiv_.ABCD[13], sizeof(double) * size);
            // DxDz
            memcpy(source_ + 50 * size, libderiv_.ABCD[14], sizeof(double) * size);
            // DyDy
            memcpy(source_ + 51 * size, libderiv_.ABCD[25], sizeof(double) * size);
            // DyDz
            memcpy(source_ + 52 * size, libderiv_.ABCD[26], sizeof(double) * size);
            // DzDz
            memcpy(source_ + 53 * size, libderiv_.ABCD[38], sizeof(double) * size);
            break;
        case DCBA:
            // Ax
            memcpy(source_ + 0 * size, libderiv_.ABCD[9], sizeof(double) * size);
            // Ay
            memcpy(source_ + 1 * size, libderiv_.ABCD[10], sizeof(double) * size);
            // Az
            memcpy(source_ + 2 * size, libderiv_.ABCD[11], sizeof(double) * size);
            // Cx
            C_DAXPY(size, -1.0, libderiv_.ABCD[0], 1, source_ + 3 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[6], 1, source_ + 3 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[9], 1, source_ + 3 * size, 1);
            // Cy
            C_DAXPY(size, -1.0, libderiv_.ABCD[1], 1, source_ + 4 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[7], 1, source_ + 4 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[10], 1, source_ + 4 * size, 1);
            // Cz
            C_DAXPY(size, -1.0, libderiv_.ABCD[2], 1, source_ + 5 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[8], 1, source_ + 5 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[11], 1, source_ + 5 * size, 1);
            // Dx
            memcpy(source_ + 6 * size, libderiv_.ABCD[0], sizeof(double) * size);
            // Dy
            memcpy(source_ + 7 * size, libderiv_.ABCD[1], sizeof(double) * size);
            // Dz
            memcpy(source_ + 8 * size, libderiv_.ABCD[2], sizeof(double) * size);
            // AxAx
            memcpy(source_ + 9 * size, libderiv_.ABCD[129], sizeof(double) * size);
            // AxAy
            memcpy(source_ + 10 * size, libderiv_.ABCD[130], sizeof(double) * size);
            // AxAz
            memcpy(source_ + 11 * size, libderiv_.ABCD[131], sizeof(double) * size);
            // AxCx
            C_DAXPY(size, -1.0, libderiv_.ABCD[21], 1, source_ + 12 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[93], 1, source_ + 12 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[129], 1, source_ + 12 * size, 1);
            // AxCy
            C_DAXPY(size, -1.0, libderiv_.ABCD[33], 1, source_ + 13 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[105], 1, source_ + 13 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[130], 1, source_ + 13 * size, 1);
            // AxCz
            C_DAXPY(size, -1.0, libderiv_.ABCD[45], 1, source_ + 14 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[117], 1, source_ + 14 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[131], 1, source_ + 14 * size, 1);
            // AxDx
            memcpy(source_ + 15 * size, libderiv_.ABCD[21], sizeof(double) * size);
            // AxDy
            memcpy(source_ + 16 * size, libderiv_.ABCD[33], sizeof(double) * size);
            // AxDz
            memcpy(source_ + 17 * size, libderiv_.ABCD[45], sizeof(double) * size);
            // AyAy
            memcpy(source_ + 18 * size, libderiv_.ABCD[142], sizeof(double) * size);
            // AyAz
            memcpy(source_ + 19 * size, libderiv_.ABCD[143], sizeof(double) * size);
            // AyCx
            C_DAXPY(size, -1.0, libderiv_.ABCD[22], 1, source_ + 20 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[94], 1, source_ + 20 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[130], 1, source_ + 20 * size, 1);
            // AyCy
            C_DAXPY(size, -1.0, libderiv_.ABCD[34], 1, source_ + 21 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[106], 1, source_ + 21 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[142], 1, source_ + 21 * size, 1);
            // AyCz
            C_DAXPY(size, -1.0, libderiv_.ABCD[46], 1, source_ + 22 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[118], 1, source_ + 22 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[143], 1, source_ + 22 * size, 1);
            // AyDx
            memcpy(source_ + 23 * size, libderiv_.ABCD[22], sizeof(double) * size);
            // AyDy
            memcpy(source_ + 24 * size, libderiv_.ABCD[34], sizeof(double) * size);
            // AyDz
            memcpy(source_ + 25 * size, libderiv_.ABCD[46], sizeof(double) * size);
            // AzAz
            memcpy(source_ + 26 * size, libderiv_.ABCD[155], sizeof(double) * size);
            // AzCx
            C_DAXPY(size, -1.0, libderiv_.ABCD[23], 1, source_ + 27 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[95], 1, source_ + 27 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[131], 1, source_ + 27 * size, 1);
            // AzCy
            C_DAXPY(size, -1.0, libderiv_.ABCD[35], 1, source_ + 28 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[107], 1, source_ + 28 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[143], 1, source_ + 28 * size, 1);
            // AzCz
            C_DAXPY(size, -1.0, libderiv_.ABCD[47], 1, source_ + 29 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[119], 1, source_ + 29 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[155], 1, source_ + 29 * size, 1);
            // AzDx
            memcpy(source_ + 30 * size, libderiv_.ABCD[23], sizeof(double) * size);
            // AzDy
            memcpy(source_ + 31 * size, libderiv_.ABCD[35], sizeof(double) * size);
            // AzDz
            memcpy(source_ + 32 * size, libderiv_.ABCD[47], sizeof(double) * size);
            // CxCx
            C_DAXPY(size, 1.0, libderiv_.ABCD[12], 1, source_ + 33 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[18], 1, source_ + 33 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[90], 1, source_ + 33 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[21], 1, source_ + 33 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[93], 1, source_ + 33 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[129], 1, source_ + 33 * size, 1);
            // CxCy
            C_DAXPY(size, 1.0, libderiv_.ABCD[13], 1, source_ + 34 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[19], 1, source_ + 34 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[22], 1, source_ + 34 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[30], 1, source_ + 34 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[91], 1, source_ + 34 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[94], 1, source_ + 34 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[33], 1, source_ + 34 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[105], 1, source_ + 34 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[130], 1, source_ + 34 * size, 1);
            // CxCz
            C_DAXPY(size, 1.0, libderiv_.ABCD[14], 1, source_ + 35 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[20], 1, source_ + 35 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[23], 1, source_ + 35 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[42], 1, source_ + 35 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[92], 1, source_ + 35 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[95], 1, source_ + 35 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[45], 1, source_ + 35 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[117], 1, source_ + 35 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[131], 1, source_ + 35 * size, 1);
            // CxDx
            C_DAXPY(size, -1.0, libderiv_.ABCD[12], 1, source_ + 36 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[18], 1, source_ + 36 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[21], 1, source_ + 36 * size, 1);
            // CxDy
            C_DAXPY(size, -1.0, libderiv_.ABCD[13], 1, source_ + 37 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[30], 1, source_ + 37 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[33], 1, source_ + 37 * size, 1);
            // CxDz
            C_DAXPY(size, -1.0, libderiv_.ABCD[14], 1, source_ + 38 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[42], 1, source_ + 38 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[45], 1, source_ + 38 * size, 1);
            // CyCy
            C_DAXPY(size, 1.0, libderiv_.ABCD[25], 1, source_ + 39 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[31], 1, source_ + 39 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[103], 1, source_ + 39 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[34], 1, source_ + 39 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[106], 1, source_ + 39 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[142], 1, source_ + 39 * size, 1);
            // CyCz
            C_DAXPY(size, 1.0, libderiv_.ABCD[26], 1, source_ + 40 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[32], 1, source_ + 40 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[35], 1, source_ + 40 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[43], 1, source_ + 40 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[104], 1, source_ + 40 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[107], 1, source_ + 40 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[46], 1, source_ + 40 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[118], 1, source_ + 40 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[143], 1, source_ + 40 * size, 1);
            // CyDx
            C_DAXPY(size, -1.0, libderiv_.ABCD[13], 1, source_ + 41 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[19], 1, source_ + 41 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[22], 1, source_ + 41 * size, 1);
            // CyDy
            C_DAXPY(size, -1.0, libderiv_.ABCD[25], 1, source_ + 42 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[31], 1, source_ + 42 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[34], 1, source_ + 42 * size, 1);
            // CyDz
            C_DAXPY(size, -1.0, libderiv_.ABCD[26], 1, source_ + 43 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[43], 1, source_ + 43 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[46], 1, source_ + 43 * size, 1);
            // CzCz
            C_DAXPY(size, 1.0, libderiv_.ABCD[38], 1, source_ + 44 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[44], 1, source_ + 44 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[116], 1, source_ + 44 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[47], 1, source_ + 44 * size, 1);
            C_DAXPY(size, 2.0, libderiv_.ABCD[119], 1, source_ + 44 * size, 1);
            C_DAXPY(size, 1.0, libderiv_.ABCD[155], 1, source_ + 44 * size, 1);
            // CzDx
            C_DAXPY(size, -1.0, libderiv_.ABCD[14], 1, source_ + 45 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[20], 1, source_ + 45 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[23], 1, source_ + 45 * size, 1);
            // CzDy
            C_DAXPY(size, -1.0, libderiv_.ABCD[26], 1, source_ + 46 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[32], 1, source_ + 46 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[35], 1, source_ + 46 * size, 1);
            // CzDz
            C_DAXPY(size, -1.0, libderiv_.ABCD[38], 1, source_ + 47 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[44], 1, source_ + 47 * size, 1);
            C_DAXPY(size, -1.0, libderiv_.ABCD[47], 1, source_ + 47 * size, 1);
            // DxDx
            memcpy(source_ + 48 * size, libderiv_.ABCD[12], sizeof(double) * size);
            // DxDy
            memcpy(source_ + 49 * size, libderiv_.ABCD[13], sizeof(double) * size);
            // DxDz
            memcpy(source_ + 50 * size, libderiv_.ABCD[14], sizeof(double) * size);
            // DyDy
            memcpy(source_ + 51 * size, libderiv_.ABCD[25], sizeof(double) * size);
            // DyDz
            memcpy(source_ + 52 * size, libderiv_.ABCD[26], sizeof(double) * size);
            // DzDz
            memcpy(source_ + 53 * size, libderiv_.ABCD[38], sizeof(double) * size);
            break;
    }
}

/**
 * @brief Fills the primitive data structure used by libint/libderiv with information from the ShellPairs
 * @param PrimQuartet The structure to hold the data.
 * @param fjt Object used to compute the fundamental integrals.
 * @param p12 ShellPair data structure for the left
 * @param p34 ShellPair data structure for the right
 * @param am Total angular momentum of this quartet
 * @param sh1eqsh2 Is the shell on center 1 identical to that on center 2?
 * @param sh3eqsh4 Is the shell on center 3 identical to that on center 4?
 * @param deriv_lvl Derivitive level of the integral
 * @return The total number of primitive combinations found. This is passed to libint/libderiv.
 */
static size_t fill_primitive_data(prim_data *PrimQuartet, Fjt *fjt, const L1ShellPair &sp12, const L1ShellPair &sp34,
                                  int am, bool sh1eqsh2, bool sh3eqsh4, int deriv_lvl) {
    double zeta, eta, ooze, rho, poz, coef1, PQx, PQy, PQz, PQ2, Wx, Wy, Wz, o12, o34, T, *F;
    double a1, a2, a3, a4;
    int p12, p34, i;
    size_t nprim = 0L;
    for (p12 = 0; p12 < sp12.nonzeroPrimPairs.size(); ++p12) {
        const PrimPair &pp12 = sp12.nonzeroPrimPairs[p12];
        a1 = pp12.ai;
        a2 = pp12.aj;
        o12 = pp12.overlap;
        zeta = pp12.gamma;
        double PAx = pp12.PA[0];
        double PAy = pp12.PA[1];
        double PAz = pp12.PA[2];
        double PBx = pp12.PB[0];
        double PBy = pp12.PB[1];
        double PBz = pp12.PB[2];
        double PABx = pp12.P[0];
        double PABy = pp12.P[1];
        double PABz = pp12.P[2];
        for (p34 = 0; p34 < sp34.nonzeroPrimPairs.size(); ++p34) {
            const PrimPair &pp34 = sp34.nonzeroPrimPairs[p34];
            a3 = pp34.ai;
            a4 = pp34.aj;
            o34 = pp34.overlap;
            eta = pp34.gamma;
            double PCx = pp34.PA[0];
            double PCy = pp34.PA[1];
            double PCz = pp34.PA[2];
            double PDx = pp34.PB[0];
            double PDy = pp34.PB[1];
            double PDz = pp34.PB[2];
            double PCDx = pp34.P[0];
            double PCDy = pp34.P[1];
            double PCDz = pp34.P[2];

            ooze = 1.0 / (zeta + eta);
            poz = eta * ooze;
            rho = zeta * poz;
            coef1 = 2.0 * sqrt(rho * M_1_PI) * o12 * o34;

            PrimQuartet[nprim].poz = poz;
            PrimQuartet[nprim].oo2zn = 0.5 * ooze;
            PrimQuartet[nprim].pon = zeta * ooze;
            PrimQuartet[nprim].oo2z = 0.5 / zeta;
            PrimQuartet[nprim].oo2n = 0.5 / eta;
            PrimQuartet[nprim].twozeta_a = 2.0 * a1;
            PrimQuartet[nprim].twozeta_b = 2.0 * a2;
            PrimQuartet[nprim].twozeta_c = 2.0 * a3;
            PrimQuartet[nprim].twozeta_d = 2.0 * a4;

            PQx = PABx - PCDx;
            PQy = PABy - PCDy;
            PQz = PABz - PCDz;
            PQ2 = PQx * PQx + PQy * PQy + PQz * PQz;

            Wx = (PABx * zeta + PCDx * eta) * ooze;
            Wy = (PABy * zeta + PCDy * eta) * ooze;
            Wz = (PABz * zeta + PCDz * eta) * ooze;

            // PA
            PrimQuartet[nprim].U[0][0] = PAx;
            PrimQuartet[nprim].U[0][1] = PAy;
            PrimQuartet[nprim].U[0][2] = PAz;
            // PB
            PrimQuartet[nprim].U[1][0] = PBx;
            PrimQuartet[nprim].U[1][1] = PBy;
            PrimQuartet[nprim].U[1][2] = PBz;
            // QC
            PrimQuartet[nprim].U[2][0] = PCx;
            PrimQuartet[nprim].U[2][1] = PCy;
            PrimQuartet[nprim].U[2][2] = PCz;
            // QD
            PrimQuartet[nprim].U[3][0] = PDx;
            PrimQuartet[nprim].U[3][1] = PDy;
            PrimQuartet[nprim].U[3][2] = PDz;
            // WP
            PrimQuartet[nprim].U[4][0] = Wx - PABx;
            PrimQuartet[nprim].U[4][1] = Wy - PABy;
            PrimQuartet[nprim].U[4][2] = Wz - PABz;
            // WQ
            PrimQuartet[nprim].U[5][0] = Wx - PCDx;
            PrimQuartet[nprim].U[5][1] = Wy - PCDy;
            PrimQuartet[nprim].U[5][2] = Wz - PCDz;

            T = rho * PQ2;
            fjt->set_rho(rho);
            F = fjt->values(am + deriv_lvl, T);

            for (i = 0; i <= am + deriv_lvl; ++i) PrimQuartet[nprim].F[i] = F[i] * coef1;

            nprim++;
        }
    }
    return nprim;
}
#endif  // ENABLE_Libint1t

}  // end namespace

#ifdef ENABLE_Libint1t
TwoElectronInt::TwoElectronInt(const IntegralFactory *integral, int deriv, bool use_shell_pairs)
    : TwoBodyAOInt(integral, deriv), use_shell_pairs_(use_shell_pairs) {
    // Initialize libint static data
    init_libint_base();
    if (deriv_) init_libderiv_base();

    // Figure out some information to initialize libint/libderiv with
    // 1. Maximum angular momentum
    int max_am = MAX(MAX(basis1()->max_am(), basis2()->max_am()), MAX(basis3()->max_am(), basis4()->max_am()));
    // 2. Maximum number of primitive combinations
    int max_nprim = basis1()->max_nprimitive() * basis2()->max_nprimitive() * basis3()->max_nprimitive() *
                    basis4()->max_nprimitive();
    // 3. Maximum Cartesian class size
    max_cart_ = ioff[basis1()->max_am() + 1] * ioff[basis2()->max_am() + 1] * ioff[basis3()->max_am() + 1] *
                ioff[basis4()->max_am() + 1];

    // Make sure libint is compiled to handle our max AM
    if (max_am >= LIBINT_MAX_AM) {
        outfile->Printf(
            "ERROR: ERI - Libint cannot handle angular momentum this high (%d).\n"
            "       Rebuild Libint with MAX_AM_ERI at least %d.\n",
            max_am, max_am);
        throw LimitExceeded<int>(
            "ERI - Libint cannot handle angular momentum this high.\n"
            "Rebuild Libint with MAX_AM_ERI at least (actual).\n",
            LIBINT_MAX_AM - 1, max_am, __FILE__, __LINE__);
    } else if (deriv_ == 1 && max_am >= LIBDERIV_MAX_AM1) {
        outfile->Printf(
            "ERROR: ERI - Libint cannot handle angular momentum this high (%d) for first derivatives.\n"
            "     Rebuild Libint with MAX_AM_ERI at least %d.\n",
            max_am, max_am + 1);
        throw LimitExceeded<int>(
            "ERI - Libint cannot handle angular momentum this high.\n"
            "Rebuild Libint with MAX_AM_ERI at least (actual + 1).\n",
            LIBDERIV_MAX_AM1 - 1, max_am, __FILE__, __LINE__);
    } else if (deriv_ == 2 && max_am >= LIBDERIV_MAX_AM12) {
        outfile->Printf(
            "ERROR: ERI - Libint cannot handle angular momentum this high (%d) for second derivatives.\n"
            "       Reconfigure Libint with MAX_AM_ERI at least %d\n",
            max_am, max_am + 2);
        throw LimitExceeded<int>(
            "ERI - Libint cannot handle angular momentum this high.\n"
            "Rebuild Libint with MAX_AM_ERI at least (actual + 2).\n",
            LIBDERIV_MAX_AM12 - 1, max_am, __FILE__, __LINE__);
    } else if (deriv_ > 2) {
        outfile->Printf("ERROR: ERI - Cannot compute higher than second derivatives.");
        throw PSIEXCEPTION("ERI - Cannot compute higher than second derivatives.");
    }

    try {
        // Initialize libint
        init_libint(&libint_, max_am, max_nprim);
        // and libderiv, if needed
        if (deriv_ == 1)
            init_libderiv1(&libderiv_, max_am, max_nprim, max_cart_);
        else if (deriv_ == 2)
            init_libderiv12(&libderiv_, max_am, max_nprim, max_cart_);
    } catch (std::bad_alloc &e) {
        outfile->Printf("Error allocating memory for libint/libderiv.\n");
        exit(EXIT_FAILURE);
    }
    size_t size = INT_NCART(basis1()->max_am()) * INT_NCART(basis2()->max_am()) * INT_NCART(basis3()->max_am()) *
                  INT_NCART(basis4()->max_am());

    // Used in pure_transform
    try {
        tformbuf_ = new double[size];
    } catch (std::bad_alloc &e) {
        outfile->Printf("Error allocating tformbuf_.\n%s\n", e.what());
        exit(EXIT_FAILURE);
    }
    memset(tformbuf_, 0, sizeof(double) * size);

    // ntypes is the number of integral types provided by libint/libderiv.
    size *= ntypes[deriv_];

    try {
        target_full_ = new double[size];
        target_ = target_full_;
    } catch (std::bad_alloc &e) {
        outfile->Printf("Error allocating target_.\n%s\n", e.what());
        exit(EXIT_FAILURE);
    }
    memset(target_, 0, sizeof(double) * size);

    try {
        source_full_ = new double[size];
        source_ = source_full_;
    } catch (std::bad_alloc &e) {
        outfile->Printf("Error allocating source_.\n%s\n", e.what());
        exit(EXIT_FAILURE);
    }
    memset(source_, 0, sizeof(double) * size);

    if (basis1() != basis2() || basis1() != basis3() || basis2() != basis4()) {
        use_shell_pairs_ = false;
    }

    if (use_shell_pairs_) {
        // Precompute a bunch of information
        init_shell_pairs12();
        // If basis3 and basis4 equals basis1 and basis2, then the following function will do nothing,
        // except assign pairs34_ to pairs12_
        init_shell_pairs34();
    }
}

TwoElectronInt::~TwoElectronInt() {
    delete[] tformbuf_;
    delete[] target_full_;
    delete[] source_full_;
    free_libint(&libint_);
    if (deriv_) free_libderiv(&libderiv_);
}

void TwoElectronInt::init_shell_pairs12() {
    Vector3 P, PA, PB, AB, A, B;
    int i, j, si, sj, np_i, np_j, np_sig;
    double a1, a2, ab2, gam, c1, c2, overlap;

    // Make a shared pointer to a vector of vectors of ShellPair
    pairs12_ = std::make_shared<std::vector<std::vector<L1ShellPair>>>(
        std::vector<std::vector<L1ShellPair>>(basis2()->nshell(), std::vector<L1ShellPair>(basis2()->nshell())));

    // Counts of primitive pairs (over all shell pairs) [FOR PRINING]
    int prim_pairs_total = 0;
    int prim_pairs_keep = 0;

    // Criteria for determining if the overlap between primitives is negligible
    // #TODO make this a user defined parameter
    double overlap_cutoff = 1.0e-12;

    // Loop over all shell pairs (si, sj) and create primitive pairs pairs
    for (si = 0; si < basis1()->nshell(); ++si) {
        A = basis1()->shell(si).center();

        for (sj = 0; sj < basis2()->nshell(); ++sj) {
            B = basis2()->shell(sj).center();

            AB = A - B;
            ab2 = AB.dot(AB);

            // Make and populate fields of screened shell pair
            L1ShellPair sp;
            sp.i = si;
            sp.j = sj;
            sp.AB[0] = AB[0];
            sp.AB[1] = AB[1];
            sp.AB[2] = AB[2];

            np_i = basis1()->shell(si).nprimitive();
            np_j = basis2()->shell(sj).nprimitive();

            // Reserve some memory for the primitives pairs of this shell pair
            sp.nonzeroPrimPairs.resize(np_i * np_j);

            // the number of significant primitive pairs just in this shell pair (not screened)
            np_sig = 0;

            // Pre-compute all data that we can:
            for (i = 0; i < np_i; ++i) {
                a1 = basis1()->shell(si).exp(i);
                c1 = basis1()->shell(si).coef(i);

                for (j = 0; j < np_j; ++j) {
                    a2 = basis2()->shell(sj).exp(j);
                    c2 = basis2()->shell(sj).coef(j);
                    gam = a1 + a2;
                    overlap = pow(M_PI / gam, 3.0 / 2.0) * exp(-a1 * a2 * ab2 / gam) * c1 * c2;

                    // Check overlap for screening condition
                    if (fabs(overlap) >= overlap_cutoff) {
                        // Compute Gaussian product and component distances
                        P = (A * a1 + B * a2) / gam;
                        PA = P - A;
                        PB = P - B;

                        // Make and populate fields of screened primitive pair
                        PrimPair pp;

                        // Copy data into pairs array
                        pp.P[0] = P[0];
                        pp.P[1] = P[1];
                        pp.P[2] = P[2];
                        pp.PA[0] = PA[0];
                        pp.PA[1] = PA[1];
                        pp.PA[2] = PA[2];
                        pp.PB[0] = PB[0];
                        pp.PB[1] = PB[1];
                        pp.PB[2] = PB[2];

                        // Save some information
                        pp.ai = a1;
                        pp.aj = a2;
                        pp.gamma = gam;
                        pp.ci = c1;
                        pp.cj = c2;
                        pp.overlap = overlap;

                        // Store this primitive pair in the shell pair
                        sp.nonzeroPrimPairs[np_sig] = pp;
                        ++np_sig;
                    } else {
                        // do nothing
                    }
                }
            }

            // Update total counts
            prim_pairs_total += np_i * np_j;
            prim_pairs_keep += np_sig;

            sp.nonzeroPrimPairs.resize(np_sig);
            sp.nonzeroPrimPairs.shrink_to_fit();

            (*pairs12_)[si][sj] = sp;
        }
    }

    double percent_keep = prim_pairs_keep * 100.0 / prim_pairs_total;
    outfile->Printf("\n  ShellPair Screening for Disk ERIs: Overlap Cutoff of %.12f \n", overlap_cutoff);
    outfile->Printf("  Using %d of %d gaussian primitive products (%.2f%%)\n\n", prim_pairs_keep, prim_pairs_total,
                    percent_keep);
}

void TwoElectronInt::init_shell_pairs34() {
    // If basis1 == basis3 && basis2 == basis4, then we don't need to do anything except use the pointer
    // of pairs12_.
    if (use_shell_pairs_ == true) {
        // This assumes init_shell_pairs12 was called and precomputed the values.
        pairs34_ = pairs12_;
        return;
    }
}

size_t TwoElectronInt::memory_to_store_shell_pairs(const std::shared_ptr<BasisSet> &bs1,
                                                   const std::shared_ptr<BasisSet> &bs2) {
    int i, j, np_i, np_j;
    size_t mem = 0;

    for (i = 0; i < bs1->nshell(); ++i) {
        np_i = bs1->shell(i).nprimitive();
        for (j = 0; j < bs2->nshell(); ++j) {
            np_j = bs2->shell(j).nprimitive();
            mem += (2 * (np_i + np_j) + 11 * np_i * np_j);
        }
    }
    return mem;
}

size_t TwoElectronInt::compute_shell(const AOShellCombinationsIterator &shellIter) {
    return compute_shell(shellIter.p(), shellIter.q(), shellIter.r(), shellIter.s());
}

size_t TwoElectronInt::compute_shell_for_sieve(const std::shared_ptr<BasisSet> bs, int sh1, int sh2, int sh3, int sh4,
                                               bool is_bra) {
    return compute_shell(sh1, sh2, sh1, sh2);
}

size_t TwoElectronInt::compute_shell(int sh1, int sh2, int sh3, int sh4) {
#ifdef MINTS_TIMER
    timer_on("ERI::compute_shell");
#endif
// Need to ensure the ordering asked by the user is valid for libint
// compute_quartet does NOT check this. SEGFAULTS should occur if order
// is not guaranteed.
#ifdef MINTS_TIMER
    timer_on("reorder");
#endif

    int s1, s2, s3, s4;
    int am1, am2, am3, am4, temp;
    std::shared_ptr<BasisSet> bs_temp;

    p13p24_ = false;
    p12_ = false;
    p34_ = false;

    // AM used for ordering
    am1 = original_bs1_->shell(sh1).am();
    am2 = original_bs2_->shell(sh2).am();
    am3 = original_bs3_->shell(sh3).am();
    am4 = original_bs4_->shell(sh4).am();
    temp = am1 + am2 + am3 + am4;

    // c1 = original_bs1_->shell(sh1).ncenter();
    // c2 = original_bs1_->shell(sh2).ncenter();
    // c3 = original_bs1_->shell(sh3).ncenter();
    // c4 = original_bs1_->shell(sh4).ncenter();

    // TODO: Check this!
    //	if (c1 == c2 && c1 == c3 && c1 && c4 && temp % 2 != 0) {
    //#ifdef MINTS_TIMER
    //		timer_off("reorder");
    //		timer_off("ERI::compute_shell");
    //#endif
    //		return 0;
    //	}

    int n1, n2, n3, n4;

    n1 = original_bs1_->shell(sh1).nfunction();
    n2 = original_bs2_->shell(sh2).nfunction();
    n3 = original_bs3_->shell(sh3).nfunction();
    n4 = original_bs4_->shell(sh4).nfunction();
    curr_buff_size_ = n1 * n2 * n3 * n4;

    // Save the original requested shell ordering. The pre-computed shell pair information
    // requires the original ordering.
    osh1_ = sh1;
    osh2_ = sh2;
    osh3_ = sh3;
    osh4_ = sh4;

    // l(a) >= l(b), l(c) >= l(d), and l(c) + l(d) >= l(a) + l(b).
    if (am1 >= am2) {
        s1 = sh1;
        s2 = sh2;

        bs1_ = original_bs1_;
        bs2_ = original_bs2_;
    } else {
        s1 = sh2;
        s2 = sh1;

        bs1_ = original_bs2_;
        bs2_ = original_bs1_;

        p12_ = true;
    }

    if (am3 >= am4) {
        s3 = sh3;
        s4 = sh4;

        bs3_ = original_bs3_;
        bs4_ = original_bs4_;

    } else {
        s3 = sh4;
        s4 = sh3;

        bs3_ = original_bs4_;
        bs4_ = original_bs3_;

        p34_ = true;
    }

    if ((am1 + am2) > (am3 + am4)) {
        // Swap s1 and s2 with s3 and s4
        temp = s1;
        s1 = s3;
        s3 = temp;

        temp = s2;
        s2 = s4;
        s4 = temp;

        bs_temp = bs1_;
        bs1_ = bs3_;
        bs3_ = bs_temp;

        bs_temp = bs2_;
        bs2_ = bs4_;
        bs4_ = bs_temp;

        p13p24_ = true;
    }
#ifdef MINTS_TIMER
    timer_off("reorder");
#endif

    // s1, s2, s3, s4 contain the shells to do in libint order
    size_t ncomputed = compute_quartet(s1, s2, s3, s4);
    if (ncomputed) {
        // Only do the following if we did any work.

        // Permute integrals back, if needed
        if (p12_ || p34_ || p13p24_) {
#ifdef MINTS_TIMER
            timer_on("permute_target");
#endif
            permute_target(source_, target_, s1, s2, s3, s4, p12_, p34_, p13p24_);
#ifdef MINTS_TIMER
            timer_off("permute_target");
#endif
        } else {
#ifdef MINTS_TIMER
            timer_on("memcpy - no resort");
#endif
            // copy the integrals to the target_
            memcpy(target_, source_, n1 * n2 * n3 * n4 * sizeof(double));
#ifdef MINTS_TIMER
            timer_off("memcpy - no resort");
#endif
        }
    }

#ifdef MINTS_TIMER
    timer_off("ERI::compute_shell");
#endif
    return ncomputed;
}

size_t TwoElectronInt::compute_quartet(int sh1, int sh2, int sh3, int sh4) {
#ifdef MINTS_TIMER
    timer_on("setup");
#endif
    const GaussianShell &s1 = bs1_->shell(sh1);
    const GaussianShell &s2 = bs2_->shell(sh2);
    const GaussianShell &s3 = bs3_->shell(sh3);
    const GaussianShell &s4 = bs4_->shell(sh4);

    int am1 = s1.am();
    int am2 = s2.am();
    int am3 = s3.am();
    int am4 = s4.am();
    int am = am1 + am2 + am3 + am4;  // total am
    int nprim1;
    int nprim2;
    int nprim3;
    int nprim4;
    double A[3], B[3], C[3], D[3];

    A[0] = s1.center()[0];
    A[1] = s1.center()[1];
    A[2] = s1.center()[2];
    B[0] = s2.center()[0];
    B[1] = s2.center()[1];
    B[2] = s2.center()[2];
    C[0] = s3.center()[0];
    C[1] = s3.center()[1];
    C[2] = s3.center()[2];
    D[0] = s4.center()[0];
    D[1] = s4.center()[1];
    D[2] = s4.center()[2];

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);
    double CD2 = 0.0;
    CD2 += (C[0] - D[0]) * (C[0] - D[0]);
    CD2 += (C[1] - D[1]) * (C[1] - D[1]);
    CD2 += (C[2] - D[2]) * (C[2] - D[2]);

    libint_.AB[0] = A[0] - B[0];
    libint_.AB[1] = A[1] - B[1];
    libint_.AB[2] = A[2] - B[2];
    libint_.CD[0] = C[0] - D[0];
    libint_.CD[1] = C[1] - D[1];
    libint_.CD[2] = C[2] - D[2];

#ifdef MINTS_TIMER
    timer_off("setup");
#endif

#ifdef MINTS_TIMER
    timer_on("Primitive setup");
#endif

    // Prepare all the data needed by libint
    size_t nprim = 0;
    nprim1 = s1.nprimitive();
    nprim2 = s2.nprimitive();
    nprim3 = s3.nprimitive();
    nprim4 = s4.nprimitive();

    // If we can, use the precomputed values found in ShellPair.
    if (use_shell_pairs_) {
        const L1ShellPair &sp12 = (*pairs12_)[sh1][sh2];
        const L1ShellPair &sp34 = (*pairs34_)[sh3][sh4];

        nprim = fill_primitive_data(libint_.PrimQuartet, fjt_, sp12, sp34, am, sh1 == sh2, sh3 == sh4, 0);

    } else {
        const double *a1s = s1.exps();
        const double *a2s = s2.exps();
        const double *a3s = s3.exps();
        const double *a4s = s4.exps();
        const double *c1s = s1.coefs();
        const double *c2s = s2.coefs();
        const double *c3s = s3.coefs();
        const double *c4s = s4.coefs();

        // Old version - without ShellPair - STILL USED BY RI CODES
        for (int p1 = 0; p1 < nprim1; ++p1) {
            double a1 = a1s[p1];
            double c1 = c1s[p1];
            for (int p2 = 0; p2 < nprim2; ++p2) {
                double a2 = a2s[p2];
                double c2 = c2s[p2];
                double zeta = a1 + a2;
                double ooz = 1.0 / zeta;
                double oo2z = 1.0 / (2.0 * zeta);

                double PA[3], PB[3];
                double P[3];

                P[0] = (a1 * A[0] + a2 * B[0]) * ooz;
                P[1] = (a1 * A[1] + a2 * B[1]) * ooz;
                P[2] = (a1 * A[2] + a2 * B[2]) * ooz;
                PA[0] = P[0] - A[0];
                PA[1] = P[1] - A[1];
                PA[2] = P[2] - A[2];
                PB[0] = P[0] - B[0];
                PB[1] = P[1] - B[1];
                PB[2] = P[2] - B[2];

                double Sab = pow(M_PI * ooz, 3.0 / 2.0) * exp(-a1 * a2 * ooz * AB2) * c1 * c2;

                for (int p3 = 0; p3 < nprim3; ++p3) {
                    double a3 = a3s[p3];
                    double c3 = c3s[p3];
                    for (int p4 = 0; p4 < nprim4; ++p4) {
                        double a4 = a4s[p4];
                        double c4 = c4s[p4];
                        double nu = a3 + a4;
                        double oon = 1.0 / nu;
                        double oo2n = 1.0 / (2.0 * nu);
                        double oo2zn = 1.0 / (2.0 * (zeta + nu));
                        double rho = (zeta * nu) / (zeta + nu);
                        double oo2rho = 1.0 / (2.0 * rho);

                        double QC[3], QD[3], WP[3], WQ[3], PQ[3];
                        double Q[3], W[3], a3C[3], a4D[3];

                        a3C[0] = a3 * C[0];
                        a3C[1] = a3 * C[1];
                        a3C[2] = a3 * C[2];

                        a4D[0] = a4 * D[0];
                        a4D[1] = a4 * D[1];
                        a4D[2] = a4 * D[2];

                        Q[0] = (a3C[0] + a4D[0]) * oon;
                        Q[1] = (a3C[1] + a4D[1]) * oon;
                        Q[2] = (a3C[2] + a4D[2]) * oon;

                        QC[0] = Q[0] - C[0];
                        QC[1] = Q[1] - C[1];
                        QC[2] = Q[2] - C[2];
                        QD[0] = Q[0] - D[0];
                        QD[1] = Q[1] - D[1];
                        QD[2] = Q[2] - D[2];
                        PQ[0] = P[0] - Q[0];
                        PQ[1] = P[1] - Q[1];
                        PQ[2] = P[2] - Q[2];

                        double PQ2 = 0.0;
                        PQ2 += (P[0] - Q[0]) * (P[0] - Q[0]);
                        PQ2 += (P[1] - Q[1]) * (P[1] - Q[1]);
                        PQ2 += (P[2] - Q[2]) * (P[2] - Q[2]);

                        W[0] = (zeta * P[0] + nu * Q[0]) / (zeta + nu);
                        W[1] = (zeta * P[1] + nu * Q[1]) / (zeta + nu);
                        W[2] = (zeta * P[2] + nu * Q[2]) / (zeta + nu);
                        WP[0] = W[0] - P[0];
                        WP[1] = W[1] - P[1];
                        WP[2] = W[2] - P[2];
                        WQ[0] = W[0] - Q[0];
                        WQ[1] = W[1] - Q[1];
                        WQ[2] = W[2] - Q[2];

                        for (int i = 0; i < 3; ++i) {
                            libint_.PrimQuartet[nprim].U[0][i] = PA[i];
                            libint_.PrimQuartet[nprim].U[2][i] = QC[i];
                            libint_.PrimQuartet[nprim].U[4][i] = WP[i];
                            libint_.PrimQuartet[nprim].U[5][i] = WQ[i];
                        }
                        libint_.PrimQuartet[nprim].oo2z = oo2z;
                        libint_.PrimQuartet[nprim].oo2n = oo2n;
                        libint_.PrimQuartet[nprim].oo2zn = oo2zn;
                        libint_.PrimQuartet[nprim].poz = rho * ooz;
                        libint_.PrimQuartet[nprim].pon = rho * oon;
                        libint_.PrimQuartet[nprim].oo2p = oo2rho;

                        double T = rho * PQ2;
                        fjt_->set_rho(rho);
                        double *F = fjt_->values(am, T);

                        // Modify F to include overlap of ab and cd, eqs 14, 15, 16 of libint manual
                        double Scd = pow(M_PI * oon, 3.0 / 2.0) * exp(-a3 * a4 * oon * CD2) * c3 * c4;
                        double val = 2.0 * sqrt(rho * M_1_PI) * Sab * Scd;
                        for (int i = 0; i <= am; ++i) {
                            libint_.PrimQuartet[nprim].F[i] = F[i] * val;
                        }
                        nprim++;
                    }
                }
            }
        }
    }
#ifdef MINTS_TIMER
    timer_off("Primitive setup");
#endif

    // How many are there?
    size_t size = INT_NCART(am1) * INT_NCART(am2) * INT_NCART(am3) * INT_NCART(am4);

#ifdef MINTS_TIMER
    timer_on("libint overhead");
#endif

    // Compute the integral
    if (am) {
        double *target_ints;

        target_ints = build_eri[am1][am2][am3][am4](&libint_, nprim);

        memcpy(source_, target_ints, sizeof(double) * size);
    } else {
        // Handle (ss|ss)
        double temp = 0.0;
        for (size_t i = 0; i < nprim; ++i) temp += (double)libint_.PrimQuartet[i].F[0];
        source_[0] = temp;
        //        outfile->Printf( "s-functions = %8.5f\n", temp);
    }

#ifdef MINTS_TIMER
    timer_off("libint overhead");
#endif

    // The following two functions time themselves.

    // Normalize the integrals for angular momentum
    // normalize_am(s1, s2, s3, s4);

    // Transform the integrals into pure angular momentum
    pure_transform(sh1, sh2, sh3, sh4, 1);

    // Results are in source_
    return size;
}

size_t TwoElectronInt::compute_shell_deriv1(int sh1, int sh2, int sh3, int sh4) {
    if (deriv_ < 1) {
        outfile->Printf("ERROR - ERI: ERI object not initialized to handle derivatives.\n");
        abort();
    }
    // Need to ensure the ordering asked by the user is valid for libint
    // compute_quartet does NOT check this. SEGFAULTS should occur if order
    // is not guaranteed.

    int s1, s2, s3, s4;
    int am1, am2, am3, am4, temp;
    std::shared_ptr<BasisSet> bs_temp;
    bool p13p24 = false, p12 = false, p34 = false;

    // AM used for ordering
    am1 = original_bs1_->shell(sh1).am();
    am2 = original_bs2_->shell(sh2).am();
    am3 = original_bs3_->shell(sh3).am();
    am4 = original_bs4_->shell(sh4).am();

    int n1, n2, n3, n4;
    n1 = original_bs1_->shell(sh1).ncartesian();
    n2 = original_bs2_->shell(sh2).ncartesian();
    n3 = original_bs3_->shell(sh3).ncartesian();
    n4 = original_bs4_->shell(sh4).ncartesian();

    // l(a) >= l(b), l(c) >= l(d), and l(c) + l(d) >= l(a) + l(b).
    if (am1 >= am2) {
        s1 = sh1;
        s2 = sh2;

        bs1_ = original_bs1_;
        bs2_ = original_bs2_;
    } else {
        s1 = sh2;
        s2 = sh1;

        bs1_ = original_bs2_;
        bs2_ = original_bs1_;

        p12 = true;
    }

    if (am3 >= am4) {
        s3 = sh3;
        s4 = sh4;

        bs3_ = original_bs3_;
        bs4_ = original_bs4_;

    } else {
        s3 = sh4;
        s4 = sh3;

        bs3_ = original_bs4_;
        bs4_ = original_bs3_;

        p34 = true;
    }

    if ((am1 + am2) > (am3 + am4)) {
        // Swap s1 and s2 with s3 and s4
        temp = s1;
        s1 = s3;
        s3 = temp;

        temp = s2;
        s2 = s4;
        s4 = temp;

        bs_temp = bs1_;
        bs1_ = bs3_;
        bs3_ = bs_temp;

        bs_temp = bs2_;
        bs2_ = bs4_;
        bs4_ = bs_temp;

        p13p24 = true;
    }

    if (p12) {
        if (p34) {
            if (p13p24) {
                // (AB|CD) -> (DC|BA)
                permuted_order_ = DCBA;
            } else {
                // (AB|CD) -> (BA|DC)
                permuted_order_ = BADC;
            }
        } else {
            if (p13p24) {
                // (AB|CD) -> (CD|BA)
                permuted_order_ = CDBA;
            } else {
                // (AB|CD) -> (BA|CD)
                permuted_order_ = BACD;
            }
        }
    } else {
        if (p34) {
            if (p13p24) {
                // (AB|CD) -> (DC|AB)
                permuted_order_ = DCAB;
            } else {
                // (AB|CD) -> (AB|DC)
                permuted_order_ = ABDC;
            }
        } else {
            if (p13p24) {
                // (AB|CD) -> (CD|AB)
                permuted_order_ = CDAB;
            } else {
                // (AB|CD) -> (AB|CD)
                permuted_order_ = ABCD;
            }
        }
    }

    // s1, s2, s3, s4 contain the shells to do in libderiv order
    compute_quartet_deriv1(s1, s2, s3, s4);  // compute 9 sets of integral derivatives

    // Need both sizes because source_ is in cartesians and target_ might be in pure am
    size_t size = n1 * n2 * n3 * n4;
    // Permute integrals back, if needed
    if (p12 || p34 || p13p24) {
        // ERI_1DER_NTYPE of them
        for (int i = 0; i < ERI_1DER_NTYPE; ++i)
            permute_target(source_ + (i * size), target_ + (i * size), s1, s2, s3, s4, p12, p34, p13p24);
    } else {
        // copy the integrals to the target_, 3n of them
        memcpy(target_, source_, ERI_1DER_NTYPE * size * sizeof(double));
    }
    return size;
}

size_t TwoElectronInt::compute_quartet_deriv1(int sh1, int sh2, int sh3, int sh4) {
    const GaussianShell &s1 = bs1_->shell(sh1);
    const GaussianShell &s2 = bs2_->shell(sh2);
    const GaussianShell &s3 = bs3_->shell(sh3);
    const GaussianShell &s4 = bs4_->shell(sh4);

    int am1 = s1.am();
    int am2 = s2.am();
    int am3 = s3.am();
    int am4 = s4.am();

    int am = am1 + am2 + am3 + am4;  // total am

    int nprim1 = s1.nprimitive();
    int nprim2 = s2.nprimitive();
    int nprim3 = s3.nprimitive();
    int nprim4 = s4.nprimitive();
    size_t nprim;

    double A[3], B[3], C[3], D[3];
    A[0] = s1.center()[0];
    A[1] = s1.center()[1];
    A[2] = s1.center()[2];

    B[0] = s2.center()[0];
    B[1] = s2.center()[1];
    B[2] = s2.center()[2];

    C[0] = s3.center()[0];
    C[1] = s3.center()[1];
    C[2] = s3.center()[2];

    D[0] = s4.center()[0];
    D[1] = s4.center()[1];
    D[2] = s4.center()[2];

    // Prefactor
    double prefactor = 1.0;

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    double CD2 = 0.0;
    CD2 += (C[0] - D[0]) * (C[0] - D[0]);
    CD2 += (C[1] - D[1]) * (C[1] - D[1]);
    CD2 += (C[2] - D[2]) * (C[2] - D[2]);

    libderiv_.AB[0] = A[0] - B[0];
    libderiv_.AB[1] = A[1] - B[1];
    libderiv_.AB[2] = A[2] - B[2];
    libderiv_.CD[0] = C[0] - D[0];
    libderiv_.CD[1] = C[1] - D[1];
    libderiv_.CD[2] = C[2] - D[2];

    // Prepare all the data needed by libderiv
    nprim = 0;

    if (use_shell_pairs_) {
        const L1ShellPair &sp12 = (*pairs12_)[sh1][sh2];
        const L1ShellPair &sp34 = (*pairs34_)[sh3][sh4];

        nprim = fill_primitive_data(libderiv_.PrimQuartet, fjt_, sp12, sp34, am, sh1 == sh2, sh3 == sh4, 1);
    } else {
        for (int p1 = 0; p1 < nprim1; ++p1) {
            double a1 = s1.exp(p1);
            double c1 = s1.coef(p1);
            for (int p2 = 0; p2 < nprim2; ++p2) {
                double a2 = s2.exp(p2);
                double c2 = s2.coef(p2);
                double zeta = a1 + a2;
                double ooz = 1.0 / zeta;
                double oo2z = 1.0 / (2.0 * zeta);

                double PA[3], PB[3];
                double P[3];

                P[0] = (a1 * A[0] + a2 * B[0]) * ooz;
                P[1] = (a1 * A[1] + a2 * B[1]) * ooz;
                P[2] = (a1 * A[2] + a2 * B[2]) * ooz;
                PA[0] = P[0] - A[0];
                PA[1] = P[1] - A[1];
                PA[2] = P[2] - A[2];
                PB[0] = P[0] - B[0];
                PB[1] = P[1] - B[1];
                PB[2] = P[2] - B[2];

                double Sab = pow(M_PI * ooz, 3.0 / 2.0) * exp(-a1 * a2 * ooz * AB2) * c1 * c2;

                for (int p3 = 0; p3 < nprim3; ++p3) {
                    double a3 = s3.exp(p3);
                    double c3 = s3.coef(p3);
                    for (int p4 = 0; p4 < nprim4; ++p4) {
                        double a4 = s4.exp(p4);
                        double c4 = s4.coef(p4);
                        double eta = a3 + a4;
                        double ooze = 1.0 / (zeta + eta);
                        double poz = eta * ooze;
                        double rho = zeta * poz;
                        double oon = 1.0 / eta;

                        double Scd = pow(M_PI * oon, 3.0 / 2.0) * exp(-a3 * a4 * oon * CD2) * c3 * c4;

                        double coef1 = 2.0 * sqrt(rho * M_1_PI) * Sab * Scd;

                        libderiv_.PrimQuartet[nprim].poz = poz;
                        libderiv_.PrimQuartet[nprim].oo2zn = 0.5 * ooze;
                        libderiv_.PrimQuartet[nprim].pon = zeta * ooze;
                        libderiv_.PrimQuartet[nprim].oo2z = 0.5 / zeta;
                        libderiv_.PrimQuartet[nprim].oo2n = 0.5 / eta;
                        libderiv_.PrimQuartet[nprim].twozeta_a = 2.0 * a1;
                        libderiv_.PrimQuartet[nprim].twozeta_b = 2.0 * a2;
                        libderiv_.PrimQuartet[nprim].twozeta_c = 2.0 * a3;
                        libderiv_.PrimQuartet[nprim].twozeta_d = 2.0 * a4;

                        // double oo2n = 1.0 / (2.0 * eta);
                        // double oo2zn = 1.0 / (2.0 * (zeta + eta));
                        // double rho = (zeta * eta) / (zeta + eta);
                        // double oo2rho = 1.0 / (2.0 * rho);

                        double QC[3], QD[3], WP[3], WQ[3], PQ[3];
                        double Q[3], W[3];

                        Q[0] = (a3 * C[0] + a4 * D[0]) * oon;
                        Q[1] = (a3 * C[1] + a4 * D[1]) * oon;
                        Q[2] = (a3 * C[2] + a4 * D[2]) * oon;
                        QC[0] = Q[0] - C[0];
                        QC[1] = Q[1] - C[1];
                        QC[2] = Q[2] - C[2];
                        QD[0] = Q[0] - D[0];
                        QD[1] = Q[1] - D[1];
                        QD[2] = Q[2] - D[2];
                        PQ[0] = P[0] - Q[0];
                        PQ[1] = P[1] - Q[1];
                        PQ[2] = P[2] - Q[2];

                        double PQ2 = 0.0;
                        PQ2 += (P[0] - Q[0]) * (P[0] - Q[0]);
                        PQ2 += (P[1] - Q[1]) * (P[1] - Q[1]);
                        PQ2 += (P[2] - Q[2]) * (P[2] - Q[2]);

                        W[0] = (zeta * P[0] + eta * Q[0]) * ooze;
                        W[1] = (zeta * P[1] + eta * Q[1]) * ooze;
                        W[2] = (zeta * P[2] + eta * Q[2]) * ooze;
                        WP[0] = W[0] - P[0];
                        WP[1] = W[1] - P[1];
                        WP[2] = W[2] - P[2];
                        WQ[0] = W[0] - Q[0];
                        WQ[1] = W[1] - Q[1];
                        WQ[2] = W[2] - Q[2];

                        for (int i = 0; i < 3; ++i) {
                            libderiv_.PrimQuartet[nprim].U[0][i] = PA[i];
                            libderiv_.PrimQuartet[nprim].U[1][i] = PB[i];
                            libderiv_.PrimQuartet[nprim].U[2][i] = QC[i];
                            libderiv_.PrimQuartet[nprim].U[3][i] = QD[i];
                            libderiv_.PrimQuartet[nprim].U[4][i] = WP[i];
                            libderiv_.PrimQuartet[nprim].U[5][i] = WQ[i];
                        }

                        double T = rho * PQ2;
                        fjt_->set_rho(rho);
                        double *F = fjt_->values(am + 1, T);

                        // Modify F to include overlap of ab and cd, eqs 14, 15, 16 of libint manual
                        // double val = 2.0 * sqrt(rho * M_1_PI) * Sab * Scd * prefactor;

                        for (int i = 0; i <= am + 1; ++i) {
                            libderiv_.PrimQuartet[nprim].F[i] = F[i] * coef1;
                        }

                        nprim++;
                    }
                }
            }
        }
    }

    // How many are there?
    size_t size = INT_NCART(am1) * INT_NCART(am2) * INT_NCART(am3) * INT_NCART(am4);

    // Compute the integral
    build_deriv1_eri[am1][am2][am3][am4](&libderiv_, nprim);

    // Zero out memory
    memset(source_, 0, sizeof(double) * size * ERI_1DER_NTYPE);

    // Copy results from libderiv into source_ (note libderiv only gives 3 of the centers).
    // The libmints array returns the following integral derivatives:
    //   0 -> A_x
    //   1 -> A_y
    //   2 -> A_z
    //   3 -> C_x
    //   4 -> C_y
    //   5 -> C_z
    //   6 -> D_x
    //   7 -> D_y
    //   8 -> D_z
    // Center B can be determined by:
    //   B_x = -(A_x + C_x + D_x)
    //   B_y = -(A_y + C_y + D_y)
    //   B_z = -(A_z + C_z + D_z)

    handle_reordering1(permuted_order_, libderiv_, source_, size);
    // Transform the integrals to the spherical basis
    pure_transform(sh1, sh2, sh3, sh4, ERI_1DER_NTYPE);

    // Results are in source_
    return size;
}

size_t TwoElectronInt::compute_shell_deriv2(int sh1, int sh2, int sh3, int sh4) {
    if (deriv_ < 2) throw PSIEXCEPTION("ERROR - ERI: ERI object not initialized to handle second derivatives.\n");

    // Need to ensure the ordering asked by the user is valid for libderiv.
    // compute_quartet_deriv2 does NOT check this. SEGFAULTS will likely occur
    // if order is not guarantee.
    int s1, s2, s3, s4;

    int am1, am2, am3, am4, temp;
    std::shared_ptr<BasisSet> bs_temp;
    bool p13p24 = false, p12 = false, p34 = false;

    // AM used for ordering
    am1 = original_bs1_->shell(sh1).am();
    am2 = original_bs2_->shell(sh2).am();
    am3 = original_bs3_->shell(sh3).am();
    am4 = original_bs4_->shell(sh4).am();

    int n1, n2, n3, n4;
    n1 = original_bs1_->shell(sh1).ncartesian();
    n2 = original_bs2_->shell(sh2).ncartesian();
    n3 = original_bs3_->shell(sh3).ncartesian();
    n4 = original_bs4_->shell(sh4).ncartesian();

    // am1 >= am2, am3 >= am4, and am3 + am3 >= am1 + am2
    if (am1 >= am2) {
        s1 = sh1;
        s2 = sh2;

        bs1_ = original_bs1_;
        bs2_ = original_bs2_;
    } else {
        s1 = sh2;
        s2 = sh1;

        bs1_ = original_bs2_;
        bs2_ = original_bs1_;

        p12 = true;
    }

    if (am3 >= am4) {
        s3 = sh3;
        s4 = sh4;

        bs3_ = original_bs3_;
        bs4_ = original_bs4_;
    } else {
        s3 = sh4;
        s4 = sh3;

        bs3_ = original_bs4_;
        bs4_ = original_bs3_;

        p34 = true;
    }

    if ((am1 + am2) > (am3 + am4)) {
        // swap s1 and s2 with s3 and s4.
        temp = s1;
        s1 = s3;
        s3 = temp;

        temp = s2;
        s2 = s4;
        s4 = temp;

        bs_temp = bs1_;
        bs1_ = bs3_;
        bs3_ = bs_temp;

        bs_temp = bs2_;
        bs2_ = bs4_;
        bs4_ = bs_temp;

        p13p24 = true;
    }

    if (p12) {
        if (p34) {
            if (p13p24) {
                // (AB|CD) -> (DC|BA)
                permuted_order_ = DCBA;
            } else {
                // (AB|CD) -> (BA|DC)
                permuted_order_ = BADC;
            }
        } else {
            if (p13p24) {
                // (AB|CD) -> (CD|BA)
                permuted_order_ = CDBA;
            } else {
                // (AB|CD) -> (BA|CD)
                permuted_order_ = BACD;
            }
        }
    } else {
        if (p34) {
            if (p13p24) {
                // (AB|CD) -> (DC|AB)
                permuted_order_ = DCAB;
            } else {
                // (AB|CD) -> (AB|DC)
                permuted_order_ = ABDC;
            }
        } else {
            if (p13p24) {
                // (AB|CD) -> (CD|AB)
                permuted_order_ = CDAB;
            } else {
                // (AB|CD) -> (AB|CD)
                permuted_order_ = ABCD;
            }
        }
    }

    // This handles computing the
    compute_quartet_deriv2(s1, s2, s3, s4);

    size_t size = n1 * n2 * n3 * n4;
    if (p12 || p34 || p13p24) {
        for (int i = 0; i < ERI_2DER_NTYPE; ++i)
            permute_target(source_ + (i * size), target_ + (i * size), s1, s2, s3, s4, p12, p34, p13p24);
    } else {
        memcpy(target_, source_, ERI_2DER_NTYPE * size * sizeof(double));
    }
    return size;
}

size_t TwoElectronInt::compute_quartet_deriv2(int sh1, int sh2, int sh3, int sh4) {
    const GaussianShell &s1 = bs1_->shell(sh1);
    const GaussianShell &s2 = bs2_->shell(sh2);
    const GaussianShell &s3 = bs3_->shell(sh3);
    const GaussianShell &s4 = bs4_->shell(sh4);

    int am1 = s1.am();
    int am2 = s2.am();
    int am3 = s3.am();
    int am4 = s4.am();
    int am = am1 + am2 + am3 + am4;

    int nprim1 = s1.nprimitive();
    int nprim2 = s2.nprimitive();
    int nprim3 = s3.nprimitive();
    int nprim4 = s4.nprimitive();
    size_t nprim = 0;

    double A[3], B[3], C[3], D[3];
    A[0] = s1.center()[0];
    A[1] = s1.center()[1];
    A[2] = s1.center()[2];

    B[0] = s2.center()[0];
    B[1] = s2.center()[1];
    B[2] = s2.center()[2];

    C[0] = s3.center()[0];
    C[1] = s3.center()[1];
    C[2] = s3.center()[2];

    D[0] = s4.center()[0];
    D[1] = s4.center()[1];
    D[2] = s4.center()[2];

    // prefactor
    double prefactor = 1.0;

    // compute intermediates
    double AB2 = 0.0;
    AB2 += (A[0] - B[0]) * (A[0] - B[0]);
    AB2 += (A[1] - B[1]) * (A[1] - B[1]);
    AB2 += (A[2] - B[2]) * (A[2] - B[2]);

    double CD2 = 0.0;
    CD2 += (C[0] - D[0]) * (C[0] - D[0]);
    CD2 += (C[1] - D[1]) * (C[1] - D[1]);
    CD2 += (C[2] - D[2]) * (C[2] - D[2]);

    libderiv_.AB[0] = A[0] - B[0];
    libderiv_.AB[1] = A[1] - B[1];
    libderiv_.AB[2] = A[2] - B[2];
    libderiv_.CD[0] = C[0] - D[0];
    libderiv_.CD[1] = C[1] - D[1];
    libderiv_.CD[2] = C[2] - D[2];

    // prepare all the data needed for libderiv
    if (use_shell_pairs_) {
        // ShellPair *p12, *p34;
        // p12 = &(pairs12_[sh1][sh2]);
        // p34 = &(pairs34_[sh3][sh4]);
        const L1ShellPair &sp12 = (*pairs12_)[sh1][sh2];
        const L1ShellPair &sp34 = (*pairs34_)[sh3][sh4];

        nprim = fill_primitive_data(libderiv_.PrimQuartet, fjt_, sp12, sp34, am, sh1 == sh2, sh3 == sh4, 2);
    } else {
        for (int p1 = 0; p1 < nprim1; ++p1) {
            double a1 = s1.exp(p1);
            double c1 = s1.coef(p1);
            for (int p2 = 0; p2 < nprim2; ++p2) {
                double a2 = s2.exp(p2);
                double c2 = s2.coef(p2);
                double zeta = a1 + a2;
                double ooz = 1.0 / zeta;
                double oo2z = 1.0 / (2.0 * zeta);

                double PA[3], PB[3];
                double P[3];

                P[0] = (a1 * A[0] + a2 * B[0]) * ooz;
                P[1] = (a1 * A[1] + a2 * B[1]) * ooz;
                P[2] = (a1 * A[2] + a2 * B[2]) * ooz;
                PA[0] = P[0] - A[0];
                PA[1] = P[1] - A[1];
                PA[2] = P[2] - A[2];
                PB[0] = P[0] - B[0];
                PB[1] = P[1] - B[1];
                PB[2] = P[2] - B[2];

                double Sab = pow(M_PI * ooz, 3.0 / 2.0) * exp(-a1 * a2 * ooz * AB2) * c1 * c2;

                for (int p3 = 0; p3 < nprim3; ++p3) {
                    double a3 = s3.exp(p3);
                    double c3 = s3.coef(p3);
                    for (int p4 = 0; p4 < nprim4; ++p4) {
                        double a4 = s4.exp(p4);
                        double c4 = s4.coef(p4);
                        double nu = a3 + a4;
                        double oon = 1.0 / nu;
                        double oo2n = 1.0 / (2.0 * nu);
                        double oo2zn = 1.0 / (2.0 * (zeta + nu));
                        double rho = (zeta * nu) / (zeta + nu);

                        double QC[3], QD[3], WP[3], WQ[3], PQ[3];
                        double Q[3], W[3];

                        Q[0] = (a3 * C[0] + a4 * D[0]) * oon;
                        Q[1] = (a3 * C[1] + a4 * D[1]) * oon;
                        Q[2] = (a3 * C[2] + a4 * D[2]) * oon;
                        QC[0] = Q[0] - C[0];
                        QC[1] = Q[1] - C[1];
                        QC[2] = Q[2] - C[2];
                        QD[0] = Q[0] - D[0];
                        QD[1] = Q[1] - D[1];
                        QD[2] = Q[2] - D[2];
                        PQ[0] = P[0] - Q[0];
                        PQ[1] = P[1] - Q[1];
                        PQ[2] = P[2] - Q[2];

                        double PQ2 = 0.0;
                        PQ2 += (P[0] - Q[0]) * (P[0] - Q[0]);
                        PQ2 += (P[1] - Q[1]) * (P[1] - Q[1]);
                        PQ2 += (P[2] - Q[2]) * (P[2] - Q[2]);

                        W[0] = (zeta * P[0] + nu * Q[0]) / (zeta + nu);
                        W[1] = (zeta * P[1] + nu * Q[1]) / (zeta + nu);
                        W[2] = (zeta * P[2] + nu * Q[2]) / (zeta + nu);
                        WP[0] = W[0] - P[0];
                        WP[1] = W[1] - P[1];
                        WP[2] = W[2] - P[2];
                        WQ[0] = W[0] - Q[0];
                        WQ[1] = W[1] - Q[1];
                        WQ[2] = W[2] - Q[2];

                        for (int i = 0; i < 3; ++i) {
                            libderiv_.PrimQuartet[nprim].U[0][i] = PA[i];
                            libderiv_.PrimQuartet[nprim].U[1][i] = PB[i];
                            libderiv_.PrimQuartet[nprim].U[2][i] = QC[i];
                            libderiv_.PrimQuartet[nprim].U[3][i] = QD[i];
                            libderiv_.PrimQuartet[nprim].U[4][i] = WP[i];
                            libderiv_.PrimQuartet[nprim].U[5][i] = WQ[i];
                        }
                        libderiv_.PrimQuartet[nprim].oo2z = oo2z;
                        libderiv_.PrimQuartet[nprim].oo2n = oo2n;
                        libderiv_.PrimQuartet[nprim].oo2zn = oo2zn;
                        libderiv_.PrimQuartet[nprim].poz = rho * ooz;
                        libderiv_.PrimQuartet[nprim].pon = rho * oon;
                        // libderiv_.PrimQuartet[nprim].oo2p = oo2rho;   // NOT SET IN CINTS
                        libderiv_.PrimQuartet[nprim].twozeta_a = 2.0 * a1;
                        libderiv_.PrimQuartet[nprim].twozeta_b = 2.0 * a2;
                        libderiv_.PrimQuartet[nprim].twozeta_c = 2.0 * a3;
                        libderiv_.PrimQuartet[nprim].twozeta_d = 2.0 * a4;

                        double T = rho * PQ2;
                        fjt_->set_rho(rho);
                        double *F = fjt_->values(am + 2, T);

                        // Modify F to include overlap of ab and cd, eqs 14, 15, 16 of libint manual
                        double Scd = pow(M_PI * oon, 3.0 / 2.0) * exp(-a3 * a4 * oon * CD2) * c3 * c4;
                        double val = 2.0 * sqrt(rho * M_1_PI) * Sab * Scd * prefactor;

                        for (int i = 0; i <= am + 2; ++i) {
                            libderiv_.PrimQuartet[nprim].F[i] = F[i] * val;
                        }

                        nprim++;
                    }
                }
            }
        }
    }

    size_t size = INT_NCART(am1) * INT_NCART(am2) * INT_NCART(am3) * INT_NCART(am4);
    build_deriv12_eri[am1][am2][am3][am4](&libderiv_, nprim);

    // zero out the memory
    memset(source_, 0, sizeof(double) * size * ERI_2DER_NTYPE);

    // Copy results from libderiv into source_ (note libderiv only gives 3 of the centers)
    handle_reordering12(permuted_order_, libderiv_, source_, size);

    // Transform the integrals to the spherical basis
    pure_transform(sh1, sh2, sh3, sh4, ERI_2DER_NTYPE);

    // Results are in source_
    return size;
}
#endif  // ENABLE_Libint1t

///// Libint2 implementation

Libint2TwoElectronInt::Libint2TwoElectronInt(const IntegralFactory *integral, int deriv, double screening_threshold,
                                             bool use_shell_pairs, bool needs_exchange)
    : TwoBodyAOInt(integral, deriv), use_shell_pairs_(use_shell_pairs) {
    // Initialize libint static data
#if LIBINT2_SHGAUSS_ORDERING == LIBINT_SHGSHELL_ORDERING_STANDARD
    libint2::initialize(libint2::SHGShellOrdering_Standard);
#elif LIBINT2_SHGAUSS_ORDERING == LIBINT_SHGSHELL_ORDERING_GAUSSIAN
    libint2::initialize(libint2::SHGShellOrdering_Gaussian);
#else
#  error "unknown value of macro LIBINT2_SHGAUSS_ORDERING"
#endif

    // Make sure there's enough space for the sieve generation.  This array is used to return an array of
    // zeros back to the caller if libint2 gave us nullptr, so the caller doesn't have to check.
    size_t sieve_size =  std::max(
                           basis1()->max_function_per_shell() * basis2()->max_function_per_shell() *
                           basis1()->max_function_per_shell() * basis2()->max_function_per_shell(),
                           basis3()->max_function_per_shell() * basis4()->max_function_per_shell() *
                           basis3()->max_function_per_shell() * basis4()->max_function_per_shell());
    size_t size = std::max((size_t) basis1()->max_function_per_shell() * basis2()->max_function_per_shell() *
                           basis3()->max_function_per_shell() * basis4()->max_function_per_shell(), sieve_size);
    zero_vec_ = std::vector<double>(size, 0.0);
}

Libint2TwoElectronInt::Libint2TwoElectronInt(const Libint2TwoElectronInt &rhs)
    : TwoBodyAOInt(rhs), schwarz_engine_(rhs.schwarz_engine_), braket_(rhs.braket_), use_shell_pairs_(rhs.use_shell_pairs_)
{
    pairs12_ = rhs.pairs12_;
    pairs34_ = rhs.pairs34_;
    zero_vec_ = rhs.zero_vec_;
    for (const auto &e : rhs.engines_) engines_.emplace_back(e);
}

void Libint2TwoElectronInt::common_init() {
    bool dummy1 = basis1()->l2_shell(0) == libint2::Shell::unit();
    bool dummy2 = basis2()->l2_shell(0) == libint2::Shell::unit();
    bool dummy3 = basis3()->l2_shell(0) == libint2::Shell::unit();
    bool dummy4 = basis4()->l2_shell(0) == libint2::Shell::unit();

    if (!dummy1 && !dummy2 && !dummy3 && !dummy4) {
        braket_ = libint2::BraKet::xx_xx;
    } else if (!dummy1 && dummy2 && !dummy3 && !dummy4) {
        braket_ = libint2::BraKet::xs_xx;
    } else if (!dummy1 && !dummy2 && !dummy3 && dummy4) {
        braket_ = libint2::BraKet::xx_xs;
    } else if (!dummy1 && dummy2 && !dummy3 && dummy4) {
        braket_ = libint2::BraKet::xs_xs;
    } else {
        throw PSIEXCEPTION("Bad BraKet type in Libint2TwoElectronInt");
    }

    for (auto &engine : engines_) engine.set(braket_);

    int num_chunks;
    switch (deriv_) {
        case 0:
            num_chunks = 1;
            break;
        case 1:
            num_chunks = 12;
            break;
        case 2:
            num_chunks = 78;
            break;
        default:
            throw PSIEXCEPTION("Libint2 engine only supports up to second derivatives currently.");
    }
    buffers_.resize(num_chunks);

    target_full_ = const_cast<double *>(engines_[0].results()[0]);
    target_ = target_full_;

    // Make sure the engine can handle the type of integral used to build a sieve
    setup_sieve();
    // Reset the engine type back to the general case needed
    create_blocks();
    const auto max_engine_precision = std::numeric_limits<double>::epsilon() * screening_threshold_;

    size_t npairs = shell_pairs_bra_.size();
    pairs12_.resize(npairs);
//#pragma omp parallel for
    for (int pair = 0; pair < npairs; ++pair) {
        auto s1 = shell_pairs_bra_[pair].first;
        auto s2 = shell_pairs_bra_[pair].second;
        pairs12_[pair] = std::make_shared<libint2::ShellPair>(basis1()->l2_shell(s1), basis2()->l2_shell(s2),
                                                              std::log(max_engine_precision));
    }
    npairs = shell_pairs_ket_.size();
    pairs34_.resize(npairs);
//#pragma omp parallel for
    for (int pair = 0; pair < npairs; ++pair) {
        auto s3 = shell_pairs_ket_[pair].first;
        auto s4 = shell_pairs_ket_[pair].second;
        pairs34_[pair] = std::make_shared<libint2::ShellPair>(basis3()->l2_shell(s3), basis4()->l2_shell(s4),
                                                              std::log(max_engine_precision));
    }
}

Libint2TwoElectronInt::~Libint2TwoElectronInt() { libint2::finalize(); }

size_t Libint2TwoElectronInt::compute_shell(const AOShellCombinationsIterator &shellIter) {
    return compute_shell(shellIter.p(), shellIter.q(), shellIter.r(), shellIter.s());
}

size_t Libint2TwoElectronInt::compute_shell_for_sieve(const std::shared_ptr<BasisSet> bs, int s1, int s2, int s3,
                                                      int s4, bool is_bra) {
#ifdef MINTS_TIMER
    timer_on("Libint2ERI::compute_shell_for_sieve");
#endif
    const auto &sh1 = bs->l2_shell(s1);
    const auto &sh2 = bs->l2_shell(s2);
    const auto &sh3 = bs->l2_shell(s3);
    const auto &sh4 = bs->l2_shell(s4);

    schwarz_engine_.compute(sh1, sh2, sh3, sh4);

    size_t ntot = sh1.size() * sh2.size() * sh3.size() * sh4.size();
    buffers_[0] = target_full_ = const_cast<double *>(schwarz_engine_.results()[0]);
    if (target_full_ == nullptr) {
        // The caller will try to read the buffer if there isn't a check on the number of ints computed
        // so we point to a valid array of zeros here to prevent memory bugs in the calling routine.
        buffers_[0] = target_full_ = zero_vec_.data();
        ntot = 0;
    }

#ifdef MINTS_TIMER
    timer_off("Libint2ERI::compute_shell_for_sieve");
#endif
    return ntot;
}

size_t Libint2TwoElectronInt::compute_shell(int s1, int s2, int s3, int s4) {
#ifdef MINTS_TIMER
    timer_on("Libint2ERI::compute_shell");
#endif

    const auto &sh1 = bs1_->l2_shell(s1);
    const auto &sh2 = bs2_->l2_shell(s2);
    const auto &sh3 = bs3_->l2_shell(s3);
    const auto &sh4 = bs4_->l2_shell(s4);

    libint2_wrapper0(sh1, sh2, sh3, sh4);

    size_t ntot = sh1.size() * sh2.size() * sh3.size() * sh4.size();

    buffers_[0] = target_full_ = const_cast<double *>(engines_[0].results()[0]);
    if (target_full_ == nullptr) {
        // The caller will try to read the buffer if there isn't a check on the number of ints computed
        // so we point to a valid array of zeros here to prevent memory bugs in the calling routine.
        buffers_[0] = target_full_ = zero_vec_.data();
        ntot = 0;
    }

#ifdef MINTS_TIMER
    timer_off("Libint2ERI::compute_shell");
#endif
    return ntot;
}

size_t Libint2TwoElectronInt::compute_shell_deriv1(int s1, int s2, int s3, int s4) {
#ifdef MINTS_TIMER
    timer_on("Libint2ERI::compute_shell_deriv1");
#endif

    const auto &sh1 = bs1_->l2_shell(s1);
    const auto &sh2 = bs2_->l2_shell(s2);
    const auto &sh3 = bs3_->l2_shell(s3);
    const auto &sh4 = bs4_->l2_shell(s4);

    libint2_wrapper1(sh1, sh2, sh3, sh4);


    size_t ntot = 0;
    bool none_computed = engines_[1].results()[0] == nullptr;
    if (none_computed) {
        for (int i = 0; i < 12; ++i) {
            buffers_[i] = zero_vec_.data();
        }
    } else {
        for (int i = 0; i < 12; ++i) {
            buffers_[i] = engines_[1].results()[i];
        }
        ntot = 12 * sh1.size() * sh2.size() * sh3.size() * sh4.size();
    }

#ifdef MINTS_TIMER
    timer_off("Libint2ERI::compute_shell_deriv1");
#endif
    return ntot;
}

size_t Libint2TwoElectronInt::compute_shell_deriv2(int s1, int s2, int s3, int s4) {
#ifdef MINTS_TIMER
    timer_on("Libint2ERI::compute_shell_deriv2");
#endif

    const auto &sh1 = bs1_->l2_shell(s1);
    const auto &sh2 = bs2_->l2_shell(s2);
    const auto &sh3 = bs3_->l2_shell(s3);
    const auto &sh4 = bs4_->l2_shell(s4);

    libint2_wrapper2(sh1, sh2, sh3, sh4);

    size_t ntot = 0;
    bool none_computed = engines_[2].results()[0] == nullptr;
    if (none_computed) {
        for (int i = 0; i < 78; ++i) {
            buffers_[i] = zero_vec_.data();
        }
    } else {
        for (int i = 0; i < 78; ++i) {
            buffers_[i] = engines_[2].results()[i];
        }
        ntot = 78 * sh1.size() * sh2.size() * sh3.size() * sh4.size();
    }

#ifdef MINTS_TIMER
    timer_off("Libint2ERI::compute_shell_deriv2");
#endif
    return ntot;
}

void Libint2TwoElectronInt::compute_shell_blocks(int shellpair12, int shellpair34, int npair12, int npair34) {
    if (npair12 != -1 || npair34 != -1)
        throw PSIEXCEPTION("npair12 and npair34 arguments are not supported by the Libint2 engine.");
#ifdef MINTS_TIMER
    timer_on("Libint2ERI::compute_shell_blocks");
#endif
    // This engine doesn't block shells, so each "block" is just 1 shell
    int s1 = blocks12_[shellpair12][0].first;
    int s2 = blocks12_[shellpair12][0].second;
    int s3 = blocks34_[shellpair34][0].first;
    int s4 = blocks34_[shellpair34][0].second;

    const auto &sh1 = bs1_->l2_shell(s1);
    const auto &sh2 = bs2_->l2_shell(s2);
    const auto &sh3 = bs3_->l2_shell(s3);
    const auto &sh4 = bs4_->l2_shell(s4);

    size_t ntot = 0;

    const auto *sp12 = pairs12_[shellpair12].get();
    const auto *sp34 = pairs34_[shellpair34].get();
    libint2_wrapper0(sh1, sh2, sh3, sh4, sp12, sp34);

    target_full_ = const_cast<double *>(engines_[0].results()[0]);
    if (target_full_) {
        buffers_[0] = engines_[0].results()[0];
        ntot = sh1.size() * sh2.size() * sh3.size() * sh4.size();
    } else {
        target_full_ = zero_vec_.data();
    }

#ifdef MINTS_TIMER
    timer_off("Libint2ERI::compute_shell_blocks");
#endif
}
