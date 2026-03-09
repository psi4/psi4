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

#include "integraltransform.h"

#include <algorithm>
#include <numeric>
#include <string>
#include <vector>

#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"

using namespace psi;

/**
 * Transforms a packed symmetric matrix.
 *
 * @param m - input matrix row dimension
 * @param n - output matrix row dimension
 * @param input - pointer to input integrals (the lower-triangle of a symmetric matrix)
 * @param pointer to output integrals (the lower-triangle of a symmetric matrix)
 * @param C transformation matrix (rectangular, m X n)
 * @param offset - the point in the full list of SOs where we want to start.  This is
 *                 useful for transforming integrals one irrep at a time and in this
 *                 case the offset would correspond to the position of the first
 *                 orbital in the current irrep.
 * @param order - a reordering array to change the order of the output
 * @param backtransform - whether this is a forward or backwards transformation
 * @param scale - the amount of the existing output buffer to mix into the result
 */

void IntegralTransform::trans_one(int m, int n, double *input, double *output, double **C, int offset, int *order,
                                  bool backtransform, double scale) {
    // TODO the order argument is actually not used right now.  I don't know that anybody will need it
    // so I haven't bothered so far...
    int dim = std::max(m, n);
    auto TMP0 = block_matrix(dim, dim);
    auto TMP1 = block_matrix(dim, dim);

    for (int p = 0; p < m; ++p) {
        for (int q = 0; q <= p; ++q) {
            size_t pq = INDEX((p + offset), (q + offset));
            TMP0[p][q] = TMP0[q][p] = input[pq];
        }
    }
    int nc;
    if (backtransform) {
        nc = m;
        if (m && n) {
            C_DGEMM('n', 't', m, n, m, 1.0, TMP0[0], dim, C[0], nc, 0.0, TMP1[0], dim);
            C_DGEMM('n', 'n', n, n, m, 1.0, C[0], nc, TMP1[0], dim, 0.0, TMP0[0], dim);
        }
    } else {
        nc = n;
        if (m && n) {
            C_DGEMM('n', 'n', m, n, m, 1.0, TMP0[0], dim, C[0], nc, 0.0, TMP1[0], dim);
            C_DGEMM('t', 'n', n, n, m, 1.0, C[0], nc, TMP1[0], dim, 0.0, TMP0[0], dim);
        }
    }

    for (int p = 0; p < nc; ++p) {
        for (int q = 0; q <= p; ++q) {
            size_t P = order[p];
            size_t Q = order[q];
            size_t PQ = INDEX(P, Q);
            output[PQ] = scale * output[PQ] + TMP0[p][q];
        }
    }

    free_block(TMP0);
    free_block(TMP1);

    return;
}
