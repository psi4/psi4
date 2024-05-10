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

#include <cmath>
#include <algorithm>

#include "psi4/libmoinfo/libmoinfo.h"
#include "psi4/libpsi4util/libpsi4util.h"
#include "psi4/libmints/mintshelper.h"
#include "psi4/libmints/matrix.h"

#include "algebra_interface.h"
#include "blas.h"
#include "transform.h"
#include "matrix.h"

#define CCTRANSFORM_USE_BLAS

#define MAX(i, j) ((i > j) ? i : j)
#define MIN(i, j) ((i > j) ? j : i)
#define INDEX(i, j) ((i > j) ? (ioff[(i)] + (j)) : (ioff[(j)] + (i)))
#define four(i, j, k, l) INDEX(INDEX(i, j), INDEX(k, l))

#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libiwl/iwl.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"

namespace psi {
namespace psimrcc {

/*!
    \fn CCTransform::read_oei_integrals()
 */
void CCTransform::read_oei_so_integrals() {
    int nso = wfn_->moinfo()->get_nso();
    // Read all the (frozen + non-frozen) OEI in Pitzer order
    allocate_oei_so();

    // Read the kinetic energy integrals
    auto T = wfn_->mintshelper()->so_kinetic()->clone();
    auto V = wfn_->mintshelper()->so_potential()->clone();
    T->add(V);

    double** Tbmat = T->to_block_matrix();

    for (int i = 0; i < nso; i++)
        for (int j = 0; j < nso; j++) oei_so[i][j] += Tbmat[i][j];
    free_block(Tbmat);
}

}  // namespace psimrcc
}  // namespace psi
