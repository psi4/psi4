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

/*! \file
    \ingroup CCDENSITY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include "psi4/libiwl/iwl.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccdensity {

/* add_ref_UHF(): This function adds the reference contributions to
** the one- and two-particle density matrices.  These contributions
** are simply the prefactors in front of the one- and two-electron
** intgegrals, respectively, in the UHF energy expression.  In the
** case of the two-pdm, however, care must be taken that only the
** permutationally unique elements be written to disk.
**
** In spin-orbitals with Mulliken-order two-electron integrals, the
** three spin contributions to the SCF energy are:
**
** E_AA(SCF) = sum_I h_II + 1/2 sum_IJ ([II|JJ] - [IJ|IJ])
** E_BB(SCF) = sum_i h_ii + 1/2 sum_ij ([ii|jj] - [ij|ij])
** E_AB(SCF) = 1/2 sum_Ij [II|jj]
**
** I use QT-standard ordering for the indices in these expressions.
*/

void add_ref_UHF(struct iwlbuf *AA, struct iwlbuf *BB, struct iwlbuf *AB) {
    int mo_offset = 0;
    for (int h = 0; h < moinfo.nirreps; h++) {
        auto clsd_h = moinfo.frdocc[h] + moinfo.clsdpi[h];
        // Alpha
        for (int i = 0; i < moinfo.frdocc[h] + moinfo.aoccpi[h]; i++) {
            // One electron alpha
            moinfo.opdm_a.add(h, i, i, 1);
            // Two electron alpha-alpha
            auto qt_i = moinfo.pitzer2qt[i + mo_offset];
            for (int qt_j = 0; qt_j < qt_i; qt_j++) {
                iwl_buf_wrt_val(AA, qt_i, qt_i, qt_j, qt_j, 0.5, 0, "outfile", 0);
                iwl_buf_wrt_val(AA, qt_i, qt_j, qt_i, qt_j, -0.25, 0, "outfile", 0);
                iwl_buf_wrt_val(AA, qt_j, qt_i, qt_j, qt_i, -0.25, 0, "outfile", 0);
                iwl_buf_wrt_val(AA, qt_i, qt_j, qt_j, qt_i, -0.25, 0, "outfile", 0);
            }
        }
        // Beta
        for (int i = 0; i < moinfo.frdocc[h] + moinfo.boccpi[h]; i++) {
            // One electron beta
            moinfo.opdm_b.add(h, i, i, 1);
            // Two electron beta-beta
            auto qt_i = moinfo.pitzer2qt[i + mo_offset];
            for (int qt_j = 0; qt_j < qt_i; qt_j++) {
                iwl_buf_wrt_val(BB, qt_i, qt_i, qt_j, qt_j, 0.5, 0, "outfile", 0);
                iwl_buf_wrt_val(BB, qt_i, qt_j, qt_i, qt_j, -0.25, 0, "outfile", 0);
                iwl_buf_wrt_val(BB, qt_j, qt_i, qt_j, qt_i, -0.25, 0, "outfile", 0);
                iwl_buf_wrt_val(BB, qt_i, qt_j, qt_j, qt_i, -0.25, 0, "outfile", 0);
            }
        }
        mo_offset += moinfo.orbspi[h];
    }

    // Two-electron alpha-beta
    for (int i = 0; i < (moinfo.nfzc + moinfo.nclsd + moinfo.nopen); i++)
        for (int j = 0; j < (moinfo.nfzc + moinfo.nclsd); j++) iwl_buf_wrt_val(AB, i, i, j, j, 1.0, 0, "outfile", 0);
}
}  // namespace ccdensity
}  // namespace psi
