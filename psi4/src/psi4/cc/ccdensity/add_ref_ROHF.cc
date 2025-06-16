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

void add_ref_ROHF(struct iwlbuf *OutBuf) {
    int mo_offset = 0;
    for (int h = 0; h < moinfo.nirreps; h++) {
        auto clsd_h = moinfo.frdocc[h] + moinfo.clsdpi[h];
        for (int i = 0; i < moinfo.frdocc[h] + moinfo.clsdpi[h]; i++) {
            // One electron closed-shell
            moinfo.opdm.add(h, i, i, 2.0);
            moinfo.opdm_a.add(h, i, i, 1.0);
            moinfo.opdm_b.add(h, i, i, 1.0);
            // Two electron closed-shell
            auto qt_i = moinfo.pitzer2qt[i + mo_offset];
            iwl_buf_wrt_val(OutBuf, qt_i, qt_i, qt_i, qt_i, 1.0, 0, "outfile", 0);
            for (int qt_j = 0; qt_j < qt_i; qt_j++) {
                iwl_buf_wrt_val(OutBuf, qt_i, qt_i, qt_j, qt_j, 2.0, 0, "outfile", 0);
                iwl_buf_wrt_val(OutBuf, qt_i, qt_j, qt_j, qt_i, -1.0, 0, "outfile", 0);
            }
        }
        for (int i = 0; i < moinfo.openpi[h]; i++) {
            // One electron open-shell
            moinfo.opdm.add(h, i + clsd_h, i + clsd_h, 1);
            moinfo.opdm_a.add(h, i + clsd_h, i + clsd_h, 1);
            auto qt_i = moinfo.pitzer2qt[i + mo_offset + clsd_h];
            // Two electron open x closed
            for (int qt_j = 0; qt_j < moinfo.nfzc + moinfo.nclsd; qt_j++) {
                iwl_buf_wrt_val(OutBuf, qt_i, qt_i, qt_j, qt_j, 1.0, 0, "outfile", 0);
                iwl_buf_wrt_val(OutBuf, qt_i, qt_j, qt_j, qt_i, -0.5, 0, "outfile", 0);
            }
            // Two electron open x open
            for (int qt_j = (moinfo.nfzc + moinfo.nclsd); qt_j < qt_i; qt_j++) {
                iwl_buf_wrt_val(OutBuf, qt_i, qt_i, qt_j, qt_j, 0.5, 0, "outfile", 0);
                iwl_buf_wrt_val(OutBuf, qt_i, qt_j, qt_j, qt_i, -0.5, 0, "outfile", 0);
            }
        }
        mo_offset += moinfo.orbspi[h];
    }
}

}  // namespace ccdensity
}  // namespace psi
