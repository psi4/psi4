/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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
#include "psi4/libciomr/libciomr.h"
#include "psi4/libiwl/iwl.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccdensity {

void add_core_ROHF(struct iwlbuf *OutBuf) {
    const auto actpi = moinfo.occpi + moinfo.virtpi;
    int mo_offset = 0;

    for (int h = 0; h < moinfo.nirreps; h++) {
        mo_offset += moinfo.frdocc[h];
        for (int p = 0; p < actpi[h]; p++) {
            for (int q = 0; q < actpi[h]; q++) {
                double value = moinfo.opdm.get(h, p, q);
                double p_qt = moinfo.pitzer2qt[p + mo_offset];
                double q_qt = moinfo.pitzer2qt[q + mo_offset];
                for (int m = 0; m < moinfo.nfzc; m++) {
                    iwl_buf_wrt_val(OutBuf, p_qt, q_qt, m, m, value, 0, "outfile", 0);
                    iwl_buf_wrt_val(OutBuf, p_qt, m, m, q_qt, -0.5 * value, 0, "outfile", 0);
                }
            }
        }
        mo_offset += moinfo.orbspi[h] - moinfo.frdocc[h];
    }
}

}  // namespace ccdensity
}  // namespace psi
