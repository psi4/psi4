/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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
#include "psi4/libiwl/iwl_writer.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#include "globals.h"

namespace psi {
namespace ccdensity {

/* doesn't work yet */

void add_core_UHF(IWLWriter &OutBuf) {
    int p, q, m, n;
    int nmo, nfzv, nfzc;
    double value;

    nmo = moinfo.nmo;
    nfzv = moinfo.nfzv;
    nfzc = moinfo.nfzc;

    return;

    for (p = nfzc; p < (nmo - nfzv); p++) {
        for (q = nfzc; q < (nmo - nfzv); q++) {
            value = moinfo.opdm_a[p][q];
            for (m = 0; m < nfzc; m++) {
                OutBuf.write(p, q, m, m, value);
                OutBuf.write(p, m, m, q, -0.5 * value);
            }
        }
    }
}

}  // namespace ccdensity
}  // namespace psi
