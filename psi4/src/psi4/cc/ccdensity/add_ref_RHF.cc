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
#include "psi4/libiwl/iwl.h"
#include "MOInfo.h"
#include "Params.h"
#include "Frozen.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccdensity {

void add_ref_RHF(struct iwlbuf *OutBuf) {
    int i, j;
    int nfzc, nclsd, nopen;

    nfzc = moinfo.nfzc;
    nclsd = moinfo.nclsd;
    nopen = moinfo.nopen;

    /*** One-electron component ***/

    for (i = 0; i < (nfzc + nclsd); i++) moinfo.opdm[i][i] += 2.0;

    for (i = nfzc + nclsd; i < (nfzc + nclsd + nopen); i++) moinfo.opdm[i][i] += 1.0;

    /*** Two-electron component ***/

    /* docc-docc */
    for (i = 0; i < (nfzc + nclsd); i++) {
        iwl_buf_wrt_val(OutBuf, i, i, i, i, 1.0, 0, "outfile", 0);
        for (j = 0; j < i; j++) {
            iwl_buf_wrt_val(OutBuf, i, i, j, j, 2.0, 0, "outfile", 0);
            iwl_buf_wrt_val(OutBuf, i, j, j, i, -1.0, 0, "outfile", 0);
        }
    }
}

}  // namespace ccdensity
}  // namespace psi
