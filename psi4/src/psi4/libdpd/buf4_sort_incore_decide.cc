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
    \ingroup DPD
    \brief Decision logic for incore vs. out-of-core sorting
*/

#include "dpd.h"
#include "psi4/libpsi4util/PsiOutStream.h"

namespace psi {
bool DPD::buf4_sort_incore_decide(const dpdbuf4 &InBuf, const int32_t nirreps, const int32_t my_irrep) {
    bool incore = true;
    int64_t core_total = 0;
    for (int32_t h = 0; h < nirreps; h++) {
        const int64_t coltot = InBuf.params->coltot[h ^ my_irrep];
        const int64_t maxrows = [coltot, this] {
            int64_t mxrw;
            if (coltot != 0) {
                mxrw = DPD_BIGNUM / coltot;
                if (mxrw < 1) {
                    outfile->Printf("\nLIBDPD Error: too many rows in buf4_sort_axpy.\n");
                    dpd_error("buf4_sort_axpy", "outfile");
                }
            } else {
                mxrw = DPD_BIGNUM;
            }
            return mxrw;
        }();

        int64_t rowtot = InBuf.params->rowtot[h];
        for (; rowtot > maxrows; rowtot -= maxrows) {
            if (core_total > (core_total + 2 * maxrows * coltot)) {
                incore = false;
            } else {
                core_total += 2 * maxrows * coltot;
            }
        }
        if (core_total > (core_total + 2 * rowtot * coltot)) incore = false;
        core_total += 2 * rowtot * coltot;
    }
    if (core_total > dpd_memfree()) incore = false;
    return incore;
}
}  // namespace psi