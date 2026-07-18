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

/*!
  \file
  \ingroup IWL
*/
#include <cmath>
#include <memory>
#include "iwl.h"
#include "iwl.hpp"
#include "iwl_impl.h"

namespace psi {

using iwl_impl::check_label_fits;

/*!
** iwl_buf_wrt()
**
** Write to an Integrals With Labels formatted buffer.
** The buffer must have been initialized with iwl_buf_init().  Don't
** forget to call iwl_buf_flush() when finished with all writes to the
** buffer to ensure that all contents are written to disk.
** David Sherrill, March 1995
**
** Revised 6/27/96 by CDS for new format
** \ingroup IWL
*/
void IWL::write(int p, int q, int pq, int pqsym, double *arr, int rmax, int *ioff, int *orbsym, int *firsti, int *lasti,
                int printflag, std::string out) {
    std::shared_ptr<PsiOutStream> printer = printflag ? ((out == "outfile") ? outfile : std::make_shared<PsiOutStream>(out))
                                                      : nullptr;

    for (int r = 0; r < rmax; r++) {
        int rsym = orbsym[r];
        int ssym = pqsym ^ rsym;
        int smax = (rsym == ssym) ? r : lasti[ssym];

        for (int s = firsti[ssym]; s <= smax; s++) {
            int rs = ioff[r] + s;
            double value = arr[rs];

            if (std::fabs(value) <= cutoff_) continue;
            check_label_fits(p, q, r, s);

            int idx = 4 * idx_;
            labels_[idx + 0] = static_cast<Label>(p);
            labels_[idx + 1] = static_cast<Label>(q);
            labels_[idx + 2] = static_cast<Label>(r);
            labels_[idx + 3] = static_cast<Label>(s);
            values_[idx_] = static_cast<Value>(value);

            idx_++;

            if (idx_ == ints_per_buf_) {
                inbuf_ = idx_;
                lastbuf_ = 0;
                put();
                idx_ = 0;
            }

            if (printflag) printer->Printf("<%d %d %d %d [%d] [%d] = %20.10f\n", p, q, r, s, pq, rs, value);
        }
    }
}
}
