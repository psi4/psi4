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
#include <algorithm>
#include "iwl.h"
#include "iwl.hpp"
#include "iwl_impl.h"

namespace psi {

using iwl_impl::check_label_fits;

namespace {

inline int packed_index(int i, int j, const int *ioff) {
    return (i > j) ? (ioff[i] + j) : (ioff[j] + i);
}

}  // namespace

/*!
** iwl_buf_wrt_mat()
**
** Write to an Integrals With Labels formatted buffer. The buffer
** must have been initialized with iwl_buf_init().  Don't forget to
** call iwl_buf_flush() when finished with all writes to the buffer
** to ensure that all contents are written to disk.
**
** This version takes as input a matrix, as might be handy for a
** matrix formulation of an integral transformation.  It assumes that
** all rs are available for a given pq. r and s are allowed to range
** from rfirst/sfirst to rlast/slast (with s<=r), and this maps to a
** matrix addressing of 0 to (rlast-rfirst) and 0 to (slast-sfirst).
**
** This routine is also compatible with a reordering of the orbitals
** before output.  We assume that p and q are already reordered,
** but r and s are not (yet).  The reordered address of an orbital
** is computed according to rtr = reorder[r] - reorder_offset.
**
** I have further modified the routine to spit out integrals whose
** reordered indices are canonical ij >= kl.  The routine does not
** care whether ptr >= qtr, but it insists that the transformed pq
** and kl indices (called ij and kl in the routine) satisfy ij >= kl.
** David Sherrill, October 1995
**
** Revised 6/27/96 by CDS for new format
** \ingroup IWL
*/
void IWL::write_matrix(int ptr, int qtr, double **mat, int rfirst, int rlast, int sfirst, int slast, int *reorder,
                       int reorder_offset, int printflag, int *ioff, std::string out) {
    std::shared_ptr<PsiOutStream> printer = printflag ? ((out == "outfile") ? outfile : std::make_shared<PsiOutStream>(out))
                                                      : nullptr;

    const int ij = packed_index(ptr, qtr, ioff);

    for (int r = rfirst, R = 0; r <= rlast; r++, R++) {
        const int rtr = reorder[r] - reorder_offset;

        for (int s = sfirst, S = 0; s <= slast && s <= r; s++, S++) {
            const int str = reorder[s] - reorder_offset;
            const int kl = packed_index(rtr, str, ioff);
            const double value = mat[R][S];

            if (ij < kl) continue;
            if (std::fabs(value) <= cutoff_) continue;

            const int p_out = std::max(ptr, qtr);
            const int q_out = std::min(ptr, qtr);
            const int r_out = std::max(rtr, str);
            const int s_out = std::min(rtr, str);
            check_label_fits(p_out, q_out, r_out, s_out);

            int idx = 4 * idx_;
            labels_[idx + 0] = static_cast<Label>(p_out);
            labels_[idx + 1] = static_cast<Label>(q_out);
            labels_[idx + 2] = static_cast<Label>(r_out);
            labels_[idx + 3] = static_cast<Label>(s_out);
            values_[idx_] = static_cast<Value>(value);

            idx_++;

            if (idx_ == ints_per_buf_) {
                lastbuf_ = 0;
                inbuf_ = idx_;
                put();
                idx_ = 0;
            }

            if (printflag) printer->Printf(">%d %d %d %d [%d] [%d] = %20.10f\n", ptr, qtr, rtr, str, ij, kl, value);
        }
    }
}
}
