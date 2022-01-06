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

/*!
  \file
  \ingroup IWL
*/
#include <cstdio>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "iwl.h"
#include "iwl.hpp"
#include "psi4/libpsi4util/PsiOutStream.h"
namespace psi {

#define MAX0(a, b) (((a) > (b)) ? (a) : (b))
#define MIN0(a, b) (((a) < (b)) ? (a) : (b))
#define INDEX(i, j) ((i > j) ? (ioff[(i)] + (j)) : (ioff[(j)] + (i)))

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
    std::shared_ptr<psi::PsiOutStream> printer = (out == "outfile" ? outfile : std::make_shared<PsiOutStream>(out));
    int idx, r, s, R, S, rtr, str;
    int ij, kl;
    double value;
    Label *lblptr;
    Value *valptr;

    lblptr = labels_;
    valptr = values_;

    ij = INDEX(ptr, qtr);

    for (r = rfirst, R = 0; r <= rlast; r++, R++) {
        rtr = reorder[r] - reorder_offset;

        for (s = sfirst, S = 0; s <= slast && s <= r; s++, S++) {
            str = reorder[s] - reorder_offset;

            kl = INDEX(rtr, str);

            value = mat[R][S];

            if (ij >= kl && std::fabs(value) > cutoff_) {
                idx = 4 * idx_;
                lblptr[idx++] = (Label)MAX0(ptr, qtr);
                lblptr[idx++] = (Label)MIN0(ptr, qtr);
                lblptr[idx++] = (Label)MAX0(rtr, str);
                lblptr[idx++] = (Label)MIN0(rtr, str);
                valptr[idx_] = (Value)value;

                idx_++;

                if (idx_ == ints_per_buf_) {
                    lastbuf_ = 0;
                    inbuf_ = idx_;
                    put();
                    idx_ = 0;
                }

                if (printflag) printer->Printf(">%d %d %d %d [%d] [%d] = %20.10f\n", ptr, qtr, rtr, str, ij, kl, value);

            } /* end if (std::fabs(value) > Buf->cutoff) ... */
        }     /* end loop over s */
    }         /* end loop over r */
}
}
