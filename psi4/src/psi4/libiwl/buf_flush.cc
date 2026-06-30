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
#include "iwl.h"
#include "iwl.hpp"

namespace psi {

namespace {

// Zero-fill the tail of a partially-filled bucket so that on-disk consumers
// always read a well-defined bucket layout, then mark it written.
void zero_fill_tail(Label *labels, Value *values, int &idx, int ints_per_buf, int lastbuf_flag,
                    int &inbuf_out, int &lastbuf_out) {
    inbuf_out = idx;
    for (int i = idx; i < ints_per_buf; ++i) {
        labels[4 * i + 0] = 0;
        labels[4 * i + 1] = 0;
        labels[4 * i + 2] = 0;
        labels[4 * i + 3] = 0;
        values[i] = 0.0;
    }
    idx = ints_per_buf;
    lastbuf_out = lastbuf_flag ? 1 : 0;
}

}  // namespace

void IWL::flush(int lastbuf) {
    zero_fill_tail(labels_, values_, idx_, ints_per_buf_, lastbuf, inbuf_, lastbuf_);
    put();
    idx_ = 0;
}

/*!
** iwl_buf_flush()
**
**	\param Buf     To be flushed buffer
**	\param lastbuf Flag for the last buffer
**
** Flush an Integrals With Labels Buffer
** All flushing should be done through this routine!
** David Sherrill, March 1995
** \ingroup IWL
*/
void iwl_buf_flush(struct iwlbuf *Buf, int lastbuf) {
    zero_fill_tail(Buf->labels, Buf->values, Buf->idx, Buf->ints_per_buf, lastbuf, Buf->inbuf, Buf->lastbuf);
    iwl_buf_put(Buf);
    Buf->idx = 0;
}
}
