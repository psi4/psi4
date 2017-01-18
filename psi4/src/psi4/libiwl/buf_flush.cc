/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*!
  \file
  \ingroup IWL
*/
#include <cstdio>
#include "psi4/libciomr/libciomr.h"
#include "iwl.h"
#include "iwl.hpp"

namespace psi {

void IWL::flush(int lastbuf)
{
    int idx;
    Label *lblptr;
    Value *valptr;

    inbuf_ = idx_;
    lblptr = labels_;
    valptr = values_;

    idx = 4 * idx_;

    while (idx_ < ints_per_buf_) {
        lblptr[idx++] = 0;
        lblptr[idx++] = 0;
        lblptr[idx++] = 0;
        lblptr[idx++] = 0;
        valptr[idx_] = 0.0;
        idx_++;
    }

    if (lastbuf) lastbuf_ = 1;
    else lastbuf_ = 0;

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
void iwl_buf_flush(struct iwlbuf *Buf, int lastbuf)
{
  int idx;
  Label *lblptr;
  Value *valptr;

  Buf->inbuf = Buf->idx;
  lblptr = Buf->labels;
  valptr = Buf->values;

  idx = 4 * Buf->idx;

  while (Buf->idx < Buf->ints_per_buf) {
    lblptr[idx++] = 0;
    lblptr[idx++] = 0;
    lblptr[idx++] = 0;
    lblptr[idx++] = 0;
    valptr[Buf->idx] = 0.0;
    Buf->idx++;
  }

  if (lastbuf) Buf->lastbuf = 1;
  else Buf->lastbuf = 0;

  iwl_buf_put(Buf);
  Buf->idx = 0;
}

}
