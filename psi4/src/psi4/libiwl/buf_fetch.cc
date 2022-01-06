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
#include "psi4/libpsio/psio.h"
#include "iwl.h"
#include "iwl.hpp"

namespace psi {

void IWL::fetch() {
    psio_->read(itap_, IWL_KEY_BUF, (char *)&(lastbuf_), sizeof(int), bufpos_, &bufpos_);
    psio_->read(itap_, IWL_KEY_BUF, (char *)&(inbuf_), sizeof(int), bufpos_, &bufpos_);
    psio_->read(itap_, IWL_KEY_BUF, (char *)labels_, ints_per_buf_ * 4 * sizeof(Label), bufpos_, &bufpos_);
    psio_->read(itap_, IWL_KEY_BUF, (char *)values_, ints_per_buf_ * sizeof(Value), bufpos_, &bufpos_);
    idx_ = 0;
}

/*!
** iwl_buf_fetch()
**
** Fetch an IWL buffer from disk
** David Sherrill, 26 June 1996
** \ingroup IWL
*/
void PSI_API iwl_buf_fetch(struct iwlbuf *Buf) {
    psio_read(Buf->itap, IWL_KEY_BUF, (char *)&(Buf->lastbuf), sizeof(int), Buf->bufpos, &Buf->bufpos);
    psio_read(Buf->itap, IWL_KEY_BUF, (char *)&(Buf->inbuf), sizeof(int), Buf->bufpos, &Buf->bufpos);
    psio_read(Buf->itap, IWL_KEY_BUF, (char *)Buf->labels, Buf->ints_per_buf * 4 * sizeof(Label), Buf->bufpos,
              &Buf->bufpos);
    psio_read(Buf->itap, IWL_KEY_BUF, (char *)Buf->values, Buf->ints_per_buf * sizeof(Value), Buf->bufpos,
              &Buf->bufpos);
    Buf->idx = 0;
}
}
