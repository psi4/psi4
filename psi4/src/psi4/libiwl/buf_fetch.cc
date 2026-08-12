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
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include "iwl.h"
#include "iwl.hpp"

namespace psi {

namespace {

// Single source of truth for the bucket read sequence. Both the C-style and
// C++ APIs go through here.
void fetch_bucket(PSIO *psio, int itap, psio_address &bufpos, int &lastbuf, int &inbuf, Label *labels, Value *values,
                  int ints_per_buf, int &idx) {
    psio->read(itap, IWL_KEY_BUF, reinterpret_cast<char *>(&lastbuf), sizeof(int), bufpos, &bufpos);
    psio->read(itap, IWL_KEY_BUF, reinterpret_cast<char *>(&inbuf), sizeof(int), bufpos, &bufpos);
    psio->read(itap, IWL_KEY_BUF, reinterpret_cast<char *>(labels), ints_per_buf * 4 * sizeof(Label), bufpos, &bufpos);
    psio->read(itap, IWL_KEY_BUF, reinterpret_cast<char *>(values), ints_per_buf * sizeof(Value), bufpos, &bufpos);
    idx = 0;
}

}  // namespace

void IWL::fetch() { fetch_bucket(psio_, itap_, bufpos_, lastbuf_, inbuf_, labels_, values_, ints_per_buf_, idx_); }

/*!
** iwl_buf_fetch()
**
** Fetch an IWL buffer from disk
** David Sherrill, 26 June 1996
** \ingroup IWL
*/
void PSI_API iwl_buf_fetch(struct iwlbuf *Buf) {
    fetch_bucket(_default_psio_lib_.get(), Buf->itap, Buf->bufpos, Buf->lastbuf, Buf->inbuf, Buf->labels, Buf->values,
                 Buf->ints_per_buf, Buf->idx);
}
}
