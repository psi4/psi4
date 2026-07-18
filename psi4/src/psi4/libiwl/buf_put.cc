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

void put_bucket(PSIO *psio, int itap, psio_address &bufpos, int lastbuf, int inbuf, const Label *labels,
                const Value *values, int ints_per_buf) {
    psio->write(itap, IWL_KEY_BUF, reinterpret_cast<char *>(const_cast<int *>(&lastbuf)), sizeof(int), bufpos, &bufpos);
    psio->write(itap, IWL_KEY_BUF, reinterpret_cast<char *>(const_cast<int *>(&inbuf)), sizeof(int), bufpos, &bufpos);
    psio->write(itap, IWL_KEY_BUF, reinterpret_cast<char *>(const_cast<Label *>(labels)),
                ints_per_buf * 4 * sizeof(Label), bufpos, &bufpos);
    psio->write(itap, IWL_KEY_BUF, reinterpret_cast<char *>(const_cast<Value *>(values)), ints_per_buf * sizeof(Value),
                bufpos, &bufpos);
}

}  // namespace

void IWL::put() {
    put_bucket(psio_, itap_, bufpos_, lastbuf_, inbuf_, labels_, values_, ints_per_buf_);
}

/*!
** iwl_buf_put(struct iwlbuf *Buf)
**
** Put an IWL buffer to disk
** David Sherrill, 26 June 1996
** \ingroup IWL
*/
void iwl_buf_put(struct iwlbuf *Buf) {
    put_bucket(_default_psio_lib_.get(), Buf->itap, Buf->bufpos, Buf->lastbuf, Buf->inbuf, Buf->labels, Buf->values,
               Buf->ints_per_buf);
}
}
