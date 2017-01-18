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
#include "psi4/libpsio/psio.h"
#include "iwl.h"
#include "iwl.hpp"

namespace psi {

void IWL::put()
{
    psio_->write(itap_, IWL_KEY_BUF, (char *) &(lastbuf_), sizeof(int),
        bufpos_, &(bufpos_));
    psio_->write(itap_, IWL_KEY_BUF, (char *) &(inbuf_), sizeof(int),
        bufpos_, &(bufpos_));
    psio_->write(itap_, IWL_KEY_BUF, (char *) labels_, ints_per_buf_ *
        4 * sizeof(Label), bufpos_, &(bufpos_));
    psio_->write(itap_, IWL_KEY_BUF, (char *) values_, ints_per_buf_ *
        sizeof(Value), bufpos_, &(bufpos_));
}

/*!
** iwl_buf_put(struct iwlbuf *Buf)
**
** Put an IWL buffer to disk
** David Sherrill, 26 June 1996
** \ingroup IWL
*/
void iwl_buf_put(struct iwlbuf *Buf)
{
  psio_write(Buf->itap, IWL_KEY_BUF, (char *) &(Buf->lastbuf), sizeof(int),
	     Buf->bufpos, &(Buf->bufpos));
  psio_write(Buf->itap, IWL_KEY_BUF, (char *) &(Buf->inbuf), sizeof(int),
	     Buf->bufpos, &(Buf->bufpos));
  psio_write(Buf->itap, IWL_KEY_BUF, (char *) Buf->labels, Buf->ints_per_buf *
	     4 * sizeof(Label), Buf->bufpos, &(Buf->bufpos));
  psio_write(Buf->itap, IWL_KEY_BUF, (char *) Buf->values, Buf->ints_per_buf *
	     sizeof(Value), Buf->bufpos, &(Buf->bufpos));
}

}
