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
#include <cstdlib>
#include "psi4/libpsio/psio.h"
#include "iwl.h"
#include "iwl.hpp"
#include "psi4/psi4-dec.h" // need outfile;

namespace psi {

void IWL::to_end()
{
    psio_tocentry *this_entry;
    ULI entry_length;

    this_entry = psio_->tocscan(itap_, IWL_KEY_BUF);
    if (this_entry == NULL) {
        outfile->Printf(
            "iwl_buf_toend: Can't find IWL buffer entry in file %d\n", itap_);
        set_keep_flag(1);
        close();
        return;
    }

    /* set up buffer pointer */
    entry_length = psio_get_length(this_entry->sadd,this_entry->eadd);
    bufpos_ = psio_get_address(PSIO_ZERO,entry_length);
}

/*!
** iwl_buf_toend()
**
** Set IWL Buffer's pointer to the end. Useful when want to append to
** an already existing file
**
** Edward Valeev, January 2001
** \ingroup IWL
*/
void iwl_buf_toend(struct iwlbuf *Buf)
{
  psio_tocentry *this_entry;
  ULI entry_length;

  this_entry = psio_tocscan(Buf->itap, IWL_KEY_BUF);
  if (this_entry == NULL) {
    outfile->Printf(
        "iwl_buf_toend: Can't find IWL buffer entry in file %d\n", Buf->itap);
    iwl_buf_close(Buf,1);
    return;
  }

  /* set up buffer pointer */
  entry_length = psio_get_length(this_entry->sadd,this_entry->eadd);
  Buf->bufpos = psio_get_address(PSIO_ZERO,entry_length);

  return;
}

}
