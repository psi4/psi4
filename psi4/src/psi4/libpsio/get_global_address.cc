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
 \ingroup PSIO
 */

#include "psi4/libpsio/psio.h"

namespace psi {
  /*!
   ** PSIO_GET_GLOBAL_ADDRESS(): Given the global starting address for a
   ** TOC entry and a relative offset within the entry, compute the global
   ** address for the offset.
   **
   ** \ingroup PSIO
   */

psio_address psio_get_global_address(psio_address entry_start,
                                     psio_address rel_address) {
  psio_address address;

  address.page = entry_start.page + rel_address.page;
  address.offset = entry_start.offset + rel_address.offset;
  if ((entry_start.offset + rel_address.offset) >= PSIO_PAGELEN) {
    address.offset -= PSIO_PAGELEN;
    address.page++;
  }

  return (address);
}

}
