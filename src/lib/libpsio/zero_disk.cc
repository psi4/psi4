/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
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

#include <boost/shared_ptr.hpp>
#include <libpsio/psio.h>
#include <libpsio/psio.hpp>
#include <string.h>

namespace psi {

void PSIO::zero_disk(unsigned int unit, const char *key, ULI rows, ULI cols) {

      double* buf = new double[cols];
      ::memset(static_cast<void*>(buf),'\0',cols*sizeof(double));

      psio_address next_psio = PSIO_ZERO;
      for (int i=0; i<rows; i++) {
          PSIO::write(unit,key,(char *) (buf),
          sizeof(double)*cols,next_psio,&next_psio);
      }

      delete[] buf;
}

  /*!
   ** PSIO_ZERO_DISK()
   **
   ** \ingroup PSIO
   */

  int psio_zero_disk(unsigned int unit, const char *key, ULI rows, ULI cols) {
    _default_psio_lib_->zero_disk(unit, key, rows, cols);
    return 1;
  }

}