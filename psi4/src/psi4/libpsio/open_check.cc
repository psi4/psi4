/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2023 The Psi4 Developers.
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
 \ingroup PSIO
 */

#include "psi4/pragma.h"
#include <memory>
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"

namespace psi {

int PSIO::open_check(size_t unit) {
    psio_ud *this_unit;

    this_unit = &(psio_unit[unit]);

    if (this_unit->vol[0].stream != -1)
        return 1;
    else
        return 0;
}

/*!
 ** PSIO_OPEN_CHECK(): Check to see if a given PSI direct access file
 ** is already open.
 **
 ** \param unit = the PSI unit number.
 **
 ** \ingroup PSIO
 */

int psio_open_check(size_t unit) { return _default_psio_lib_->open_check(unit); }

}  // namespace psi
