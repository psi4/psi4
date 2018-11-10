/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
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
PRAGMA_WARNING_PUSH
PRAGMA_WARNING_IGNORE_DEPRECATED_DECLARATIONS
#include <memory>
PRAGMA_WARNING_POP
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"
#include <string.h>

namespace psi {

void PSIO::zero_disk(size_t unit, const char *key, size_t rows, size_t cols) {
    double *buf = new double[cols];
    ::memset(static_cast<void *>(buf), '\0', cols * sizeof(double));

    psio_address next_psio = PSIO_ZERO;
    for (int i = 0; i < rows; i++) {
        PSIO::write(unit, key, (char *)(buf), sizeof(double) * cols, next_psio, &next_psio);
    }

    delete[] buf;
}

/*!
 ** PSIO_ZERO_DISK()
 **
 ** \ingroup PSIO
 */

int psio_zero_disk(size_t unit, const char *key, size_t rows, size_t cols) {
    _default_psio_lib_->zero_disk(unit, key, rows, cols);
    return 1;
}

}  // namespace psi
