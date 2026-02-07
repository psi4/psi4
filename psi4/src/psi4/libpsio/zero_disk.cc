/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2025 The Psi4 Developers.
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
#include <string.h>

namespace psi {

/*!
 ** PSIO_ZERO_DISK()
 **
 ** \ingroup PSIO
 */
void PSIO::zero_disk(size_t unit, const char *key, size_t rows, size_t cols) {
    double *buf = new double[cols];
    ::memset(static_cast<void *>(buf), '\0', cols * sizeof(double));

    psio_address next_psio = PSIO_ZERO;
    for (int i = 0; i < rows; i++) {
        PSIO::write(unit, key, (char *)(buf), sizeof(double) * cols, next_psio, &next_psio);
    }

    delete[] buf;
}

void PSIO::zero_disk(size_t unit, const std::string& key, size_t rows, size_t cols) { zero_disk(unit, key.c_str(), rows, cols); };

}  // namespace psi
