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
 ** \file
 ** \ingroup PSIO
 */

#include <cstdlib>
#include "psi4/pragma.h"
#include <memory>
#include "psi4/libpsio/psio.h"
#include "psi4/libpsio/psio.hpp"

#ifdef PSIO_STATS
#include <ctime>
#endif

namespace psi {

PSIO::~PSIO() {
#ifdef PSIO_STATS
    int i;
    size_t total_read = 0, total_write = 0;
    FILE *io_out;
    std::time_t my_time;
    my_time = std::time(nullptr);
    io_out = fopen("io.dat", "a+");
    fprintf(io_out, "\nLIBPSIO Read/Write Statistics\n\n");
    fprintf(io_out, "Run at: %s\n", ctime(&my_time));
    fprintf(io_out, "Unit      Read(kB)    Write(kB)\n");
    fprintf(io_out, "-------------------------------\n");
    for (i = 0; i < PSIO_MAXUNIT; i++) {
        total_read += psio_readlen[i];
        total_write += psio_writlen[i];

        if (psio_readlen[i] || psio_writlen[i])
            fprintf(io_out, "%3d   %10.1f   %10.1f\n", i, ((double)psio_readlen[i]) / ((double)1024),
                    ((double)psio_writlen[i]) / ((double)1024));
    }
    fprintf(io_out, "-------------------------------\n");
    fprintf(io_out, "Total %10.1f   %10.1f\n", ((double)total_read) / ((double)1024),
            ((double)total_write) / ((double)1024));
    fclose(io_out);
    free(psio_readlen);
    free(psio_writlen);
#endif

    free(psio_unit);
    state_ = 0;
    files_keywords_.clear();
}
}  // namespace psi
