/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

/*! \file
    \ingroup DPD
    \brief Enter brief description of file here
*/
#include "psi4/libpsio/psio.h"
#include "dpd.h"

namespace psi {

int DPD::file4_mat_irrep_wrt(dpdfile4 *File, int irrep) {
    int rowtot, coltot, my_irrep;
    psio_address irrep_ptr, next_address;
    long int size;

    if (File->incore) {
        file4_cache_dirty(File); /* Flag this cache entry for writing */
        return 0;                /* We're keeping this data in core */
    }

    my_irrep = File->my_irrep;
    irrep_ptr = File->lfiles[irrep];
    rowtot = File->params->rowtot[irrep];
    coltot = File->params->coltot[irrep ^ my_irrep];
    size = ((long)rowtot) * ((long)coltot);

    if (rowtot && coltot)
        psio_write(File->filenum, File->label, (char *)File->matrix[irrep][0], size * ((long)sizeof(double)), irrep_ptr,
                   &next_address);

    return 0;
}

}  // namespace psi
