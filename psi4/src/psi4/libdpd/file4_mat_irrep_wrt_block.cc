/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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
#include <cstdio>
#include "psi4/libpsio/psio.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "dpd.h"

namespace psi {

int DPD::file4_mat_irrep_wrt_block(dpdfile4 *File, int irrep, int start_pq, int num_pq) {
    int rowtot, coltot, my_irrep;
    int seek_block;
    psio_address irrep_ptr, next_address;
    long int size;

    if (File->incore) {
        file4_cache_dirty(File); /* Flag this cache entry for writing */
        return 0;                /* We're keeping this data in core */
    }

    my_irrep = File->my_irrep;
    irrep_ptr = File->lfiles[irrep];
    rowtot = num_pq;
    coltot = File->params->coltot[irrep ^ my_irrep];
    size = ((long)rowtot) * ((long)coltot);

    /* Advance file pointer to current row */
    if (coltot) {
        seek_block = DPD_BIGNUM / (coltot * sizeof(double)); /* no. of rows for which we can compute the address */
        if (seek_block < 1) {
            outfile->Printf("\nLIBDPD Error: each row of %s is too long to compute an address.\n", File->label);
            dpd_error("dpd_file4_mat_irrep_rd_block", "outfile");
        }
        for (; start_pq > seek_block; start_pq -= seek_block)
            irrep_ptr = psio_get_address(irrep_ptr, sizeof(double) * seek_block * coltot);
        irrep_ptr = psio_get_address(irrep_ptr, sizeof(double) * start_pq * coltot);
    }

    if (rowtot && coltot)
        psio_write(File->filenum, File->label, (char *)File->matrix[irrep][0], size * ((long)sizeof(double)), irrep_ptr,
                   &next_address);

    return 0;
}

}  // namespace psi
