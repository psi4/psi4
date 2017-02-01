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

/*! \file
    \ingroup DPD
    \brief Enter brief description of file here
*/
#include "psi4/libpsio/psio.h"
#include "dpd.h"

namespace psi {

int DPD::file2_mat_wrt(dpdfile2 *File)
{
    int h, my_irrep, rowtot, coltot;
    psio_address irrep_ptr, next_address;

    my_irrep = File->my_irrep;

    if(File->incore) {
        file2_cache_dirty(File); /* Flag this cache entry for writing */
        return 0;  /* We're keeping this data in core */
    }

    for(h=0; h < File->params->nirreps; h++) {
        irrep_ptr = File->lfiles[h];
        rowtot = File->params->rowtot[h];
        coltot = File->params->coltot[h^my_irrep];

        if(rowtot && coltot)
            psio_write(File->filenum, File->label, (char *) File->matrix[h][0],
                    rowtot*coltot*sizeof(double), irrep_ptr, &next_address);
    }

    return 0;
}

}
