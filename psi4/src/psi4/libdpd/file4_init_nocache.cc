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
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "psi4/libpsio/psio.h"
#include "dpd.h"

namespace psi {

/** This file4_init() does not add the file4 to the cache **/

/* dpd_file4_init(): Prepares a dpd four-index file on disk for
** reading/writing.
**
** Arguments:
**   dpdfile4 *File: A pointer to the dpdfile4 to be initialized.
**   int filenum: The PSI unit number for this file.
**   int irrep: The irrep of this quantity.
**   int pqnum: The index combination for the bra indices for the
**              data as it will be stored on disk.
**   int rsnum: The index combination for the ket indices for the
**              data as it will be stored on disk.
**   char *label: A string labelling for this buffer.
*/

int DPD::file4_init_nocache(dpdfile4 *File, int filenum, int irrep, int pqnum,
                            int rsnum,  const char *label)
{
    int i;
    int maxrows, rowtot, coltot;
    dpd_file4_cache_entry *this_entry;
    psio_address irrep_ptr;

    File->dpdnum = dpd_default;
    File->params = &(dpd_list[dpd_default]->params4[pqnum][rsnum]);

    strcpy(File->label,label);
    File->filenum = filenum;
    File->my_irrep = irrep;

    this_entry = file4_cache_scan(filenum, irrep, pqnum, rsnum, label, dpd_default);
    if(this_entry != NULL) {
        File->incore = 1;
        File->matrix = this_entry->matrix;
    }
    else {
        File->incore = 0;
        File->matrix = (double ***) malloc(File->params->nirreps*sizeof(double **));
    }

    /* Construct logical subfile pointers */
    File->lfiles = (psio_address *) malloc(File->params->nirreps *
                                           sizeof(psio_address));
    File->lfiles[0] = PSIO_ZERO;
    for(i=1; i < File->params->nirreps; i++) {

        rowtot = File->params->rowtot[i-1];
        coltot = File->params->coltot[(i-1)^irrep];

        if(coltot) {
            /* number of rows for which we can compute the address offset directly */
            maxrows = DPD_BIGNUM/(coltot*sizeof(double));
            if(maxrows < 1) {
                outfile->Printf( "\nLIBDPD Error: each row of %s is too long to compute an address.\n",
                        File->label);
                dpd_error("dpd_file4_init_nocache", "outfile");
            }
        }
        else maxrows = DPD_BIGNUM;

        /* compute the file offset by increments */
        irrep_ptr = File->lfiles[i-1];
        for(; rowtot > maxrows; rowtot -= maxrows)
            irrep_ptr = psio_get_address(irrep_ptr, maxrows*coltot*sizeof(double));
        irrep_ptr = psio_get_address(irrep_ptr, rowtot*coltot*sizeof(double));

        File->lfiles[i] = irrep_ptr;
    }

    return 0;
}


}
