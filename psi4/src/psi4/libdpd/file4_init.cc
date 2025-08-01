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

/*! \file
    \ingroup DPD
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "psi4/libpsio/psio.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "dpd.h"

namespace psi {

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

int DPD::file4_init(dpdfile4 *File, int filenum, int irrep, int pqnum, int rsnum, const char *label) {
    int i;
    int maxrows, rowtot, coltot;
    size_t priority;
    dpd_file4_cache_entry *this_entry;
    psio_address irrep_ptr;

    File->dpdnum = dpd_default;
    File->params = &(dpd_list[dpd_default]->params4[pqnum][rsnum]);

    strcpy(File->label, label);
    File->filenum = filenum;
    File->my_irrep = irrep;

    this_entry = file4_cache_scan(filenum, irrep, pqnum, rsnum, label, dpd_default);
    if (this_entry != nullptr) {
        File->incore = true;
        File->matrix = this_entry->matrix;
    } else {
        File->incore = false;
        File->matrix = (double ***)malloc(File->params->nirreps * sizeof(double **));
    }

    /* Construct logical subfile pointers */
    File->lfiles = (psio_address *)malloc(File->params->nirreps * sizeof(psio_address));
    File->lfiles[0] = PSIO_ZERO;
    for (i = 1; i < File->params->nirreps; i++) {
        rowtot = File->params->rowtot[i - 1];
        coltot = File->params->coltot[(i - 1) ^ irrep];

        if (coltot) {
            /* number of rows for which we can compute the address offset directly */
            maxrows = DPD_BIGNUM / (coltot * sizeof(double));
            if (maxrows < 1) {
                outfile->Printf("\nLIBDPD Error: each row of %s is too long to compute an address.\n", File->label);
                dpd_error("dpd_file4_init", "outfile");
            }
        } else
            maxrows = DPD_BIGNUM;

        /* compute the file offset by increments */
        irrep_ptr = File->lfiles[i - 1];
        for (; rowtot > maxrows; rowtot -= maxrows)
            irrep_ptr = psio_get_address(irrep_ptr, sizeof(double) * maxrows * coltot);
        irrep_ptr = psio_get_address(irrep_ptr, sizeof(double) * rowtot * coltot);

        File->lfiles[i] = irrep_ptr;
    }

    /* Put this file4 into cache if requested */
    if (dpd_main.cachefiles[filenum] && dpd_main.cachelist[pqnum][rsnum]) {
        /* Get the file4's cache priority */
        if (dpd_main.cachetype == 1)
            priority = file4_cache_get_priority(File);
        else
            priority = 0;

        file4_cache_add(File, priority);

        /* Make sure this cache entry can't be deleted until we're done */
        file4_cache_lock(File);
    }

    return 0;
}

}  // namespace psi
