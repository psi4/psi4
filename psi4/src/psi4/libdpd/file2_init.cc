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

/* dpd_file2_init(): Initializes a dpd two-index file for reading
** or writing data.
**
** Arguments:
**   dpdfile2 *File: A pointer to the two-index dpdfile.
**   int filenum: The PSI unit number for this file.
**   int irrep: The symmetry of the data (=0 for totally-symmetric)
**   int pnum: The orbital subspace number for the left index [see
**             dpd_init()].
**   int qnum: The orbital subspace number for the right index [see
**             dpd_init()].
**   char *label: A string labelling for this buffer.
**   Note: Make sure that you use the correct label and inputfile combination.
**      If you intend to read from or write to an existing quantity on disk be sure
**      that the label string/file number point to that quantity. If you intend to
**      create and populate a new quantity on disk, ensure that the label is not
**      already used in the file.  PSIO::tocprint(int filenum) can be used to print
**      the labels currently used in in filenum and is quite useful for debugging.
*/

int DPD::file2_init(dpdfile2 *File, int filenum, int irrep, int pnum,
                    int qnum, const char *label)
{
    int i, q, rs, nirreps;
    dpd_file2_cache_entry *this_entry;

    File->dpdnum = dpd_default;
    File->params = &(dpd_list[dpd_default]->params2[pnum][qnum]);
    strcpy(File->label,label);
    File->filenum = filenum;
    File->my_irrep = irrep;

    nirreps = File->params->nirreps;

    this_entry = file2_cache_scan(filenum, irrep, pnum, qnum, label, dpd_default);
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
    for(i=1; i < File->params->nirreps; i++)
        File->lfiles[i] = psio_get_address(File->lfiles[i-1],
                (File->params->rowtot[i-1] *
                File->params->coltot[(i-1)^irrep] *
                sizeof(double)));

    /* Force all two-index files into cache */
    /*  dpd_file2_cache_add(File); */

    return 0;
}

}
