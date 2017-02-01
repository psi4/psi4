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
#include <cmath>
#include "psi4/libqt/qt.h"
#include "dpd.h"

namespace psi {

double DPD::buf4_dot(dpdbuf4 *BufA, dpdbuf4 *BufB)
{
    int h, nirreps, n, my_irrep;
    double dot;
    int incore, nbuckets;
    long int memoryd, rows_per_bucket, rows_left;

    nirreps = BufA->params->nirreps;
    my_irrep = BufA->file.my_irrep;

    dot = 0.0;

    for(h=0; h < nirreps; h++) {

        memoryd = dpd_memfree();

        if(BufA->params->rowtot[h] && BufA->params->coltot[h^my_irrep]) {

            /* Compute the memory for one row of A/B */
            if(BufA->params->coltot[h^my_irrep])
                /* NB: we need at least one row of both A and B */
                rows_per_bucket = memoryd/(2 * BufA->params->coltot[h^my_irrep]);
            else rows_per_bucket = -1;

            if(rows_per_bucket > BufA->params->rowtot[h])
                rows_per_bucket = BufA->params->rowtot[h];

            if(!rows_per_bucket)
                dpd_error("buf4_dot: Not enough memory for one row!", "outfile");

            nbuckets = (int) ceil((double) BufA->params->rowtot[h]/
                                  (double) rows_per_bucket);

            rows_left = BufA->params->rowtot[h] % rows_per_bucket;

            incore = 1;
            if(nbuckets > 1) incore = 0;

        }
        else incore = 1;

        if(incore) {

            buf4_mat_irrep_init(BufA, h);
            buf4_mat_irrep_init(BufB, h);
            buf4_mat_irrep_rd(BufA, h);
            buf4_mat_irrep_rd(BufB, h);

            dot += dot_block(BufA->matrix[h], BufB->matrix[h],
                             BufA->params->rowtot[h],
                             BufA->params->coltot[h^my_irrep], 1.0);

            buf4_mat_irrep_close(BufA, h);
            buf4_mat_irrep_close(BufB, h);
        }
        else {
            buf4_mat_irrep_init_block(BufA, h, rows_per_bucket);
            buf4_mat_irrep_init_block(BufB, h, rows_per_bucket);

            for(n=0; n < (rows_left ? nbuckets-1: nbuckets); n++) {

                buf4_mat_irrep_rd_block(BufA, h, n*rows_per_bucket,
                                        rows_per_bucket);
                buf4_mat_irrep_rd_block(BufB, h, n*rows_per_bucket,
                                        rows_per_bucket);

                dot += dot_block(BufA->matrix[h], BufB->matrix[h],
                                 rows_per_bucket,
                                 BufA->params->coltot[h^my_irrep], 1.0);
            }

            if(rows_left) {

                buf4_mat_irrep_rd_block(BufA, h, n*rows_per_bucket,
                                        rows_left);
                buf4_mat_irrep_rd_block(BufB, h, n*rows_per_bucket,
                                        rows_left);

                dot += dot_block(BufA->matrix[h], BufB->matrix[h],
                                 rows_left,
                                 BufA->params->coltot[h^my_irrep], 1.0);
            }

            buf4_mat_irrep_close_block(BufA, h, rows_per_bucket);
            buf4_mat_irrep_close_block(BufB, h, rows_per_bucket);

        }
    }

    return dot;
}

}
