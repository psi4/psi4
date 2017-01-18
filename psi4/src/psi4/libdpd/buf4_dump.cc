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
#include "psi4/libiwl/iwl.h"
#include "dpd.h"

namespace psi {

int DPD::buf4_dump(dpdbuf4 *DPDBuf, struct iwlbuf *IWLBuf,
                   int *prel, int *qrel, int *rrel, int *srel,
                   int bk_pack, int swap23)
{
    int h, row, col, p, q, r, s, P, Q, R, S, my_irrep;
    double value;

    my_irrep = DPDBuf->file.my_irrep;

    for(h=0; h < DPDBuf->params->nirreps; h++) {
        buf4_mat_irrep_init(DPDBuf, h);
        buf4_mat_irrep_rd(DPDBuf, h);
        for(row=0; row < DPDBuf->params->rowtot[h]; row++) {
            p = DPDBuf->params->roworb[h][row][0]; P = prel[p];
            q = DPDBuf->params->roworb[h][row][1]; Q = qrel[q];
            if(bk_pack) {
                for(col=0; col <= row; col++) {
                    r = DPDBuf->params->colorb[h^my_irrep][col][0]; R = rrel[r];
                    s = DPDBuf->params->colorb[h^my_irrep][col][1]; S = srel[s];

                    value = DPDBuf->matrix[h][row][col];

                    if(swap23)
                        iwl_buf_wrt_val(IWLBuf, P, R, Q, S, value, 0,
                                        "NULL", 0);
                    else
                        iwl_buf_wrt_val(IWLBuf, P, Q, R, S, value, 0,
                                        "NULL", 0);
                }
            }
            else {
                for(col=0; col < DPDBuf->params->coltot[h^my_irrep]; col++) {
                    r = DPDBuf->params->colorb[h^my_irrep][col][0]; R = rrel[r];
                    s = DPDBuf->params->colorb[h^my_irrep][col][1]; S = srel[s];

                    value = DPDBuf->matrix[h][row][col];

                    if(swap23)
                        iwl_buf_wrt_val(IWLBuf, P, R, Q, S, value, 0,
                                        "NULL", 0);
                    else
                        iwl_buf_wrt_val(IWLBuf, P, Q, R, S, value, 0,
                                        "NULL", 0);
                }
            }
        }
        buf4_mat_irrep_close(DPDBuf, h);
    }

    return 0;
}

}
