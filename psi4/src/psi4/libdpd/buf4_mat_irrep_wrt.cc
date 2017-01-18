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
#include "psi4/libqt/qt.h"
#include "dpd.h"
#include "psi4/psi4-dec.h"
namespace psi {

/* dpd_buf4_mat_irrep_wrt(): Writes an entire irrep from disk into a dpd
** four-index buffer using the "rules" specified when the buffer was
** initialized by dpd_buf4_init().
**
** Arguments:
**   dpdbuf4 *Buf: A pointer to the dpdbuf where the data is
**                       stored.
**   int irrep: The irrep number to be written.
**
** The case numbers are the same as those from buf4_mat_irrep_rd().
** Unpacking indices on the write has not been tested, and should be currently
** disabled.
**
** T. Daniel Crawford
** December 1996
**
** Modified for caching.
** TDC
** September 1999
*/

int DPD::buf4_mat_irrep_wrt(dpdbuf4 *Buf, int irrep)
{
    int method, filerow, all_buf_irrep;
    int pq, rs;  /* dpdfile row and column indices */
    int p, q, r, s;  /* orbital indices */
    int bufpq, bufrs;  /* Input dpdbuf row and column indices */
    int rowtot, coltot;  /* dpdfile row and column dimensions */
    int b_perm_pq, b_perm_rs, b_peq, b_res;
    int f_perm_pq, f_perm_rs, f_peq, f_res;
    int permute;
    double value;
    long int size;

    all_buf_irrep = Buf->file.my_irrep;

    /* Row and column dimensions in the DPD file */
    rowtot = Buf->file.params->rowtot[irrep];
    coltot = Buf->file.params->coltot[irrep^all_buf_irrep];
    size = ((long) rowtot) * ((long) coltot);

    /* Index packing information */
    b_perm_pq = Buf->params->perm_pq; b_perm_rs = Buf->params->perm_rs;
    f_perm_pq = Buf->file.params->perm_pq; f_perm_rs = Buf->file.params->perm_rs;
    b_peq = Buf->params->peq; b_res = Buf->params->res;
    f_peq = Buf->file.params->peq; f_res = Buf->file.params->res;

    /* Exit if buffer is antisymmetrized */
    if(Buf->anti) {
        printf("\n\tCannot write antisymmetrized buffer\n");
        printf(  "\tback to original DPD file!\n");
        exit(PSI_RETURN_FAILURE);
    }

    if((b_perm_pq == f_perm_pq) && (b_perm_rs == f_perm_rs) &&
            (b_peq == f_peq) && (b_res == f_res))   method = 12;
    else if((b_perm_pq != f_perm_pq) && (b_perm_rs == f_perm_rs) &&
            (b_res == f_res)) {
        if(f_perm_pq && !b_perm_pq) method = 21;
        else if(!f_perm_pq && b_perm_pq) method = 23;
        else {
            printf("\n\tInvalid second-level method!\n");
            exit(PSI_RETURN_FAILURE);
        }
    }
    else if((b_perm_pq == f_perm_pq) && (b_perm_rs != f_perm_rs) &&
            (b_peq == f_peq)) {
        if(f_perm_rs && !b_perm_rs) method = 31;
        else if(!f_perm_rs && b_perm_rs) method = 33;
        else {
            printf("\n\tInvalid third-level method!\n");
            exit(PSI_RETURN_FAILURE);
        }
    }
    else if((b_perm_pq != f_perm_pq) && (b_perm_rs != f_perm_rs)) {
        if(f_perm_pq && !b_perm_pq) {
            if(f_perm_rs && !b_perm_rs) method = 41;
            else if(!f_perm_rs && b_perm_rs) method = 42;
        }
        else if(!f_perm_pq && b_perm_pq) {
            if(f_perm_rs && !b_perm_rs) method = 43;
            else if(!f_perm_rs && b_perm_rs) method = 45;
        }
        else {
            printf("\n\tInvalid fourth-level method!\n");
            exit(PSI_RETURN_FAILURE);
        }
    }
    else {
        printf("\n\tInvalid method in dpd_buf_mat_irrep_rd!\n");
        exit(PSI_RETURN_FAILURE);
    }


    switch(method) {
    case 12: /* No change in pq or rs */

#ifdef DPD_TIMER
        timer_on("buf_wrt_12");
#endif

        if(Buf->file.incore && size) {
            file4_cache_dirty(&(Buf->file));
        }
        /*
      memcpy((void *) &(Buf->file.matrix[irrep][0][0]),
         (const void *) &(Buf->matrix[irrep][0][0]),
         sizeof(double)*rowtot*coltot);
*/
        else {
            Buf->file.matrix[irrep] = Buf->matrix[irrep];
            file4_mat_irrep_wrt(&(Buf->file), irrep);
        }

#ifdef DPD_TIMER
        timer_off("buf_wrt_12");
#endif

        break;
    case 21: /* Pack pq; no change in rs */
        /* Prepare the output buffer for the output DPD file */
        file4_mat_irrep_row_init(&(Buf->file), irrep);

        /* Loop over rows in the dpdfile */
        for(pq=0; pq < rowtot; pq++) {
            p = Buf->file.params->roworb[irrep][pq][0];
            q = Buf->file.params->roworb[irrep][pq][1];
            bufpq = Buf->params->rowidx[p][q];
            filerow = Buf->file.incore ? pq : 0;

            /* Loop over the columns in the dpdbuf */
            for(rs=0; rs < coltot; rs++) {
                bufrs = rs;

                value = Buf->matrix[irrep][bufpq][bufrs];

                /* Assign the value */
                Buf->file.matrix[irrep][filerow][rs] = value;
            }

            /* Write out the row */
            file4_mat_irrep_row_wrt(&(Buf->file), irrep, pq);
        }

        /* Close the input buffer */
        file4_mat_irrep_row_close(&(Buf->file), irrep);

        break;
    case 23: /* Unpack pq; no change in rs */
        /* I don't know if I'll ever use this, so I'll avoid it for now */
        printf("\n\tShould you be using method %d?\n", method);
        exit(PSI_RETURN_FAILURE);
        /* Prepare the output buffer for the output DPD file */
        file4_mat_irrep_row_init(&(Buf->file), irrep);

        /* Loop over rows in the dpdfile */
        for(pq=0; pq < rowtot; pq++) {
            p = Buf->file.params->roworb[irrep][pq][0];
            q = Buf->file.params->roworb[irrep][pq][1];
            bufpq = Buf->params->rowidx[p][q];

            filerow = Buf->file.incore ? pq : 0;

            /* Set the permutation operator's value */
            permute = ((p < q) && (b_perm_pq < 0) ? -1 : 1);

            /* Loop over the columns in the dpdbuf */
            for(rs=0; rs < coltot; rs++) {
                bufrs = rs;

                value = Buf->matrix[irrep][bufpq][bufrs];

                /* Assign the value */
                Buf->file.matrix[irrep][filerow][rs] = permute*value;
            }

            /* Write out the row */
            file4_mat_irrep_row_wrt(&(Buf->file), irrep, pq);
        }

        /* Close the input buffer */
        file4_mat_irrep_row_close(&(Buf->file), irrep);

        break;
    case 31: /* No change in pq; pack rs */
        /* Prepare the output buffer for the output DPD file */
        file4_mat_irrep_row_init(&(Buf->file), irrep);

        /* Loop over rows in the dpdbuf/dpdfile */
        for(pq=0; pq < rowtot; pq++) {
            bufpq = pq;

            filerow = Buf->file.incore ? pq : 0;

            /* Loop over the columns in the dpdfile */
            for(rs=0; rs < coltot; rs++) {
                r = Buf->file.params->colorb[irrep^all_buf_irrep][rs][0];
                s = Buf->file.params->colorb[irrep^all_buf_irrep][rs][1];
                bufrs = Buf->params->colidx[r][s];

                value = Buf->matrix[irrep][bufpq][bufrs];

                /* Assign the value */
                Buf->file.matrix[irrep][filerow][rs] = value;
            }

            /* Write out the row */
            file4_mat_irrep_row_wrt(&(Buf->file), irrep, pq);
        }

        /* Close the input buffer */
        file4_mat_irrep_row_close(&(Buf->file), irrep);

        break;
    case 33: /* No change in pq; unpack rs */
        /* I'm not sure if I'll ever need this, so I'm removing it for now */
        printf("\n\tShould you be using method %d?\n", method);
        exit(PSI_RETURN_FAILURE);
        /* Prepare the output buffer for the output DPD file */
        file4_mat_irrep_row_init(&(Buf->file), irrep);

        /* Loop over rows in the dpdbuf/dpdfile */
        for(pq=0; pq < rowtot; pq++) {
            bufpq = pq;

            filerow = Buf->file.incore ? pq : 0;

            /* Loop over the columns in the dpdfile */
            for(rs=0; rs < coltot; rs++) {
                r = Buf->file.params->colorb[irrep^all_buf_irrep][rs][0];
                s = Buf->file.params->colorb[irrep^all_buf_irrep][rs][1];
                bufrs = Buf->params->colidx[r][s];

                value = Buf->matrix[irrep][bufpq][bufrs];

                /* Assign the value */
                Buf->file.matrix[irrep][filerow][rs] = value;
            }

            /* Write out the row */
            file4_mat_irrep_row_wrt(&(Buf->file), irrep, pq);
        }

        /* Close the input buffer */
        file4_mat_irrep_row_close(&(Buf->file), irrep);

        break;
    case 41: /* Pack pq and rs */
        /* Prepare the output buffer for the output DPD file */
        file4_mat_irrep_row_init(&(Buf->file), irrep);

        /* Loop over rows in the dpdfile */
        for(pq=0; pq < rowtot; pq++) {
            p = Buf->file.params->roworb[irrep][pq][0];
            q = Buf->file.params->roworb[irrep][pq][1];
            bufpq = Buf->params->rowidx[p][q];

            filerow = Buf->file.incore ? pq : 0;

            /* Loop over the columns in the dpdfile */
            for(rs=0; rs < coltot; rs++) {
                r = Buf->file.params->colorb[irrep^all_buf_irrep][rs][0];
                s = Buf->file.params->colorb[irrep^all_buf_irrep][rs][1];
                bufrs = Buf->params->colidx[r][s];

                value = Buf->matrix[irrep][bufpq][bufrs];

                /* Assign the value */
                Buf->file.matrix[irrep][filerow][rs] = value;
            }

            /* Write out the row */
            file4_mat_irrep_row_wrt(&(Buf->file), irrep, pq);
        }

        /* Close the input buffer */
        file4_mat_irrep_row_close(&(Buf->file), irrep);

        break;
    case 42: /* Pack pq; unpack rs */
        printf("\n\tHaven't programmed method 42 yet!\n");
        exit(PSI_RETURN_FAILURE);

        break;
    case 43: /* Unpack pq; pack rs */
        printf("\n\tHaven't programmed method 43 yet!\n");
        exit(PSI_RETURN_FAILURE);

        break;
    case 45: /* Unpack pq and rs */
        /* I'm not sure if I'll ever need this, so I'm removing it for now */
        printf("\n\tShould you be using method %d?\n", method);
        exit(PSI_RETURN_FAILURE);
        /* Prepare the output buffer for the output DPD file */
        file4_mat_irrep_row_init(&(Buf->file), irrep);

        /* Loop over rows in the dpdfile */
        for(pq=0; pq < rowtot; pq++) {
            p = Buf->file.params->roworb[irrep][pq][0];
            q = Buf->file.params->roworb[irrep][pq][1];
            bufpq = Buf->params->rowidx[p][q];

            filerow = Buf->file.incore ? pq : 0;

            /* Loop over the columns in the dpdfile */
            for(rs=0; rs < coltot; rs++) {
                r = Buf->file.params->colorb[irrep^all_buf_irrep][rs][0];
                s = Buf->file.params->colorb[irrep^all_buf_irrep][rs][1];
                bufrs = Buf->params->colidx[r][s];

                value = Buf->matrix[irrep][bufpq][bufrs];

                /* Assign the value */
                Buf->file.matrix[irrep][filerow][rs] = value;
            }

            /* Write out the row */
            file4_mat_irrep_row_wrt(&(Buf->file), irrep, pq);
        }

        /* Close the input buffer */
        file4_mat_irrep_row_close(&(Buf->file), irrep);

        break;
    default:  /* Error trapping */
        printf("\n\tInvalid switch case in dpd_buf_mat_irrep_rd!\n");
        exit(PSI_RETURN_FAILURE);
        break;
    }

    return 0;

}

}
