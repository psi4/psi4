/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2024 The Psi4 Developers.
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

#include "dpd.h"

#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libpsi4util/exception.h"

#include <cstdio>
#include <cstdlib>

namespace psi {

int DPD::buf4_mat_irrep_row_rd(dpdbuf4 *Buf, int irrep, int pq) {
    int method = 0;
    int filerow, all_buf_irrep;
    int rs;                     /* dpdbuf row and column indices */
    int p, q, r, s;             /* orbital indices */
    int filepq, filers, filesr; /* Input dpdfile row and column indices */
    int rowtot, coltot;         /* dpdbuf row and column dimensions */
    int pq_permute, permute;
    double value;

    all_buf_irrep = Buf->file.my_irrep;
    rowtot = Buf->params->rowtot[irrep];
    coltot = Buf->params->coltot[irrep ^ all_buf_irrep];

    const auto& b_perm_pq = Buf->params->perm_pq;
    const auto& b_perm_rs = Buf->params->perm_rs;
    const auto& f_perm_pq = Buf->file.params->perm_pq;
    const auto& f_perm_rs = Buf->file.params->perm_rs;
    const auto& b_peq = Buf->params->peq;
    const auto& b_res = Buf->params->res;
    const auto& f_peq = Buf->file.params->peq;
    const auto& f_res = Buf->file.params->res;

    const auto& AllPolicy = dpdparams4::DiagPolicy::All;

    if ((b_perm_pq == f_perm_pq) && (b_perm_rs == f_perm_rs) && (b_peq == f_peq) && (b_res == f_res)) {
        if (Buf->anti)
            method = 11;
        else
            method = 12;
    } else if ((b_perm_pq != f_perm_pq) && (b_perm_rs == f_perm_rs) && (b_res == f_res)) {
        if (b_perm_pq == AllPolicy) {
            if (Buf->anti) {
                outfile->Printf("\n\tUnpack pq and antisymmetrize?\n");
                throw PSIEXCEPTION("Unpack pq and antisymmetrize?");
            }
            method = 21;
        } else if (f_perm_pq == AllPolicy) {
            if (Buf->anti)
                method = 22;
            else
                method = 23;
        } else {
            outfile->Printf("\n\tInvalid second-level method!\n");
            throw PSIEXCEPTION("Invalid second-level method!");
        }
    } else if ((b_perm_pq == f_perm_pq) && (b_perm_rs != f_perm_rs) && (b_peq == f_peq)) {
        if (b_perm_rs == AllPolicy) {
            if (Buf->anti) {
                outfile->Printf("\n\tUnpack rs and antisymmetrize?\n");
                throw PSIEXCEPTION("Unpack rs and antisymmetrize?");
            }
            method = 31;
        } else if (f_perm_rs == AllPolicy) {
            if (Buf->anti)
                method = 32;
            else
                method = 33;
        } else {
            outfile->Printf("\n\tInvalid third-level method!\n");
            throw PSIEXCEPTION("Invalid third-level method!");
        }
    } else if ((b_perm_pq != f_perm_pq) && (b_perm_rs != f_perm_rs)) {
        if (b_perm_pq == AllPolicy) {
            if (b_perm_rs == AllPolicy) {
                if (Buf->anti) {
                    outfile->Printf("\n\tUnpack pq and rs and antisymmetrize?\n");
                    throw PSIEXCEPTION("Unpack pq and rs and antisymmetrize?");
                } else
                    method = 41;
            } else if (f_perm_rs == AllPolicy) {
                if (Buf->anti) {
                    outfile->Printf("\n\tUnpack pq and antisymmetrize?\n");
                    throw PSIEXCEPTION("Unpack pq and antisymmetrize?");
                } else
                    method = 42;
            }
        } else if (f_perm_pq == AllPolicy) {
            if (b_perm_rs == AllPolicy) {
                if (Buf->anti) {
                    outfile->Printf("\n\tUnpack rs and antisymmetrize?\n");
                    throw PSIEXCEPTION("Unpack rs and antisymmetrize?");
                } else
                    method = 43;
            } else if (f_perm_rs == AllPolicy) {
                if (Buf->anti)
                    method = 44;
                else
                    method = 45;
            }
        } else {
            outfile->Printf("\n\tInvalid fourth-level method!\n");
            throw PSIEXCEPTION("Invalid fourth-level method!");
        }
    } else {
        outfile->Printf("\n\tInvalid method in dpd_buf_mat_irrep_rd!\n");
        throw PSIEXCEPTION("Invalid method in dpd_buf_mat_irrep_rd!");
    }

    switch (method) {
        case 11: /* No change in pq or rs; antisymmetrize */

            filerow = Buf->file.incore ? pq : 0;

            /* Prepare the input buffer from the input file */
            file4_mat_irrep_row_init(&(Buf->file), irrep);

            /* Fill the buffer */
            file4_mat_irrep_row_rd(&(Buf->file), irrep, pq);

            /* Loop over the columns in the dpdbuf */
            for (rs = 0; rs < coltot; rs++) {
                r = Buf->params->colorb[irrep ^ all_buf_irrep][rs][0];
                s = Buf->params->colorb[irrep ^ all_buf_irrep][rs][1];

                /* Column indices in the dpdfile */
                filers = rs;
                filesr = Buf->file.params->colidx[s][r];

                value = Buf->file.matrix[irrep][filerow][filers];
                value -= Buf->file.matrix[irrep][filerow][filesr];

                /* Assign the value */
                Buf->matrix[irrep][0][rs] = value;
            }

            /* Close the input buffer */
            file4_mat_irrep_row_close(&(Buf->file), irrep);

            break;

        case 12: /* No change in pq or rs */

            if (Buf->file.incore)
                for (rs = 0; rs < coltot; rs++) Buf->matrix[irrep][0][rs] = Buf->file.matrix[irrep][pq][rs];
            else {
                Buf->file.matrix[irrep] = Buf->matrix[irrep];
                file4_mat_irrep_row_rd(&(Buf->file), irrep, pq);
            }

            break;
        case 21: /* Unpack pq; no change in rs */
            /* Prepare the input buffer from the input file */
            file4_mat_irrep_row_init(&(Buf->file), irrep);

            p = Buf->params->roworb[irrep][pq][0];
            q = Buf->params->roworb[irrep][pq][1];
            filepq = Buf->file.params->rowidx[p][q];

            filerow = Buf->file.incore ? filepq : 0;

            /* Set the permutation operator's value */
            permute = ((p < q) && (f_perm_pq == dpdparams4::DiagPolicy::AntiSymm) ? -1 : 1);

            /* Fill the buffer */
            if (filepq >= 0)
                file4_mat_irrep_row_rd(&(Buf->file), irrep, filepq);
            else
                file4_mat_irrep_row_zero(&(Buf->file), irrep, filepq);

            /* Loop over the columns in the dpdbuf */
            for (rs = 0; rs < coltot; rs++) {
                filers = rs;

                value = Buf->file.matrix[irrep][filerow][filers];

                /* Assign the value, keeping track of the sign */
                Buf->matrix[irrep][0][rs] = permute * value;
            }

            /* Close the input buffer */
            file4_mat_irrep_row_close(&(Buf->file), irrep);

            break;
        case 22: /* Pack pq; no change in rs; antisymmetrize */
            /* Prepare the input buffer from the input file */
            file4_mat_irrep_row_init(&(Buf->file), irrep);

            p = Buf->params->roworb[irrep][pq][0];
            q = Buf->params->roworb[irrep][pq][1];
            filepq = Buf->file.params->rowidx[p][q];

            filerow = Buf->file.incore ? filepq : 0;

            /* Fill the buffer */
            file4_mat_irrep_row_rd(&(Buf->file), irrep, filepq);

            /* Loop over the columns in the dpdbuf */
            for (rs = 0; rs < coltot; rs++) {
                r = Buf->params->colorb[irrep ^ all_buf_irrep][rs][0];
                s = Buf->params->colorb[irrep ^ all_buf_irrep][rs][1];

                /* Column indices in the dpdfile */
                filers = rs;
                filesr = Buf->file.params->colidx[s][r];

                value = Buf->file.matrix[irrep][filerow][filers];

                value -= Buf->file.matrix[irrep][filerow][filesr];

                /* Assign the value */
                Buf->matrix[irrep][0][rs] = value;
            }

            /* Close the input buffer */
            file4_mat_irrep_row_close(&(Buf->file), irrep);

            break;
        case 23: /* Pack pq; no change in rs */
            /* Prepare the input buffer from the input file */
            file4_mat_irrep_row_init(&(Buf->file), irrep);

            p = Buf->params->roworb[irrep][pq][0];
            q = Buf->params->roworb[irrep][pq][1];
            filepq = Buf->file.params->rowidx[p][q];

            filerow = Buf->file.incore ? filepq : 0;

            file4_mat_irrep_row_rd(&(Buf->file), irrep, filepq);

            /* Loop over the columns in the dpdbuf */
            for (rs = 0; rs < coltot; rs++) {
                filers = rs;

                value = Buf->file.matrix[irrep][filerow][filers];

                /* Assign the value */
                Buf->matrix[irrep][0][rs] = value;
            }

            /* Close the input buffer */
            file4_mat_irrep_row_close(&(Buf->file), irrep);

            break;
        case 31: /* No change in pq; unpack rs */
            /* Prepare the input buffer from the input file */
            file4_mat_irrep_row_init(&(Buf->file), irrep);

            filepq = pq;

            filerow = Buf->file.incore ? filepq : 0;

            /* Fill the buffer */
            file4_mat_irrep_row_rd(&(Buf->file), irrep, filepq);

            /* Loop over the columns in the dpdbuf */
            for (rs = 0; rs < coltot; rs++) {
                r = Buf->params->colorb[irrep ^ all_buf_irrep][rs][0];
                s = Buf->params->colorb[irrep ^ all_buf_irrep][rs][1];
                filers = Buf->file.params->colidx[r][s];

                /* rs permutation operator */
                permute = ((r < s) && (f_perm_rs == dpdparams4::DiagPolicy::AntiSymm) ? -1 : 1);

                /* Is this fast enough? */
                value = ((filers < 0) ? 0 : Buf->file.matrix[irrep][filerow][filers]);

                /* Assign the value */
                Buf->matrix[irrep][0][rs] = permute * value;
            }

            /* Close the input buffer */
            file4_mat_irrep_row_close(&(Buf->file), irrep);

            break;
        case 32: /* No change in pq; pack rs; antisymmetrize */
            /* Prepare the input buffer from the input file */
            file4_mat_irrep_row_init(&(Buf->file), irrep);

            filepq = pq;

            filerow = Buf->file.incore ? filepq : 0;

            /* Fill the buffer */
            file4_mat_irrep_row_rd(&(Buf->file), irrep, filepq);

            /* Loop over the columns in the dpdbuf */
            for (rs = 0; rs < coltot; rs++) {
                r = Buf->params->colorb[irrep ^ all_buf_irrep][rs][0];
                s = Buf->params->colorb[irrep ^ all_buf_irrep][rs][1];

                /* Column indices in the dpdfile */
                filers = Buf->file.params->colidx[r][s];
                filesr = Buf->file.params->colidx[s][r];

                value = Buf->file.matrix[irrep][filerow][filers];
                value -= Buf->file.matrix[irrep][filerow][filesr];

                /* Assign the value */
                Buf->matrix[irrep][0][rs] = value;
            }

            /* Close the input buffer */
            file4_mat_irrep_row_close(&(Buf->file), irrep);

            break;
        case 33: /* No change in pq; pack rs */
            /* Prepare the input buffer from the input file */
            file4_mat_irrep_row_init(&(Buf->file), irrep);

            filepq = pq;
            filerow = Buf->file.incore ? filepq : 0;

            /* Fill the buffer */
            file4_mat_irrep_row_rd(&(Buf->file), irrep, filepq);

            /* Loop over the columns in the dpdbuf */
            for (rs = 0; rs < coltot; rs++) {
                r = Buf->params->colorb[irrep ^ all_buf_irrep][rs][0];
                s = Buf->params->colorb[irrep ^ all_buf_irrep][rs][1];
                filers = Buf->file.params->colidx[r][s];

                value = Buf->file.matrix[irrep][filerow][filers];

                /* Assign the value */
                Buf->matrix[irrep][0][rs] = value;
            }

            /* Close the input buffer */
            file4_mat_irrep_row_close(&(Buf->file), irrep);

            break;
        case 41: /* Unpack pq and rs */
            /* Prepare the input buffer from the input file */
            file4_mat_irrep_row_init(&(Buf->file), irrep);

            p = Buf->params->roworb[irrep][pq][0];
            q = Buf->params->roworb[irrep][pq][1];
            filepq = Buf->file.params->rowidx[p][q];

            filerow = Buf->file.incore ? filepq : 0;

            /* Set the value of the pq permutation operator */
            pq_permute = ((p < q) && (f_perm_pq == dpdparams4::DiagPolicy::AntiSymm) ? -1 : 1);

            /* Fill the buffer */
            if (filepq >= 0)
                file4_mat_irrep_row_rd(&(Buf->file), irrep, filepq);
            else
                file4_mat_irrep_row_zero(&(Buf->file), irrep, filepq);

            /* Loop over the columns in the dpdbuf */
            for (rs = 0; rs < coltot; rs++) {
                r = Buf->params->colorb[irrep ^ all_buf_irrep][rs][0];
                s = Buf->params->colorb[irrep ^ all_buf_irrep][rs][1];
                filers = Buf->file.params->colidx[r][s];

                /* Set the value of the pqrs permutation operator */
                permute = ((r < s) && (f_perm_rs == dpdparams4::DiagPolicy::AntiSymm) ? -1 : 1) * pq_permute;

                value = 0;

                if (filers >= 0) value = Buf->file.matrix[irrep][filerow][filers];

                /* Assign the value */
                Buf->matrix[irrep][0][rs] = permute * value;
            }

            /* Close the input buffer */
            file4_mat_irrep_row_close(&(Buf->file), irrep);

            break;
        case 42: /* Pack pq; unpack rs */
            outfile->Printf("\n\tHaven't programmed method 42 yet!\n");
            throw PSIEXCEPTION("Haven't programmed method 42 yet!");

            break;
        case 43: /* Unpack pq; pack rs */
            outfile->Printf("\n\tHaven't programmed method 43 yet!\n");
            throw PSIEXCEPTION("Haven't programmed method 43 yet!");

            break;
        case 44: /* Pack pq; pack rs; antisymmetrize */
            /* Prepare the input buffer from the input file */
            file4_mat_irrep_row_init(&(Buf->file), irrep);

            p = Buf->params->roworb[irrep][pq][0];
            q = Buf->params->roworb[irrep][pq][1];
            filepq = Buf->file.params->rowidx[p][q];

            filerow = Buf->file.incore ? filepq : 0;

            /* Fill the buffer */
            file4_mat_irrep_row_rd(&(Buf->file), irrep, filepq);

            /* Loop over the columns in the dpdbuf */
            for (rs = 0; rs < coltot; rs++) {
                r = Buf->params->colorb[irrep ^ all_buf_irrep][rs][0];
                s = Buf->params->colorb[irrep ^ all_buf_irrep][rs][1];

                /* Column indices in the dpdfile */
                filers = Buf->file.params->colidx[r][s];
                filesr = Buf->file.params->colidx[s][r];

                value = Buf->file.matrix[irrep][filerow][filers];
                value -= Buf->file.matrix[irrep][filerow][filesr];

                /* Assign the value */
                Buf->matrix[irrep][0][rs] = value;
            }

            /* Close the input buffer */
            file4_mat_irrep_row_close(&(Buf->file), irrep);

            break;
        case 45: /* Pack pq and rs */
            /* Prepare the input buffer from the input file */
            file4_mat_irrep_row_init(&(Buf->file), irrep);

            p = Buf->params->roworb[irrep][pq][0];
            q = Buf->params->roworb[irrep][pq][1];
            filepq = Buf->file.params->rowidx[p][q];

            filerow = Buf->file.incore ? filepq : 0;

            file4_mat_irrep_row_rd(&(Buf->file), irrep, filepq);

            /* Loop over the columns in the dpdbuf */
            for (rs = 0; rs < coltot; rs++) {
                r = Buf->params->colorb[irrep ^ all_buf_irrep][rs][0];
                s = Buf->params->colorb[irrep ^ all_buf_irrep][rs][1];
                filers = Buf->file.params->colidx[r][s];

                if (filers < 0) {
                    outfile->Printf("\n\tNegative colidx in method 44?\n");
                    throw PSIEXCEPTION("Negative colidx in method 44?");
                }

                value = Buf->file.matrix[irrep][filerow][filers];

                /* Assign the value */
                Buf->matrix[irrep][0][rs] = value;
            }

            /* Close the input buffer */
            file4_mat_irrep_row_close(&(Buf->file), irrep);

            break;
        default: /* Error trapping */
            outfile->Printf("\n\tInvalid switch case in dpd_buf_mat_irrep_rd!\n");
            throw PSIEXCEPTION("Invalid switch case in dpd_buf_mat_irrep_rd!");
            break;
    }

    return 0;
}

}  // namespace psi
