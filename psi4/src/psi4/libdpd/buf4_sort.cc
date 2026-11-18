/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2026 The Psi4 Developers.
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

#include "psi4/libqt/qt.h"
#include "psi4/psi4-dec.h"
#include "psi4/libpsi4util/PsiOutStream.h"

#include <algorithm>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <cmath>

using std::string;
namespace psi {

/*
** dpd_buf4_sort(): A general DPD buffer sorting function that will
** (eventually) handle all 24 possible permutations of four-index
** buffers.  This code assumes that all symmetry-blocks of both source and
** target buffers can be stored in core.
**
** Arguments:
**   dpdbuf4 *InBuf: A pointer to the alread-initialized input
**     buffer.
**   int outfilenum: The PSI unit number for the target data.
**   enum indices index: The desired sorting pattern (see dpd.h).
**   int pqnum: The index combination for the bra indices for the new
**     dpd file4.  N.B. this is NOT error checked for consistency
**     with the input buffer.
**   int rsnum: The index combination for the ket indices for the new
**     dpd file4.  N.B. this is NOT error checked for consistency
**     with the input buffer.
**   std::string &label: A string labelling for this buffer.
**
** Note that the notation for some of these multi-index sorts may be
** confusing.  When the caller requests ordering "psqr" for example,
** the desired arrangement returned is indeed Out[ps][qr] =
** In[pq][rs].  However, the index notation *within the code below*
** will appear as Out[pq][rs] = In[pr][sq]. To compute this easily,
** take the desired ordering psqr, translate its indices to p->p,
** s->q, q->r, and r->s, and finally write the source indices as
** prsq.  (This "problem" arises because I want to use the same index
** notation Out[pq][rs] for the target for every rearrangement.  For
** all pair permutations, the notation is straightforward. )
**
** Note that for pqsr, qprs, and rspq (others?), we assume that
** the InBuf and OutBuf have identical row/column orderings for all
** unswapped indices.  For example, for pqsr, we assume that the
** pq-ordering of both buffers is identical.
**
** -Daniel, December 1998
**
** Modified for new naming conventions and non-totally-symmetric data.
** TDC
** September 1999
**
** NB: Timing tests have indicated that sorts which mix bra and ket
** indices are *substantially* more expensive that those which either
** transpose bra and ket or simply mix bra or mix ket indices separately.
** The source of the problem is the multiple buf4_mat_irrep_rd() calls
** required for the b-k mixing sorts (cf. pqsr vs. prqs).  If possible, one
** should arrange contractions to use fewer b-k mixing sorts.
** TDC
** May 2000
**
** Added fully out-of-core (multipass) sorting algorithm to qpsr
** sorting case. More cases will follow as I need them.
**
** -TDC, April 2005
**
** the enum-argument labelling is used in this list
** IC=in-core capable; OOC=out-of-core capable
** pqrs: error  ** pqsr: IC/OOC
** prqs: IC/OOC ** prsq: IC/OOC
** psqr: IC     ** psrq: IC
** qprs: IC/OOC ** qpsr: IC/OOC
** qrps: IC     ** qrsp: IC
** qspr: IC     ** qsrp: IC
** rqps: IC     ** rqsp: IC
** rpqs: IC     ** rpsq: IC
** rsqp: IC     ** rspq: IC/OOC
** sqrp: IC     ** sqpr: none
** srqp: IC     ** srpq: IC
** spqr: IC     ** sprq: IC
** -RAK, Nov. 2005*/

int DPD::buf4_sort(dpdbuf4 *InBuf, int outfilenum, enum indices index, int pqnum, int rsnum, const std::string& label) {
    int h, nirreps, row, col, my_irrep, r_irrep;
    int p, q, r, s, P, Q, R, S, pq, rs, sr, pr, qs, qp, rq, qr, ps, sp, rp, sq;
    int PQ, RS;
    int Gp, Gq, Gr, Gs, Gpq, Grs, Gpr, Gqs, Grq, Gqr, Gps, Gsp, Grp, Gsq;
    dpdbuf4 OutBuf;
    long int rowtot, coltot, core_total, maxrows;
    int Grow, Gcol;
    int64_t out_rows_per_bucket, in_rows_per_bucket, rows_per_bucket;
    int out_nbuckets, out_rows_left, out_row_start, n;
    int in_nbuckets, in_rows_left, in_row_start, m;
    int nbuckets, rows_left;

    nirreps = InBuf->params->nirreps;
    my_irrep = InBuf->file.my_irrep;

#ifdef DPD_TIMER
    timer_on("buf4_sort");
#endif

    buf4_init(&OutBuf, outfilenum, my_irrep, pqnum, rsnum, pqnum, rsnum, 0, label);

    /* select in-core vs. out-of-core algorithms */
    bool incore = true;
    core_total = 0;
    for (h = 0; h < nirreps; h++) {
        coltot = InBuf->params->coltot[h ^ my_irrep];
        if (coltot) {
            maxrows = DPD_BIGNUM / coltot;
            if (maxrows < 1) {
                outfile->Printf("\nLIBDPD Error: too many rows in buf4_sort_axpy.\n");
                dpd_error("buf4_sort_axpy", "outfile");
            }
        } else
            maxrows = DPD_BIGNUM;
        rowtot = InBuf->params->rowtot[h];
        for (; rowtot > maxrows; rowtot -= maxrows) {
            if (core_total > (core_total + 2 * maxrows * coltot))
                incore = false;
            else
                core_total += 2 * maxrows * coltot;
        }
        if (core_total > (core_total + 2 * rowtot * coltot)) incore = false;
        core_total += 2 * rowtot * coltot;
    }
    if (core_total > dpd_memfree()) incore = false;

#ifdef DPD_DEBUG
    if (incore == false) {
        switch (index) {
            case (pqsr):
                printf("Doing out-of-core pqsr sort.\n");
                break;
            case (prqs):
                printf("Doing out-of-core prqs sort.\n");
                break;
            case (prsq):
                printf("Doing out-of-core prsq sort.\n");
                break;
            case (qprs):
                printf("Doing out-of-core qprs sort.\n");
                break;
            case (qpsr):
                printf("Doing out-of-core qpsr sort.\n");
                break;
            case (sqpr):
                printf("Doing out-of-core sqpr sort.\n");
                break;
            case (rspq):
                printf("Doing out-of-core rspq sort.\n");
                break;
        }
    }
#endif

#ifdef ALL_BUF4_SORT_OOC
    switch (index) {
        case (pqsr):
        case (prqs):
        case (prsq):
        case (qprs):
        case (qpsr):
        case (sqpr):
        case (rspq):
            incore = false;
            break;
    }
#endif

    if (!incore) note_ooc_sort("buf4_sort", index);

    /* Init input and output buffers and read in all blocks of the input */
    if (incore) {
        for (h = 0; h < nirreps; h++) {
            buf4_mat_irrep_init(&OutBuf, h);
            buf4_mat_irrep_init(InBuf, h);
            buf4_mat_irrep_rd(InBuf, h);
        }
    }

    switch (index) {
        case pqrs:
            outfile->Printf("\nDPD sort error: invalid index ordering.\n");
            dpd_error("buf_sort", "outfile");
            break;

        case pqsr:

#ifdef DPD_TIMER
            timer_on("pqsr");
#endif

            /* p->p; q->q; s->r; r->s = pqsr */
            if (incore) {
                for (h = 0; h < nirreps; h++) {
                    r_irrep = h ^ my_irrep;

                    for (pq = 0; pq < OutBuf.params->rowtot[h]; pq++) {
                        p = OutBuf.params->roworb[h][pq][0];
                        q = OutBuf.params->roworb[h][pq][1];

                        row = InBuf->params->rowidx[p][q];

                        for (rs = 0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
                            r = OutBuf.params->colorb[r_irrep][rs][0];
                            s = OutBuf.params->colorb[r_irrep][rs][1];

                            sr = InBuf->params->colidx[s][r];

                            OutBuf.matrix[h][pq][rs] = InBuf->matrix[h][row][sr];
                        }
                    }
                }
            } else {
                buf4_sort_ooc_pairs(InBuf, &OutBuf, SortPair::PQ, SortPair::SR, false, 1.0,
                                    "buf4_sort (pqsr): not enough memory for one row");
            }

#ifdef DPD_TIMER
            timer_off("pqsr");
#endif
            break;

        case prqs:

#ifdef DPD_TIMER
            timer_on("prqs");
#endif

            /* p->p; r->q; q->r; s->s = prqs */
            if (incore) {
                for (h = 0; h < nirreps; h++) {
                    r_irrep = h ^ my_irrep;

                    for (Gp = 0; Gp < nirreps; Gp++) {
                        Gq = Gp ^ h;
                        for (Gr = 0; Gr < nirreps; Gr++) {
                            Gs = Gr ^ r_irrep;

                            /* Irreps on the source */
                            Gpr = Gp ^ Gr;
                            Gqs = Gq ^ Gs;

                            for (p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                                P = OutBuf.params->poff[Gp] + p;
                                for (q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                                    Q = OutBuf.params->qoff[Gq] + q;
                                    pq = OutBuf.params->rowidx[P][Q];

                                    for (r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                                        R = OutBuf.params->roff[Gr] + r;
                                        pr = InBuf->params->rowidx[P][R];

                                        for (s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                            S = OutBuf.params->soff[Gs] + s;
                                            rs = OutBuf.params->colidx[R][S];
                                            qs = InBuf->params->colidx[Q][S];

                                            OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gpr][pr][qs];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } else { /* pqrs <- prqs */

                for (Gpq = 0; Gpq < nirreps; Gpq++) {
                    Grs = Gpq ^ my_irrep;
                    /* an empty or zero-width block has no elements to sort (and a null block matrix) */
                    if (!OutBuf.params->rowtot[Gpq] || !OutBuf.params->coltot[Grs]) continue;

                    out_rows_per_bucket = dpd_rows_per_bucket_or_die(dpd_memfree(), 2LL * OutBuf.params->coltot[Grs],
                                                                     OutBuf.params->rowtot[Gpq],
                                                                     "buf4_sort (prqs): not enough memory for one row");
                    out_nbuckets = (int)ceil((double)OutBuf.params->rowtot[Gpq] / (double)out_rows_per_bucket);
                    if (out_nbuckets == 1)
                        out_rows_left = out_rows_per_bucket;
                    else
                        out_rows_left = OutBuf.params->rowtot[Gpq] % out_rows_per_bucket;

                    /* allocate space for the bucket of rows */
                    buf4_mat_irrep_init_block(&OutBuf, Gpq, out_rows_per_bucket);

                    for (n = 0; n < (out_rows_left ? out_nbuckets - 1 : out_nbuckets); n++) {
                        out_row_start = n * out_rows_per_bucket;

                        for (Grow = 0; Grow < nirreps; Grow++) { /*Grow = Gpr*/
                            Gcol = Grow ^ my_irrep;              /*Gcol = Gqs*/

                            /* determine how many rows of InBuf we can store in the other half of the core */
                            in_rows_per_bucket =
                                dpd_rows_per_bucket_or_die(dpd_memfree(), 2LL * InBuf->params->coltot[Gcol],
                                                           InBuf->params->rowtot[Grow],
                                                           "buf4_sort (prqs): not enough memory for one row");
                            in_nbuckets = (int)ceil((double)InBuf->params->rowtot[Grow] / (double)in_rows_per_bucket);
                            if (in_nbuckets == 1)
                                in_rows_left = in_rows_per_bucket;
                            else
                                in_rows_left = InBuf->params->rowtot[Grow] % in_rows_per_bucket;

                            /* allocate space for the bucket of rows */
                            buf4_mat_irrep_init_block(InBuf, Grow, in_rows_per_bucket);

                            /* pqrs <- prqs */
                            for (m = 0; m < (in_rows_left ? in_nbuckets - 1 : in_nbuckets); m++) {
                                in_row_start = m * in_rows_per_bucket;
                                buf4_mat_irrep_rd_block(InBuf, Grow, in_row_start, in_rows_per_bucket);

                                for (pq = 0; pq < out_rows_per_bucket; pq++) {
                                    p = OutBuf.params->roworb[Gpq][pq + out_row_start][0];
                                    q = OutBuf.params->roworb[Gpq][pq + out_row_start][1];
                                    Gp = OutBuf.params->psym[p];
                                    Gq = Gpq ^ Gp;
                                    for (rs = 0; rs < OutBuf.params->coltot[Grs]; rs++) {
                                        r = OutBuf.params->colorb[Grs][rs][0];
                                        s = OutBuf.params->colorb[Grs][rs][1];
                                        Gr = OutBuf.params->rsym[r];
                                        Gs = Grs ^ Gr;
                                        Gpr = Gp ^ Gr;

                                        if (Gpr == Grow) {
                                            pr = InBuf->params->rowidx[p][r] - in_row_start;
                                            /* check if the current value is in the current in_bucket or not */
                                            if (pr >= 0 && pr < in_rows_per_bucket) {
                                                qs = InBuf->params->colidx[q][s];
                                                OutBuf.matrix[Gpq][pq][rs] = InBuf->matrix[Grow][pr][qs];
                                            }
                                        }
                                    }
                                }
                            }
                            if (in_rows_left) {
                                in_row_start = m * in_rows_per_bucket;
                                buf4_mat_irrep_rd_block(InBuf, Grow, in_row_start, in_rows_left);

                                /* pqrs <- prqs */
                                for (pq = 0; pq < out_rows_per_bucket; pq++) {
                                    p = OutBuf.params->roworb[Gpq][pq + out_row_start][0];
                                    q = OutBuf.params->roworb[Gpq][pq + out_row_start][1];
                                    Gp = OutBuf.params->psym[p];
                                    Gq = Gpq ^ Gp;
                                    for (rs = 0; rs < OutBuf.params->coltot[Grs]; rs++) {
                                        r = OutBuf.params->colorb[Grs][rs][0];
                                        s = OutBuf.params->colorb[Grs][rs][1];
                                        Gr = OutBuf.params->rsym[r];
                                        Gs = Grs ^ Gr;
                                        Gpr = Gp ^ Gr;

                                        if (Gpr == Grow) {
                                            pr = InBuf->params->rowidx[p][r] - in_row_start;
                                            /* check if the current value is in core or not */
                                            if (pr >= 0 && pr < in_rows_left) {
                                                qs = InBuf->params->colidx[q][s];
                                                OutBuf.matrix[Gpq][pq][rs] = InBuf->matrix[Grow][pr][qs];
                                            }
                                        }
                                    }
                                }
                            }
                            buf4_mat_irrep_close_block(InBuf, Grow, in_rows_per_bucket);
                        }
                        buf4_mat_irrep_wrt_block(&OutBuf, Gpq, out_row_start, out_rows_per_bucket);
                    }
                    if (out_rows_left) {
                        out_row_start = n * out_rows_per_bucket;

                        for (Grow = 0; Grow < nirreps; Grow++) {
                            Gcol = Grow ^ my_irrep;

                            /* determine how many rows of InBuf we can store in the other half of the core */
                            in_rows_per_bucket =
                                dpd_rows_per_bucket_or_die(dpd_memfree(), 2LL * InBuf->params->coltot[Gcol],
                                                           InBuf->params->rowtot[Grow],
                                                           "buf4_sort (prqs): not enough memory for one row");
                            in_nbuckets = (int)ceil((double)InBuf->params->rowtot[Grow] / (double)in_rows_per_bucket);
                            if (in_nbuckets == 1)
                                in_rows_left = in_rows_per_bucket;
                            else
                                in_rows_left = InBuf->params->rowtot[Grow] % in_rows_per_bucket;

                            /* allocate space for the bucket of rows */
                            buf4_mat_irrep_init_block(InBuf, Grow, in_rows_per_bucket);

                            /* pqrs <- prqs */
                            for (m = 0; m < (in_rows_left ? in_nbuckets - 1 : in_nbuckets); m++) {
                                in_row_start = m * in_rows_per_bucket;
                                buf4_mat_irrep_rd_block(InBuf, Grow, in_row_start, in_rows_per_bucket);

                                for (pq = 0; pq < out_rows_left; pq++) {
                                    p = OutBuf.params->roworb[Gpq][pq + out_row_start][0];
                                    q = OutBuf.params->roworb[Gpq][pq + out_row_start][1];
                                    Gp = OutBuf.params->psym[p];
                                    Gq = Gpq ^ Gp;
                                    for (rs = 0; rs < OutBuf.params->coltot[Grs]; rs++) {
                                        r = OutBuf.params->colorb[Grs][rs][0];
                                        s = OutBuf.params->colorb[Grs][rs][1];
                                        Gr = OutBuf.params->rsym[r];
                                        Gs = Grs ^ Gr;
                                        Gpr = Gp ^ Gr;

                                        if (Gpr == Grow) {
                                            pr = InBuf->params->rowidx[p][r] - in_row_start;
                                            /* check if the current value is in the current in_bucket or not */
                                            if (pr >= 0 && pr < in_rows_per_bucket) {
                                                qs = InBuf->params->colidx[q][s];
                                                OutBuf.matrix[Gpq][pq][rs] = InBuf->matrix[Grow][pr][qs];
                                            }
                                        }
                                    }
                                }
                            }
                            if (in_rows_left) {
                                in_row_start = m * in_rows_per_bucket;
                                buf4_mat_irrep_rd_block(InBuf, Grow, in_row_start, in_rows_left);

                                for (pq = 0; pq < out_rows_left; pq++) {
                                    p = OutBuf.params->roworb[Gpq][pq + out_row_start][0];
                                    q = OutBuf.params->roworb[Gpq][pq + out_row_start][1];
                                    Gp = OutBuf.params->psym[p];
                                    Gq = Gpq ^ Gp;
                                    for (rs = 0; rs < OutBuf.params->coltot[Grs]; rs++) {
                                        r = OutBuf.params->colorb[Grs][rs][0];
                                        s = OutBuf.params->colorb[Grs][rs][1];
                                        Gr = OutBuf.params->rsym[r];
                                        Gs = Grs ^ Gr;
                                        Gpr = Gp ^ Gr;

                                        if (Gpr == Grow) {
                                            pr = InBuf->params->rowidx[p][r] - in_row_start;
                                            /* check if the current value is in core or not */
                                            if (pr >= 0 && pr < in_rows_left) {
                                                qs = InBuf->params->colidx[q][s];
                                                OutBuf.matrix[Gpq][pq][rs] = InBuf->matrix[Grow][pr][qs];
                                            }
                                        }
                                    }
                                }
                            }

                            buf4_mat_irrep_close_block(InBuf, Grow, in_rows_per_bucket);
                        }

                        buf4_mat_irrep_wrt_block(&OutBuf, Gpq, out_row_start, out_rows_left);
                    }

                    buf4_mat_irrep_close_block(&OutBuf, Gpq, out_rows_per_bucket);
                }
            }

#ifdef DPD_TIMER
            timer_off("prqs");
#endif
            break;

        case prsq:

#ifdef DPD_TIMER
            timer_on("prsq");
#endif

            /* p->p; r->q; s->r; q->s = psqr */

            if (incore) {
                for (h = 0; h < nirreps; h++) {
                    r_irrep = h ^ my_irrep;

                    for (Gp = 0; Gp < nirreps; Gp++) {
                        Gq = Gp ^ h;
                        for (Gr = 0; Gr < nirreps; Gr++) {
                            Gs = Gr ^ r_irrep;

                            Gps = Gp ^ Gs;
                            Gqr = Gq ^ Gr;

                            for (p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                                P = OutBuf.params->poff[Gp] + p;
                                for (q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                                    Q = OutBuf.params->qoff[Gq] + q;
                                    pq = OutBuf.params->rowidx[P][Q];

                                    for (r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                                        R = OutBuf.params->roff[Gr] + r;
                                        qr = InBuf->params->colidx[Q][R];

                                        for (s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                            S = OutBuf.params->soff[Gs] + s;
                                            rs = OutBuf.params->colidx[R][S];
                                            ps = InBuf->params->rowidx[P][S];

                                            OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gps][ps][qr];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                for (Gpq = 0; Gpq < nirreps; Gpq++) {
                    Grs = Gpq ^ my_irrep;
                    /* an empty or zero-width block has no elements to sort (and a null block matrix) */
                    if (!OutBuf.params->rowtot[Gpq] || !OutBuf.params->coltot[Grs]) continue;

                    /* determine how many rows of OutBuf we can store in half of the core */
                    out_rows_per_bucket = dpd_rows_per_bucket_or_die(dpd_memfree(), 2LL * OutBuf.params->coltot[Grs],
                                                                     OutBuf.params->rowtot[Gpq],
                                                                     "buf4_sort (prsq): not enough memory for one row");
                    out_nbuckets = (int)ceil((double)OutBuf.params->rowtot[Gpq] / (double)out_rows_per_bucket);
                    if (out_nbuckets == 1)
                        out_rows_left = out_rows_per_bucket;
                    else
                        out_rows_left = OutBuf.params->rowtot[Gpq] % out_rows_per_bucket;

                    /* allocate space for the bucket of rows */
                    buf4_mat_irrep_init_block(&OutBuf, Gpq, out_rows_per_bucket);

                    for (n = 0; n < (out_rows_left ? out_nbuckets - 1 : out_nbuckets); n++) {
                        out_row_start = n * out_rows_per_bucket;

                        for (Grow = 0; Grow < nirreps; Grow++) {
                            Gcol = Grow ^ my_irrep;

                            /* determine how many rows of InBuf we can store in the other half of the core */
                            in_rows_per_bucket =
                                dpd_rows_per_bucket_or_die(dpd_memfree(), 2LL * InBuf->params->coltot[Gcol],
                                                           InBuf->params->rowtot[Grow],
                                                           "buf4_sort (prsq): not enough memory for one row");
                            in_nbuckets = (int)ceil((double)InBuf->params->rowtot[Grow] / (double)in_rows_per_bucket);
                            if (in_nbuckets == 1)
                                in_rows_left = in_rows_per_bucket;
                            else
                                in_rows_left = InBuf->params->rowtot[Grow] % in_rows_per_bucket;

                            /* allocate space for the bucket of rows */
                            buf4_mat_irrep_init_block(InBuf, Grow, in_rows_per_bucket);

                            for (m = 0; m < (in_rows_left ? in_nbuckets - 1 : in_nbuckets); m++) {
                                in_row_start = m * in_rows_per_bucket;
                                buf4_mat_irrep_rd_block(InBuf, Grow, in_row_start, in_rows_per_bucket);

                                for (pq = 0; pq < out_rows_per_bucket; pq++) {
                                    p = OutBuf.params->roworb[Gpq][pq + out_row_start][0];
                                    q = OutBuf.params->roworb[Gpq][pq + out_row_start][1];
                                    Gp = OutBuf.params->psym[p];
                                    Gq = Gpq ^ Gp;
                                    for (rs = 0; rs < OutBuf.params->coltot[Grs]; rs++) {
                                        r = OutBuf.params->colorb[Grs][rs][0];
                                        s = OutBuf.params->colorb[Grs][rs][1];
                                        Gr = OutBuf.params->rsym[r];
                                        Gs = Grs ^ Gr;

                                        Gps = Gp ^ Gs;

                                        if (Gps == Grow) {
                                            ps = InBuf->params->rowidx[p][s] - in_row_start;
                                            /* check if the current value is in the current in_bucket or not */
                                            if (ps >= 0 && ps < in_rows_per_bucket) {
                                                qr = InBuf->params->colidx[q][r];
                                                OutBuf.matrix[Gpq][pq][rs] = InBuf->matrix[Grow][ps][qr];
                                            }
                                        }
                                    }
                                }
                            }
                            if (in_rows_left) {
                                in_row_start = m * in_rows_per_bucket;
                                buf4_mat_irrep_rd_block(InBuf, Grow, in_row_start, in_rows_left);

                                for (pq = 0; pq < out_rows_per_bucket; pq++) {
                                    p = OutBuf.params->roworb[Gpq][pq + out_row_start][0];
                                    q = OutBuf.params->roworb[Gpq][pq + out_row_start][1];
                                    Gp = OutBuf.params->psym[p];
                                    Gq = Gpq ^ Gp;
                                    for (rs = 0; rs < OutBuf.params->coltot[Grs]; rs++) {
                                        r = OutBuf.params->colorb[Grs][rs][0];
                                        s = OutBuf.params->colorb[Grs][rs][1];
                                        Gr = OutBuf.params->rsym[r];
                                        Gs = Grs ^ Gr;

                                        Gps = Gp ^ Gs;

                                        if (Gps == Grow) {
                                            ps = InBuf->params->rowidx[p][s] - in_row_start;
                                            /* check if the current value is in core or not */
                                            if (ps >= 0 && ps < in_rows_left) {
                                                qr = InBuf->params->colidx[q][r];
                                                OutBuf.matrix[Gpq][pq][rs] = InBuf->matrix[Grow][ps][qr];
                                            }
                                        }
                                    }
                                }
                            }
                            buf4_mat_irrep_close_block(InBuf, Grow, in_rows_per_bucket);
                        }
                        buf4_mat_irrep_wrt_block(&OutBuf, Gpq, out_row_start, out_rows_per_bucket);
                    }
                    if (out_rows_left) {
                        out_row_start = n * out_rows_per_bucket;

                        for (Grow = 0; Grow < nirreps; Grow++) {
                            Gcol = Grow ^ my_irrep;

                            /* determine how many rows of InBuf we can store in the other half of the core */
                            in_rows_per_bucket =
                                dpd_rows_per_bucket_or_die(dpd_memfree(), 2LL * InBuf->params->coltot[Gcol],
                                                           InBuf->params->rowtot[Grow],
                                                           "buf4_sort (prsq): not enough memory for one row");
                            in_nbuckets = (int)ceil((double)InBuf->params->rowtot[Grow] / (double)in_rows_per_bucket);
                            if (in_nbuckets == 1)
                                in_rows_left = in_rows_per_bucket;
                            else
                                in_rows_left = InBuf->params->rowtot[Grow] % in_rows_per_bucket;

                            /* allocate space for the bucket of rows */
                            buf4_mat_irrep_init_block(InBuf, Grow, in_rows_per_bucket);

                            for (m = 0; m < (in_rows_left ? in_nbuckets - 1 : in_nbuckets); m++) {
                                in_row_start = m * in_rows_per_bucket;
                                buf4_mat_irrep_rd_block(InBuf, Grow, in_row_start, in_rows_per_bucket);

                                for (pq = 0; pq < out_rows_left; pq++) {
                                    p = OutBuf.params->roworb[Gpq][pq + out_row_start][0];
                                    q = OutBuf.params->roworb[Gpq][pq + out_row_start][1];
                                    Gp = OutBuf.params->psym[p];
                                    Gq = Gpq ^ Gp;
                                    for (rs = 0; rs < OutBuf.params->coltot[Grs]; rs++) {
                                        r = OutBuf.params->colorb[Grs][rs][0];
                                        s = OutBuf.params->colorb[Grs][rs][1];
                                        Gr = OutBuf.params->rsym[r];
                                        Gs = Grs ^ Gr;

                                        Gps = Gp ^ Gs;

                                        if (Gps == Grow) {
                                            ps = InBuf->params->rowidx[p][s] - in_row_start;
                                            /* check if the current value is in the current in_bucket or not */
                                            if (ps >= 0 && ps < in_rows_per_bucket) {
                                                qr = InBuf->params->colidx[q][r];
                                                OutBuf.matrix[Gpq][pq][rs] = InBuf->matrix[Grow][ps][qr];
                                            }
                                        }
                                    }
                                }
                            }
                            if (in_rows_left) {
                                in_row_start = m * in_rows_per_bucket;
                                buf4_mat_irrep_rd_block(InBuf, Grow, in_row_start, in_rows_left);

                                for (pq = 0; pq < out_rows_left; pq++) {
                                    p = OutBuf.params->roworb[Gpq][pq + out_row_start][0];
                                    q = OutBuf.params->roworb[Gpq][pq + out_row_start][1];
                                    Gp = OutBuf.params->psym[p];
                                    Gq = Gpq ^ Gp;
                                    for (rs = 0; rs < OutBuf.params->coltot[Grs]; rs++) {
                                        r = OutBuf.params->colorb[Grs][rs][0];
                                        s = OutBuf.params->colorb[Grs][rs][1];
                                        Gr = OutBuf.params->rsym[r];
                                        Gs = Grs ^ Gr;

                                        Gps = Gp ^ Gs;

                                        if (Gps == Grow) {
                                            ps = InBuf->params->rowidx[p][s] - in_row_start;
                                            /* check if the current value is in core or not */
                                            if (ps >= 0 && ps < in_rows_left) {
                                                qr = InBuf->params->colidx[q][r];
                                                OutBuf.matrix[Gpq][pq][rs] = InBuf->matrix[Grow][ps][qr];
                                            }
                                        }
                                    }
                                }
                            }

                            buf4_mat_irrep_close_block(InBuf, Grow, in_rows_per_bucket);
                        }

                        buf4_mat_irrep_wrt_block(&OutBuf, Gpq, out_row_start, out_rows_left);
                    }

                    buf4_mat_irrep_close_block(&OutBuf, Gpq, out_rows_per_bucket);
                }
            }

#ifdef DPD_TIMER
            timer_off("prsq");
#endif
            break;

        case psqr:

#ifdef DPD_TIMER
            timer_on("psqr");
#endif

            /* p->p; s->q; q->r; r->s = prsq */

            if (incore) {
                for (h = 0; h < nirreps; h++) {
                    r_irrep = h ^ my_irrep;

                    for (Gp = 0; Gp < nirreps; Gp++) {
                        Gq = Gp ^ h;
                        for (Gr = 0; Gr < nirreps; Gr++) {
                            Gs = Gr ^ r_irrep;

                            Gpr = Gp ^ Gr;
                            Gsq = Gs ^ Gq;

                            for (p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                                P = OutBuf.params->poff[Gp] + p;
                                for (q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                                    Q = OutBuf.params->qoff[Gq] + q;
                                    pq = OutBuf.params->rowidx[P][Q];

                                    for (r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                                        R = OutBuf.params->roff[Gr] + r;
                                        pr = InBuf->params->rowidx[P][R];

                                        for (s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                            S = OutBuf.params->soff[Gs] + s;
                                            rs = OutBuf.params->colidx[R][S];
                                            sq = InBuf->params->colidx[S][Q];

                                            OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gpr][pr][sq];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                outfile->Printf("LIBDPD: Out-of-core algorithm not yet coded for psqr sort.\n");
                dpd_error("buf4_sort", "outfile");
            }

#ifdef DPD_TIMER
            timer_off("psqr");
#endif
            break;

        case psrq:

#ifdef DPD_TIMER
            timer_on("psrq");
#endif

            /* p->p; s->q; r->r; q->s = psrq */

            if (incore) {
                for (h = 0; h < nirreps; h++) {
                    r_irrep = h ^ my_irrep;

                    for (Gp = 0; Gp < nirreps; Gp++) {
                        Gq = Gp ^ h;
                        for (Gr = 0; Gr < nirreps; Gr++) {
                            Gs = Gr ^ r_irrep;

                            Gps = Gp ^ Gs;
                            Grq = Gr ^ Gq;

                            for (p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                                P = OutBuf.params->poff[Gp] + p;
                                for (q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                                    Q = OutBuf.params->qoff[Gq] + q;
                                    pq = OutBuf.params->rowidx[P][Q];

                                    for (r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                                        R = OutBuf.params->roff[Gr] + r;
                                        rq = InBuf->params->colidx[R][Q];

                                        for (s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                            S = OutBuf.params->soff[Gs] + s;
                                            rs = OutBuf.params->colidx[R][S];
                                            ps = InBuf->params->rowidx[P][S];

                                            OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gps][ps][rq];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                outfile->Printf("LIBDPD: Out-of-core algorithm not yet coded for psrq sort.\n");
                dpd_error("buf4_sort", "outfile");
            }

#ifdef DPD_TIMER
            timer_off("psrq");
#endif
            break;

        case qprs:

#ifdef DPD_TIMER
            timer_on("qprs");
#endif

            /* q->p; p->q; r->r; s->s = qprs */

            if (incore) {
                for (h = 0; h < nirreps; h++) {
                    r_irrep = h ^ my_irrep;

                    for (pq = 0; pq < OutBuf.params->rowtot[h]; pq++) {
                        p = OutBuf.params->roworb[h][pq][0];
                        q = OutBuf.params->roworb[h][pq][1];
                        qp = InBuf->params->rowidx[q][p];

                        for (rs = 0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
                            r = OutBuf.params->colorb[r_irrep][rs][0];
                            s = OutBuf.params->colorb[r_irrep][rs][1];

                            col = InBuf->params->colidx[r][s];

                            OutBuf.matrix[h][pq][rs] = InBuf->matrix[h][qp][col];
                        }
                    }
                }
            } else {
                buf4_sort_ooc_pairs(InBuf, &OutBuf, SortPair::QP, SortPair::RS, false, 1.0,
                                    "buf4_sort (qprs): not enough memory for one row");
            }
#ifdef DPD_TIMER
            timer_off("qprs");
#endif
            break;

        case qpsr:

#ifdef DPD_TIMER
            timer_on("qpsr");
#endif

            /* q->p; p->q; s->r; r->s = qpsr */

            if (incore) {
                for (h = 0; h < nirreps; h++) {
                    r_irrep = h ^ my_irrep;

                    for (pq = 0; pq < OutBuf.params->rowtot[h]; pq++) {
                        p = OutBuf.params->roworb[h][pq][0];
                        q = OutBuf.params->roworb[h][pq][1];
                        qp = InBuf->params->rowidx[q][p];

                        for (rs = 0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
                            r = OutBuf.params->colorb[r_irrep][rs][0];
                            s = OutBuf.params->colorb[r_irrep][rs][1];
                            sr = InBuf->params->colidx[s][r];

                            OutBuf.matrix[h][pq][rs] = InBuf->matrix[h][qp][sr];
                        }
                    }
                }
            } else {
                buf4_sort_ooc_pairs(InBuf, &OutBuf, SortPair::QP, SortPair::SR, false, 1.0,
                                    "buf4_sort (qpsr): not enough memory for one row");
            }

#ifdef DPD_TIMER
            timer_off("qpsr");
#endif
            break;

        case qrps:
#ifdef DPD_TIMER
            timer_on("qrps");
#endif

            /* q->p; r->q; p->r; s->s = rpqs */

            if (incore) {
                for (h = 0; h < nirreps; h++) {
                    r_irrep = h ^ my_irrep;

                    for (Gp = 0; Gp < nirreps; Gp++) {
                        Gq = Gp ^ h;
                        for (Gr = 0; Gr < nirreps; Gr++) {
                            Gs = Gr ^ r_irrep;

                            Grp = Gr ^ Gp;
                            Gqs = Gq ^ Gs;

                            for (p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                                P = OutBuf.params->poff[Gp] + p;
                                for (q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                                    Q = OutBuf.params->qoff[Gq] + q;
                                    pq = OutBuf.params->rowidx[P][Q];

                                    for (r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                                        R = OutBuf.params->roff[Gr] + r;
                                        rp = InBuf->params->rowidx[R][P];

                                        for (s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                            S = OutBuf.params->soff[Gs] + s;
                                            rs = OutBuf.params->colidx[R][S];
                                            qs = InBuf->params->colidx[Q][S];

                                            OutBuf.matrix[h][pq][rs] = InBuf->matrix[Grp][rp][qs];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                outfile->Printf("LIBDPD: Out-of-core algorithm not yet coded for qrps sort.\n");
                dpd_error("buf4_sort", "outfile");
            }

#ifdef DPD_TIMER
            timer_off("qrps");
#endif
            break;

        case qrsp:

#ifdef DPD_TIMER
            timer_on("qrsp");
#endif

            /* q->p; r->q; s->r; p->s = spqr */

            if (incore) {
                for (h = 0; h < nirreps; h++) {
                    r_irrep = h ^ my_irrep;

                    for (Gp = 0; Gp < nirreps; Gp++) {
                        Gq = Gp ^ h;
                        for (Gr = 0; Gr < nirreps; Gr++) {
                            Gs = Gr ^ r_irrep;

                            Gsp = Gs ^ Gp;
                            Gqr = Gq ^ Gr;

                            for (p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                                P = OutBuf.params->poff[Gp] + p;
                                for (q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                                    Q = OutBuf.params->qoff[Gq] + q;
                                    pq = OutBuf.params->rowidx[P][Q];

                                    for (r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                                        R = OutBuf.params->roff[Gr] + r;
                                        qr = InBuf->params->colidx[Q][R];

                                        for (s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                            S = OutBuf.params->soff[Gs] + s;
                                            rs = OutBuf.params->colidx[R][S];
                                            sp = InBuf->params->rowidx[S][P];

                                            OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gsp][sp][qr];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                outfile->Printf("LIBDPD: Out-of-core algorithm not yet coded for qrsp sort.\n");
                dpd_error("buf4_sort", "outfile");
            }

#ifdef DPD_TIMER
            timer_off("qrsp");
#endif
            break;

        case qspr:
#ifdef DPD_TIMER
            timer_on("qspr");
#endif

            /* q->p; s->q; p->r; r->s = rpsq */

            if (incore) {
                for (h = 0; h < nirreps; h++) {
                    r_irrep = h ^ my_irrep;

                    for (Gp = 0; Gp < nirreps; Gp++) {
                        Gq = Gp ^ h;
                        for (Gr = 0; Gr < nirreps; Gr++) {
                            Gs = Gr ^ r_irrep;

                            Grp = Gr ^ Gp;
                            Gsq = Gs ^ Gq;

                            for (p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                                P = OutBuf.params->poff[Gp] + p;
                                for (q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                                    Q = OutBuf.params->qoff[Gq] + q;
                                    pq = OutBuf.params->rowidx[P][Q];

                                    for (r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                                        R = OutBuf.params->roff[Gr] + r;
                                        rp = InBuf->params->rowidx[R][P];

                                        for (s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                            S = OutBuf.params->soff[Gs] + s;
                                            rs = OutBuf.params->colidx[R][S];
                                            sq = InBuf->params->colidx[S][Q];

                                            OutBuf.matrix[h][pq][rs] = InBuf->matrix[Grp][rp][sq];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                outfile->Printf("LIBDPD: Out-of-core algorithm not yet coded for qspr sort.\n");
                dpd_error("buf4_sort", "outfile");
            }

#ifdef DPD_TIMER
            timer_off("qspr");
#endif
            break;

        case qsrp:

#ifdef DPD_TIMER
            timer_on("qsrp");
#endif

            /* q->p; s->q; r->r; p->s = sprq */
            if (incore) {
                for (h = 0; h < nirreps; h++) {
                    r_irrep = h ^ my_irrep;

                    for (Gp = 0; Gp < nirreps; Gp++) {
                        Gq = Gp ^ h;
                        for (Gr = 0; Gr < nirreps; Gr++) {
                            Gs = Gr ^ r_irrep;

                            Gsp = Gs ^ Gp;
                            Grq = Gr ^ Gq;

                            for (p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                                P = OutBuf.params->poff[Gp] + p;
                                for (q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                                    Q = OutBuf.params->qoff[Gq] + q;
                                    pq = OutBuf.params->rowidx[P][Q];

                                    for (r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                                        R = OutBuf.params->roff[Gr] + r;
                                        rq = InBuf->params->colidx[R][Q];

                                        for (s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                            S = OutBuf.params->soff[Gs] + s;
                                            rs = OutBuf.params->colidx[R][S];
                                            sp = InBuf->params->rowidx[S][P];

                                            OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gsp][sp][rq];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                outfile->Printf("LIBDPD: Out-of-core algorithm not yet coded for qsrp sort.\n");
                dpd_error("buf4_sort", "outfile");
            }

#ifdef DPD_TIMER
            timer_off("qsrp");
#endif
            break;

        case rqps:

#ifdef DPD_TIMER
            timer_on("rqps");
#endif

            /* r->p; q->q; p->r; s->s = rqps */

            if (incore) {
                for (h = 0; h < nirreps; h++) {
                    r_irrep = h ^ my_irrep;

                    for (Gp = 0; Gp < nirreps; Gp++) {
                        Gq = Gp ^ h;
                        for (Gr = 0; Gr < nirreps; Gr++) {
                            Gs = Gr ^ r_irrep;

                            Grq = Gr ^ Gq;
                            Gps = Gp ^ Gs;

                            for (p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                                P = OutBuf.params->poff[Gp] + p;
                                for (q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                                    Q = OutBuf.params->qoff[Gq] + q;
                                    pq = OutBuf.params->rowidx[P][Q];

                                    for (r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                                        R = OutBuf.params->roff[Gr] + r;
                                        rq = InBuf->params->rowidx[R][Q];

                                        for (s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                            S = OutBuf.params->soff[Gs] + s;
                                            rs = OutBuf.params->colidx[R][S];
                                            ps = InBuf->params->colidx[P][S];

                                            OutBuf.matrix[h][pq][rs] = InBuf->matrix[Grq][rq][ps];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                outfile->Printf("LIBDPD: Out-of-core algorithm not yet coded for rqps sort.\n");
                dpd_error("buf4_sort", "outfile");
            }

#ifdef DPD_TIMER
            timer_off("rqps");
#endif
            break;

        case rqsp:

#ifdef DPD_TIMER
            timer_on("rqsp");
#endif

            /* r->p; q->q; s->r; p->s = sqpr */

            if (incore) {
                for (h = 0; h < nirreps; h++) {
                    r_irrep = h ^ my_irrep;

                    for (Gp = 0; Gp < nirreps; Gp++) {
                        Gq = Gp ^ h;
                        for (Gr = 0; Gr < nirreps; Gr++) {
                            Gs = Gr ^ r_irrep;

                            Gsq = Gs ^ Gq;
                            Gpr = Gp ^ Gr;

                            for (p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                                P = OutBuf.params->poff[Gp] + p;
                                for (q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                                    Q = OutBuf.params->qoff[Gq] + q;
                                    pq = OutBuf.params->rowidx[P][Q];

                                    for (r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                                        R = OutBuf.params->roff[Gr] + r;
                                        pr = InBuf->params->colidx[P][R];

                                        for (s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                            S = OutBuf.params->soff[Gs] + s;
                                            rs = OutBuf.params->colidx[R][S];
                                            sq = InBuf->params->rowidx[S][Q];

                                            OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gsq][sq][pr];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                outfile->Printf("LIBDPD: Out-of-core algorithm not yet coded for rqsp sort.\n");
                dpd_error("buf4_sort", "outfile");
            }

#ifdef DPD_TIMER
            timer_off("rqsp");
#endif
            break;

        case rpqs:

#ifdef DPD_TIMER
            timer_on("rpqs");
#endif

            /* r->p; p->q; q->r; s->s = qrps */

            if (incore) {
                for (h = 0; h < nirreps; h++) {
                    r_irrep = h ^ my_irrep;

                    for (Gp = 0; Gp < nirreps; Gp++) {
                        Gq = Gp ^ h;
                        for (Gr = 0; Gr < nirreps; Gr++) {
                            Gs = Gr ^ r_irrep;

                            Gqr = Gq ^ Gr;
                            Gps = Gp ^ Gs;

                            for (p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                                P = OutBuf.params->poff[Gp] + p;
                                for (q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                                    Q = OutBuf.params->qoff[Gq] + q;
                                    pq = OutBuf.params->rowidx[P][Q];

                                    for (r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                                        R = OutBuf.params->roff[Gr] + r;
                                        qr = InBuf->params->rowidx[Q][R];

                                        for (s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                            S = OutBuf.params->soff[Gs] + s;
                                            rs = OutBuf.params->colidx[R][S];
                                            ps = InBuf->params->colidx[P][S];

                                            OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gqr][qr][ps];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                outfile->Printf("LIBDPD: Out-of-core algorithm not yet coded for rpqs sort.\n");
                dpd_error("buf4_sort", "outfile");
            }

#ifdef DPD_TIMER
            timer_off("rpqs");
#endif
            break;

        case rpsq:

#ifdef DPD_TIMER
            timer_on("rpsq");
#endif

            /* r->p; p->q; s->r; q->s = qspr */

            if (incore) {
                for (h = 0; h < nirreps; h++) {
                    r_irrep = h ^ my_irrep;

                    for (Gp = 0; Gp < nirreps; Gp++) {
                        Gq = Gp ^ h;
                        for (Gr = 0; Gr < nirreps; Gr++) {
                            Gs = Gr ^ r_irrep;

                            Gqs = Gq ^ Gs;
                            Gpr = Gp ^ Gr;

                            for (p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                                P = OutBuf.params->poff[Gp] + p;
                                for (q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                                    Q = OutBuf.params->qoff[Gq] + q;
                                    pq = OutBuf.params->rowidx[P][Q];

                                    for (r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                                        R = OutBuf.params->roff[Gr] + r;
                                        pr = InBuf->params->colidx[P][R];

                                        for (s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                            S = OutBuf.params->soff[Gs] + s;
                                            rs = OutBuf.params->colidx[R][S];
                                            qs = InBuf->params->rowidx[Q][S];

                                            OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gqs][qs][pr];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                outfile->Printf("LIBDPD: Out-of-core algorithm not yet coded for rpsq sort.\n");
                dpd_error("buf4_sort", "outfile");
            }

#ifdef DPD_TIMER
            timer_off("rpsq");
#endif
            break;

        case rsqp:

#ifdef DPD_TIMER
            timer_on("rsqp");
#endif

            /* r->p; s->q; q->r; p->s = srpq */

            if (incore) {
                for (h = 0; h < nirreps; h++) {
                    r_irrep = h ^ my_irrep;

                    for (pq = 0; pq < OutBuf.params->rowtot[h]; pq++) {
                        p = OutBuf.params->roworb[h][pq][0];
                        q = OutBuf.params->roworb[h][pq][1];

                        col = InBuf->params->colidx[p][q];

                        for (rs = 0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
                            r = OutBuf.params->colorb[r_irrep][rs][0];
                            s = OutBuf.params->colorb[r_irrep][rs][1];

                            row = InBuf->params->rowidx[s][r];

                            OutBuf.matrix[h][pq][rs] = InBuf->matrix[h][row][col];
                        }
                    }
                }
            } else {
                outfile->Printf("LIBDPD: Out-of-core algorithm not yet coded for rsqp sort.\n");
                dpd_error("buf4_sort", "outfile");
            }

#ifdef DPD_TIMER
            timer_off("rsqp");
#endif
            break;

        case rspq:

#ifdef DPD_TIMER
            timer_on("rspq");
#endif

            /* r->p; s->q; p->r; q->s = rspq */

            if (incore) {
                for (h = 0; h < nirreps; h++) {
                    r_irrep = h ^ my_irrep;

                    for (pq = 0; pq < OutBuf.params->rowtot[h]; pq++) {
                        p = OutBuf.params->roworb[h][pq][0];
                        q = OutBuf.params->roworb[h][pq][1];

                        col = InBuf->params->colidx[p][q];

                        for (rs = 0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
                            r = OutBuf.params->colorb[r_irrep][rs][0];
                            s = OutBuf.params->colorb[r_irrep][rs][1];

                            row = InBuf->params->rowidx[r][s];

                            OutBuf.matrix[h][pq][rs] = InBuf->matrix[r_irrep][row][col];
                        }
                    }
                }
            } else {
                buf4_sort_ooc_pairs(InBuf, &OutBuf, SortPair::RS, SortPair::PQ, false, 1.0,
                                    "buf4_sort (rspq): not enough memory for one row");
            }

#ifdef DPD_TIMER
            timer_off("rspq");
#endif
            break;

        case sqrp:

#ifdef DPD_TIMER
            timer_on("sqrp");
#endif

            /* s->p; q->q; r->r; p->s = sqrp */

            if (incore) {
                for (h = 0; h < nirreps; h++) {
                    r_irrep = h ^ my_irrep;

                    for (Gp = 0; Gp < nirreps; Gp++) {
                        Gq = Gp ^ h;
                        for (Gr = 0; Gr < nirreps; Gr++) {
                            Gs = Gr ^ r_irrep;

                            Gsq = Gs ^ Gq;
                            Grp = Gr ^ Gp;

                            for (p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                                P = OutBuf.params->poff[Gp] + p;
                                for (q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                                    Q = OutBuf.params->qoff[Gq] + q;
                                    pq = OutBuf.params->rowidx[P][Q];

                                    for (r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                                        R = OutBuf.params->roff[Gr] + r;
                                        rp = InBuf->params->colidx[R][P];

                                        for (s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                            S = OutBuf.params->soff[Gs] + s;
                                            rs = OutBuf.params->colidx[R][S];
                                            sq = InBuf->params->rowidx[S][Q];

                                            OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gsq][sq][rp];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                outfile->Printf("LIBDPD: Out-of-core algorithm not yet coded for sqrp sort.\n");
                dpd_error("buf4_sort", "outfile");
            }

#ifdef DPD_TIMER
            timer_off("sqrp");
#endif
            break;

        case sqpr:
            outfile->Printf("\nDPD sort error: sqpr index ordering not yet coded.\n");
            dpd_error("buf_sort", "outfile");
            break;

        case srqp:

#ifdef DPD_TIMER
            timer_on("srqp");
#endif

            /* s->p; r->q; q->r; p->s = srqp */

            if (incore) {
                for (h = 0; h < nirreps; h++) {
                    r_irrep = h ^ my_irrep;

                    for (pq = 0; pq < OutBuf.params->rowtot[h]; pq++) {
                        p = OutBuf.params->roworb[h][pq][0];
                        q = OutBuf.params->roworb[h][pq][1];

                        col = InBuf->params->colidx[q][p];

                        for (rs = 0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
                            r = OutBuf.params->colorb[r_irrep][rs][0];
                            s = OutBuf.params->colorb[r_irrep][rs][1];

                            row = InBuf->params->rowidx[s][r];

                            OutBuf.matrix[h][pq][rs] = InBuf->matrix[r_irrep][row][col];
                        }
                    }
                }
            } else {
                outfile->Printf("LIBDPD: Out-of-core algorithm not yet coded for srqp sort.\n");
                dpd_error("buf4_sort", "outfile");
            }

#ifdef DPD_TIMER
            timer_off("srqp");
#endif
            break;

        case srpq:
#ifdef DPD_TIMER
            timer_on("srpq");
#endif

            /* s->p; r->q; p->r; q->s = rsqp */

            if (incore) {
                for (h = 0; h < nirreps; h++) {
                    r_irrep = h ^ my_irrep;

                    for (pq = 0; pq < OutBuf.params->rowtot[h]; pq++) {
                        p = OutBuf.params->roworb[h][pq][0];
                        q = OutBuf.params->roworb[h][pq][1];

                        col = InBuf->params->colidx[q][p];

                        for (rs = 0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
                            r = OutBuf.params->colorb[r_irrep][rs][0];
                            s = OutBuf.params->colorb[r_irrep][rs][1];

                            row = InBuf->params->rowidx[r][s];

                            OutBuf.matrix[h][pq][rs] = InBuf->matrix[h][row][col];
                        }
                    }
                }
            } else {
                outfile->Printf("LIBDPD: Out-of-core algorithm not yet coded for srpq sort.\n");
                dpd_error("buf4_sort", "outfile");
            }

#ifdef DPD_TIMER
            timer_off("srpq");
#endif
            break;

        case spqr:

#ifdef DPD_TIMER
            timer_on("spqr");
#endif

            /* s->p; p->q; q->r; r->s = qrsp */

            if (incore) {
                for (h = 0; h < nirreps; h++) {
                    r_irrep = h ^ my_irrep;

                    for (Gp = 0; Gp < nirreps; Gp++) {
                        Gq = Gp ^ h;
                        for (Gr = 0; Gr < nirreps; Gr++) {
                            Gs = Gr ^ r_irrep;

                            Gqr = Gq ^ Gr;
                            Gsp = Gs ^ Gp;

                            for (p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                                P = OutBuf.params->poff[Gp] + p;
                                for (q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                                    Q = OutBuf.params->qoff[Gq] + q;
                                    pq = OutBuf.params->rowidx[P][Q];

                                    for (r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                                        R = OutBuf.params->roff[Gr] + r;
                                        qr = InBuf->params->rowidx[Q][R];

                                        for (s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                            S = OutBuf.params->soff[Gs] + s;
                                            rs = OutBuf.params->colidx[R][S];
                                            sp = InBuf->params->colidx[S][P];

                                            OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gqr][qr][sp];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                outfile->Printf("LIBDPD: Out-of-core algorithm not yet coded for spqr sort.\n");
                dpd_error("buf4_sort", "outfile");
            }

#ifdef DPD_TIMER
            timer_off("spqr");
#endif
            break;

        case sprq:

#ifdef DPD_TIMER
            timer_on("sprq");
#endif

            /* s->p; p->q; r->r; q->s = qsrp */

            if (incore) {
                for (h = 0; h < nirreps; h++) {
                    r_irrep = h ^ my_irrep;

                    for (Gp = 0; Gp < nirreps; Gp++) {
                        Gq = Gp ^ h;
                        for (Gr = 0; Gr < nirreps; Gr++) {
                            Gs = Gr ^ r_irrep;

                            Gqs = Gq ^ Gs;
                            Grp = Gr ^ Gp;

                            for (p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                                P = OutBuf.params->poff[Gp] + p;
                                for (q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                                    Q = OutBuf.params->qoff[Gq] + q;
                                    pq = OutBuf.params->rowidx[P][Q];

                                    for (r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                                        R = OutBuf.params->roff[Gr] + r;
                                        rp = InBuf->params->colidx[R][P];

                                        for (s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                            S = OutBuf.params->soff[Gs] + s;
                                            rs = OutBuf.params->colidx[R][S];
                                            qs = InBuf->params->rowidx[Q][S];

                                            OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gqs][qs][rp];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                outfile->Printf("LIBDPD: Out-of-core algorithm not yet coded for sprq sort.\n");
                dpd_error("buf4_sort", "outfile");
            }

#ifdef DPD_TIMER
            timer_off("sprq");
#endif
            break;
    }

    if (incore) {
        for (h = 0; h < nirreps; h++) {
            buf4_mat_irrep_wrt(&OutBuf, h);
            buf4_mat_irrep_close(&OutBuf, h);
            buf4_mat_irrep_close(InBuf, h);
        }
    }

    buf4_close(&OutBuf);

#ifdef DPD_TIMER
    timer_off("buf4_sort");
#endif

    return 0;
}

void DPD::note_ooc_sort(const char *function, enum indices index) {
    static const char *names[] = {"pqrs", "pqsr", "prqs", "prsq", "psqr", "psrq", "qprs", "qpsr",
                                  "qrps", "qrsp", "qspr", "qsrp", "rqps", "rqsp", "rpqs", "rpsq",
                                  "rsqp", "rspq", "sqrp", "sqpr", "srqp", "srpq", "spqr", "sprq"};
    if (ooc_sort_noted_) return;
    ooc_sort_noted_ = true;
    outfile->Printf("\n  DPD: not enough memory for an in-core %s (%s); using out-of-core sorts.\n\n", function,
                    names[index]);
}

namespace {
inline int sort_pair_index(int **idx, DPD::SortPair which, int p, int q, int r, int s) {
    switch (which) {
        case DPD::SortPair::PQ:
            return idx[p][q];
        case DPD::SortPair::QP:
            return idx[q][p];
        case DPD::SortPair::RS:
            return idx[r][s];
        default:
            return idx[s][r];
    }
}
}  // namespace

int DPD::buf4_sort_ooc_pairs(dpdbuf4 *InBuf, dpdbuf4 *OutBuf, SortPair row_from, SortPair col_from, bool accumulate,
                             double alpha, const char *caller) {
    const int nirreps = OutBuf->params->nirreps;
    const int my_irrep = OutBuf->file.my_irrep;
    const bool in_block_is_pq = (row_from == SortPair::PQ || row_from == SortPair::QP);

    for (int Gpq = 0; Gpq < nirreps; Gpq++) {
        const int Grs = Gpq ^ my_irrep;
        const int out_rows = OutBuf->params->rowtot[Gpq];
        const int out_cols = OutBuf->params->coltot[Grs];
        /* an empty or zero-width output block has no elements (and a null block matrix) */
        if (!out_rows || !out_cols) continue;

        const int Gin = in_block_is_pq ? Gpq : Grs;
        const int in_rows = InBuf->params->rowtot[Gin];
        const int in_cols = InBuf->params->coltot[Gin ^ InBuf->file.my_irrep];

        /* half the free core for a bucket of output rows, the rest for a bucket of input rows */
        const int out_rows_per_bucket = dpd_rows_per_bucket_or_die(dpd_memfree() / 2, out_cols, out_rows, caller);
        buf4_mat_irrep_init_block(OutBuf, Gpq, out_rows_per_bucket);
        const int in_rows_per_bucket = dpd_rows_per_bucket_or_die(dpd_memfree(), in_cols, in_rows, caller);
        buf4_mat_irrep_init_block(InBuf, Gin, in_rows_per_bucket);

        for (int out_start = 0; out_start < out_rows; out_start += out_rows_per_bucket) {
            const int nout = std::min(out_rows_per_bucket, out_rows - out_start);
            if (accumulate)
                buf4_mat_irrep_rd_block(OutBuf, Gpq, out_start, nout);
            else
                std::memset(OutBuf->matrix[Gpq][0], 0, sizeof(double) * nout * static_cast<size_t>(out_cols));

            for (int in_start = 0; in_start < in_rows; in_start += in_rows_per_bucket) {
                const int nin = std::min(in_rows_per_bucket, in_rows - in_start);
                buf4_mat_irrep_rd_block(InBuf, Gin, in_start, nin);

                for (int pq = 0; pq < nout; pq++) {
                    const int p = OutBuf->params->roworb[Gpq][out_start + pq][0];
                    const int q = OutBuf->params->roworb[Gpq][out_start + pq][1];
                    for (int rs = 0; rs < out_cols; rs++) {
                        const int r = OutBuf->params->colorb[Grs][rs][0];
                        const int s = OutBuf->params->colorb[Grs][rs][1];
                        /* is this element's input row in the current input bucket? */
                        const int row = sort_pair_index(InBuf->params->rowidx, row_from, p, q, r, s) - in_start;
                        if (row < 0 || row >= nin) continue;
                        const int col = sort_pair_index(InBuf->params->colidx, col_from, p, q, r, s);
                        if (col < 0) continue;
                        const double value = InBuf->matrix[Gin][row][col];
                        if (accumulate)
                            OutBuf->matrix[Gpq][pq][rs] += alpha * value;
                        else
                            OutBuf->matrix[Gpq][pq][rs] = value;
                    }
                }
            }
            buf4_mat_irrep_wrt_block(OutBuf, Gpq, out_start, nout);
        }

        buf4_mat_irrep_close_block(InBuf, Gin, in_rows_per_bucket);
        buf4_mat_irrep_close_block(OutBuf, Gpq, out_rows_per_bucket);
    }

    return 0;
}

int DPD::buf4_sort(dpdbuf4 *InBuf, int outfilenum, enum indices index, string pq, string rs, const std::string& label) {
    return buf4_sort(InBuf, outfilenum, index, pairnum(pq), pairnum(rs), label);
}

}  // namespace psi
