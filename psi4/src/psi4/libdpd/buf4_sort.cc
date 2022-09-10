/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

#include <cstdio>
#include <cstdlib>
#include <cmath>

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

int DPD::buf4_sort(dpdbuf4* InBuf, const int outfilenum, const enum indices index, const int pqnum, const int rsnum,
                   const std::string& label) {
    dpdbuf4 OutBuf;
    int in_nbuckets, in_rows_left, in_row_start, m;
    int rows_per_bucket, nbuckets, rows_left;

    const int nirreps = InBuf->params->nirreps;
    const int my_irrep = InBuf->file.my_irrep;

#ifdef DPD_TIMER
    timer_on("buf4_sort");
#endif

    buf4_init(&OutBuf, outfilenum, my_irrep, pqnum, rsnum, pqnum, rsnum, 0, label);

    /* select in-core vs. out-of-core algorithms */
    int incore = 1;
    {
        long int core_total = 0;
        for (int h = 0; h < nirreps; h++) {
            const long int coltot = InBuf->params->coltot[h ^ my_irrep];
            long int maxrows;
            if (coltot) {
                maxrows = DPD_BIGNUM / coltot;
                if (maxrows < 1) {
                    outfile->Printf("\nLIBDPD Error: too many rows in buf4_sort_axpy.\n");
                    dpd_error("buf4_sort_axpy", "outfile");
                }
            } else
                maxrows = DPD_BIGNUM;
            long int rowtot = InBuf->params->rowtot[h];
            for (; rowtot > maxrows; rowtot -= maxrows) {
                if (core_total > (core_total + 2 * maxrows * coltot))
                    incore = 0;
                else
                    core_total += 2 * maxrows * coltot;
            }
            if (core_total > (core_total + 2 * rowtot * coltot)) incore = 0;
            core_total += 2 * rowtot * coltot;
        }
        if (core_total > dpd_memfree()) incore = 0;
    }

#ifdef DPD_DEBUG
    if (incore == 0) {
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
            incore = 0;
            break;
        case (prqs):
            incore = 0;
            break;
        case (prsq):
            incore = 0;
            break;
        case (qprs):
            incore = 0;
            break;
        case (qpsr):
            incore = 0;
            break;
        case (sqpr):
            incore = 0;
            break;
        case (rspq):
            incore = 0;
            break;
    }
#endif

    /* Init input and output buffers and read in all blocks of the input */
    if (incore) {
        for (int h = 0; h < nirreps; h++) {
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
                for (int h = 0; h < nirreps; h++) {
                    const int r_irrep = h ^ my_irrep;

                    for (int pq = 0; pq < OutBuf.params->rowtot[h]; pq++) {
                        const int p = OutBuf.params->roworb[h][pq][0];
                        const int q = OutBuf.params->roworb[h][pq][1];

                        const int row = InBuf->params->rowidx[p][q];

                        for (int rs = 0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
                            const int r = OutBuf.params->colorb[r_irrep][rs][0];
                            const int s = OutBuf.params->colorb[r_irrep][rs][1];

                            const int sr = InBuf->params->colidx[s][r];

                            OutBuf.matrix[h][pq][rs] = InBuf->matrix[h][row][sr];
                        }
                    }
                }
            } else { /* out-of-core pqsr -> pqrs */

                for (int Gpq = 0; Gpq < nirreps; Gpq++) {
                    const int Grs = Gpq ^ my_irrep;
                    rows_per_bucket = dpd_memfree() / 2 / InBuf->params->coltot[Grs];

                    if (rows_per_bucket > InBuf->params->rowtot[Gpq]) rows_per_bucket = InBuf->params->rowtot[Gpq];
                    if (!rows_per_bucket) dpd_error("buf4_sort_pqsr: Not enough memory for one row!", "outfile");

                    nbuckets = (int)ceil(((double)InBuf->params->rowtot[Gpq]) / ((double)rows_per_bucket));
                    if (nbuckets == 1)
                        rows_left = rows_per_bucket;
                    else
                        rows_left = InBuf->params->rowtot[Gpq] % rows_per_bucket;

                    buf4_mat_irrep_init_block(InBuf, Gpq, rows_per_bucket);
                    buf4_mat_irrep_init_block(&OutBuf, Gpq, rows_per_bucket);

                    int n;
                    for (n = 0; n < (rows_left ? nbuckets - 1 : nbuckets); n++) {
                        buf4_mat_irrep_rd_block(InBuf, Gpq, n * rows_per_bucket, rows_per_bucket);

                        for (int pq = 0; pq < rows_per_bucket; pq++) {
                            for (int rs = 0; rs < OutBuf.params->coltot[Grs]; rs++) {
                                const int r = OutBuf.params->colorb[Grs][rs][0];
                                const int s = OutBuf.params->colorb[Grs][rs][1];

                                const int sr = InBuf->params->colidx[s][r];

                                OutBuf.matrix[Gpq][pq][rs] = InBuf->matrix[Gpq][pq][sr];
                            }
                        }

                        buf4_mat_irrep_wrt_block(&OutBuf, Gpq, n * rows_per_bucket, rows_per_bucket);
                    }
                    if (rows_left) {
                        buf4_mat_irrep_rd_block(InBuf, Gpq, n * rows_per_bucket, rows_left);

                        for (int pq = 0; pq < rows_left; pq++) {
                            for (int rs = 0; rs < OutBuf.params->coltot[Grs]; rs++) {
                                const int r = OutBuf.params->colorb[Grs][rs][0];
                                const int s = OutBuf.params->colorb[Grs][rs][1];

                                const int sr = InBuf->params->colidx[s][r];

                                OutBuf.matrix[Gpq][pq][rs] = InBuf->matrix[Gpq][pq][sr];
                            }
                        }

                        buf4_mat_irrep_wrt_block(&OutBuf, Gpq, n * rows_per_bucket, rows_left);
                    }

                    buf4_mat_irrep_close_block(InBuf, Gpq, rows_per_bucket);
                    buf4_mat_irrep_close_block(&OutBuf, Gpq, rows_per_bucket);
                }
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
                for (int h = 0; h < nirreps; h++) {
                    const int r_irrep = h ^ my_irrep;

                    for (int Gp = 0; Gp < nirreps; Gp++) {
                        const int Gq = Gp ^ h;
                        for (int Gr = 0; Gr < nirreps; Gr++) {
                            const int Gs = Gr ^ r_irrep;

                            /* Irreps on the source */
                            const int Gpr = Gp ^ Gr;
                            const int Gqs = Gq ^ Gs;

                            for (int p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                                const int P = OutBuf.params->poff[Gp] + p;
                                for (int q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                                    const int Q = OutBuf.params->qoff[Gq] + q;
                                    const int pq = OutBuf.params->rowidx[P][Q];

                                    for (int r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                                        const int R = OutBuf.params->roff[Gr] + r;
                                        const int pr = InBuf->params->rowidx[P][R];

                                        for (int s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                            const int S = OutBuf.params->soff[Gs] + s;
                                            const int rs = OutBuf.params->colidx[R][S];
                                            const int qs = InBuf->params->colidx[Q][S];

                                            OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gpr][pr][qs];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } else { /* pqrs <- prqs */

                for (int Gpq = 0; Gpq < nirreps; Gpq++) {
                    const int Grs = Gpq ^ my_irrep;

                    const int out_rows_per_bucket = [&OutBuf, Grs, Gpq] {
                        int orpb = dpd_memfree() / (2 * OutBuf.params->coltot[Grs]);
                        if (orpb > OutBuf.params->rowtot[Gpq]) orpb = OutBuf.params->rowtot[Gpq];
                        return orpb;
                    }();

                    const int out_nbuckets =
                        (int)ceil((double)OutBuf.params->rowtot[Gpq] / (double)out_rows_per_bucket);
                    const int out_rows_left =
                        (out_nbuckets == 1) ? out_rows_per_bucket : OutBuf.params->rowtot[Gpq] % out_rows_per_bucket;

                    /* allocate space for the bucket of rows */
                    buf4_mat_irrep_init_block(&OutBuf, Gpq, out_rows_per_bucket);

                    int n;
                    for (n = 0; n < (out_rows_left ? out_nbuckets - 1 : out_nbuckets); n++) {
                        const int out_row_start = n * out_rows_per_bucket;

                        for (int Grow = 0; Grow < nirreps; Grow++) { /*Grow = Gpr*/
                            const int Gcol = Grow ^ my_irrep;        /*Gcol = Gqs*/

                            /* determine how many rows of InBuf we can store in the other half of the core */
                            const int in_rows_per_bucket = [&InBuf, Gcol, Grow] {
                                int irpb = dpd_memfree() / (2 * InBuf->params->coltot[Gcol]);
                                if (irpb > InBuf->params->rowtot[Grow]) irpb = InBuf->params->rowtot[Grow];
                                return irpb;
                            }();

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

                                for (int pq = 0; pq < out_rows_per_bucket; pq++) {
                                    const int p = OutBuf.params->roworb[Gpq][pq + out_row_start][0];
                                    const int q = OutBuf.params->roworb[Gpq][pq + out_row_start][1];
                                    const int Gp = OutBuf.params->psym[p];
                                    const int Gq = Gpq ^ Gp;
                                    for (int rs = 0; rs < OutBuf.params->coltot[Grs]; rs++) {
                                        const int r = OutBuf.params->colorb[Grs][rs][0];
                                        const int s = OutBuf.params->colorb[Grs][rs][1];
                                        const int Gr = OutBuf.params->rsym[r];
                                        const int Gs = Grs ^ Gr;
                                        const int Gpr = Gp ^ Gr;

                                        if (Gpr == Grow) {
                                            const int pr = InBuf->params->rowidx[p][r] - in_row_start;
                                            /* check if the current value is in the current in_bucket or not */
                                            if (pr >= 0 && pr < in_rows_per_bucket) {
                                                const int qs = InBuf->params->colidx[q][s];
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
                                for (int pq = 0; pq < out_rows_per_bucket; pq++) {
                                    const int p = OutBuf.params->roworb[Gpq][pq + out_row_start][0];
                                    const int q = OutBuf.params->roworb[Gpq][pq + out_row_start][1];
                                    const int Gp = OutBuf.params->psym[p];
                                    const int Gq = Gpq ^ Gp;
                                    for (int rs = 0; rs < OutBuf.params->coltot[Grs]; rs++) {
                                        const int r = OutBuf.params->colorb[Grs][rs][0];
                                        const int s = OutBuf.params->colorb[Grs][rs][1];
                                        const int Gr = OutBuf.params->rsym[r];
                                        const int Gs = Grs ^ Gr;
                                        const int Gpr = Gp ^ Gr;

                                        if (Gpr == Grow) {
                                            const int pr = InBuf->params->rowidx[p][r] - in_row_start;
                                            /* check if the current value is in core or not */
                                            if (pr >= 0 && pr < in_rows_left) {
                                                const int qs = InBuf->params->colidx[q][s];
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
                        const int out_row_start = n * out_rows_per_bucket;

                        for (int Grow = 0; Grow < nirreps; Grow++) {
                            const int Gcol = Grow ^ my_irrep;

                            /* determine how many rows of InBuf we can store in the other half of the core */
                            const int in_rows_per_bucket = [&InBuf, Gcol, Grow] {
                                int irpb = dpd_memfree() / (2 * InBuf->params->coltot[Gcol]);
                                if (irpb > InBuf->params->rowtot[Grow]) irpb = InBuf->params->rowtot[Grow];
                                return irpb;
                            }();

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

                                for (int pq = 0; pq < out_rows_left; pq++) {
                                    const int p = OutBuf.params->roworb[Gpq][pq + out_row_start][0];
                                    const int q = OutBuf.params->roworb[Gpq][pq + out_row_start][1];
                                    const int Gp = OutBuf.params->psym[p];
                                    const int Gq = Gpq ^ Gp;
                                    for (int rs = 0; rs < OutBuf.params->coltot[Grs]; rs++) {
                                        const int r = OutBuf.params->colorb[Grs][rs][0];
                                        const int s = OutBuf.params->colorb[Grs][rs][1];
                                        const int Gr = OutBuf.params->rsym[r];
                                        const int Gs = Grs ^ Gr;
                                        const int Gpr = Gp ^ Gr;

                                        if (Gpr == Grow) {
                                            const int pr = InBuf->params->rowidx[p][r] - in_row_start;
                                            /* check if the current value is in the current in_bucket or not */
                                            if (pr >= 0 && pr < in_rows_per_bucket) {
                                                const int qs = InBuf->params->colidx[q][s];
                                                OutBuf.matrix[Gpq][pq][rs] = InBuf->matrix[Grow][pr][qs];
                                            }
                                        }
                                    }
                                }
                            }
                            if (in_rows_left) {
                                in_row_start = m * in_rows_per_bucket;
                                buf4_mat_irrep_rd_block(InBuf, Grow, in_row_start, in_rows_left);

                                for (int pq = 0; pq < out_rows_left; pq++) {
                                    const int p = OutBuf.params->roworb[Gpq][pq + out_row_start][0];
                                    const int q = OutBuf.params->roworb[Gpq][pq + out_row_start][1];
                                    const int Gp = OutBuf.params->psym[p];
                                    const int Gq = Gpq ^ Gp;
                                    for (int rs = 0; rs < OutBuf.params->coltot[Grs]; rs++) {
                                        const int r = OutBuf.params->colorb[Grs][rs][0];
                                        const int s = OutBuf.params->colorb[Grs][rs][1];
                                        const int Gr = OutBuf.params->rsym[r];
                                        const int Gs = Grs ^ Gr;
                                        const int Gpr = Gp ^ Gr;

                                        if (Gpr == Grow) {
                                            const int pr = InBuf->params->rowidx[p][r] - in_row_start;
                                            /* check if the current value is in core or not */
                                            if (pr >= 0 && pr < in_rows_left) {
                                                const int qs = InBuf->params->colidx[q][s];
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
                for (int h = 0; h < nirreps; h++) {
                    const int r_irrep = h ^ my_irrep;

                    for (int Gp = 0; Gp < nirreps; Gp++) {
                        const int Gq = Gp ^ h;
                        for (int Gr = 0; Gr < nirreps; Gr++) {
                            const int Gs = Gr ^ r_irrep;

                            const int Gps = Gp ^ Gs;
                            const int Gqr = Gq ^ Gr;

                            for (int p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                                const int P = OutBuf.params->poff[Gp] + p;
                                for (int q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                                    const int Q = OutBuf.params->qoff[Gq] + q;
                                    const int pq = OutBuf.params->rowidx[P][Q];

                                    for (int r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                                        const int R = OutBuf.params->roff[Gr] + r;
                                        const int qr = InBuf->params->colidx[Q][R];

                                        for (int s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                            const int S = OutBuf.params->soff[Gs] + s;
                                            const int rs = OutBuf.params->colidx[R][S];
                                            const int ps = InBuf->params->rowidx[P][S];

                                            OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gps][ps][qr];
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } else {
                for (int Gpq = 0; Gpq < nirreps; Gpq++) {
                    const int Grs = Gpq ^ my_irrep;

                    /* determine how many rows of OutBuf we can store in half of the core */
                    const int out_rows_per_bucket = [&OutBuf, Grs, Gpq] {
                        int orpb = dpd_memfree() / (2 * OutBuf.params->coltot[Grs]);
                        if (orpb > OutBuf.params->rowtot[Gpq]) orpb = OutBuf.params->rowtot[Gpq];
                        return orpb;
                    }();

                    const int out_nbuckets =
                        (int)ceil((double)OutBuf.params->rowtot[Gpq] / (double)out_rows_per_bucket);
                    const int out_rows_left =
                        (out_nbuckets == 1) ? out_rows_per_bucket : OutBuf.params->rowtot[Gpq] % out_rows_per_bucket;

                    /* allocate space for the bucket of rows */
                    buf4_mat_irrep_init_block(&OutBuf, Gpq, out_rows_per_bucket);

                    int n;
                    for (n = 0; n < (out_rows_left ? out_nbuckets - 1 : out_nbuckets); n++) {
                        const int out_row_start = n * out_rows_per_bucket;

                        for (int Grow = 0; Grow < nirreps; Grow++) {
                            const int Gcol = Grow ^ my_irrep;

                            /* determine how many rows of InBuf we can store in the other half of the core */
                            const int in_rows_per_bucket = [&InBuf, Gcol, Grow] {
                                int irpb = dpd_memfree() / (2 * InBuf->params->coltot[Gcol]);
                                if (irpb > InBuf->params->rowtot[Grow]) irpb = InBuf->params->rowtot[Grow];
                                return irpb;
                            }();

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

                                for (int pq = 0; pq < out_rows_per_bucket; pq++) {
                                    const int p = OutBuf.params->roworb[Gpq][pq + out_row_start][0];
                                    const int q = OutBuf.params->roworb[Gpq][pq + out_row_start][1];
                                    const int Gp = OutBuf.params->psym[p];
                                    const int Gq = Gpq ^ Gp;
                                    for (int rs = 0; rs < OutBuf.params->coltot[Grs]; rs++) {
                                        const int r = OutBuf.params->colorb[Grs][rs][0];
                                        const int s = OutBuf.params->colorb[Grs][rs][1];
                                        const int Gr = OutBuf.params->rsym[r];
                                        const int Gs = Grs ^ Gr;
                                        const int Gps = Gp ^ Gs;

                                        if (Gps == Grow) {
                                            const int ps = InBuf->params->rowidx[p][s] - in_row_start;
                                            /* check if the current value is in the current in_bucket or not */
                                            if (ps >= 0 && ps < in_rows_per_bucket) {
                                                const int qr = InBuf->params->colidx[q][r];
                                                OutBuf.matrix[Gpq][pq][rs] = InBuf->matrix[Grow][ps][qr];
                                            }
                                        }
                                    }
                                }
                            }
                            if (in_rows_left) {
                                in_row_start = m * in_rows_per_bucket;
                                buf4_mat_irrep_rd_block(InBuf, Grow, in_row_start, in_rows_left);

                                for (int pq = 0; pq < out_rows_per_bucket; pq++) {
                                    const int p = OutBuf.params->roworb[Gpq][pq + out_row_start][0];
                                    const int q = OutBuf.params->roworb[Gpq][pq + out_row_start][1];
                                    const int Gp = OutBuf.params->psym[p];
                                    const int Gq = Gpq ^ Gp;
                                    for (int rs = 0; rs < OutBuf.params->coltot[Grs]; rs++) {
                                        const int r = OutBuf.params->colorb[Grs][rs][0];
                                        const int s = OutBuf.params->colorb[Grs][rs][1];
                                        const int Gr = OutBuf.params->rsym[r];
                                        const int Gs = Grs ^ Gr;
                                        const int Gps = Gp ^ Gs;

                                        if (Gps == Grow) {
                                            const int ps = InBuf->params->rowidx[p][s] - in_row_start;
                                            /* check if the current value is in core or not */
                                            if (ps >= 0 && ps < in_rows_left) {
                                                const int qr = InBuf->params->colidx[q][r];
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
                        const int out_row_start = n * out_rows_per_bucket;

                        for (int Grow = 0; Grow < nirreps; Grow++) {
                            const int Gcol = Grow ^ my_irrep;

                            /* determine how many rows of InBuf we can store in the other half of the core */
                            const int in_rows_per_bucket = [&InBuf, Gcol, Grow] {
                                int irpb = dpd_memfree() / (2 * InBuf->params->coltot[Gcol]);
                                if (irpb > InBuf->params->rowtot[Grow]) irpb = InBuf->params->rowtot[Grow];
                                return irpb;
                            }();

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

                                for (int pq = 0; pq < out_rows_left; pq++) {
                                    const int p = OutBuf.params->roworb[Gpq][pq + out_row_start][0];
                                    const int q = OutBuf.params->roworb[Gpq][pq + out_row_start][1];
                                    const int Gp = OutBuf.params->psym[p];
                                    const int Gq = Gpq ^ Gp;
                                    for (int rs = 0; rs < OutBuf.params->coltot[Grs]; rs++) {
                                        const int r = OutBuf.params->colorb[Grs][rs][0];
                                        const int s = OutBuf.params->colorb[Grs][rs][1];
                                        const int Gr = OutBuf.params->rsym[r];
                                        const int Gs = Grs ^ Gr;
                                        const int Gps = Gp ^ Gs;

                                        if (Gps == Grow) {
                                            const int ps = InBuf->params->rowidx[p][s] - in_row_start;
                                            /* check if the current value is in the current in_bucket or not */
                                            if (ps >= 0 && ps < in_rows_per_bucket) {
                                                const int qr = InBuf->params->colidx[q][r];
                                                OutBuf.matrix[Gpq][pq][rs] = InBuf->matrix[Grow][ps][qr];
                                            }
                                        }
                                    }
                                }
                            }
                            if (in_rows_left) {
                                in_row_start = m * in_rows_per_bucket;
                                buf4_mat_irrep_rd_block(InBuf, Grow, in_row_start, in_rows_left);

                                for (int pq = 0; pq < out_rows_left; pq++) {
                                    const int p = OutBuf.params->roworb[Gpq][pq + out_row_start][0];
                                    const int q = OutBuf.params->roworb[Gpq][pq + out_row_start][1];
                                    const int Gp = OutBuf.params->psym[p];
                                    const int Gq = Gpq ^ Gp;
                                    for (int rs = 0; rs < OutBuf.params->coltot[Grs]; rs++) {
                                        const int r = OutBuf.params->colorb[Grs][rs][0];
                                        const int s = OutBuf.params->colorb[Grs][rs][1];
                                        const int Gr = OutBuf.params->rsym[r];
                                        const int Gs = Grs ^ Gr;
                                        const int Gps = Gp ^ Gs;

                                        if (Gps == Grow) {
                                            const int ps = InBuf->params->rowidx[p][s] - in_row_start;
                                            /* check if the current value is in core or not */
                                            if (ps >= 0 && ps < in_rows_left) {
                                                const int qr = InBuf->params->colidx[q][r];
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
                for (int h = 0; h < nirreps; h++) {
                    const int r_irrep = h ^ my_irrep;

                    for (int Gp = 0; Gp < nirreps; Gp++) {
                        const int Gq = Gp ^ h;
                        for (int Gr = 0; Gr < nirreps; Gr++) {
                            const int Gs = Gr ^ r_irrep;

                            const int Gpr = Gp ^ Gr;
                            const int Gsq = Gs ^ Gq;

                            for (int p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                                const int P = OutBuf.params->poff[Gp] + p;
                                for (int q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                                    const int Q = OutBuf.params->qoff[Gq] + q;
                                    const int pq = OutBuf.params->rowidx[P][Q];

                                    for (int r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                                        const int R = OutBuf.params->roff[Gr] + r;
                                        const int pr = InBuf->params->rowidx[P][R];

                                        for (int s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                            const int S = OutBuf.params->soff[Gs] + s;
                                            const int rs = OutBuf.params->colidx[R][S];
                                            const int sq = InBuf->params->colidx[S][Q];

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
                for (int h = 0; h < nirreps; h++) {
                    const int r_irrep = h ^ my_irrep;

                    for (int Gp = 0; Gp < nirreps; Gp++) {
                        const int Gq = Gp ^ h;
                        for (int Gr = 0; Gr < nirreps; Gr++) {
                            const int Gs = Gr ^ r_irrep;

                            const int Gps = Gp ^ Gs;
                            const int Grq = Gr ^ Gq;

                            for (int p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                                const int P = OutBuf.params->poff[Gp] + p;
                                for (int q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                                    const int Q = OutBuf.params->qoff[Gq] + q;
                                    const int pq = OutBuf.params->rowidx[P][Q];

                                    for (int r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                                        const int R = OutBuf.params->roff[Gr] + r;
                                        const int rq = InBuf->params->colidx[R][Q];

                                        for (int s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                            const int S = OutBuf.params->soff[Gs] + s;
                                            const int rs = OutBuf.params->colidx[R][S];
                                            const int ps = InBuf->params->rowidx[P][S];

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
                for (int h = 0; h < nirreps; h++) {
                    const int r_irrep = h ^ my_irrep;

                    for (int pq = 0; pq < OutBuf.params->rowtot[h]; pq++) {
                        const int p = OutBuf.params->roworb[h][pq][0];
                        const int q = OutBuf.params->roworb[h][pq][1];
                        const int qp = InBuf->params->rowidx[q][p];

                        for (int rs = 0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
                            const int r = OutBuf.params->colorb[r_irrep][rs][0];
                            const int s = OutBuf.params->colorb[r_irrep][rs][1];

                            const int col = InBuf->params->colidx[r][s];

                            OutBuf.matrix[h][pq][rs] = InBuf->matrix[h][qp][col];
                        }
                    }
                }
            } else {
                for (int Gpq = 0; Gpq < nirreps; Gpq++) {
                    const int Grs = Gpq ^ my_irrep;

                    /* determine how many rows of OutBuf/InBuf we can store in half the core */
                    rows_per_bucket = dpd_memfree() / (2 * OutBuf.params->coltot[Grs]);
                    if (rows_per_bucket > OutBuf.params->rowtot[Gpq]) rows_per_bucket = OutBuf.params->rowtot[Gpq];
                    nbuckets = (int)ceil((double)OutBuf.params->rowtot[Gpq] / (double)rows_per_bucket);
                    if (nbuckets == 1)
                        rows_left = rows_per_bucket;
                    else
                        rows_left = OutBuf.params->rowtot[Gpq] % rows_per_bucket;

                    /* allocate space for the bucket of rows */
                    buf4_mat_irrep_init_block(&OutBuf, Gpq, rows_per_bucket);
                    buf4_mat_irrep_init_block(InBuf, Gpq, rows_per_bucket);

                    int out_row_start, n;
                    for (n = 0; n < (rows_left ? nbuckets - 1 : nbuckets); n++) {
                        out_row_start = n * rows_per_bucket;

                        for (m = 0; m < (rows_left ? nbuckets - 1 : nbuckets); m++) {
                            in_row_start = m * rows_per_bucket;
                            buf4_mat_irrep_rd_block(InBuf, Gpq, in_row_start, rows_per_bucket);
                            for (int pq = 0; pq < rows_per_bucket; pq++) {
                                /* check to see if this row is contained in the current input-bucket */
                                const int p = OutBuf.params->roworb[Gpq][pq + out_row_start][0];
                                const int q = OutBuf.params->roworb[Gpq][pq + out_row_start][1];
                                const int qp = InBuf->params->rowidx[q][p] - in_row_start;
                                if (qp >= 0 && qp < rows_per_bucket) {
                                    C_DCOPY(OutBuf.params->coltot[Grs], InBuf->matrix[Gpq][qp], 1,
                                            OutBuf.matrix[Gpq][pq], 1);
                                }
                            }
                        }

                        if (rows_left) {
                            in_row_start = m * rows_per_bucket;
                            buf4_mat_irrep_rd_block(InBuf, Gpq, in_row_start, rows_left);
                            for (int pq = 0; pq < rows_per_bucket; pq++) {
                                /* check to see if this row is contained in the current input-bucket */
                                const int p = OutBuf.params->roworb[Gpq][pq + out_row_start][0];
                                const int q = OutBuf.params->roworb[Gpq][pq + out_row_start][1];
                                const int qp = InBuf->params->rowidx[q][p] - in_row_start;
                                if (qp >= 0 && qp < rows_left) {
                                    C_DCOPY(OutBuf.params->coltot[Grs], InBuf->matrix[Gpq][qp], 1,
                                            OutBuf.matrix[Gpq][pq], 1);
                                }
                            }
                        }

                        buf4_mat_irrep_wrt_block(&OutBuf, Gpq, out_row_start, rows_per_bucket);

                    } /* n */
                    if (rows_left) {
                        out_row_start = n * rows_per_bucket;

                        for (m = 0; m < (rows_left ? nbuckets - 1 : nbuckets); m++) {
                            in_row_start = m * rows_per_bucket;
                            buf4_mat_irrep_rd_block(InBuf, Gpq, in_row_start, rows_per_bucket);
                            for (int pq = 0; pq < rows_left; pq++) {
                                /* check to see if this row is contained in the current input-bucket */
                                const int p = OutBuf.params->roworb[Gpq][pq + out_row_start][0];
                                const int q = OutBuf.params->roworb[Gpq][pq + out_row_start][1];
                                const int qp = InBuf->params->rowidx[q][p] - in_row_start;
                                if (qp >= 0 && qp < rows_per_bucket) {
                                    C_DCOPY(OutBuf.params->coltot[Grs], InBuf->matrix[Gpq][qp], 1,
                                            OutBuf.matrix[Gpq][pq], 1);
                                }
                            }
                        }

                        if (rows_left) {
                            in_row_start = m * rows_per_bucket;
                            buf4_mat_irrep_rd_block(InBuf, Gpq, in_row_start, rows_left);
                            for (int pq = 0; pq < rows_left; pq++) {
                                /* check to see if this row is contained in the current input-bucket */
                                const int p = OutBuf.params->roworb[Gpq][pq + out_row_start][0];
                                const int q = OutBuf.params->roworb[Gpq][pq + out_row_start][1];
                                const int qp = InBuf->params->rowidx[q][p] - in_row_start;
                                if (qp >= 0 && qp < rows_left) {
                                    C_DCOPY(OutBuf.params->coltot[Grs], InBuf->matrix[Gpq][qp], 1,
                                            OutBuf.matrix[Gpq][pq], 1);
                                }
                            }
                        }
                    }

                    buf4_mat_irrep_wrt_block(&OutBuf, Gpq, out_row_start, rows_left);
                    buf4_mat_irrep_close_block(InBuf, Gpq, rows_per_bucket);
                    buf4_mat_irrep_close_block(&OutBuf, Gpq, rows_per_bucket);
                } /* Gpq */
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
                for (int h = 0; h < nirreps; h++) {
                    const int r_irrep = h ^ my_irrep;

                    for (int pq = 0; pq < OutBuf.params->rowtot[h]; pq++) {
                        const int p = OutBuf.params->roworb[h][pq][0];
                        const int q = OutBuf.params->roworb[h][pq][1];
                        const int qp = InBuf->params->rowidx[q][p];

                        for (int rs = 0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
                            const int r = OutBuf.params->colorb[r_irrep][rs][0];
                            const int s = OutBuf.params->colorb[r_irrep][rs][1];
                            const int sr = InBuf->params->colidx[s][r];

                            OutBuf.matrix[h][pq][rs] = InBuf->matrix[h][qp][sr];
                        }
                    }
                }
            } else {
                for (int Gpq = 0; Gpq < nirreps; Gpq++) {
                    const int Grs = Gpq ^ my_irrep;

                    /* determine how many rows of OutBuf/InBuf we can store in half the core */
                    rows_per_bucket = dpd_memfree() / (2 * OutBuf.params->coltot[Grs]);
                    if (rows_per_bucket > OutBuf.params->rowtot[Gpq]) rows_per_bucket = OutBuf.params->rowtot[Gpq];
                    nbuckets = (int)ceil((double)OutBuf.params->rowtot[Gpq] / (double)rows_per_bucket);
                    if (nbuckets == 1)
                        rows_left = rows_per_bucket;
                    else
                        rows_left = OutBuf.params->rowtot[Gpq] % rows_per_bucket;

                    /* allocate space for the bucket of rows */
                    buf4_mat_irrep_init_block(&OutBuf, Gpq, rows_per_bucket);
                    buf4_mat_irrep_init_block(InBuf, Gpq, rows_per_bucket);

                    int out_row_start, n;
                    for (n = 0; n < (rows_left ? nbuckets - 1 : nbuckets); n++) {
                        out_row_start = n * rows_per_bucket;

                        for (m = 0; m < (rows_left ? nbuckets - 1 : nbuckets); m++) {
                            in_row_start = m * rows_per_bucket;
                            buf4_mat_irrep_rd_block(InBuf, Gpq, in_row_start, rows_per_bucket);

                            for (int pq = 0; pq < rows_per_bucket; pq++) {
                                /* check to see if this row is contained in the current input-bucket */
                                const int p = OutBuf.params->roworb[Gpq][pq + out_row_start][0];
                                const int q = OutBuf.params->roworb[Gpq][pq + out_row_start][1];
                                const int qp = InBuf->params->rowidx[q][p] - in_row_start;
                                if (qp >= 0 && qp < rows_per_bucket) {
                                    for (int rs = 0; rs < OutBuf.params->coltot[Grs]; rs++) {
                                        const int r = OutBuf.params->colorb[Grs][rs][0];
                                        const int s = OutBuf.params->colorb[Grs][rs][1];
                                        const int sr = InBuf->params->colidx[s][r];
                                        OutBuf.matrix[Gpq][pq][rs] = InBuf->matrix[Gpq][qp][sr];
                                    }
                                }
                            }
                        }
                        if (rows_left) {
                            in_row_start = m * rows_per_bucket;
                            buf4_mat_irrep_rd_block(InBuf, Gpq, in_row_start, rows_left);

                            for (int pq = 0; pq < rows_per_bucket; pq++) {
                                /* check to see if this row is contained in the current input-bucket */
                                const int p = OutBuf.params->roworb[Gpq][pq + out_row_start][0];
                                const int q = OutBuf.params->roworb[Gpq][pq + out_row_start][1];
                                const int qp = InBuf->params->rowidx[q][p] - in_row_start;
                                if (qp >= 0 && qp < rows_left) {
                                    for (int rs = 0; rs < OutBuf.params->coltot[Grs]; rs++) {
                                        const int r = OutBuf.params->colorb[Grs][rs][0];
                                        const int s = OutBuf.params->colorb[Grs][rs][1];
                                        const int sr = InBuf->params->colidx[s][r];
                                        OutBuf.matrix[Gpq][pq][rs] = InBuf->matrix[Gpq][qp][sr];
                                    }
                                }
                            }
                        }

                        buf4_mat_irrep_wrt_block(&OutBuf, Gpq, out_row_start, rows_per_bucket);

                    } /* n */
                    if (rows_left) {
                        out_row_start = n * rows_per_bucket;

                        for (m = 0; m < (rows_left ? nbuckets - 1 : nbuckets); m++) {
                            in_row_start = m * rows_per_bucket;
                            buf4_mat_irrep_rd_block(InBuf, Gpq, in_row_start, rows_per_bucket);

                            for (int pq = 0; pq < rows_left; pq++) {
                                /* check to see if this row is contained in the current input-bucket */
                                const int p = OutBuf.params->roworb[Gpq][pq + out_row_start][0];
                                const int q = OutBuf.params->roworb[Gpq][pq + out_row_start][1];
                                const int qp = InBuf->params->rowidx[q][p] - in_row_start;
                                if (qp >= 0 && qp < rows_per_bucket) {
                                    for (int rs = 0; rs < OutBuf.params->coltot[Grs]; rs++) {
                                        const int r = OutBuf.params->colorb[Grs][rs][0];
                                        const int s = OutBuf.params->colorb[Grs][rs][1];
                                        const int sr = InBuf->params->colidx[s][r];
                                        OutBuf.matrix[Gpq][pq][rs] = InBuf->matrix[Gpq][qp][sr];
                                    }
                                }
                            }
                        }
                        if (rows_left) {
                            in_row_start = m * rows_per_bucket;
                            buf4_mat_irrep_rd_block(InBuf, Gpq, in_row_start, rows_left);

                            for (int pq = 0; pq < rows_left; pq++) {
                                /* check to see if this row is contained in the current input-bucket */
                                const int p = OutBuf.params->roworb[Gpq][pq + out_row_start][0];
                                const int q = OutBuf.params->roworb[Gpq][pq + out_row_start][1];
                                const int qp = InBuf->params->rowidx[q][p] - in_row_start;
                                if (qp >= 0 && qp < rows_left) {
                                    for (int rs = 0; rs < OutBuf.params->coltot[Grs]; rs++) {
                                        const int r = OutBuf.params->colorb[Grs][rs][0];
                                        const int s = OutBuf.params->colorb[Grs][rs][1];
                                        const int sr = InBuf->params->colidx[s][r];
                                        OutBuf.matrix[Gpq][pq][rs] = InBuf->matrix[Gpq][qp][sr];
                                    }
                                }
                            }
                        }
                    }

                    buf4_mat_irrep_wrt_block(&OutBuf, Gpq, out_row_start, rows_left);
                    buf4_mat_irrep_close_block(InBuf, Gpq, rows_per_bucket);
                    buf4_mat_irrep_close_block(&OutBuf, Gpq, rows_per_bucket);

                } /* Gpq */
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
                for (int h = 0; h < nirreps; h++) {
                    const int r_irrep = h ^ my_irrep;

                    for (int Gp = 0; Gp < nirreps; Gp++) {
                        const int Gq = Gp ^ h;
                        for (int Gr = 0; Gr < nirreps; Gr++) {
                            const int Gs = Gr ^ r_irrep;

                            const int Grp = Gr ^ Gp;
                            const int Gqs = Gq ^ Gs;

                            for (int p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                                const int P = OutBuf.params->poff[Gp] + p;
                                for (int q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                                    const int Q = OutBuf.params->qoff[Gq] + q;
                                    const int pq = OutBuf.params->rowidx[P][Q];

                                    for (int r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                                        const int R = OutBuf.params->roff[Gr] + r;
                                        const int rp = InBuf->params->rowidx[R][P];

                                        for (int s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                            const int S = OutBuf.params->soff[Gs] + s;
                                            const int rs = OutBuf.params->colidx[R][S];
                                            const int qs = InBuf->params->colidx[Q][S];

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
                for (int h = 0; h < nirreps; h++) {
                    const int r_irrep = h ^ my_irrep;

                    for (int Gp = 0; Gp < nirreps; Gp++) {
                        const int Gq = Gp ^ h;
                        for (int Gr = 0; Gr < nirreps; Gr++) {
                            const int Gs = Gr ^ r_irrep;

                            const int Gsp = Gs ^ Gp;
                            const int Gqr = Gq ^ Gr;

                            for (int p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                                const int P = OutBuf.params->poff[Gp] + p;
                                for (int q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                                    const int Q = OutBuf.params->qoff[Gq] + q;
                                    const int pq = OutBuf.params->rowidx[P][Q];

                                    for (int r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                                        const int R = OutBuf.params->roff[Gr] + r;
                                        const int qr = InBuf->params->colidx[Q][R];

                                        for (int s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                            const int S = OutBuf.params->soff[Gs] + s;
                                            const int rs = OutBuf.params->colidx[R][S];
                                            const int sp = InBuf->params->rowidx[S][P];

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
                for (int h = 0; h < nirreps; h++) {
                    const int r_irrep = h ^ my_irrep;

                    for (int Gp = 0; Gp < nirreps; Gp++) {
                        const int Gq = Gp ^ h;
                        for (int Gr = 0; Gr < nirreps; Gr++) {
                            const int Gs = Gr ^ r_irrep;

                            const int Grp = Gr ^ Gp;
                            const int Gsq = Gs ^ Gq;

                            for (int p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                                const int P = OutBuf.params->poff[Gp] + p;
                                for (int q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                                    const int Q = OutBuf.params->qoff[Gq] + q;
                                    const int pq = OutBuf.params->rowidx[P][Q];

                                    for (int r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                                        const int R = OutBuf.params->roff[Gr] + r;
                                        const int rp = InBuf->params->rowidx[R][P];

                                        for (int s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                            const int S = OutBuf.params->soff[Gs] + s;
                                            const int rs = OutBuf.params->colidx[R][S];
                                            const int sq = InBuf->params->colidx[S][Q];

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
                for (int h = 0; h < nirreps; h++) {
                    const int r_irrep = h ^ my_irrep;

                    for (int Gp = 0; Gp < nirreps; Gp++) {
                        const int Gq = Gp ^ h;
                        for (int Gr = 0; Gr < nirreps; Gr++) {
                            const int Gs = Gr ^ r_irrep;

                            const int Gsp = Gs ^ Gp;
                            const int Grq = Gr ^ Gq;

                            for (int p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                                const int P = OutBuf.params->poff[Gp] + p;
                                for (int q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                                    const int Q = OutBuf.params->qoff[Gq] + q;
                                    const int pq = OutBuf.params->rowidx[P][Q];

                                    for (int r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                                        const int R = OutBuf.params->roff[Gr] + r;
                                        const int rq = InBuf->params->colidx[R][Q];

                                        for (int s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                            const int S = OutBuf.params->soff[Gs] + s;
                                            const int rs = OutBuf.params->colidx[R][S];
                                            const int sp = InBuf->params->rowidx[S][P];

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
                for (int h = 0; h < nirreps; h++) {
                    const int r_irrep = h ^ my_irrep;

                    for (int Gp = 0; Gp < nirreps; Gp++) {
                        const int Gq = Gp ^ h;
                        for (int Gr = 0; Gr < nirreps; Gr++) {
                            const int Gs = Gr ^ r_irrep;

                            const int Grq = Gr ^ Gq;
                            const int Gps = Gp ^ Gs;

                            for (int p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                                const int P = OutBuf.params->poff[Gp] + p;
                                for (int q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                                    const int Q = OutBuf.params->qoff[Gq] + q;
                                    const int pq = OutBuf.params->rowidx[P][Q];

                                    for (int r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                                        const int R = OutBuf.params->roff[Gr] + r;
                                        const int rq = InBuf->params->rowidx[R][Q];

                                        for (int s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                            const int S = OutBuf.params->soff[Gs] + s;
                                            const int rs = OutBuf.params->colidx[R][S];
                                            const int ps = InBuf->params->colidx[P][S];

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
                for (int h = 0; h < nirreps; h++) {
                    const int r_irrep = h ^ my_irrep;

                    for (int Gp = 0; Gp < nirreps; Gp++) {
                        const int Gq = Gp ^ h;
                        for (int Gr = 0; Gr < nirreps; Gr++) {
                            const int Gs = Gr ^ r_irrep;

                            const int Gsq = Gs ^ Gq;
                            const int Gpr = Gp ^ Gr;

                            for (int p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                                const int P = OutBuf.params->poff[Gp] + p;
                                for (int q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                                    const int Q = OutBuf.params->qoff[Gq] + q;
                                    const int pq = OutBuf.params->rowidx[P][Q];

                                    for (int r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                                        const int R = OutBuf.params->roff[Gr] + r;
                                        const int pr = InBuf->params->colidx[P][R];

                                        for (int s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                            const int S = OutBuf.params->soff[Gs] + s;
                                            const int rs = OutBuf.params->colidx[R][S];
                                            const int sq = InBuf->params->rowidx[S][Q];

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
                for (int h = 0; h < nirreps; h++) {
                    const int r_irrep = h ^ my_irrep;

                    for (int Gp = 0; Gp < nirreps; Gp++) {
                        const int Gq = Gp ^ h;
                        for (int Gr = 0; Gr < nirreps; Gr++) {
                            const int Gs = Gr ^ r_irrep;

                            const int Gqr = Gq ^ Gr;
                            const int Gps = Gp ^ Gs;

                            for (int p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                                const int P = OutBuf.params->poff[Gp] + p;
                                for (int q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                                    const int Q = OutBuf.params->qoff[Gq] + q;
                                    const int pq = OutBuf.params->rowidx[P][Q];

                                    for (int r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                                        const int R = OutBuf.params->roff[Gr] + r;
                                        const int qr = InBuf->params->rowidx[Q][R];

                                        for (int s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                            const int S = OutBuf.params->soff[Gs] + s;
                                            const int rs = OutBuf.params->colidx[R][S];
                                            const int ps = InBuf->params->colidx[P][S];

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
                for (int h = 0; h < nirreps; h++) {
                    const int r_irrep = h ^ my_irrep;

                    for (int Gp = 0; Gp < nirreps; Gp++) {
                        const int Gq = Gp ^ h;
                        for (int Gr = 0; Gr < nirreps; Gr++) {
                            const int Gs = Gr ^ r_irrep;

                            const int Gqs = Gq ^ Gs;
                            const int Gpr = Gp ^ Gr;

                            for (int p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                                const int P = OutBuf.params->poff[Gp] + p;
                                for (int q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                                    const int Q = OutBuf.params->qoff[Gq] + q;
                                    const int pq = OutBuf.params->rowidx[P][Q];

                                    for (int r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                                        const int R = OutBuf.params->roff[Gr] + r;
                                        const int pr = InBuf->params->colidx[P][R];

                                        for (int s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                            const int S = OutBuf.params->soff[Gs] + s;
                                            const int rs = OutBuf.params->colidx[R][S];
                                            const int qs = InBuf->params->rowidx[Q][S];

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
                for (int h = 0; h < nirreps; h++) {
                    const int r_irrep = h ^ my_irrep;

                    for (int pq = 0; pq < OutBuf.params->rowtot[h]; pq++) {
                        const int p = OutBuf.params->roworb[h][pq][0];
                        const int q = OutBuf.params->roworb[h][pq][1];

                        const int col = InBuf->params->colidx[p][q];

                        for (int rs = 0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
                            const int r = OutBuf.params->colorb[r_irrep][rs][0];
                            const int s = OutBuf.params->colorb[r_irrep][rs][1];

                            const int row = InBuf->params->rowidx[s][r];

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
                for (int h = 0; h < nirreps; h++) {
                    const int r_irrep = h ^ my_irrep;

                    for (int pq = 0; pq < OutBuf.params->rowtot[h]; pq++) {
                        const int p = OutBuf.params->roworb[h][pq][0];
                        const int q = OutBuf.params->roworb[h][pq][1];

                        const int col = InBuf->params->colidx[p][q];

                        for (int rs = 0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
                            const int r = OutBuf.params->colorb[r_irrep][rs][0];
                            const int s = OutBuf.params->colorb[r_irrep][rs][1];

                            const int row = InBuf->params->rowidx[r][s];

                            OutBuf.matrix[h][pq][rs] = InBuf->matrix[r_irrep][row][col];
                        }
                    }
                }
            } else {
                for (int Gpq = 0; Gpq < nirreps; Gpq++) {
                    const int Grs = Gpq ^ my_irrep;

                    const int out_rows_per_bucket = [&OutBuf, Grs, Gpq] {
                        int orpb = (dpd_memfree() - OutBuf.params->coltot[Grs]) / (2 * OutBuf.params->coltot[Grs]);
                        if (orpb > OutBuf.params->rowtot[Gpq]) orpb = OutBuf.params->rowtot[Gpq];
                        return orpb;
                    }();

                    const int out_nbuckets =
                        (int)ceil((double)OutBuf.params->rowtot[Gpq] / (double)out_rows_per_bucket);
                    const int out_rows_left =
                        (out_nbuckets == 1) ? out_rows_per_bucket : OutBuf.params->rowtot[Gpq] % out_rows_per_bucket;

                    const int in_rows_per_bucket = [&InBuf, Gpq, Grs] {
                        int irpb = (dpd_memfree() - InBuf->params->coltot[Gpq]) / (2 * InBuf->params->coltot[Gpq]);
                        if (irpb > InBuf->params->rowtot[Grs]) irpb = InBuf->params->rowtot[Grs];
                        return irpb;
                    }();

                    in_nbuckets = (int)ceil((double)InBuf->params->rowtot[Grs] / (double)in_rows_per_bucket);
                    if (in_nbuckets == 1)
                        in_rows_left = in_rows_per_bucket;
                    else
                        in_rows_left = InBuf->params->rowtot[Grs] % in_rows_per_bucket;

#ifdef DPD_DEBUG
                    outfile->Printf(stdout, "Gpq = %d\n", Gpq);
                    outfile->Printf(stdout, "OutBuf.rowtot[Gpq]  = %d\n", OutBuf.params->rowtot[Gpq]);
                    outfile->Printf(stdout, "OutBuf.coltot[Grs]  = %d\n", OutBuf.params->coltot[Grs]);
                    outfile->Printf(stdout, "out_nbuckets        = %d\n", out_nbuckets);
                    outfile->Printf(stdout, "out_rows_per_bucket = %d\n", out_rows_per_bucket);
                    outfile->Printf(stdout, "out_rows_left       = %d\n", out_rows_left);
                    outfile->Printf(stdout, "InBuf.rowtot[Grs]   = %d\n", InBuf->params->rowtot[Grs]);
                    outfile->Printf(stdout, "InBuf.coltot[Gpq]   = %d\n", InBuf->params->coltot[Gpq]);
                    outfile->Printf(stdout, "in_nbuckets         = %d\n", in_nbuckets);
                    outfile->Printf(stdout, "in_rows_per_bucket  = %d\n", in_rows_per_bucket);
                    outfile->Printf(stdout, "in_rows_left        = %d\n", in_rows_left);
                    fflush(stdout);
#endif

                    buf4_mat_irrep_init_block(&OutBuf, Gpq, out_rows_per_bucket);
                    buf4_mat_irrep_init_block(InBuf, Grs, in_rows_per_bucket);

                    for (int n = 0; n < out_nbuckets; n++) {
                        const int out_row_start = n * out_rows_per_bucket;

                        for (m = 0; m < in_nbuckets; m++) {
                            in_row_start = m * in_rows_per_bucket;
                            buf4_mat_irrep_rd_block(InBuf, Grs, in_row_start,
                                                    (m == in_nbuckets - 1 ? in_rows_left : in_rows_per_bucket));

                            for (int pq = 0; pq < (n == out_nbuckets - 1 ? out_rows_left : out_rows_per_bucket); pq++) {
                                const int PQ = pq + n * out_rows_per_bucket;
                                for (int RS = 0; RS < (m == in_nbuckets - 1 ? in_rows_left : in_rows_per_bucket);
                                     RS++) {
                                    const int rs = RS + m * in_rows_per_bucket;
                                    OutBuf.matrix[Gpq][pq][rs] = InBuf->matrix[Grs][RS][PQ];
                                }
                            }
                        }

                        buf4_mat_irrep_wrt_block(&OutBuf, Gpq, out_row_start,
                                                 (n == out_nbuckets - 1 ? out_rows_left : out_rows_per_bucket));
                    }

                    buf4_mat_irrep_close_block(&OutBuf, Gpq, out_rows_per_bucket);
                    buf4_mat_irrep_close_block(InBuf, Grs, in_rows_per_bucket);

                } /* Gpq */
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
                for (int h = 0; h < nirreps; h++) {
                    const int r_irrep = h ^ my_irrep;

                    for (int Gp = 0; Gp < nirreps; Gp++) {
                        const int Gq = Gp ^ h;
                        for (int Gr = 0; Gr < nirreps; Gr++) {
                            const int Gs = Gr ^ r_irrep;

                            const int Gsq = Gs ^ Gq;
                            const int Grp = Gr ^ Gp;

                            for (int p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                                const int P = OutBuf.params->poff[Gp] + p;
                                for (int q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                                    const int Q = OutBuf.params->qoff[Gq] + q;
                                    const int pq = OutBuf.params->rowidx[P][Q];

                                    for (int r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                                        const int R = OutBuf.params->roff[Gr] + r;
                                        const int rp = InBuf->params->colidx[R][P];

                                        for (int s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                            const int S = OutBuf.params->soff[Gs] + s;
                                            const int rs = OutBuf.params->colidx[R][S];
                                            const int sq = InBuf->params->rowidx[S][Q];

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
                for (int h = 0; h < nirreps; h++) {
                    const int r_irrep = h ^ my_irrep;

                    for (int pq = 0; pq < OutBuf.params->rowtot[h]; pq++) {
                        const int p = OutBuf.params->roworb[h][pq][0];
                        const int q = OutBuf.params->roworb[h][pq][1];

                        const int col = InBuf->params->colidx[q][p];

                        for (int rs = 0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
                            const int r = OutBuf.params->colorb[r_irrep][rs][0];
                            const int s = OutBuf.params->colorb[r_irrep][rs][1];

                            const int row = InBuf->params->rowidx[s][r];

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
                for (int h = 0; h < nirreps; h++) {
                    const int r_irrep = h ^ my_irrep;

                    for (int pq = 0; pq < OutBuf.params->rowtot[h]; pq++) {
                        const int p = OutBuf.params->roworb[h][pq][0];
                        const int q = OutBuf.params->roworb[h][pq][1];

                        const int col = InBuf->params->colidx[q][p];

                        for (int rs = 0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
                            const int r = OutBuf.params->colorb[r_irrep][rs][0];
                            const int s = OutBuf.params->colorb[r_irrep][rs][1];

                            const int row = InBuf->params->rowidx[r][s];

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
                for (int h = 0; h < nirreps; h++) {
                    const int r_irrep = h ^ my_irrep;

                    for (int Gp = 0; Gp < nirreps; Gp++) {
                        const int Gq = Gp ^ h;
                        for (int Gr = 0; Gr < nirreps; Gr++) {
                            const int Gs = Gr ^ r_irrep;

                            const int Gqr = Gq ^ Gr;
                            const int Gsp = Gs ^ Gp;

                            for (int p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                                const int P = OutBuf.params->poff[Gp] + p;
                                for (int q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                                    const int Q = OutBuf.params->qoff[Gq] + q;
                                    const int pq = OutBuf.params->rowidx[P][Q];

                                    for (int r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                                        const int R = OutBuf.params->roff[Gr] + r;
                                        const int qr = InBuf->params->rowidx[Q][R];

                                        for (int s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                            const int S = OutBuf.params->soff[Gs] + s;
                                            const int rs = OutBuf.params->colidx[R][S];
                                            const int sp = InBuf->params->colidx[S][P];

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
                for (int h = 0; h < nirreps; h++) {
                    const int r_irrep = h ^ my_irrep;

                    for (int Gp = 0; Gp < nirreps; Gp++) {
                        const int Gq = Gp ^ h;
                        for (int Gr = 0; Gr < nirreps; Gr++) {
                            const int Gs = Gr ^ r_irrep;

                            const int Gqs = Gq ^ Gs;
                            const int Grp = Gr ^ Gp;

                            for (int p = 0; p < OutBuf.params->ppi[Gp]; p++) {
                                const int P = OutBuf.params->poff[Gp] + p;
                                for (int q = 0; q < OutBuf.params->qpi[Gq]; q++) {
                                    const int Q = OutBuf.params->qoff[Gq] + q;
                                    const int pq = OutBuf.params->rowidx[P][Q];

                                    for (int r = 0; r < OutBuf.params->rpi[Gr]; r++) {
                                        const int R = OutBuf.params->roff[Gr] + r;
                                        const int rp = InBuf->params->colidx[R][P];

                                        for (int s = 0; s < OutBuf.params->spi[Gs]; s++) {
                                            const int S = OutBuf.params->soff[Gs] + s;
                                            const int rs = OutBuf.params->colidx[R][S];
                                            const int qs = InBuf->params->rowidx[Q][S];

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
        for (int h = 0; h < nirreps; h++) {
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

int DPD::buf4_sort(dpdbuf4* InBuf, const int outfilenum, const enum indices index, const std::string pq,
                   const std::string rs, const std::string& label) {
    return buf4_sort(InBuf, outfilenum, index, pairnum(pq), pairnum(rs), label);
}

}  // namespace psi
