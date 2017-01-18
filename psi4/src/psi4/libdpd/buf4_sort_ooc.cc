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
#include <cmath>
#include "psi4/libqt/qt.h"
#include "dpd.h"
#include "psi4/psi4-dec.h"
namespace psi {

/*
** dpd_buf4_sort_ooc(): A general DPD buffer sorting function that will
** (eventually) handle all 24 possible permutations of four-index
** buffers.  This version uses an out-of-core algorithm that should only be
** applied to large cases.  See the comments in dpd_buf4_sort() for
** argument details.
**
** TDC
** May 2000
*/

int DPD::buf4_sort_ooc(dpdbuf4 *InBuf, int outfilenum, enum indices index,
                       int pqnum, int rsnum, const char *label)
{
    int h,nirreps, row, col, all_buf_irrep, r_irrep;
    int p, q, r, s, P, Q, R, S, pq, rs, sr, pr, qs, qp, rq, qr, ps, sp, rp, sq;
    int Gp, Gq, Gr, Gs, Gpq, Grs, Gpr, Gqs, Grq, Gqr, Gps, Gsp, Grp, Gsq;
    int memoryd, rows_per_bucket, nbuckets, rows_left, incore, n;
    dpdbuf4 OutBuf;

    nirreps = InBuf->params->nirreps;
    all_buf_irrep = InBuf->file.my_irrep;

#ifdef DPD_TIMER
    timer_on("buf4_sort");
#endif

    buf4_init(&OutBuf, outfilenum, all_buf_irrep, pqnum, rsnum,
              pqnum, rsnum, 0, label);

    for(h=0; h < nirreps; h++) {

        r_irrep = h^all_buf_irrep;

        switch(index) {
        case pqrs:
            outfile->Printf( "\nDPD sort error: invalid index ordering.\n");
            dpd_error("buf_sort", "outfile");
            break;

        case pqsr:

#ifdef DPD_TIMER
            timer_on("pqsr");
#endif

            /* p->p; q->q; s->r; r->s = pqsr */

            /* select algorithm for certain simple cases */
            memoryd = dpd_memfree()/2; /* use half the memory for each buf4 in the sort */
            if(InBuf->params->rowtot[h] && InBuf->params->coltot[h^all_buf_irrep]) {

                rows_per_bucket = memoryd/InBuf->params->coltot[h^all_buf_irrep];

                /* enough memory for the whole matrix? */
                if(rows_per_bucket > InBuf->params->rowtot[h])
                    rows_per_bucket = InBuf->params->rowtot[h];

                if(!rows_per_bucket) dpd_error("buf4_sort_pqsr: Not enough memory for one row!", "outfile");

                nbuckets = (int) ceil(((double) InBuf->params->rowtot[h])/((double) rows_per_bucket));

                rows_left = InBuf->params->rowtot[h] % rows_per_bucket;

                incore = 1;
                if(nbuckets > 1) {
                    incore = 0;
#if DPD_DEBUG
                    outfile->Printf( "buf4_sort_pqsr: memory information.\n");
                    outfile->Printf( "buf4_sort_pqsr: rowtot[%d] = %d\n", h, InBuf->params->rowtot[h]);
                    outfile->Printf( "buf4_sort_pqsr: nbuckets = %d\n", nbuckets);
                    outfile->Printf( "buf4_sort_pqsr: rows_per_bucket = %d\n", rows_per_bucket);
                    outfile->Printf( "buf4_sort_pqsr: rows_left = %d\n", rows_left);
                    outfile->Printf( "buf4_sort_pqsr: out-of-core algorithm used\n");
#endif
                }

            }
            else incore = 1;

            if(incore) {

                buf4_mat_irrep_init(&OutBuf, h);

                buf4_mat_irrep_init(InBuf, h);
                buf4_mat_irrep_rd(InBuf, h);

                for(pq=0; pq < OutBuf.params->rowtot[h]; pq++) {
                    p = OutBuf.params->roworb[h][pq][0];
                    q = OutBuf.params->roworb[h][pq][1];

                    row = InBuf->params->rowidx[p][q];

                    for(rs=0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
                        r = OutBuf.params->colorb[r_irrep][rs][0];
                        s = OutBuf.params->colorb[r_irrep][rs][1];

                        sr = InBuf->params->colidx[s][r];

                        OutBuf.matrix[h][pq][rs] = InBuf->matrix[h][row][sr];
                    }
                }

                buf4_mat_irrep_close(InBuf, h);
                buf4_mat_irrep_wrt(&OutBuf, h);

                buf4_mat_irrep_close(&OutBuf, h);
            }
            else {  /* out-of-core sort option */

                buf4_mat_irrep_init_block(InBuf, h, rows_per_bucket);
                buf4_mat_irrep_init_block(&OutBuf, h, rows_per_bucket);

                for(n=0; n < (rows_left ? nbuckets-1 : nbuckets); n++) {

                    buf4_mat_irrep_rd_block(InBuf, h, n*rows_per_bucket, rows_per_bucket);

                    for(pq=0; pq < rows_per_bucket; pq++) {
                        for(rs=0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
                            r = OutBuf.params->colorb[r_irrep][rs][0];
                            s = OutBuf.params->colorb[r_irrep][rs][1];

                            sr = InBuf->params->colidx[s][r];

                            OutBuf.matrix[h][pq][rs] = InBuf->matrix[h][pq][sr];
                        }
                    }

                    buf4_mat_irrep_wrt_block(&OutBuf, h, n*rows_per_bucket, rows_per_bucket);
                }
                if(rows_left) {

                    buf4_mat_irrep_rd_block(InBuf, h, n*rows_per_bucket, rows_left);

                    for(pq=0; pq < rows_left; pq++) {
                        for(rs=0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
                            r = OutBuf.params->colorb[r_irrep][rs][0];
                            s = OutBuf.params->colorb[r_irrep][rs][1];

                            sr = InBuf->params->colidx[s][r];

                            OutBuf.matrix[h][pq][rs] = InBuf->matrix[h][pq][sr];
                        }
                    }

                    buf4_mat_irrep_wrt_block(&OutBuf, h, n*rows_per_bucket, rows_left);
                }

                buf4_mat_irrep_close_block(InBuf, h, rows_per_bucket);
                buf4_mat_irrep_close_block(&OutBuf, h, rows_per_bucket);

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
            buf4_mat_irrep_init(&OutBuf, h);

            for(Gp=0; Gp < nirreps; Gp++) {
                Gq = Gp^h;
                for(Gr=0; Gr < nirreps; Gr++) {
                    Gs = Gr^r_irrep;

                    /* Irreps on the source */
                    Gpr = Gp^Gr;  Gqs = Gq^Gs;

                    buf4_mat_irrep_init(InBuf, Gpr);
                    buf4_mat_irrep_rd(InBuf, Gpr);

                    for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
                        P = OutBuf.params->poff[Gp] + p;
                        for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
                            Q = OutBuf.params->qoff[Gq] + q;
                            pq = OutBuf.params->rowidx[P][Q];

                            for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
                                R = OutBuf.params->roff[Gr] + r;
                                pr = InBuf->params->rowidx[P][R];

                                for(s=0; s < OutBuf.params->spi[Gs]; s++) {
                                    S = OutBuf.params->soff[Gs] + s;
                                    rs = OutBuf.params->colidx[R][S];
                                    qs = InBuf->params->colidx[Q][S];

                                    OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gpr][pr][qs];

                                }
                            }
                        }
                    }

                    buf4_mat_irrep_close(InBuf, Gpr);
                }
            }

            buf4_mat_irrep_wrt(&OutBuf, h);
            buf4_mat_irrep_close(&OutBuf, h);

#ifdef DPD_TIMER
            timer_off("prqs");
#endif

            break;

        case prsq:

#ifdef DPD_TIMER
            timer_on("prsq");
#endif

            /* p->p; r->q; s->r; q->s = psqr */
            buf4_mat_irrep_init(&OutBuf, h);

            for(Gp=0; Gp < nirreps; Gp++) {
                Gq = Gp^h;
                for(Gr=0; Gr < nirreps; Gr++) {
                    Gs = Gr^r_irrep;

                    Gps = Gp^Gs;  Gqr = Gq^Gr;

                    buf4_mat_irrep_init(InBuf, Gps);
                    buf4_mat_irrep_rd(InBuf, Gps);

                    for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
                        P = OutBuf.params->poff[Gp] + p;
                        for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
                            Q = OutBuf.params->qoff[Gq] + q;
                            pq = OutBuf.params->rowidx[P][Q];

                            for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
                                R = OutBuf.params->roff[Gr] + r;
                                qr = InBuf->params->colidx[Q][R];

                                for(s=0; s < OutBuf.params->spi[Gs]; s++) {
                                    S = OutBuf.params->soff[Gs] + s;
                                    rs = OutBuf.params->colidx[R][S];
                                    ps = InBuf->params->rowidx[P][S];

                                    OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gps][ps][qr];

                                }
                            }
                        }
                    }

                    buf4_mat_irrep_close(InBuf, Gps);
                }
            }

            buf4_mat_irrep_wrt(&OutBuf, h);
            buf4_mat_irrep_close(&OutBuf, h);

#ifdef DPD_TIMER
            timer_off("prsq");
#endif

            break;

        case psqr:

#ifdef DPD_TIMER
            timer_on("psqr");
#endif

            /* p->p; s->q; q->r; r->s = prsq */
            buf4_mat_irrep_init(&OutBuf, h);

            for(Gp=0; Gp < nirreps; Gp++) {
                Gq = Gp^h;
                for(Gr=0; Gr < nirreps; Gr++) {
                    Gs = Gr^r_irrep;

                    Gpr = Gp^Gr;  Gsq = Gs^Gq;

                    buf4_mat_irrep_init(InBuf, Gpr);
                    buf4_mat_irrep_rd(InBuf, Gpr);

                    for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
                        P = OutBuf.params->poff[Gp] + p;
                        for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
                            Q = OutBuf.params->qoff[Gq] + q;
                            pq = OutBuf.params->rowidx[P][Q];

                            for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
                                R = OutBuf.params->roff[Gr] + r;
                                pr = InBuf->params->rowidx[P][R];

                                for(s=0; s < OutBuf.params->spi[Gs]; s++) {
                                    S = OutBuf.params->soff[Gs] + s;
                                    rs = OutBuf.params->colidx[R][S];
                                    sq = InBuf->params->colidx[S][Q];

                                    OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gpr][pr][sq];

                                }
                            }
                        }
                    }

                    buf4_mat_irrep_close(InBuf, Gpr);
                }
            }
            buf4_mat_irrep_wrt(&OutBuf, h);
            buf4_mat_irrep_close(&OutBuf, h);

#ifdef DPD_TIMER
            timer_off("psqr");
#endif
            break;

        case psrq:

#ifdef DPD_TIMER
            timer_on("psrq");
#endif

            /* p->p; s->q; r->r; q->s = psrq */
            buf4_mat_irrep_init(&OutBuf, h);

            for(Gp=0; Gp < nirreps; Gp++) {
                Gq = Gp^h;
                for(Gr=0; Gr < nirreps; Gr++) {
                    Gs = Gr^r_irrep;

                    Gps = Gp^Gs;  Grq = Gr^Gq;

                    buf4_mat_irrep_init(InBuf, Gps);
                    buf4_mat_irrep_rd(InBuf, Gps);

                    for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
                        P = OutBuf.params->poff[Gp] + p;
                        for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
                            Q = OutBuf.params->qoff[Gq] + q;
                            pq = OutBuf.params->rowidx[P][Q];

                            for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
                                R = OutBuf.params->roff[Gr] + r;
                                rq = InBuf->params->colidx[R][Q];

                                for(s=0; s < OutBuf.params->spi[Gs]; s++) {
                                    S = OutBuf.params->soff[Gs] + s;
                                    rs = OutBuf.params->colidx[R][S];
                                    ps = InBuf->params->rowidx[P][S];

                                    OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gps][ps][rq];

                                }
                            }
                        }
                    }

                    buf4_mat_irrep_close(InBuf, Gps);
                }
            }
            buf4_mat_irrep_wrt(&OutBuf, h);
            buf4_mat_irrep_close(&OutBuf, h);

#ifdef DPD_TIMER
            timer_off("psrq");
#endif
            break;

        case qprs:

#ifdef DPD_TIMER
            timer_on("qprs");
#endif

            /* q->p; p->q; r->r; s->s = qprs */
            buf4_mat_irrep_init(&OutBuf, h);

            buf4_mat_irrep_init(InBuf, h);
            buf4_mat_irrep_rd(InBuf, h);

            for(pq=0; pq < OutBuf.params->rowtot[h]; pq++) {
                p = OutBuf.params->roworb[h][pq][0];
                q = OutBuf.params->roworb[h][pq][1];
                qp = InBuf->params->rowidx[q][p];

                for(rs=0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
                    r = OutBuf.params->colorb[r_irrep][rs][0];
                    s = OutBuf.params->colorb[r_irrep][rs][1];

                    col = InBuf->params->colidx[r][s];

                    OutBuf.matrix[h][pq][rs] = InBuf->matrix[h][qp][col];
                }
            }

            buf4_mat_irrep_close(InBuf, h);
            buf4_mat_irrep_wrt(&OutBuf, h);
            buf4_mat_irrep_close(&OutBuf, h);

#ifdef DPD_TIMER
            timer_off("qprs");
#endif
            break;

        case qpsr:

#ifdef DPD_TIMER
            timer_on("qpsr");
#endif

            /* q->p; p->q; s->r; r->s = qpsr */
            buf4_mat_irrep_init(&OutBuf, h);

            buf4_mat_irrep_init(InBuf, h);
            buf4_mat_irrep_rd(InBuf, h);

            for(pq=0; pq < OutBuf.params->rowtot[h]; pq++) {
                p = OutBuf.params->roworb[h][pq][0];
                q = OutBuf.params->roworb[h][pq][1];
                qp = InBuf->params->rowidx[q][p];

                for(rs=0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
                    r = OutBuf.params->colorb[r_irrep][rs][0];
                    s = OutBuf.params->colorb[r_irrep][rs][1];
                    sr = InBuf->params->colidx[s][r];

                    OutBuf.matrix[h][pq][rs] = InBuf->matrix[h][qp][sr];
                }
            }

            buf4_mat_irrep_close(InBuf, h);
            buf4_mat_irrep_wrt(&OutBuf, h);
            buf4_mat_irrep_close(&OutBuf, h);

#ifdef DPD_TIMER
            timer_off("qpsr");
#endif
            break;

        case qrps:
#ifdef DPD_TIMER
            timer_on("qrps");
#endif

            /* q->p; r->q; p->r; s->s = rpqs */
            buf4_mat_irrep_init(&OutBuf, h);

            for(Gp=0; Gp < nirreps; Gp++) {
                Gq = Gp^h;
                for(Gr=0; Gr < nirreps; Gr++) {
                    Gs = Gr^r_irrep;

                    Grp = Gr^Gp; Gqs = Gq^Gs;

                    buf4_mat_irrep_init(InBuf, Grp);
                    buf4_mat_irrep_rd(InBuf, Grp);

                    for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
                        P = OutBuf.params->poff[Gp] + p;
                        for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
                            Q = OutBuf.params->qoff[Gq] + q;
                            pq = OutBuf.params->rowidx[P][Q];

                            for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
                                R = OutBuf.params->roff[Gr] + r;
                                rp = InBuf->params->rowidx[R][P];

                                for(s=0; s < OutBuf.params->spi[Gs]; s++) {
                                    S = OutBuf.params->soff[Gs] + s;
                                    rs = OutBuf.params->colidx[R][S];
                                    qs = InBuf->params->colidx[Q][S];

                                    OutBuf.matrix[h][pq][rs] = InBuf->matrix[Grp][rp][qs];

                                }
                            }
                        }
                    }

                    buf4_mat_irrep_close(InBuf, Grp);
                }
            }

            buf4_mat_irrep_wrt(&OutBuf, h);
            buf4_mat_irrep_close(&OutBuf, h);

#ifdef DPD_TIMER
            timer_off("qrps");
#endif
            break;

        case qrsp:

#ifdef DPD_TIMER
            timer_on("qrsp");
#endif

            /* q->p; r->q; s->r; p->s = spqr */
            buf4_mat_irrep_init(&OutBuf, h);

            for(Gp=0; Gp < nirreps; Gp++) {
                Gq = Gp^h;
                for(Gr=0; Gr < nirreps; Gr++) {
                    Gs = Gr^r_irrep;

                    Gsp = Gs^Gp; Gqr = Gq^Gr;

                    buf4_mat_irrep_init(InBuf, Gsp);
                    buf4_mat_irrep_rd(InBuf, Gsp);

                    for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
                        P = OutBuf.params->poff[Gp] + p;
                        for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
                            Q = OutBuf.params->qoff[Gq] + q;
                            pq = OutBuf.params->rowidx[P][Q];

                            for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
                                R = OutBuf.params->roff[Gr] + r;
                                qr = InBuf->params->colidx[Q][R];

                                for(s=0; s < OutBuf.params->spi[Gs]; s++) {
                                    S = OutBuf.params->soff[Gs] + s;
                                    rs = OutBuf.params->colidx[R][S];
                                    sp = InBuf->params->rowidx[S][P];

                                    OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gsp][sp][qr];

                                }
                            }
                        }
                    }

                    buf4_mat_irrep_close(InBuf, Gsp);
                }
            }

            buf4_mat_irrep_wrt(&OutBuf, h);
            buf4_mat_irrep_close(&OutBuf, h);

#ifdef DPD_TIMER
            timer_off("qrsp");
#endif
            break;

        case qspr:
            outfile->Printf("\nDPD sort error: index ordering not yet coded.\n");
            dpd_error("buf_sort", "outfile");
            break;

        case qsrp:
            outfile->Printf("\nDPD sort error: index ordering not yet coded.\n");
            dpd_error("buf_sort", "outfile");
            break;

        case rqps:

#ifdef DPD_TIMER
            timer_on("rqps");
#endif

            /* r->p; q->q; p->r; s->s = rqps */
            buf4_mat_irrep_init(&OutBuf, h);

            for(Gp=0; Gp < nirreps; Gp++) {
                Gq = Gp^h;
                for(Gr=0; Gr < nirreps; Gr++) {
                    Gs = Gr^r_irrep;

                    Grq = Gr^Gq; Gps = Gp^Gs;

                    buf4_mat_irrep_init(InBuf, Grq);
                    buf4_mat_irrep_rd(InBuf, Grq);

                    for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
                        P = OutBuf.params->poff[Gp] + p;
                        for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
                            Q = OutBuf.params->qoff[Gq] + q;
                            pq = OutBuf.params->rowidx[P][Q];

                            for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
                                R = OutBuf.params->roff[Gr] + r;
                                rq = InBuf->params->rowidx[R][Q];

                                for(s=0; s < OutBuf.params->spi[Gs]; s++) {
                                    S = OutBuf.params->soff[Gs] + s;
                                    rs = OutBuf.params->colidx[R][S];
                                    ps = InBuf->params->colidx[P][S];

                                    OutBuf.matrix[h][pq][rs] = InBuf->matrix[Grq][rq][ps];

                                }
                            }
                        }
                    }

                    buf4_mat_irrep_close(InBuf, Grq);
                }
            }
            buf4_mat_irrep_wrt(&OutBuf, h);
            buf4_mat_irrep_close(&OutBuf, h);

#ifdef DPD_TIMER
            timer_off("rqps");
#endif
            break;

        case rqsp:

#ifdef DPD_TIMER
            timer_on("rqsp");
#endif

            /* r->p; q->q; s->r; p->s = sqpr */
            buf4_mat_irrep_init(&OutBuf, h);

            for(Gp=0; Gp < nirreps; Gp++) {
                Gq = Gp^h;
                for(Gr=0; Gr < nirreps; Gr++) {
                    Gs = Gr^r_irrep;

                    Gsq = Gs^Gq;  Gpr = Gp^Gr;

                    buf4_mat_irrep_init(InBuf, Gsq);
                    buf4_mat_irrep_rd(InBuf, Gsq);

                    for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
                        P = OutBuf.params->poff[Gp] + p;
                        for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
                            Q = OutBuf.params->qoff[Gq] + q;
                            pq = OutBuf.params->rowidx[P][Q];

                            for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
                                R = OutBuf.params->roff[Gr] + r;
                                pr = InBuf->params->colidx[P][R];

                                for(s=0; s < OutBuf.params->spi[Gs]; s++) {
                                    S = OutBuf.params->soff[Gs] + s;
                                    rs = OutBuf.params->colidx[R][S];
                                    sq = InBuf->params->rowidx[S][Q];

                                    OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gsq][sq][pr];

                                }
                            }
                        }
                    }

                    buf4_mat_irrep_close(InBuf, Gsq);
                }
            }

            buf4_mat_irrep_wrt(&OutBuf, h);
            buf4_mat_irrep_close(&OutBuf, h);

#ifdef DPD_TIMER
            timer_off("rqsp");
#endif

            break;

        case rpqs:

#ifdef DPD_TIMER
            timer_on("rpqs");
#endif

            /* r->p; p->q; q->r; s->s = qrps */
            buf4_mat_irrep_init(&OutBuf, h);

            for(Gp=0; Gp < nirreps; Gp++) {
                Gq = Gp^h;
                for(Gr=0; Gr < nirreps; Gr++) {
                    Gs = Gr^r_irrep;

                    Gqr = Gq^Gr;  Gps = Gp^Gs;

                    buf4_mat_irrep_init(InBuf, Gqr);
                    buf4_mat_irrep_rd(InBuf, Gqr);

                    for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
                        P = OutBuf.params->poff[Gp] + p;
                        for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
                            Q = OutBuf.params->qoff[Gq] + q;
                            pq = OutBuf.params->rowidx[P][Q];

                            for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
                                R = OutBuf.params->roff[Gr] + r;
                                qr = InBuf->params->rowidx[Q][R];

                                for(s=0; s < OutBuf.params->spi[Gs]; s++) {
                                    S = OutBuf.params->soff[Gs] + s;
                                    rs = OutBuf.params->colidx[R][S];
                                    ps = InBuf->params->colidx[P][S];

                                    OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gqr][qr][ps];

                                }
                            }
                        }
                    }

                    buf4_mat_irrep_close(InBuf, Gqr);
                }
            }

            buf4_mat_irrep_wrt(&OutBuf, h);
            buf4_mat_irrep_close(&OutBuf, h);

#ifdef DPD_TIMER
            timer_off("rpqs");
#endif

            break;

        case rpsq:

#ifdef DPD_TIMER
            timer_on("rpsq");
#endif

            /* r->p; p->q; s->r; q->s = qspr */
            buf4_mat_irrep_init(&OutBuf, h);

            for(Gp=0; Gp < nirreps; Gp++) {
                Gq = Gp^h;
                for(Gr=0; Gr < nirreps; Gr++) {
                    Gs = Gr^r_irrep;

                    Gqs = Gq^Gs;  Gpr = Gp^Gr;

                    buf4_mat_irrep_init(InBuf, Gqs);
                    buf4_mat_irrep_rd(InBuf, Gqs);

                    for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
                        P = OutBuf.params->poff[Gp] + p;
                        for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
                            Q = OutBuf.params->qoff[Gq] + q;
                            pq = OutBuf.params->rowidx[P][Q];

                            for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
                                R = OutBuf.params->roff[Gr] + r;
                                pr = InBuf->params->colidx[P][R];

                                for(s=0; s < OutBuf.params->spi[Gs]; s++) {
                                    S = OutBuf.params->soff[Gs] + s;
                                    rs = OutBuf.params->colidx[R][S];
                                    qs = InBuf->params->rowidx[Q][S];

                                    OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gqs][qs][pr];

                                }
                            }
                        }
                    }

                    buf4_mat_irrep_close(InBuf, Gqs);
                }
            }

            buf4_mat_irrep_wrt(&OutBuf, h);
            buf4_mat_irrep_close(&OutBuf, h);

#ifdef DPD_TIMER
            timer_off("rpsq");
#endif

            break;

        case rsqp:

#ifdef DPD_TIMER
            timer_on("rsqp");
#endif

            /* r->p; s->q; q->r; p->s = srpq */
            buf4_mat_irrep_init(&OutBuf, h);

            buf4_mat_irrep_init(InBuf, h);
            buf4_mat_irrep_rd(InBuf, h);

            for(pq=0; pq < OutBuf.params->rowtot[h]; pq++) {
                p = OutBuf.params->roworb[h][pq][0];
                q = OutBuf.params->roworb[h][pq][1];

                col = InBuf->params->colidx[p][q];

                for(rs=0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
                    r = OutBuf.params->colorb[r_irrep][rs][0];
                    s = OutBuf.params->colorb[r_irrep][rs][1];

                    row = InBuf->params->rowidx[s][r];

                    OutBuf.matrix[h][pq][rs] = InBuf->matrix[h][row][col];

                }
            }

            buf4_mat_irrep_close(InBuf, h);

            buf4_mat_irrep_wrt(&OutBuf, h);
            buf4_mat_irrep_close(&OutBuf, h);

#ifdef DPD_TIMER
            timer_off("rsqp");
#endif

            break;

        case rspq:

#ifdef DPD_TIMER
            timer_on("rspq");
#endif

            /* r->p; s->q; p->r; q->s = rspq */
            buf4_mat_irrep_init(&OutBuf, h);

            buf4_mat_irrep_init(InBuf, h);
            buf4_mat_irrep_rd(InBuf, h);

            for(pq=0; pq < OutBuf.params->rowtot[h]; pq++) {
                p = OutBuf.params->roworb[h][pq][0];
                q = OutBuf.params->roworb[h][pq][1];

                col = InBuf->params->colidx[p][q];

                for(rs=0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
                    r = OutBuf.params->colorb[r_irrep][rs][0];
                    s = OutBuf.params->colorb[r_irrep][rs][1];

                    row = InBuf->params->rowidx[r][s];

                    OutBuf.matrix[h][pq][rs] = InBuf->matrix[h][row][col];

                }
            }

            buf4_mat_irrep_close(InBuf, h);

            buf4_mat_irrep_wrt(&OutBuf, h);
            buf4_mat_irrep_close(&OutBuf, h);

#ifdef DPD_TIMER
            timer_off("rspq");
#endif

            break;

        case sqrp:

#ifdef DPD_TIMER
            timer_on("sqrp");
#endif

            /* s->p; q->q; r->r; p->s = sqrp */
            buf4_mat_irrep_init(&OutBuf, h);

            for(Gp=0; Gp < nirreps; Gp++) {
                Gq = Gp^h;
                for(Gr=0; Gr < nirreps; Gr++) {
                    Gs = Gr^r_irrep;

                    Gsq = Gs^Gq;  Grp = Gr^Gp;

                    buf4_mat_irrep_init(InBuf, Gsq);
                    buf4_mat_irrep_rd(InBuf, Gsq);

                    for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
                        P = OutBuf.params->poff[Gp] + p;
                        for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
                            Q = OutBuf.params->qoff[Gq] + q;
                            pq = OutBuf.params->rowidx[P][Q];

                            for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
                                R = OutBuf.params->roff[Gr] + r;
                                rp = InBuf->params->colidx[R][P];

                                for(s=0; s < OutBuf.params->spi[Gs]; s++) {
                                    S = OutBuf.params->soff[Gs] + s;
                                    rs = OutBuf.params->colidx[R][S];
                                    sq = InBuf->params->rowidx[S][Q];

                                    OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gsq][sq][rp];

                                }
                            }
                        }
                    }
                    buf4_mat_irrep_close(InBuf, Gsq);
                }
            }
            buf4_mat_irrep_wrt(&OutBuf, h);
            buf4_mat_irrep_close(&OutBuf, h);

#ifdef DPD_TIMER
            timer_off("sqrp");
#endif
            break;

        case sqpr:
            outfile->Printf("\nDPD sort error: index ordering not yet coded.\n");
            dpd_error("buf_sort", "outfile");
            break;

        case srqp:

#ifdef DPD_TIMER
            timer_on("srqp");
#endif

            /* s->p; r->q; q->r; p->s = srqp */
            buf4_mat_irrep_init(&OutBuf, h);

            buf4_mat_irrep_init(InBuf, h);
            buf4_mat_irrep_rd(InBuf, h);

            for(pq=0; pq < OutBuf.params->rowtot[h]; pq++) {
                p = OutBuf.params->roworb[h][pq][0];
                q = OutBuf.params->roworb[h][pq][1];

                col = InBuf->params->colidx[q][p];

                for(rs=0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
                    r = OutBuf.params->colorb[r_irrep][rs][0];
                    s = OutBuf.params->colorb[r_irrep][rs][1];

                    row = InBuf->params->rowidx[s][r];

                    OutBuf.matrix[h][pq][rs] = InBuf->matrix[h][row][col];

                }
            }

            buf4_mat_irrep_close(InBuf, h);

            buf4_mat_irrep_wrt(&OutBuf, h);
            buf4_mat_irrep_close(&OutBuf, h);

#ifdef DPD_TIMER
            timer_off("srqp");
#endif
            break;

        case srpq:
#ifdef DPD_TIMER
            timer_on("srpq");
#endif

            /* s->p; r->q; p->r; q->s = rsqp */
            buf4_mat_irrep_init(&OutBuf, h);

            buf4_mat_irrep_init(InBuf, h);
            buf4_mat_irrep_rd(InBuf, h);

            for(pq=0; pq < OutBuf.params->rowtot[h]; pq++) {
                p = OutBuf.params->roworb[h][pq][0];
                q = OutBuf.params->roworb[h][pq][1];

                col = InBuf->params->colidx[q][p];

                for(rs=0; rs < OutBuf.params->coltot[r_irrep]; rs++) {
                    r = OutBuf.params->colorb[r_irrep][rs][0];
                    s = OutBuf.params->colorb[r_irrep][rs][1];

                    row = InBuf->params->rowidx[r][s];

                    OutBuf.matrix[h][pq][rs] = InBuf->matrix[h][row][col];

                }
            }

            buf4_mat_irrep_close(InBuf, h);

            buf4_mat_irrep_wrt(&OutBuf, h);
            buf4_mat_irrep_close(&OutBuf, h);

#ifdef DPD_TIMER
            timer_off("srpq");
#endif

            break;

        case spqr:

#ifdef DPD_TIMER
            timer_on("spqr");
#endif

            /* s->p; p->q; q->r; r->s = qrsp */
            buf4_mat_irrep_init(&OutBuf, h);

            for(Gp=0; Gp < nirreps; Gp++) {
                Gq = Gp^h;
                for(Gr=0; Gr < nirreps; Gr++) {
                    Gs = Gr^r_irrep;

                    Gqr = Gq^Gr;  Gsp = Gs^Gp;

                    buf4_mat_irrep_init(InBuf, Gqr);
                    buf4_mat_irrep_rd(InBuf, Gqr);

                    for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
                        P = OutBuf.params->poff[Gp] + p;
                        for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
                            Q = OutBuf.params->qoff[Gq] + q;
                            pq = OutBuf.params->rowidx[P][Q];

                            for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
                                R = OutBuf.params->roff[Gr] + r;
                                qr = InBuf->params->rowidx[Q][R];

                                for(s=0; s < OutBuf.params->spi[Gs]; s++) {
                                    S = OutBuf.params->soff[Gs] + s;
                                    rs = OutBuf.params->colidx[R][S];
                                    sp = InBuf->params->colidx[S][P];

                                    OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gqr][qr][sp];

                                }
                            }
                        }
                    }
                    buf4_mat_irrep_close(InBuf, Gqr);
                }
            }
            buf4_mat_irrep_wrt(&OutBuf, h);
            buf4_mat_irrep_close(&OutBuf, h);

#ifdef DPD_TIMER
            timer_off("spqr");
#endif
            break;

        case sprq:

#ifdef DPD_TIMER
            timer_on("sprq");
#endif

            /* s->p; p->q; r->r; q->s = qsrp */
            buf4_mat_irrep_init(&OutBuf, h);

            for(Gp=0; Gp < nirreps; Gp++) {
                Gq = Gp^h;
                for(Gr=0; Gr < nirreps; Gr++) {
                    Gs = Gr^r_irrep;

                    Gqs = Gq^Gs;  Grp = Gr^Gp;

                    buf4_mat_irrep_init(InBuf, Gqs);
                    buf4_mat_irrep_rd(InBuf, Gqs);

                    for(p=0; p < OutBuf.params->ppi[Gp]; p++) {
                        P = OutBuf.params->poff[Gp] + p;
                        for(q=0; q < OutBuf.params->qpi[Gq]; q++) {
                            Q = OutBuf.params->qoff[Gq] + q;
                            pq = OutBuf.params->rowidx[P][Q];

                            for(r=0; r < OutBuf.params->rpi[Gr]; r++) {
                                R = OutBuf.params->roff[Gr] + r;
                                rp = InBuf->params->colidx[R][P];

                                for(s=0; s < OutBuf.params->spi[Gs]; s++) {
                                    S = OutBuf.params->soff[Gs] + s;
                                    rs = OutBuf.params->colidx[R][S];
                                    qs = InBuf->params->rowidx[Q][S];

                                    OutBuf.matrix[h][pq][rs] = InBuf->matrix[Gqs][qs][rp];

                                }
                            }
                        }
                    }
                    buf4_mat_irrep_close(InBuf, Gqs);
                }
            }

            buf4_mat_irrep_wrt(&OutBuf, h);
            buf4_mat_irrep_close(&OutBuf, h);

#ifdef DPD_TIMER
            timer_off("sprq");
#endif
            break;
        }

    }

    buf4_close(&OutBuf);

#ifdef DPD_TIMER
    timer_off("buf4_sort");
#endif
    return 0;
}

}
