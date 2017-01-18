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


#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "dcft.h"
#include "psi4/libqt/qt.h"
#include "psi4/libdpd/dpd.h"

namespace psi{ namespace dcft{

void
DCFTSolver::AO_contribute(dpdbuf4 *tau1_AO, dpdbuf4 *tau2_AO, int p, int q,
        int r, int s, double value, dpdfile2 *s1, dpdfile2 *s1b, dpdfile2 *s2)
{
    int Gp, Gq, Gr, Gs, Gpr, Grp, Gps, Gsp, Gsq, Gqs, Gqr, Grq, Gpq, Gqp, Grs, Gsr;
    int prel, qrel, rrel, srel;
    int pq, rs, pr, rp, ps, sp, qr, rq, qs, sq;

    Gp = tau1_AO->params->psym[p];
    Gq = tau1_AO->params->psym[q];
    Gr = tau1_AO->params->psym[r];
    Gs = tau1_AO->params->psym[s];

    if(s1){
        prel = p - s1->params->poff[Gp];
        qrel = q - s1->params->poff[Gq];
        rrel = r - s1->params->poff[Gr];
        srel = s - s1->params->poff[Gs];
    }

    Gpr = Grp = Gp^Gr;
    Gps = Gsp = Gp^Gs;
    Gqr = Grq = Gq^Gr;
    Gqs = Gsq = Gq^Gs;
    Gpq = Gqp = Gp^Gq;
    Grs = Gsr = Gr^Gs;

    pq = tau1_AO->params->rowidx[p][q];
    rs = tau1_AO->params->rowidx[r][s];

    pr = tau1_AO->params->rowidx[p][r];
    rp = tau1_AO->params->rowidx[r][p];
    ps = tau1_AO->params->rowidx[p][s];
    sp = tau1_AO->params->rowidx[s][p];
    qr = tau1_AO->params->rowidx[q][r];
    rq = tau1_AO->params->rowidx[r][q];
    qs = tau1_AO->params->rowidx[q][s];
    sq = tau1_AO->params->rowidx[s][q];

    /* ####(pq|rs)#### */
    if(tau1_AO->params->coltot[Gpr])
        C_DAXPY(tau1_AO->params->coltot[Gpr], value, tau1_AO->matrix[Gpr][qs], 1,
                tau2_AO->matrix[Gpr][pr], 1);
    if(s1 && Gp==Gq && Gr==Gs){
        s2->matrix[Gp][prel][qrel] += value * s1->matrix[Gr][rrel][srel];
        s2->matrix[Gp][prel][qrel] += value * s1b->matrix[Gr][rrel][srel];
    }
    if(s1 && Gp==Gs && Gr==Gq)
        s2->matrix[Gp][prel][srel] -= value * s1->matrix[Gr][rrel][qrel];


    if(p!=q && r!=s && pq != rs){

        /* ####(pq|sr)#### */
        if(tau1_AO->params->coltot[Gps])
            C_DAXPY(tau1_AO->params->coltot[Gps], value, tau1_AO->matrix[Gps][qr], 1,
                    tau2_AO->matrix[Gps][ps], 1);
        if(s1 && Gp == Gq && Gs == Gr){
            s2->matrix[Gp][prel][qrel] += value * s1->matrix[Gs][srel][rrel];
            s2->matrix[Gp][prel][qrel] += value * s1b->matrix[Gs][srel][rrel];
        }
        if(s1 && Gp==Gr && Gs==Gq)
            s2->matrix[Gp][prel][rrel] -= value * s1->matrix[Gs][srel][qrel];

        /* ####(qp|rs)#### */
        if(tau1_AO->params->coltot[Gqr])
            C_DAXPY(tau1_AO->params->coltot[Gqr], value, tau1_AO->matrix[Gqr][ps], 1,
                    tau2_AO->matrix[Gqr][qr], 1);
        if(s1 && Gq==Gp && Gr==Gs){
            s2->matrix[Gq][qrel][prel] += value * s1->matrix[Gr][rrel][srel];
            s2->matrix[Gq][qrel][prel] += value * s1b->matrix[Gr][rrel][srel];
        }
        if(s1 && Gq==Gs && Gr==Gp)
            s2->matrix[Gq][qrel][srel] -= value * s1->matrix[Gr][rrel][prel];

        /* ####(qp|sr)#### */
        if(tau1_AO->params->coltot[Gqs])
            C_DAXPY(tau1_AO->params->coltot[Gqs], value, tau1_AO->matrix[Gqs][pr], 1,
                    tau2_AO->matrix[Gqs][qs], 1);
        if(s1 && Gq==Gp && Gs==Gr){
            s2->matrix[Gq][qrel][prel] += value * s1->matrix[Gs][srel][rrel];
            s2->matrix[Gq][qrel][prel] += value * s1b->matrix[Gs][srel][rrel];
        }
        if(s1 && Gq==Gr && Gs==Gp)
            s2->matrix[Gq][qrel][rrel] -= value * s1->matrix[Gs][srel][prel];

        /* ####(rs|pq)#### */
        if(tau1_AO->params->coltot[Grp])
            C_DAXPY(tau1_AO->params->coltot[Grp], value, tau1_AO->matrix[Grp][sq], 1,
                    tau2_AO->matrix[Grp][rp], 1);
        if(s1 && Gr==Gs && Gp==Gq){
            s2->matrix[Gr][rrel][srel] += value * s1->matrix[Gp][prel][qrel];
            s2->matrix[Gr][rrel][srel] += value * s1b->matrix[Gp][prel][qrel];
        }
        if(s1 && Gr==Gq && Gp==Gs)
            s2->matrix[Gr][rrel][qrel] -= value * s1->matrix[Gp][prel][srel];

        /* ####(sr|pq)#### */
        if(tau1_AO->params->coltot[Gsp])
            C_DAXPY(tau1_AO->params->coltot[Gsp], value, tau1_AO->matrix[Gsp][rq], 1,
                    tau2_AO->matrix[Gsp][sp], 1);
        if(s1 && Gs==Gr && Gp==Gq){
            s2->matrix[Gs][srel][rrel] += value * s1->matrix[Gp][prel][qrel];
            s2->matrix[Gs][srel][rrel] += value * s1b->matrix[Gp][prel][qrel];
        }
        if(s1 && Gs==Gq && Gp==Gr)
            s2->matrix[Gs][srel][qrel] -= value * s1->matrix[Gp][prel][rrel];

        /* ####(rs|qp)#### */
        if(tau1_AO->params->coltot[Grq])
            C_DAXPY(tau1_AO->params->coltot[Grq], value, tau1_AO->matrix[Grq][sp], 1,
                    tau2_AO->matrix[Grq][rq], 1);
        if(s1 && Gr==Gs && Gq==Gp){
            s2->matrix[Gr][rrel][srel] += value * s1->matrix[Gq][qrel][prel];
            s2->matrix[Gr][rrel][srel] += value * s1b->matrix[Gq][qrel][prel];
        }
        if(s1 && Gr==Gp && Gq==Gs)
            s2->matrix[Gr][rrel][prel] -= value * s1->matrix[Gq][qrel][srel];

        /* ####(sr|qp)#### */
        if(tau1_AO->params->coltot[Gsq])
            C_DAXPY(tau1_AO->params->coltot[Gsq], value, tau1_AO->matrix[Gsq][rp], 1,
                    tau2_AO->matrix[Gsq][sq],1 );
        if(s1 && Gs==Gr && Gq==Gp){
            s2->matrix[Gs][srel][rrel] += value * s1->matrix[Gq][qrel][prel];
            s2->matrix[Gs][srel][rrel] += value * s1b->matrix[Gq][qrel][prel];
        }
        if(s1 && Gs==Gp && Gq==Gr)
            s2->matrix[Gs][srel][prel] -= value * s1->matrix[Gq][qrel][rrel];

    }
    else if(p!=q && r!=s && pq==rs) {

        /* (pq|sr) */
        if(tau1_AO->params->coltot[Gps])
            C_DAXPY(tau1_AO->params->coltot[Gps], value, tau1_AO->matrix[Gps][qr], 1,
                    tau2_AO->matrix[Gps][ps], 1);
        if(s1 && Gp==Gq && Gs==Gr){
            s2->matrix[Gp][prel][qrel] += value * s1->matrix[Gs][srel][rrel];
            s2->matrix[Gp][prel][qrel] += value * s1b->matrix[Gs][srel][rrel];
        }
        if(s1 && Gp==Gr && Gs==Gq)
            s2->matrix[Gp][prel][rrel] -= value * s1->matrix[Gs][srel][qrel];

        /* (qp|rs) */
        if(tau1_AO->params->coltot[Gqr])
            C_DAXPY(tau1_AO->params->coltot[Gqr], value, tau1_AO->matrix[Gqr][ps], 1,
                    tau2_AO->matrix[Gqr][qr], 1);
        if(s1 && Gq==Gp && Gr==Gs){
            s2->matrix[Gq][qrel][prel] += value * s1->matrix[Gr][rrel][srel];
            s2->matrix[Gq][qrel][prel] += value * s1b->matrix[Gr][rrel][srel];
        }
        if(s1 && Gq==Gs && Gr==Gp)
            s2->matrix[Gq][qrel][srel] -= value * s1->matrix[Gr][rrel][prel];

        /* (qp|sr) */
        if(tau1_AO->params->coltot[Gqs])
            C_DAXPY(tau1_AO->params->coltot[Gqs], value, tau1_AO->matrix[Gqs][pr], 1,
                    tau2_AO->matrix[Gqs][qs], 1);
        if(s1 && Gq==Gp && Gs==Gr){
            s2->matrix[Gq][qrel][prel] += value * s1->matrix[Gs][srel][rrel];
            s2->matrix[Gq][qrel][prel] += value * s1b->matrix[Gs][srel][rrel];
        }
        if(s1 && Gq==Gr && Gs==Gp)
            s2->matrix[Gq][qrel][rrel] -= value * s1->matrix[Gs][srel][prel];
    }
    else if(p!=q && r==s) {

        /* (qp|rs) */
        if(tau1_AO->params->coltot[Gqr])
            C_DAXPY(tau1_AO->params->coltot[Gqr], value, tau1_AO->matrix[Gqr][ps], 1,
                    tau2_AO->matrix[Gqr][qr], 1);
        if(s1 && Gq==Gp && Gr==Gs){
            s2->matrix[Gq][qrel][prel] += value * s1->matrix[Gr][rrel][srel];
            s2->matrix[Gq][qrel][prel] += value * s1b->matrix[Gr][rrel][srel];
        }
        if(s1 && Gq==Gs && Gr==Gp)
            s2->matrix[Gq][qrel][srel] -= value * s1->matrix[Gr][rrel][prel];

        /* (rs|pq) */
        if(tau1_AO->params->coltot[Grp])
            C_DAXPY(tau1_AO->params->coltot[Grp], value, tau1_AO->matrix[Grp][sq], 1,
                    tau2_AO->matrix[Grp][rp], 1);
        if(s1 && Gr==Gs && Gp==Gq){
            s2->matrix[Gr][rrel][srel] += value * s1->matrix[Gp][prel][qrel];
            s2->matrix[Gr][rrel][srel] += value * s1b->matrix[Gp][prel][qrel];
        }
        if(s1 && Gr==Gq && Gp==Gs)
            s2->matrix[Gr][rrel][qrel] -= value * s1->matrix[Gp][prel][srel];

        /* (rs|qp) */
        if(tau1_AO->params->coltot[Grq])
            C_DAXPY(tau1_AO->params->coltot[Grq], value, tau1_AO->matrix[Grq][sp], 1,
                    tau2_AO->matrix[Grq][rq], 1);
        if(s1 && Gr==Gs && Gq==Gp){
            s2->matrix[Gr][rrel][srel] += value * s1->matrix[Gq][qrel][prel];
            s2->matrix[Gr][rrel][srel] += value * s1b->matrix[Gq][qrel][prel];
        }
        if(s1 && Gr==Gp && Gq==Gs)
            s2->matrix[Gr][rrel][prel] -= value * s1->matrix[Gq][qrel][srel];
    }

    else if(p==q && r!=s) {

        /* (pq|sr) */
        if(tau1_AO->params->coltot[Gps])
            C_DAXPY(tau1_AO->params->coltot[Gps], value, tau1_AO->matrix[Gps][qr], 1,
                    tau2_AO->matrix[Gps][ps], 1);
        if(s1 && Gp==Gq && Gs==Gr){
            s2->matrix[Gp][prel][qrel] += value * s1->matrix[Gs][srel][rrel];
            s2->matrix[Gp][prel][qrel] += value * s1b->matrix[Gs][srel][rrel];
        }
        if(s1 && Gp==Gr && Gs==Gq)
            s2->matrix[Gp][prel][rrel] -= value * s1->matrix[Gs][srel][qrel];

        /* (rs|pq) */
        if(tau1_AO->params->coltot[Grp])
            C_DAXPY(tau1_AO->params->coltot[Grp], value, tau1_AO->matrix[Grp][sq], 1,
                    tau2_AO->matrix[Grp][rp], 1);
        if(s1 && Gr==Gs && Gp==Gq){
            s2->matrix[Gr][rrel][srel] += value * s1->matrix[Gp][prel][qrel];
            s2->matrix[Gr][rrel][srel] += value * s1b->matrix[Gp][prel][qrel];
        }
        if(s1 && Gr==Gq && Gp==Gs)
            s2->matrix[Gr][rrel][qrel] -= value * s1->matrix[Gp][prel][srel];

        /* (sr|pq) */
        if(tau1_AO->params->coltot[Gsp])
            C_DAXPY(tau1_AO->params->coltot[Gsp], value, tau1_AO->matrix[Gsp][rq], 1,
                    tau2_AO->matrix[Gsp][sp], 1);
        if(s1 && Gs==Gr && Gp==Gq){
            s2->matrix[Gs][srel][rrel] += value * s1->matrix[Gp][prel][qrel];
            s2->matrix[Gs][srel][rrel] += value * s1b->matrix[Gp][prel][qrel];
        }
        if(s1 && Gs==Gq && Gp==Gr)
            s2->matrix[Gs][srel][qrel] -= value * s1->matrix[Gp][prel][rrel];

    }

    else if(p==q && r==s && pq != rs) {

        /* (rs|pq) */
        if(tau1_AO->params->coltot[Grp])
            C_DAXPY(tau1_AO->params->coltot[Grp], value, tau1_AO->matrix[Grp][sq], 1,
                    tau2_AO->matrix[Grp][rp], 1);
        if(s1 && Gr==Gs && Gp==Gq){
            s2->matrix[Gr][rrel][srel] += value * s1->matrix[Gp][prel][qrel];
            s2->matrix[Gr][rrel][srel] += value * s1b->matrix[Gp][prel][qrel];
        }
        if(s1 && Gr==Gq && Gp==Gs)
            s2->matrix[Gr][rrel][qrel] -= value * s1->matrix[Gp][prel][srel];

    }
}


}}
