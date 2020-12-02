/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2019 The Psi4 Developers.
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
    \ingroup CCRESPONSE
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libciomr/libciomr.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccresponse {

void denom1(dpdfile2 *X1, double omega);
void local_filter_T1(dpdfile2 *T1);

void Y1_homogenous_build(const char *pert, int irrep, double omega) {
    int Gim, Gi, Gm, Ga, Gam, nrows, ncols, A, a, am;
    int Gei, ei, e, i, Gef, Ge, Gf, E, I, af, fa, f;
    dpdfile2 F, z1, z2, GMI, GAE;
    dpdbuf4 Y2, W, WL, D, X2, lx_iajb;
    dpdfile2 Y1, Y1new, X1, mu1, L1, lt, lx, lx_AB, Y1inhomo;
    dpdbuf4 L2, mu2;
    char lbl[32];
    double Y1_norm, Y2_norm;
    double *X, *Y;

    // Homogenous terms
    sprintf(lbl, "New Y_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&Y1new, PSIF_CC_OEI, irrep, 0, 1, lbl);

    global_dpd_->file2_scm(&Y1new, 0); //Do I need this??

    sprintf(lbl, "Y_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&Y1, PSIF_CC_OEI, irrep, 0, 1, lbl);

    sprintf(lbl, "Inhomo Y_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&Y1inhomo, PSIF_CC_OEI, irrep, 0, 1, lbl);

    global_dpd_->file2_axpy(&Y1inhomo, &Y1new, 1, 0);
    global_dpd_->file2_close(&Y1new);

    sprintf(lbl, "New Y_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&Y1new, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->file2_axpy(&Y1, &Y1new, omega, 0);

    // Y1 RHS += Yie*Fea
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "FAE");
    global_dpd_->contract222(&Y1, &F, &Y1new, 0, 1, 1.0, 1.0);
    global_dpd_->file2_close(&F);

    // Y1 RHS += -Yma*Fim 
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 0, "FMI");
    global_dpd_->contract222(&F, &Y1, &Y1new, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&F);

    // Y1 RHS += Yme*Wieam  
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "2 W(ME,jb) + W(Me,Jb)");
    global_dpd_->contract422(&W, &Y1, &Y1new, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->file2_close(&Y1);

    // Y1 RHS += 1/2 Yimef*Wefam 
    // Y(i,a) += [ 2 Y(im,ef) - Y(im,fe) ] * W(am,ef) 
    // Note: W(am,ef) is really Wabei (ei,ab) //
    //sprintf(lbl, "Y_%s_(2IjAb-IjbA) (%5.3f)", pert, omega);
    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert, omega); 
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAbEi (Ei,Ab)");
    //global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, "2 YIjAb - YIjBa");
        //    dpd_contract442(&Y2, &W, &newYIA, 0, 0, 1.0, 1.0);
        global_dpd_->file2_mat_init(&Y1new);
        global_dpd_->file2_mat_rd(&Y1new);
        for (Gam = 0; Gam < moinfo.nirreps; Gam++) {
            Gef = Gam; // W is totally symmetric 
            Gim = Gef ^ irrep;
            global_dpd_->buf4_mat_irrep_init(&Y2, Gim);
            global_dpd_->buf4_mat_irrep_rd(&Y2, Gim);
            global_dpd_->buf4_mat_irrep_shift13(&Y2, Gim);

            for (Gi = 0; Gi < moinfo.nirreps; Gi++) {
                Ga = Gi ^ irrep;
                Gm = Ga ^ Gam;
                W.matrix[Gam] = global_dpd_->dpd_block_matrix(moinfo.occpi[Gm], W.params->coltot[Gam]);

                nrows = moinfo.occpi[Gi];
                ncols = moinfo.occpi[Gm] * W.params->coltot[Gam];

                for (A = 0; A < moinfo.virtpi[Ga]; A++) {
                    a = moinfo.vir_off[Ga] + A; 
                    am = W.row_offset[Gam][a];
          
                    global_dpd_->buf4_mat_irrep_rd_block(&W, Gam, am, moinfo.occpi[Gm]);
        
                    if (nrows && ncols && moinfo.virtpi[Ga])
                        C_DGEMV('n', nrows, ncols, 1, Y2.shift.matrix[Gim][Gi][0], ncols, W.matrix[Gam][0], 1, 1,
                                &(Y1new.matrix[Gi][0][A]), moinfo.virtpi[Ga]);
                }
    
                global_dpd_->free_dpd_block(W.matrix[Gam], moinfo.occpi[Gm], W.params->coltot[Gam]);
            }
            global_dpd_->buf4_mat_irrep_close(&Y2, Gim);
        }
        global_dpd_->file2_mat_wrt(&Y1new);
        global_dpd_->file2_mat_close(&Y1new);
        global_dpd_->buf4_close(&Y2);
        global_dpd_->buf4_close(&W);

    // Y1 RHS += -1/2 Ymnae*Wiemn 
    //sprintf(lbl, "Y_%s_(2IjAb-IjbA) (%5.3f)", pert, omega);
    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 0, 10, 0, 0, "WMbIj");
    global_dpd_->contract442(&W, &Y2, &Y1new, 0, 2, -1.0, 1.0);
    global_dpd_->buf4_close(&Y2);
    global_dpd_->buf4_close(&W);

    //Y1 RHS += -Gef*Weifa //
    //sprintf(lbl, "G_%s_AE (%5.3f)", pert, omega);
    //global_dpd_->file2_init(&GAE, PSIF_CC_OEI, irrep, 1, 1, lbl);
    //global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)"); 
    //global_dpd_->dot13(&GAE,&W,&Y1new, 0, 0, -1.0, 1.0);
    //global_dpd_->buf4_close(&W); 
     //global_dpd_->file2_close(&GAE); 
    // Above code replaced to remove disk-space and memory bottlenecks 7/26/05, -TDC /
    sprintf(lbl, "G_%s_AE (%5.3f)", pert, omega);
    global_dpd_->file2_init(&GAE, PSIF_CC_OEI, irrep, 1, 1, lbl);
    global_dpd_->file2_mat_init(&GAE);
    global_dpd_->file2_mat_rd(&GAE);
    global_dpd_->file2_mat_init(&Y1new);
    global_dpd_->file2_mat_rd(&Y1new);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");

    for (Gei = 0; Gei < moinfo.nirreps; Gei++) {
        global_dpd_->buf4_mat_irrep_row_init(&W, Gei);
        X = init_array(W.params->coltot[Gei]);
        for (ei = 0; ei < W.params->rowtot[Gei]; ei++) {
            global_dpd_->buf4_mat_irrep_row_rd(&W, Gei, ei);
            e = W.params->roworb[Gei][ei][0];
            i = W.params->roworb[Gei][ei][1];
            Ge = W.params->psym[e];
            Gf = Ge ^ irrep;
            Gi = Ge ^ Gei;
            Ga = Gi ^ irrep;
            E = e - moinfo.vir_off[Ge];
            I = i - moinfo.occ_off[Gi];

            zero_arr(X, W.params->coltot[Gei]);

            for (fa = 0; fa < W.params->coltot[Gei]; fa++) {
                f = W.params->colorb[Gei][fa][0];
                a = W.params->colorb[Gei][fa][1];
                af = W.params->colidx[a][f];
                X[fa] = 2.0 * W.matrix[Gei][0][fa] - W.matrix[Gei][0][af];
            }

            nrows = moinfo.virtpi[Gf];
            ncols = moinfo.virtpi[Ga];
            if (nrows && ncols)
                  C_DGEMV('t', nrows, ncols, -1, &X[W.col_offset[Gei][Gf]], ncols, GAE.matrix[Ge][E], 1, 1,
                         Y1new.matrix[Gi][I], 1);
        }
        global_dpd_->buf4_mat_irrep_row_close(&W, Gei);
        free(X);
    }

    global_dpd_->buf4_close(&W);
    global_dpd_->file2_mat_wrt(&Y1new);
    global_dpd_->file2_mat_close(&Y1new);
    global_dpd_->file2_mat_close(&GAE);
    global_dpd_->file2_close(&GAE);

    // Y1 RHS += -Gmn*Wmina 
    sprintf(lbl, "G_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&GMI, PSIF_CC_OEI, irrep, 0, 0, lbl);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (Mn,eI)");
    global_dpd_->dot14(&GMI, &W, &Y1new, 0, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->file2_close(&GMI);

    if (params.local && local.filter_singles)
        local_filter_T1(&Y1new);
    else
        denom1(&Y1new, omega);

    global_dpd_->file2_close(&Y1new);

    return;
}
}  // namespace ccresponse
}  // namespace psi
