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
    \ingroup ccresponse
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

void denom2(dpdbuf4 *Y2, double omega);
void local_filter_T2(dpdbuf4 *T2);

void Y2_homogenous_build(const char *pert, int irrep, double omega) {
    dpdfile2 Y1, z, F, GAE, GMI, t1;
    dpdbuf4 Y2, Y2new, Y2inhomo, W, D, Z, Z1, Z2, B, I, T2, Z_final;
    char lbl[32];
    int Gej, Gab, Gij, Ge, Gj, Gi, Ga, i, j, ij, ab, nrows, length, E, e, II;
    int Gbm, Gfe, bm, a, b, m, Gb, Gm, Gf, M, fe, f, ef, ncols;
    double *Y;
    dpdbuf4 S, A, B_s;
    int Gc, C, c, cc;
    int rows_per_bucket, nbuckets, row_start, rows_left, nlinks;
    psio_address next;
    double **Y_diag, **B_diag;
    double Y2_norm;

    //Set of homogenous terms 

    //a factor of 0.5 because teh amplitudes are symmetric
   
    sprintf(lbl, "New Y_%s_IjAb (%5.3f)", pert, omega); 
    global_dpd_->buf4_init(&Y2new, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl); 
    global_dpd_->buf4_scm(&Y2new, 0); //Do I need this??  

    //Add inhomogenous terms
    sprintf(lbl, "Inhomo Y_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&Y2inhomo, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl); 
    global_dpd_->buf4_axpy(&Y2inhomo, &Y2new, 1);
    global_dpd_->buf4_close(&Y2inhomo);

    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_axpy(&Y2, &Y2new, 0.5*omega);    //Make sure about 0.5

    sprintf(lbl, "YF_%s_ijab (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&Z, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);

    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "FME"); 
    sprintf(lbl, "Y_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&Y1, PSIF_CC_OEI, irrep, 0, 1, lbl);

    global_dpd_->file2_mat_init(&F);
    global_dpd_->file2_mat_rd(&F);
    global_dpd_->file2_mat_init(&Y1);
    global_dpd_->file2_mat_rd(&Y1);

    for (Gej = 0; Gej < moinfo.nirreps; Gej++) {
        Gab = Gej;  //Z is totally symmetric 
        Gij = Gab ^ irrep;
        global_dpd_->buf4_mat_irrep_init(&Z, Gij);
        global_dpd_->buf4_mat_irrep_shift13(&Z, Gij);
       for(Gj = 0; Gj < moinfo.nirreps; Gj++) { // irreps of A
           Ga = Gj ^ irrep;
           Gi = Gij ^ Gj;
           Gb = Gab ^ Ga;
           for(ij = 0; ij < Z.params->rowtot[Gij]; ij++) {
               i = Z.params->roworb[Gej][ij][0];
               j = Z.params->roworb[Gej][ij][1];
               Gj = Ge ^ Gej;
               Gi = Gj ^ Gij;
               for(ab = 0; ab < Z.params->coltot[Gij]; ab++) {
                   a = Z.params->colorb[Gab][ab][0];
                   b = Z.params->colorb[Gab][ab][1];
                   //Z.matrix[Gij][ij][ab]  = 2*Y1.matrix[Gi][i][a] * F.matrix[Gj][j][b]; 
                   //Z.matrix[Gij][ij][ab] -=  Y1.matrix[Gj][j][a] * F.matrix[Gi][i][b];
                   Z.matrix[Gij][ij][ab]  = Y1.matrix[Gi][i][a] * F.matrix[Gj][j][b];
               }
           }
       }
        global_dpd_->buf4_mat_irrep_wrt(&Z, Gij);
        global_dpd_->buf4_mat_irrep_close(&Z, Gij);
    }

    global_dpd_->file2_mat_close(&F);
    global_dpd_->file2_close(&F);
    global_dpd_->file2_mat_close(&Y1);
    global_dpd_->file2_close(&Y1);

    global_dpd_->buf4_axpy(&Z, &Y2new, 2);

    //sprintf(lbl, "New Y_%s_IjAb (%5.3f)", pert, omega);
    //global_dpd_->buf4_sort_axpy(&Y2new, PSIF_CC_LR, pqsr, 0, 5, lbl, -1);
    sprintf(lbl, "YF_%s_ijba (%5.3f)", pert, omega);
    global_dpd_->buf4_sort(&Z, PSIF_CC_LR, pqsr, 0, 5, lbl);
    global_dpd_->buf4_close(&Z);
 
    sprintf(lbl, "YF_%s_ijba (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&Z, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_axpy(&Z, &Y2new, -1);
    global_dpd_->buf4_close(&Z);

    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 1, 1, "FAE");
    global_dpd_->contract244(&F, &Y2, &Y2new, 0, 2, 1, 1.0, 1.0);
    global_dpd_->file2_close(&F); 

    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 0, "FMI");
    global_dpd_->contract244(&F, &Y2, &Y2new, 1, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&F);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    global_dpd_->contract444(&W, &Y2, &Y2new, 0, 1, 0.5, 1);
    global_dpd_->buf4_close(&W);     
    global_dpd_->buf4_close(&Y2);	
 

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)"); //Compute it out of core
    sprintf(lbl, "Y_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&Y1, PSIF_CC_OEI, irrep, 0, 1, lbl); 
    global_dpd_->contract244(&Y1, &W, &Y2new, 1, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&W);

    
    sprintf(lbl, "WMnIeY1_%s_ijab (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, lbl);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (nM,eI)"); 
    global_dpd_->contract424(&W, &Y1, &Z, 3, 0, 0, 1.0, 0);
    global_dpd_->buf4_axpy(&Z, &Y2new, -1);
    //global_dpd_->contract424(&W, &Y1, &Y2new, 3, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&Y1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Z);


    // r_y2 += ndot('ieam,mjeb->ijab', self.Hovvo, self.y2, prefactor=2.0)
    //sort
    //global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 11, 10, 11, 0, "WMbEj");
    //global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, prqs, 10, 11, "WMEbj");
    //global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "Z (ia|jb)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 11, 10, 11, 0, "WMbEj");
    sprintf(lbl, "Y_%s_IAjb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&W, &Y2, &Z, 0, 1, 2.0, 0);

    global_dpd_->buf4_sort(&Z, PSIF_CC_HBAR, prqs, 0, 5, "Z (ij|ab)");

    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y2);

    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "Z (ij|ab)");

    global_dpd_->buf4_axpy(&Z, &Y2new, 1.0);
    global_dpd_->buf4_close(&Z);


    //r_y2 += ndot('iema,mjeb->ijab', self.Hovov, self.y2, prefactor=-1.0)

    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "Z (ia|jb)"); 
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
    sprintf(lbl, "Y_%s_IAjb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&W, &Y2, &Z, 0, 1, 1.0, 0);

    global_dpd_->buf4_sort(&Z, PSIF_CC_HBAR, prqs, 0, 5, "Z (ij|ab)");
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y2);

    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "Z (ij|ab)");
    global_dpd_->buf4_axpy(&Z, &Y2new, 1.0);
    global_dpd_->buf4_close(&Z);

 
    //r_y2 -= ndot('mibe,jema->ijab', self.y2, self.Hovov)

    /*
    //sort
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, psrq, 10, 10, "WMebJ");
    global_dpd_->buf4_close(&W);
    */

    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "Z (ia|jb)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
    sprintf(lbl, "Y_%s_IbjA (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&Y2, &W, &Z, 1, 0, 1.0, 0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y2);     
        
    global_dpd_->buf4_sort(&Z, PSIF_CC_HBAR, prsq, 0, 5, "Z (ij|ab)");
    global_dpd_->buf4_close(&Z);
    
    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "Z (ij|ab)");
    global_dpd_->buf4_axpy(&Z, &Y2new, 1.0);
    global_dpd_->buf4_close(&Z);


    // r_y2 -= ndot('mieb,jeam->ijab', self.y2, self.Hovvo)

    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "Z (ia|jb)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 11, 10, 11, 0, "WMbEj");   
    sprintf(lbl, "Y_%s_IAjb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&Y2, &W, &Z, 1, 0, 1.0, 0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y2);

    global_dpd_->buf4_sort(&Z, PSIF_CC_HBAR, prsq, 0, 5, "Z (ij|ab)");
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "Z (ij|ab)");

    global_dpd_->buf4_axpy(&Z, &Y2new, -1.0);
    global_dpd_->buf4_close(&Z);

//------------------------------------------------------------------
    sprintf(lbl, "G_%s_AE (%5.3f)", pert, omega);
    global_dpd_->file2_init(&GAE, PSIF_CC_OEI, irrep, 1, 1, lbl); 
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->contract244(&GAE, &D, &Y2new, 1, 2, 1, 1.0, 1.0);
    global_dpd_->file2_close(&GAE);


    sprintf(lbl, "G_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&GMI, PSIF_CC_OEI, irrep, 0, 0, lbl);
    global_dpd_->contract244(&GMI, &D, &Y2new, 0, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&GMI);
    global_dpd_->buf4_close(&D);    


    //global_dpd_->buf4_close(&Y2new);

//##########################TEST################################################

    //r_y2 += ndot('ijef,efab->ijab', self.y2, self.Hvvvv, prefactor=0.5)

    sprintf(lbl, "HvvvvY2 (ij,ab) (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&Z_final, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);

    global_dpd_->buf4_scm(&Z_final, 0);

    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);

/*
    sprintf(lbl, "Z(Ab,Ij) %s", pert);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 5, 0, 5, 0, 0, lbl);
    global_dpd_->buf4_init(&I, PSIF_CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
    global_dpd_->contract444(&I, &Y2, &Z, 0, 0, 1, 0);
    global_dpd_->buf4_close(&I);


    global_dpd_->buf4_close(&Z_final); // Need to close X2new to avoid collisions /
    //sprintf(lbl, "New Y_%s_IjAb (%5.3f)", pert, omega);
    sprintf(lbl, "HvvvvY2 (ij,ab) (%5.3f)", pert, omega);  
    global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_LR, rspq, 0, 5, lbl, 1);
    global_dpd_->buf4_init(&Z_final, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl); // re-open X2new here /
    global_dpd_->buf4_close(&Z);
*/

    if (params.abcd == "OLD") {
        sprintf(lbl, "Z(Ab,Ij) %s", pert);
        global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 5, 0, 5, 0, 0, lbl);
        global_dpd_->buf4_init(&I, PSIF_CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
        global_dpd_->contract444(&I, &Y2, &Z, 0, 0, 1, 0);
        global_dpd_->buf4_close(&I);

        global_dpd_->buf4_close(&Z_final);
        //global_dpd_->buf4_close(&X2new); /* Need to close X2new to avoid collisions */
        //sprintf(lbl, "New X_%s_IjAb (%5.3f)", pert, omega);
        sprintf(lbl, "HvvvvY2 (ij,ab) (%5.3f)", pert, omega);
        global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_LR, rspq, 0, 5, lbl, 1);
        global_dpd_->buf4_init(&Z_final, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl); /* re-open X2new here */
        global_dpd_->buf4_close(&Z);


    } else if (params.abcd == "NEW") {
        timer_on("ABCD:new");

        global_dpd_->buf4_close(&Y2);

        timer_on("ABCD:S");
        sprintf(lbl, "Y_%s_(+)(ij,ab) (%5.3f)", pert, omega);
        global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep, 3, 8, 3, 8, 0, lbl);
        global_dpd_->buf4_init(&I, PSIF_CC_BINTS, 0, 8, 8, 8, 8, 0, "B(+) <ab|cd> + <ab|dc>");
        sprintf(lbl, "S_%s_(ab,ij)", pert);
        global_dpd_->buf4_init(&S, PSIF_CC_TMP0, irrep, 8, 3, 8, 3, 0, lbl);
        global_dpd_->contract444(&I, &Y2, &S, 0, 0, 0.5, 0);
        global_dpd_->buf4_close(&S);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_close(&Y2);
        timer_off("ABCD:S");

        // Y_diag(ij,c)  = 2 * Y(ij,cc)
        // NB: Gcc = 0 and B is totally symmetry, so Gab = 0 
        // But Gij = irrep ^ Gab = irrep 
        sprintf(lbl, "Y_%s_(+)(ij,ab) (%5.3f)", pert, omega);
        global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep, 3, 8, 3, 8, 0, lbl);
        global_dpd_->buf4_mat_irrep_init(&Y2, irrep);
        global_dpd_->buf4_mat_irrep_rd(&Y2, irrep);
        Y_diag = global_dpd_->dpd_block_matrix(Y2.params->rowtot[irrep], moinfo.nvirt);
        for (ij = 0; ij < Y2.params->rowtot[irrep]; ij++)
            for (Gc = 0; Gc < moinfo.nirreps; Gc++)
                for (C = 0; C < moinfo.virtpi[Gc]; C++) {
                    c = C + moinfo.vir_off[Gc];
                    cc = Y2.params->colidx[c][c];
                    Y_diag[ij][c] = Y2.matrix[irrep][ij][cc];
                }
        global_dpd_->buf4_mat_irrep_close(&Y2, irrep);

        global_dpd_->buf4_init(&B_s, PSIF_CC_BINTS, 0, 8, 8, 8, 8, 0, "B(+) <ab|cd> + <ab|dc>");
        sprintf(lbl, "S_%s_(ab,ij)", pert);
        global_dpd_->buf4_init(&S, PSIF_CC_TMP0, irrep, 8, 3, 8, 3, 0, lbl);
        global_dpd_->buf4_mat_irrep_init(&S, 0);
        global_dpd_->buf4_mat_irrep_rd(&S, 0);

        rows_per_bucket = dpd_memfree() / (B_s.params->coltot[0] + moinfo.nvirt);
        if (rows_per_bucket > B_s.params->rowtot[0]) rows_per_bucket = B_s.params->rowtot[0];
        nbuckets = (int)ceil((double)B_s.params->rowtot[0] / (double)rows_per_bucket);
        rows_left = B_s.params->rowtot[0] % rows_per_bucket;

        B_diag = global_dpd_->dpd_block_matrix(rows_per_bucket, moinfo.nvirt);
        next = PSIO_ZERO;
        ncols = Y2.params->rowtot[irrep];
        nlinks = moinfo.nvirt;
        for (m = 0; m < (rows_left ? nbuckets - 1 : nbuckets); m++) {
            row_start = m * rows_per_bucket;
            nrows = rows_per_bucket;
            if (nrows && ncols && nlinks) {
                psio_read(PSIF_CC_BINTS, "B(+) <ab|cc>", (char *)B_diag[0], sizeof(double) * nrows * nlinks, next,
                          &next);
                C_DGEMM('n', 't', nrows, ncols, nlinks, -0.25, B_diag[0], nlinks, Y_diag[0], nlinks, 1,
                        S.matrix[0][row_start], ncols);
            }
        }
        if (rows_left) {
            row_start = m * rows_per_bucket;
            nrows = rows_left;
            if (nrows && ncols && nlinks) {
                psio_read(PSIF_CC_BINTS, "B(+) <ab|cc>", (char *)B_diag[0], sizeof(double) * nrows * nlinks, next,
                          &next);
                C_DGEMM('n', 't', nrows, ncols, nlinks, -0.25, B_diag[0], nlinks, Y_diag[0], nlinks, 1,
                        S.matrix[0][row_start], ncols);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&S, 0);
        global_dpd_->buf4_mat_irrep_close(&S, 0);
        global_dpd_->buf4_close(&S);
        global_dpd_->buf4_close(&B_s);
        global_dpd_->free_dpd_block(B_diag, rows_per_bucket, moinfo.nvirt);
        global_dpd_->free_dpd_block(Y_diag, Y2.params->rowtot[irrep], moinfo.nvirt);
        global_dpd_->buf4_close(&Y2);

        timer_on("ABCD:A");
        sprintf(lbl, "Y_%s_(-)(ij,ab) (%5.3f)", pert, omega);
        global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep, 4, 9, 4, 9, 0, lbl);
        global_dpd_->buf4_init(&I, PSIF_CC_BINTS, 0, 9, 9, 9, 9, 0, "B(-) <ab|cd> - <ab|dc>");
        sprintf(lbl, "A_%s_(ab,ij)", pert);
        global_dpd_->buf4_init(&A, PSIF_CC_TMP0, irrep, 9, 4, 9, 4, 0, lbl);
        global_dpd_->contract444(&I, &Y2, &A, 0, 0, 0.5, 0);
        global_dpd_->buf4_close(&A);
        global_dpd_->buf4_close(&I);
        global_dpd_->buf4_close(&Y2);
        timer_off("ABCD:A");

        timer_on("ABCD:axpy");
        global_dpd_->buf4_close(&Z_final); // Need to close X2new to avoid collisions /
        sprintf(lbl, "S_%s_(ab,ij)", pert);
        global_dpd_->buf4_init(&S, PSIF_CC_TMP0, irrep, 5, 0, 8, 3, 0, lbl);
        //sprintf(lbl, "New Y_%s_IjAb (%5.3f)", pert, omega);
        sprintf(lbl, "HvvvvY2 (ij,ab) (%5.3f)", pert, omega);
        global_dpd_->buf4_sort_axpy(&S, PSIF_CC_LR, rspq, 0, 5, lbl, 1);
        global_dpd_->buf4_close(&S);
        sprintf(lbl, "A_%s_(ab,ij)", pert);
        global_dpd_->buf4_init(&A, PSIF_CC_TMP0, irrep, 5, 0, 9, 4, 0, lbl);
        //sprintf(lbl, "New Y_%s_IjAb (%5.3f)", pert, omega);
        sprintf(lbl, "HvvvvY2 (ij,ab) (%5.3f)", pert, omega);
        global_dpd_->buf4_sort_axpy(&A, PSIF_CC_LR, rspq, 0, 5, lbl, 1);
        global_dpd_->buf4_close(&A);
        global_dpd_->buf4_init(&Z_final, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl); // re-open X2new here /

        timer_off("ABCD:axpy");

        sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert, omega);
        global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);

        timer_off("ABCD:new");
    }


    sprintf(lbl, "Z(Ij,am) %s", pert);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 0, 11, 0, 11, 0, lbl);
    global_dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract424(&Y2, &t1, &Z, 3, 1, 0, 1, 0);
    global_dpd_->file2_close(&t1);


    sprintf(lbl, "Z(Ij,Ab) %s", pert);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&I, PSIF_CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
    global_dpd_->contract444(&Z, &I, &Z1, 0, 1, 1, 0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&Z);

    //global_dpd_->buf4_axpy(&Z1, &Y2new, -1);
    global_dpd_->buf4_axpy(&Z1, &Z_final, -1);

    //global_dpd_->buf4_close(&Y2new); // Need to close X2new to avoid collisions /
    global_dpd_->buf4_close(&Z_final);
    //sprintf(lbl, "New Y_%s_IjAb (%5.3f)", pert, omega);
    sprintf(lbl, "HvvvvY2 (ij,ab) (%5.3f)", pert, omega);
    global_dpd_->buf4_sort_axpy(&Z1, PSIF_CC_LR, qpsr, 0, 5, lbl, -1);
    //global_dpd_->buf4_init(&Y2new, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl); // re-open X2new here 
    global_dpd_->buf4_init(&Z_final, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_close(&Z1);


    sprintf(lbl, "Z(Ij,Mn) %s", pert);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 0, 0, 0, 0, 0, lbl);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    global_dpd_->contract444(&Y2, &T2, &Z, 0, 0, 1, 0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&I, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract444(&Z, &I, &Z_final, 0, 1, 1, 1);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_axpy(&Z_final, &Y2new, 0.5);
   
    global_dpd_->buf4_close(&Z_final);

//########################END TEST##############################################


    if (params.local)
        local_filter_T2(&Y2new);
    else
        denom2(&Y2new, omega);
    global_dpd_->buf4_close(&Y2new);

//global_dpd_->buf4_close(&Y2);
//global_dpd_->buf4_close(&Y2new);
}

}  // namespace ccresponse
}  // namespace psi
