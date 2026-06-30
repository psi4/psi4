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
#include <cmath>
#include "psi4/libdpd/dpd.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccresponse {


double Y1HX1X1(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
		      const char *pert_z, int irrep_z, double omega_z) {

    double result = 0.0;
    dpdfile2 X1, Y1, F, z, z1, Z_final, t1;
    dpdbuf4 W, Z, Z2, Y2, Z1, T2, YF; 
    char lbl[32];
    int i, j, a, b, ab, ij;
    int Gej, Gab, Gij, Gi, Gj, Ga, Gb, Ge;
    double Y1_norm;
    
    // *** <O|Y1(A)[[Hbar(0),X1(B),X1(C)]]|0> ***

    sprintf(lbl, "Y_%s_IA (%5.3f)", pert_x, omega_x);
    global_dpd_->file2_init(&Y1, PSIF_CC_OEI, irrep_x, 0, 1, lbl);

/*
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, irrep_x, 0, 5, 0, 5, 0, "Z2 (ij|ab)");
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "FME");

    global_dpd_->file2_mat_init(&Y1);
    global_dpd_->file2_mat_rd(&Y1);
    global_dpd_->file2_mat_init(&F);
    global_dpd_->file2_mat_rd(&F);

    for (Gej = 0; Gej < moinfo.nirreps; Gej++) {
        Gab = Gej;  //Z2 is totally symmetric 
        Gij = Gab ^ irrep_x;
        global_dpd_->buf4_mat_irrep_init(&Z2, Gij);
        global_dpd_->buf4_mat_irrep_shift13(&Z2, Gij);
       for(Gj = 0; Gj < moinfo.nirreps; Gj++) { // irreps of A
           Ga = Gj ^ irrep_x;
           Gi = Gij ^ Gj;
           Gb = Gab ^ Ga;
           for(ij = 0; ij < Z2.params->rowtot[Gij]; ij++) {
               i = Z2.params->roworb[Gej][ij][0];
               j = Z2.params->roworb[Gej][ij][1];
               Gj = Ge ^ Gej;
               Gi = Gj ^ Gij;
               for(ab = 0; ab < Z2.params->coltot[Gij]; ab++) {
                   a = Z2.params->colorb[Gab][ab][0];
                   b = Z2.params->colorb[Gab][ab][1];
                   Z2.matrix[Gij][ij][ab] -= F.matrix[Gi][i][a] * Y1.matrix[Gj][j][b]; 
                   Z2.matrix[Gij][ij][ab] -= Y1.matrix[Gj][i][a] * F.matrix[Gi][j][b];
               }
           }
       }
        global_dpd_->buf4_mat_irrep_wrt(&Z2, Gij);
        global_dpd_->buf4_mat_irrep_close(&Z2, Gij);
    }

    global_dpd_->file2_mat_close(&Y1);
    global_dpd_->file2_mat_close(&F);
    global_dpd_->file2_close(&F);
*/

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z (ij,ab)");
    global_dpd_->buf4_scm(&Z, 0);

    sprintf(lbl, "YF_%s_ijab (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&YF, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_axpy(&YF, &Z, -1); 

    sprintf(lbl, "YF_%s_jiba (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_sort(&YF, PSIF_CC_LR, qpsr, 0, 5, lbl);
    global_dpd_->buf4_close(&YF);

    sprintf(lbl, "YF_%s_jiba (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&YF, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_axpy(&YF, &Z, -1);    
    global_dpd_->buf4_close(&YF);

    //global_dpd_->buf4_axpy(&Z2, &Z, 1);
    //global_dpd_->buf4_close(&Z2);

//    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (Mn,eI)");
//    global_dpd_->contract424(&W, &Y1, &Z, 3, 0, 0, -1, 1);
//    global_dpd_->buf4_close(&W);

   /*
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z2 (ij,ab)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (nM,eI)");
    global_dpd_->contract424(&W, &Y1, &Z2, 3, 0, 0, -1, 0);   
    global_dpd_->buf4_close(&W);
    */
    sprintf(lbl, "WMnIeY1_%s_ijab (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Z2, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, lbl);

    //global_dpd_->buf4_axpy(&Z2, &Z, -1);

    sprintf(lbl, "WMnIeY1_%s_ijba (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_sort(&Z2, PSIF_CC_HBAR, pqsr, 0, 5, lbl);   //sort
    sprintf(lbl, "WMnIeY1_%s_jiab (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_sort(&Z2, PSIF_CC_HBAR, qprs, 0, 5, lbl);
    global_dpd_->buf4_close(&Z2);

    sprintf(lbl, "WMnIeY1_%s_ijba (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Z2, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, lbl); 
    global_dpd_->buf4_axpy(&Z2, &Z, -1);  
    global_dpd_->buf4_close(&Z2);

    sprintf(lbl, "WMnIeY1_%s_jiab (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Z2, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, lbl); 
    global_dpd_->buf4_axpy(&Z2, &Z, -1);  
    global_dpd_->buf4_close(&Z2);

   

//Here Compute out of core!!!!!!!!

/*
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    global_dpd_->buf4_scmcopy(&W, PSIF_CC_HBAR, "WAmEf 2(Am,Ef) - (Am,fE)", 2);
    global_dpd_->buf4_sort_axpy(&W, PSIF_CC_HBAR, pqsr, 11, 5, "WAmEf 2(Am,Ef) - (Am,fE)", -1);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, qprs, 10, 5, "WAmEf 2(mA,Ef) - (mA,fE)"); //Compute this part out of core
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, pqsr, 11, 5, "WAmEf 2(Am,Ef) - (Am,fE) (am,fe)"); //Compute this part out of core
    global_dpd_->buf4_close(&W);
*/

    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z2 (ij,ab)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 5, 10, 5, 0, "WAmEf 2(mA,Ef) - (mA,fE)");
    global_dpd_->contract424(&W, &Y1, &Z2, 1, 1, 1, 1, 0);   

    global_dpd_->buf4_axpy(&Z2, &Z, 1);
    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&W); 

    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z2 (ij,ab)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE) (am,fe)");
    global_dpd_->contract244(&Y1, &W, &Z2, 1, 0, 0, 1, 0);
    global_dpd_->buf4_close(&W);
    global_dpd_->file2_close(&Y1);

    global_dpd_->buf4_axpy(&Z2, &Z, 1);
    global_dpd_->buf4_close(&Z2);


    //-------------------------------------------------------------//

    sprintf(lbl, "Z_%s_Final (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&Z_final, PSIF_CC_OEI, irrep_y, 0, 1, lbl);

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->dot14(&X1, &Z, &Z_final, 0, 0, 1, 0);

    global_dpd_->buf4_close(&Z);

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_z, omega_z);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_z, 0, 1, lbl);

    result = global_dpd_->file2_dot(&Z_final, &X1); 

    global_dpd_->file2_close(&X1);
    global_dpd_->file2_close(&Z_final);


    return result;
}


double Y2HX1X1(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
                      const char *pert_z, int irrep_z, double omega_z) {

    double result = 0.0;
    dpdfile2 X1, Y1, GAE, GMI, z_ij, z_ia, z2_ia,F, FX, z, z1, Z_final, t1;
    dpdbuf4 D, tIjAb, W, Z, Z2, Y2, Z1, T2, I, Zfinal ;
    char lbl[32];
    int i, j, a, b, ab, ij;
    int Gej, Gab, Gij, Gi, Gj, Ga, Gb, Ge;
    double Y1_norm;

    // *** <O|Y2(A)[[Hbar(0),X1(B),X1(C)]]|0> ***

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (ia,jb)");

    sprintf(lbl, "Y_%s_IbjA (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 10, 10, 10, 10, 0, lbl);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 10, 10, 10, 0, "WMbeJ");
    global_dpd_->contract444(&W, &Y2, &Z, 0, 1, -1.0, 0.0);

    global_dpd_->contract444(&Y2, &W, &Z, 1, 0, -1.0, 1.0);
    global_dpd_->buf4_close(&Y2);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 11, 10, 11, 0, "WMbEj");
    sprintf(lbl, "Y_%s_IAjb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&W, &Y2, &Z, 0, 1, 1.0, 1.0);

    sprintf(lbl, "Y_%s_IAjb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&Y2, &W, &Z, 1, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y2);

    //sort
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj");
    global_dpd_->buf4_sort(&W, PSIF_CC_HBAR, pqsr, 0, 0, "WMnIj (Mn,jI)");
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z2 (ij,ab)");

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj (Mn,jI)");
    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract444(&W, &Y2, &Z2, 0, 1, 1.0, 0.0); 
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Y2);

    global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP0, qspr, 10, 10, "Z2 (ia,jb)"); //Are the same?
    global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP0, prqs, 10, 10, "Z2 (jb,ia)"); //Are the same?
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z2 (ia,jb)");
    global_dpd_->buf4_axpy(&Z2, &Z, -0.5);
    global_dpd_->buf4_close(&Z2);


    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z2 (jb,ia)");
    global_dpd_->buf4_axpy(&Z2, &Z, -0.5);
    global_dpd_->buf4_close(&Z2);


    //tmp  += 0.5*np.einsum('fabc,jkfa->jckb',self.Hvvvv,self.y2_A)
    //tmp  += 0.5*np.einsum('facb,kjfa->jckb',self.Hvvvv,self.y2_A)

    sprintf(lbl, "HvvvvY2 (ij,ab) (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Zfinal, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);

    global_dpd_->buf4_scm(&Zfinal, 0);

    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);

    sprintf(lbl, "Z(Ab,Ij) %s", pert_x);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, irrep_x, 5, 0, 5, 0, 0, lbl);
    global_dpd_->buf4_init(&I, PSIF_CC_BINTS, 0, 5, 5, 5, 5, 0, "B <ab|cd>");
    global_dpd_->contract444(&I, &Y2, &Z1, 0, 0, 1, 0);
    global_dpd_->buf4_close(&I);


    global_dpd_->buf4_close(&Zfinal); // Need to close X2new to avoid collisions /
    sprintf(lbl, "HvvvvY2 (ij,ab) (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_sort_axpy(&Z1, PSIF_CC_LR, rspq, 0, 5, lbl, 1);
    global_dpd_->buf4_init(&Zfinal, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl); // re-open X2new here /
    global_dpd_->buf4_close(&Z1);

    sprintf(lbl, "Z(Ij,am) %s", pert_x);
    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, irrep_x, 0, 11, 0, 11, 0, lbl);
    global_dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "tIA");
    global_dpd_->contract424(&Y2, &t1, &Z1, 3, 1, 0, 1, 0);
    global_dpd_->file2_close(&t1);

    sprintf(lbl, "Z(Ij,Ab) %s", pert_x);
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, irrep_x, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&I, PSIF_CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
    global_dpd_->contract444(&Z1, &I, &Z2, 0, 1, 1, 0);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&Z1);

    global_dpd_->buf4_axpy(&Z2, &Zfinal, -1);


    global_dpd_->buf4_close(&Zfinal); // Need to close X2new to avoid collisions //
    sprintf(lbl, "HvvvvY2 (ij,ab) (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_sort_axpy(&Z2, PSIF_CC_LR, qpsr, 0, 5, lbl, -1);
    global_dpd_->buf4_init(&Zfinal, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl); //// re-open X2new here
    global_dpd_->buf4_close(&Z2);


    sprintf(lbl, "Z(Ij,Mn) %s", pert_x);
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, irrep_x, 0, 0, 0, 0, 0, lbl);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");
    global_dpd_->contract444(&Y2, &T2, &Z2, 0, 0, 1, 0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_init(&I, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract444(&Z2, &I, &Zfinal, 0, 1, 1, 1);
    global_dpd_->buf4_close(&I);
    global_dpd_->buf4_close(&Z2);


    sprintf(lbl, "HvvvvY2 (ia,jb) (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_sort(&Zfinal, PSIF_CC_LR, psqr, 10, 10, lbl);
    sprintf(lbl, "HvvvvY2 (ja,ib) (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_sort(&Zfinal, PSIF_CC_LR, qrps, 10, 10, lbl);
    global_dpd_->buf4_close(&Zfinal);

    sprintf(lbl, "HvvvvY2 (ia,jb) (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Zfinal, PSIF_CC_LR, 0, 10, 10, 10, 10, 0, lbl);
    global_dpd_->buf4_axpy(&Zfinal, &Z, -0.5); 
    global_dpd_->buf4_close(&Zfinal);     


    sprintf(lbl, "HvvvvY2 (ja,ib) (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Zfinal, PSIF_CC_LR, 0, 10, 10, 10, 10, 0, lbl);
    global_dpd_->buf4_axpy(&Zfinal, &Z, -0.5); 
    global_dpd_->buf4_close(&Zfinal);


    global_dpd_->file2_init(&z_ia, PSIF_CC_OEI, 0, 0, 1, "z2_ia");

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->dot14(&X1, &Z, &z_ia, 0, 1, 1, 0);
    global_dpd_->buf4_close(&Z);
    global_dpd_->file2_close(&X1);


    sprintf(lbl, "X_%s_IA (%5.3f)", pert_z, omega_z);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_z, 0, 1, lbl);
    
    result -= global_dpd_->file2_dot(&z_ia, &X1);
    
    global_dpd_->file2_close(&z_ia);
    global_dpd_->file2_close(&X1);



    global_dpd_->file2_init(&z_ia, PSIF_CC_OEI, 0, 0, 1, "z_ia");

    sprintf(lbl, "G_%s_AE (%5.3f)", pert_x, omega_x);
    global_dpd_->file2_init(&GAE, PSIF_CC_OEI, irrep_x, 1, 1, lbl);
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_z, omega_z);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_z, 0, 1, lbl);
    global_dpd_->contract222(&X1, &GAE, &z_ia, 0, 1, -1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->file2_close(&GAE);

    global_dpd_->file2_init(&z2_ia, PSIF_CC_OEI, 0, 0, 1, "z2_ia");

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->dot24(&z_ia, &D, &z2_ia, 0, 0, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->file2_close(&z_ia);

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);

    result -= global_dpd_->file2_dot(&z2_ia, &X1);

    global_dpd_->file2_close(&X1);
    global_dpd_->file2_close(&z2_ia);


    // tmp = np.einsum('ijab,jlba->il',self.t2,self.y2_A)
    // tmp2 = np.einsum('kc,kicd->id',self.x1_B,self.Loovv)
    // tmp2 = np.einsum('id,ld->il',tmp2,self.x1_C)
    // self.Bcon1 -= ndot('il,il->',tmp2,tmp)


    global_dpd_->file2_init(&z_ia, PSIF_CC_OEI, 0, 0, 1, "z_ia");

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->dot13(&X1, &D, &z_ia, 0, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&D);


    global_dpd_->file2_init(&z_ij, PSIF_CC_OEI, 0, 0, 0, "z_ij");

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_z, omega_z);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_z, 0, 1, lbl);    
    global_dpd_->contract222(&z_ia, &X1, &z_ij, 0, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->file2_close(&z_ia);


   // Intermediate????
    sprintf(lbl, "G_%s_MI (%5.3f) test", pert_x, omega_x);
    global_dpd_->file2_init(&GMI, PSIF_CC_OEI, irrep_x, 0, 0, lbl);

    // Y(Mj,Ab) * [ 2 Y(Ij,Ab) - Y(Ij,Ba) ] --> G(M,I) //
    global_dpd_->buf4_init(&tIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    //sprintf(lbl, "Y_%s_(2IjAb-IjbA) (%5.3f)", pert, omega);
    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract442(&tIjAb, &Y2, &GMI, 0, 0, 1, 0);
 
    //sprintf(lbl, "G_%s_IA (%5.3f)", pert_x, omega_x);
    //global_dpd_->file2_init(&GMI, PSIF_CC_OEI, irrep_x, 0, 0, lbl);

    global_dpd_->buf4_close(&Y2);
    global_dpd_->buf4_close(&tIjAb);

    result -= global_dpd_->file2_dot(&GMI, &z_ij);

    global_dpd_->file2_close(&GMI);
    global_dpd_->file2_close(&z_ij);


    // tmp = np.einsum('ijab,jkba->ik',self.t2,self.y2_A)
    // tmp2 = np.einsum('ld,lidc->ic',self.x1_C,self.Loovv)
    // tmp2 = np.einsum('ic,kc->ik',tmp2,self.x1_B)
    // self.Bcon1 -= ndot('ik,ik->',tmp2,tmp)

    global_dpd_->file2_init(&z_ia, PSIF_CC_OEI, 0, 0, 1, "z_ia");

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_z, omega_z);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_z, 0, 1, lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->dot13(&X1, &D, &z_ia, 0, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&D);


    global_dpd_->file2_init(&z_ij, PSIF_CC_OEI, 0, 0, 0, "z_ij");

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->contract222(&z_ia, &X1, &z_ij, 0, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->file2_close(&z_ia);


   // Intermediate????
    sprintf(lbl, "G_%s_MI (%5.3f) test", pert_x, omega_x);
    global_dpd_->file2_init(&GMI, PSIF_CC_OEI, irrep_x, 0, 0, lbl);

    result -= global_dpd_->file2_dot(&GMI, &z_ij);

    global_dpd_->file2_close(&GMI);
    global_dpd_->file2_close(&z_ij);


   
    // tmp = np.einsum('ijab,ijcb->ac',self.t2,self.y2_A)
    // tmp = np.einsum('kc,ac->ka',self.x1_B,tmp)
    // tmp2 = np.einsum('ld,lkda->ka',self.x1_C,self.Loovv)
    // self.Bcon1 -= ndot('ka,ka->',tmp2,tmp)

    global_dpd_->file2_init(&z_ia, PSIF_CC_OEI, 0, 0, 1, "z_ia");

    sprintf(lbl, "G_%s_AE (%5.3f)", pert_x, omega_x);
    global_dpd_->file2_init(&GAE, PSIF_CC_OEI, irrep_x, 1, 1, lbl);


    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->contract222(&X1, &GAE, &z_ia, 0, 1, -1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->file2_close(&GAE);


    global_dpd_->file2_init(&z2_ia, PSIF_CC_OEI, 0, 0, 1, "z2_ia");

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_z, omega_z);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_z, 0, 1, lbl);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->dot24(&X1, &D, &z2_ia, 0, 0, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->file2_close(&X1);

    result -= global_dpd_->file2_dot(&z_ia, &z2_ia);

    global_dpd_->file2_close(&z_ia);
    global_dpd_->file2_close(&z2_ia);


    return result;
}



double Y2HX2X2(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
                      const char *pert_z, int irrep_z, double omega_z) {

    double result = 0.0;
    dpdfile2 X1, Y1, z_ij, z_ab, F, FX, z, z1; 
    dpdbuf4 X2,Z, D, Y2, Z_final ;
    char lbl[32];
    int i, j, a, b, ab, ij;
    int Gej, Gab, Gij, Gi, Gj, Ga, Gb, Ge;
    double Y1_norm;

    // *** <O|Y2(A)[[Hbar(0),X1(B),X1(C)]]|0> ***

    //tmp = np.einsum("klcd,ijcd->ijkl",self.x2_C,self.y2_A)    
    //tmp = np.einsum("ijkl,ijab->klab",tmp,self.x2_B)
    //self.Bcon1 += 0.5*np.einsum('klab,klab->',tmp,self.Goovv)

    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z (ij,ab)");
    global_dpd_->buf4_scm(&Z_final, 0);

    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 0, 0, 0, 0, "Z_ijkl");    
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract444(&Y2, &X2, &Z, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y2);
    global_dpd_->buf4_close(&X2);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");  
    global_dpd_->contract444(&Z, &D, &Z_final, 0, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Z);


    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);

    result = 0.5 * global_dpd_->buf4_dot(&Z_final, &X2);

    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&Z_final);
    

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 0, 0, 0, 0, "Z (ij,kl)");
    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract444(&X2, &Y2, &Z, 0, 0, 1, 0);
    global_dpd_->buf4_close(&Y2);
    global_dpd_->buf4_close(&X2);

    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z (ij,ab)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->contract444(&Z, &D, &Z_final, 1, 1, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Z);


    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 0, 5, 0, 5, 0, lbl);

    result += 0.5 * global_dpd_->buf4_dot(&Z_final, &X2);

    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&Z_final);


    //tmp = np.einsum("ijab,ikbd->jkad",self.x2_B,self.y2_A)    
    //tmp = np.einsum("jkad,klcd->jlac",tmp,self.x2_C)
    //self.Bcon1 += np.einsum('jlac,jlac->',tmp,self.Goovv) 

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (ia,jb)");
    sprintf(lbl, "Y_%s_IAjb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 10, 10, 10, 10, 0, lbl);
    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_y, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&X2, &Y2, &Z, 1, 1, 1, 0);
    global_dpd_->buf4_close(&Y2);
    global_dpd_->buf4_close(&X2);

    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (ia,jb)");
    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&Z, &X2, &Z_final, 0, 1, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_sort(&Z_final, PSIF_CC_TMP0, prqs, 0, 5, "Z (ij,ab)");   
    global_dpd_->buf4_close(&Z_final); 


    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z (ij,ab)");    
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");

    result += global_dpd_->buf4_dot(&Z_final, &D);

    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Z_final);


    //tmp = np.einsum("klcd,ikdb->licb",self.x2_C,self.y2_A)
    //tmp = np.einsum("licb,ijab->ljca",tmp,self.x2_B)
    //self.Bcon1 += np.einsum('ljca,ljac->',tmp,self.Goovv) 

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (ia,jb)");
    sprintf(lbl, "Y_%s_IbjA (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 10, 10, 10, 10, 0, lbl);
    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&X2, &Y2, &Z, 1, 0, 1, 0);
    global_dpd_->buf4_close(&Y2);
    global_dpd_->buf4_close(&X2);

    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z_final (ia,jb)");
    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_y, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&Z, &X2, &Z_final, 0, 1, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_sort(&Z_final, PSIF_CC_TMP0, prsq, 0, 5, "Z_final (ij,ab)");
    global_dpd_->buf4_close(&Z_final);

    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z_final (ij,ab)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");

    result += global_dpd_->buf4_dot(&Z_final, &D);

    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Z_final);


    //tmp = np.einsum("ijab,ijac->bc",self.x2_B,self.Loovv)  
    //tmp = np.einsum("bc,klcd->klbd",tmp,self.x2_C)

    global_dpd_->file2_init(&z_ab, PSIF_CC_OEI, 0, 1, 1, "z_ab");

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->contract442(&X2, &D, &z_ab, 3, 3, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z (ij,ab)");
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract244(&z_ab, &X2, &Z_final, 1, 2, 1, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->file2_close(&z_ab);

    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);

    result -= global_dpd_->buf4_dot(&Z_final, &Y2);

    global_dpd_->buf4_close(&Y2);
    global_dpd_->buf4_close(&Z_final);



    //tmp = np.einsum("ijab,ikab->jk",self.x2_B,self.Loovv)  
    //tmp = np.einsum("jk,klcd->jlcd",tmp,self.x2_C)
    //self.Bcon1 -= np.einsum("jlcd,jlcd->",tmp,self.y2_A)

    global_dpd_->file2_init(&z_ij, PSIF_CC_OEI, 0, 0, 0, "z_ij");

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->contract442(&X2, &D, &z_ij, 1, 1, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z (ij,ab)");
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract244(&z_ij, &X2, &Z_final, 1, 0, 0, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->file2_close(&z_ij);

    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
    result -= global_dpd_->buf4_dot(&Z_final, &Y2);

    global_dpd_->buf4_close(&Y2);
    global_dpd_->buf4_close(&Z_final);


    //tmp = np.einsum("ikbc,klcd->ilbd",self.Loovv,self.x2_C)
    //tmp = np.einsum("ilbd,ijab->jlad",tmp,self.x2_B)
    //self.Bcon1 -= np.einsum("jlad,jlad->",tmp,self.y2_A)

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (ia|jb)");
    sprintf(lbl, "X_%s_IAjb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 10, 10, 10, 10, 0, lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");
    global_dpd_->contract444(&D, &X2, &Z, 0, 1, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (ia,jb)");
    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_y, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&Z, &X2, &Z_final, 1, 1, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_sort(&Z_final, PSIF_CC_TMP0, rpsq, 0, 5, "Z (ij|ab)");
    global_dpd_->buf4_close(&Z_final);

    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z (ij|ab)");
    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);

    result -= global_dpd_->buf4_dot(&Z_final, &Y2);

    global_dpd_->buf4_close(&Z_final);
    global_dpd_->buf4_close(&Y2);



    //tmp = np.einsum("ijab,jlbc->ilac",self.x2_B,self.y2_A)
    //tmp = np.einsum("ilac,klcd->ikad",tmp,self.x2_C)
    //self.Bcon1 -= np.einsum("ikad,ikad->",tmp,self.Loovv)

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (ia|jb)");
    sprintf(lbl, "X_%s_IAjb (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_y, 10, 10, 10, 10, 0, lbl);
    sprintf(lbl, "Y_%s_IAjb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&X2, &Y2, &Z, 0, 1, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&Y2);


    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (ia,jb)");
    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&Z, &X2, &Z_final, 0, 0, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&Z);


    global_dpd_->buf4_sort(&Z_final, PSIF_CC_TMP0, prqs, 0, 5, "Z (ij|ab)");
    global_dpd_->buf4_close(&Z_final);

    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z (ij|ab)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");

    result -= global_dpd_->buf4_dot(&Z_final, &D);

    global_dpd_->buf4_close(&Z_final);
    global_dpd_->buf4_close(&D);


    //tmp = np.einsum("klca,klcd->ad",self.Loovv,self.x2_C)
    //tmp = np.einsum("ad,ijdb->ijab",tmp,self.y2_A)
    //self.Bcon1 -= np.einsum("ijab,ijab->",tmp,self.x2_B)

    global_dpd_->file2_init(&z_ab, PSIF_CC_OEI, 0, 1, 1, "z_ab");

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->contract442(&D, &X2, &z_ab, 3, 3, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z (ij,ab)");
    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract244(&z_ab, &Y2, &Z_final, 1, 2, 1, 1, 0);

    global_dpd_->buf4_close(&Y2);
    global_dpd_->file2_close(&z_ab);

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);

    result -= global_dpd_->buf4_dot(&Z_final, &X2);

    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&Z_final);


    //tmp = np.einsum("kicd,klcd->il",self.Loovv,self.x2_C)
    //tmp = np.einsum("ijab,il->ljab",self.x2_B,tmp)
    //self.Bcon1 -= np.einsum("ljab,ljab->",tmp,self.y2_A)

    global_dpd_->file2_init(&z_ij, PSIF_CC_OEI, 0, 0, 0, "z_ij");

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->contract442(&X2, &D, &z_ij, 1, 1, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&D);

    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z (ij,ab)");
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract244(&z_ij, &X2, &Z_final, 1, 0, 0, 1, 0);

    global_dpd_->buf4_close(&X2);
    global_dpd_->file2_close(&z_ij);

    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
    
    result -= global_dpd_->buf4_dot(&Z_final, &Y2);

    global_dpd_->buf4_close(&Y2);
    global_dpd_->buf4_close(&Z_final);


    //tmp = np.einsum("klcd,ikac->lida",self.x2_C,self.y2_A)
    //tmp = np.einsum("lida,jlbd->ijab",tmp,self.Loovv)
    //self.Bcon1 += 2.*np.einsum("ijab,ijab->",tmp,self.x2_B)

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (ia|jb)");
    sprintf(lbl, "X_%s_IAjb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 10, 10, 10, 10, 0, lbl);
    sprintf(lbl, "Y_%s_IAjb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&X2, &Y2, &Z, 1, 0, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&Y2);

    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (ia,jb)");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");
    global_dpd_->contract444(&Z, &D, &Z_final, 1, 0, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_sort(&Z_final, PSIF_CC_TMP0, prqs, 0, 5, "Z (ij|ab)");
    global_dpd_->buf4_close(&Z_final);

    global_dpd_->buf4_init(&Z_final, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z (ij|ab)");
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_y, 0, 5, 0, 5, 0, lbl);

    result += 2 * global_dpd_->buf4_dot(&Z_final, &X2);

    global_dpd_->buf4_close(&Z_final);
    global_dpd_->buf4_close(&X2);


    return result;
}


double Y1HX1X2(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
                      const char *pert_z, int irrep_z, double omega_z) {

    double result = 0.0;
    dpdfile2 X1, Y1, z_ij, z_jb, z_ia, z_ab;
    dpdbuf4 X2, Y2, D;
    char lbl[32];

    //tmp  = 2.*np.einsum("jkbc,kc->jb",self.x2_C,self.y1_A)
    //tmp -= np.einsum("jkcb,kc->jb",self.x2_C,self.y1_A)
    //tmp = np.einsum('ijab,jb->ia',self.Loovv,tmp)
    //self.Bcon1 += np.einsum("ia,ia->",tmp,self.x1_B)

    global_dpd_->file2_init(&z_jb, PSIF_CC_OEI, 0, 0, 1, "z_jb");

    sprintf(lbl, "X_%s_IAjb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 10, 10, 10, 10, 0, lbl);
    sprintf(lbl, "Y_%s_IA (%5.3f)", pert_x, omega_x);
    global_dpd_->file2_init(&Y1, PSIF_CC_OEI, irrep_x, 0, 1, lbl);
    global_dpd_->contract422(&X2, &Y1, &z_jb, 0, 0, 2, 0);
    global_dpd_->buf4_close(&X2);


    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract422(&X2, &Y1, &z_jb, 0, 0, -1, 1);
    global_dpd_->buf4_close(&X2);
    global_dpd_->file2_close(&Y1);


    global_dpd_->file2_init(&z_ia, PSIF_CC_OEI, 0, 0, 1, "z_ia");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");
    global_dpd_->contract422(&D, &z_jb, &z_ia, 0, 0, 1, 0);
    global_dpd_->buf4_close(&D);
    global_dpd_->file2_close(&z_jb);


    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    
    result += global_dpd_->file2_dot(&X1, &z_ia);

    global_dpd_->file2_close(&X1); 
    global_dpd_->file2_close(&z_ia);


    //tmp = np.einsum("jkbc,jkba->ca",self.x2_C,self.Loovv)
    //tmp = np.einsum("ia,ca->ic",self.x1_B,tmp)
    //self.Bcon1 -= np.einsum("ic,ic->",tmp,self.y1_A)

    global_dpd_->file2_init(&z_ab, PSIF_CC_OEI, 0, 1, 1, "z_ab");

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->contract442(&X2, &D, &z_ab, 3, 3, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&D);


    global_dpd_->file2_init(&z_ia, PSIF_CC_OEI, 0, 0, 1, "z_ia");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->contract222(&X1, &z_ab, &z_ia, 0, 0, 1, 0); 
    global_dpd_->file2_close(&X1);
    global_dpd_->file2_close(&z_ab);


    sprintf(lbl, "Y_%s_IA (%5.3f)", pert_x, omega_x);
    global_dpd_->file2_init(&Y1, PSIF_CC_OEI, irrep_x, 0, 1, lbl);


    result -= global_dpd_->file2_dot(&Y1, &z_ia);

    global_dpd_->file2_close(&Y1);
    global_dpd_->file2_close(&z_ia);


    //tmp = np.einsum("jkbc,jibc->ki",self.x2_C,self.Loovv)
    //tmp = np.einsum("ki,ia->ka",tmp,self.x1_B)
    //self.Bcon1 -= np.einsum("ka,ka->",tmp,self.y1_A) 

    global_dpd_->file2_init(&z_ij, PSIF_CC_OEI, 0, 0, 0, "z_ij");

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->contract442(&X2, &D, &z_ij, 1, 1, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&D);


    global_dpd_->file2_init(&z_ia, PSIF_CC_OEI, 0, 0, 1, "z_ia");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->contract222(&z_ij, &X1, &z_ia, 0, 1, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->file2_close(&z_ij);


    sprintf(lbl, "Y_%s_IA (%5.3f)", pert_x, omega_x);
    global_dpd_->file2_init(&Y1, PSIF_CC_OEI, irrep_x, 0, 1, lbl);


    result -= global_dpd_->file2_dot(&Y1, &z_ia);

    global_dpd_->file2_close(&Y1);
    global_dpd_->file2_close(&z_ia);
   
    //outfile->Printf("\n\tResult B1:  %20.15f\n", result);

    return result;

}


double Y2HX1X2(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
                      const char *pert_z, int irrep_z, double omega_z) {

    double result = 0.0;
    dpdfile2 X1, Y1, z_ij, z2_ij, z_ab, z2_ab, z_ia, F;
    dpdbuf4 X2, Y2, W, WX1, Z, Z2;
    char lbl[32];
    double Y1_norm;

    // *** <O|L2(A)[[Hbar(0),X1(B)],X2(C)]]|0> ***

    // tmp = np.einsum("klcd,lkdb->cb",self.x2_C,self.y2_A)
    // tmp = np.einsum("jb,cb->jc",self.x1_B,tmp)
    // self.Bcon1 -= np.einsum("jc,jc->",tmp,self.Hov)

     //***We can replace by this part of the code ****
    // tmp = np.einsum("jc,jb->cb",self.Hov,self.x1_B)
    // tmp2 = np.einsum("klcd,lkdb->cb",self.x2_C,self.y2_A)  
    // self.Bcon1 -= np.einsum("cb,cb->",tmp,tmp2)


    global_dpd_->file2_init(&z_ab, PSIF_CC_OEI, 0, 1, 1, "z_ab");

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 0, 5, 0, 5, 0, lbl);
    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract442(&X2, &Y2, &z_ab, 2, 2, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&Y2);

    global_dpd_->file2_init(&z_ia, PSIF_CC_OEI, 0, 0, 1, "z_ia");
   
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl); 
    global_dpd_->contract222(&X1, &z_ab, &z_ia, 0, 0, 1.0, 0.0);	
    global_dpd_->file2_close(&z_ab);
    global_dpd_->file2_close(&X1);

    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "FME");

    result -= global_dpd_->file2_dot(&F, &z_ia);

    global_dpd_->file2_close(&F);
    global_dpd_->file2_close(&z_ia);


    // tmp = np.einsum("klcd,ljdc->kj",self.x2_C,self.y2_A)
    // tmp = np.einsum("kj,jb->kb",tmp,self.x1_B)
    // self.Bcon1 -= np.einsum("kb,kb->",tmp,self.Hov)


    global_dpd_->file2_init(&z_ij, PSIF_CC_OEI, 0, 0, 0, "z_ij");

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 0, 5, 0, 5, 0, lbl);
    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract442(&X2, &Y2, &z_ij, 0, 0, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&Y2);

    global_dpd_->file2_init(&z_ia, PSIF_CC_OEI, 0, 0, 1, "z_ia");

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->contract222(&z_ij, &X1, &z_ia, 0, 1, 1.0, 0.0);
    global_dpd_->file2_close(&z_ij);
    global_dpd_->file2_close(&X1);

    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "FME");

    result -= global_dpd_->file2_dot(&F, &z_ia);

    global_dpd_->file2_close(&F);
    global_dpd_->file2_close(&z_ia);


    // tmp = np.einsum('lkda,klcd->ac',self.y2_A,self.x2_C)
    // tmp2 = np.einsum('jb,ajcb->ac',self.x1_B,self.Hvovv)
    // self.Bcon1 += 2.*np.einsum('ac,ac->',tmp,tmp2)

    global_dpd_->file2_init(&z_ab, PSIF_CC_OEI, 0, 1, 1, "z_ab");

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 0, 5, 0, 5, 0, lbl);
    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract442(&Y2, &X2, &z_ab, 2, 2, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&Y2);

    global_dpd_->file2_init(&z2_ab, PSIF_CC_OEI, 0, 1, 1, "z2_ab");

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf"); 
    global_dpd_->dot24(&X1, &W, &z2_ab, 0, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&W);
    
    result += 2*global_dpd_->file2_dot(&z_ab, &z2_ab);

    global_dpd_->file2_close(&z_ab);
    global_dpd_->file2_close(&z2_ab);

    // tmp = np.einsum('lkda,klcd->ac',self.y2_A,self.x2_C)
    // tmp2 = np.einsum('jb,ajbc->ac',self.x1_B,self.Hvovv)
    // self.Bcon1 -= np.einsum('ac,ac->',tmp,tmp2)

    global_dpd_->file2_init(&z_ab, PSIF_CC_OEI, 0, 1, 1, "z_ab");

/*
    // I am reusing this part
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 0, 5, 0, 5, 0, lbl);
    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract442(&Y2, &X2, &z_ab, 2, 2, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&Y2);
*/

    global_dpd_->file2_init(&z2_ab, PSIF_CC_OEI, 0, 1, 1, "z2_ab");

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    //global_dpd_->file2_init(&z_ia, PSIF_CC_OEI, 0, 0, 1, "z_ia");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    global_dpd_->dot23(&X1, &W, &z2_ab, 0, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&W);

    result -= global_dpd_->file2_dot(&z_ab, &z2_ab);


    // Hv = 2*self.Hvovv - self.Hvovv.swapaxes(2,3)
    // tmp = np.einsum('klcd,ljda->kjca',self.x2_C,self.y2_A)
    // tmp2 = np.einsum('jb,akbc->akjc',self.x1_B,Hv)
    // Bcon1 += np.einsum('akjc,kjca->',tmp2,tmp)

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (ia,jb)");

    global_dpd_->file2_init(&z_ab, PSIF_CC_OEI, 0, 1, 1, "z_ab");
    sprintf(lbl, "X_%s_IAjb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 10, 10, 10, 10, 0, lbl);
    sprintf(lbl, "Y_%s_IAjb (%5.3f)", pert_x, omega_x);

    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&X2, &Y2, &Z, 0, 1, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&Y2);

    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, sprq, 11, 10, "Z (bi,ja)");   //sort
    global_dpd_->buf4_close(&Z);

//    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 11, 10, 11, 10, 0, "Z (ak,jc)");

//MK I am here ...
/*
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)");
    global_dpd_->contract244(&X1, &W, &Z2, 1, 2, 1, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 11, 10, 11, 10, 0, "Z (bi,ja)");
*/

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 11, 10, 11, 10, 0, "Z (bi,ja)");
    sprintf(lbl, "WX_%s_fima (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&WX1, PSIF_CC_LR, irrep_y, 11, 10, 11, 10, 0, lbl);

    result += global_dpd_->buf4_dot(&Z, &WX1);
    //result += global_dpd_->buf4_dot(&Z, &Z2);

    global_dpd_->buf4_close(&Z);
    //global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&WX1);

 //It is not working!!!!

    //tmp = np.einsum('ia,fkba->fkbi',self.x1_B,self.Hvovv)
    //tmp = np.einsum('fkbi,jifc->kjbc',tmp,self.y2_A)
    //self.Bcon1 -= np.einsum('jkbc,kjbc->',self.x2_C,tmp)

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 11, 11, 11, 11, 0, "Z (ai,bj)");

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
 

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    global_dpd_->contract424(&W, &X1, &Z, 3, 1, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&W);


    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, rqsp, 11, 10, "Z (bi,ja)");
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 11, 10, 11, 10, 0, "Z (ai,jb)");
    sprintf(lbl, "Y_%s_IbjA (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 10, 10, 10, 10, 0, lbl);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 11, 10, 11, 10, 0, "Z (bi,ja)");
    global_dpd_->contract444(&Z, &Y2, &Z2, 0, 0, 1, 0);

    global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP0, rqps, 0, 5, "Z (ij,ab)");
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z (ij,ab)");
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 0, 5, 0, 5, 0, lbl);

    result -= global_dpd_->buf4_dot(&X2, &Z2);

    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&X2);



    // tmp = np.einsum('ia,fjac->fjic',self.x1_B,self.Hvovv)
    // tmp = np.einsum('fjic,ikfb->jkbc',tmp,self.y2_A)
    // self.Bcon1 -= np.einsum('jkbc,jkbc->',self.x2_C,tmp)

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 11, 10, 11, 10, 0, "Z (ai,jb)");

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    global_dpd_->contract244(&X1, &W, &Z, 1, 2, 1, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, rpqs, 10, 10, "Z (ja,ib)");   
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (ia,jb)");

    sprintf(lbl, "Y_%s_IAjb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 10, 10, 10, 10, 0, lbl);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 10, 10, 10, 10, 0, "Z (ja,ib)");
    global_dpd_->contract444(&Z, &Y2, &Z2, 1, 1, 1, 0);
    global_dpd_->buf4_close(&Y2);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP0, prsq, 0, 5, "Z (ij,ab)");
    global_dpd_->buf4_close(&Z2);


    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z (ij,ab)");
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 0, 5, 0, 5, 0, lbl);
 
    result -= global_dpd_->buf4_dot(&X2, &Z2);

    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&X2);   


   //tmp = np.einsum('ia,jkfa->jkfi',self.x1_B,self.y2_A)
   //tmp2 = np.einsum('jkbc,fibc->jkfi',self.x2_C,self.Hvovv)
   //self.Bcon1 -= np.einsum('jkfi,jkfi->',tmp2,tmp)


    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 11, 0, 11, 0, "Z (ij,ak)");

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract424(&Y2, &X1, &Z, 3, 1, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&Y2);


    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 0, 11, 0, 11, 0, "Z2 (ij,ak)");

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");   
    global_dpd_->contract444(&X2, &W, &Z2, 0, 0, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&W);


    result -= global_dpd_->buf4_dot(&Z, &Z2);

    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&Z2);



    // tmp = 2*np.einsum('jb,kjib->ki',self.x1_B,self.Hooov)
    // tmp += np.einsum('jb,jkib->ki',self.x1_B,self.Hooov)
    // tmp2 = np.einsum('klcd,ilcd->ki',self.x2_C,self.y2_A)
    // self.Bcon1 -= np.einsum('ki,ki->',tmp,tmp2)
   
/*
    global_dpd_->file2_init(&z_ij, PSIF_CC_OEI, 0, 0, 0, "z_ij");

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    //global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe"); 
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (Mn,eI)");
    //global_dpd_->dot24(&X1, &W, &z_ij, 0, 0, 1, 0);
    global_dpd_->dot23(&X1, &W, &z_ij, 0, 0, 1, 0);
    global_dpd_->file2_close(&X1);
*/

    global_dpd_->file2_init(&z2_ij, PSIF_CC_OEI, 0, 0, 0, "z2_ij");
    
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 0, 5, 0, 5, 0, lbl);
    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract442(&X2, &Y2, &z2_ij, 0, 0, 1, 0); 
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&Y2);

    sprintf(lbl, "WX_%s_IJ (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&z_ij, PSIF_CC_OEI, irrep_y, 0, 0, lbl);

    //result -= 2*global_dpd_->file2_dot(&z_ij, &z2_ij);
    result -= global_dpd_->file2_dot(&z_ij, &z2_ij);


    // tmp = np.einsum('jb,jkib->ki',self.x1_B,self.Hooov)
    // tmp2 = np.einsum('klcd,ilcd->ki',self.x2_C,self.y2_A)
    // self.Bcon1 += np.einsum('ki,ki->',tmp,tmp2)

/*
    global_dpd_->file2_init(&z_ij, PSIF_CC_OEI, 0, 0, 0, "z_ij");

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
    global_dpd_->dot14(&X1, &W, &z_ij, 0, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&W);


    global_dpd_->file2_init(&z2_ij, PSIF_CC_OEI, 0, 0, 0, "z2_ij");

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 0, 5, 0, 5, 0, lbl);
    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract442(&X2, &Y2, &z2_ij, 0, 0, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&Y2);


    result += global_dpd_->file2_dot(&z_ij, &z2_ij);

    global_dpd_->file2_close(&z_ij);
    global_dpd_->file2_close(&z2_ij);
*/

    // tmp  = 2.*np.einsum('jkic,klcd->jild',self.Hooov,self.x2_C)
    // tmp -= np.einsum('kjic,klcd->jild',self.Hooov,self.x2_C)
    // tmp  = np.einsum('jild,jb->bild',tmp,self.x1_B)
    // self.Bcon1 -= np.einsum('bild,ilbd->',tmp,self.y2_A)


    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 10, 0, 10, 0, "Z (ij,ka)");
    
    sprintf(lbl, "X_%s_IAjb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 10, 10, 10, 10, 0, lbl);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "2WMnIe - WnMIe (MI,nE)");
    global_dpd_->contract444(&W, &X2, &Z, 0, 1, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&W);


    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 11, 10, 11, 10, 0, "Z2 (bj,ka)");
   
    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->contract244(&X1, &Z, &Z2, 0, 0, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&Z);


    global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP0, qrps, 0, 5, "Z2 (ij,ab)");
    global_dpd_->buf4_close(&Z2);


    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z2 (ij,ab)"); 

    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);

    result -= global_dpd_->buf4_dot(&Y2, &Z2);

    global_dpd_->buf4_close(&Y2);
    global_dpd_->buf4_close(&Z2);


    // tmp  = np.einsum('ia,jkna->jkni',self.x1_B,self.Hooov)
    // tmp2  = np.einsum('jkbc,nibc->jkni',self.x2_C,self.y2_A)
    // self.Bcon1 += np.einsum('jkni,jkni->',tmp2,tmp)


/*
    global_dpd_->buf4_init(&Z, PSIF_CC_LR, 0, 0, 0, 0, 0, 0, "Z(ij,kl)");

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
    global_dpd_->contract424(&W, &X1, &Z, 3, 1, 0, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&W);
*/

    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 0, 0, 0, 0, 0, "Z2 (ij,kl)");

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 0, 5, 0, 5, 0, lbl);
    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 0, 5, 0, 5, 0, lbl);
    global_dpd_->contract444(&X2, &Y2, &Z2, 0, 0, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&Y2);

    sprintf(lbl, "WX_%s_ijkl (%5.3f)", pert_y, omega_y);
    global_dpd_->buf4_init(&Z, PSIF_CC_LR, 0, 0, 0, 0, 0, 0, lbl);

    result += global_dpd_->buf4_dot(&Z, &Z2);

    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&Z2);


   // Improve?
   // tmp  = np.einsum('ia,nkab->nkib',self.x1_B,self.y2_A)
   // tmp  = np.einsum('jkbc,nkib->jnic',self.x2_C,tmp)
   // self.Bcon1 += np.einsum('jnic,ijnc->',tmp,self.Hooov)

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 10, 0, 10, 0, "Z (ij,ka)");

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    sprintf(lbl, "Y_%s_IAjb (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract424(&Y2, &X1, &Z, 1, 1, 1, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&Y2);


    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 0, 10, 0, 10, 0, "Z2 (jk,ia)");

    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&Z, &X2, &Z2, 0, 0, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP0, qrps, 0, 10, "Z2 (ik,ja)"); 
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 0, 10, 0, 10, 0, "Z2 (ik,ja)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");

    result += global_dpd_->buf4_dot(&Z2, &W);

    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&W);



    //Improve?
    //tmp  = np.einsum('ia,nkba->nkbi',self.x1_B,self.y2_A)
    //tmp  = np.einsum('jkbc,nkbi->jnci',self.x2_C,tmp)
    //self.Bcon1 += np.einsum('jnci,jinc->',tmp,self.Hooov)

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 0, 10, 0, 10, 0, "Z (ij,ka)");

    sprintf(lbl, "X_%s_IA (%5.3f)", pert_y, omega_y);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep_y, 0, 1, lbl);
    sprintf(lbl, "Y_%s_IbjA (%5.3f)", pert_x, omega_x);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep_x, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract424(&Y2, &X1, &Z, 1, 1, 1, 1, 0);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&Y2);


    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 0, 10, 0, 10, 0, "Z2 (jk,ia)");

    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert_z, omega_z);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep_z, 10, 10, 10, 10, 0, lbl);
    global_dpd_->contract444(&Z, &X2, &Z2, 0, 0, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&Z);


    global_dpd_->buf4_sort(&Z2, PSIF_CC_TMP0, rqps, 0, 10, "Z2 (ik,ja)");
    global_dpd_->buf4_close(&Z2);

    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 0, 10, 0, 10, 0, "Z2 (ik,ja)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");

    result += global_dpd_->buf4_dot(&Z2, &W);

    global_dpd_->buf4_close(&Z2);
    global_dpd_->buf4_close(&W);

    return result;

}

double YHXX(const char *pert_x, int irrep_x, double omega_x, const char *pert_y, int irrep_y, double omega_y,
		     const char *pert_z, int irrep_z, double omega_z) {

    double hyper = 0.0;

    // *** <O|Y1(A)[[Hbar(0),X1(B),X1(C)]]|0> ***

    hyper += Y1HX1X1(pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z); 

   // ***  <O|Y2(A)|[[Hbar(0),X1(B)],X1(C)]|0> ***

    hyper += Y2HX1X1(pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z);

   // ***  <O|Y2(A)|[[Hbar(0),X2(B)],X2(C)]|0> ***

    hyper += Y2HX2X2(pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z);    

   // *** <O|Y1(A)[[Hbar(0),X1(B)],X2(C)]]|0> ***  
   
    hyper += Y1HX1X2(pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z);

   // *** <O|Y2(A)[[Hbar(0),X1(B)],X2(C)]]|0> ***  

    hyper += Y2HX1X2(pert_x, irrep_x, omega_x, pert_y, irrep_y, omega_y, pert_z, irrep_z, omega_z);

   // *** <O|Y1(A)[[Hbar(0),X2(B),X1(C)]]|0> ***

    hyper += Y1HX1X2(pert_x, irrep_x, omega_x, pert_z, irrep_z, omega_z, pert_y, irrep_y, omega_y);

   // *** <O|Y2(A)[[Hbar(0),X2(B)],X1(C)]]|0> ***   

    hyper += Y2HX1X2(pert_x, irrep_x, omega_x, pert_z, irrep_z, omega_z, pert_y, irrep_y, omega_y);

    return hyper;
}

}  // namespace ccresponse
}  // namespace psi
