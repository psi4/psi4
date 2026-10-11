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

void local_filter_T1(dpdfile2 *);
void lambda_residuals();
void L2HX2(const char *pert, int irrep, double omega);

void Y1_inhomogenous_build(const char *pert, int irrep, double omega) {
    int Gim, Gi, Gm, Ga, Gam, nrows, ncols, A, a, am;
    int Gei, ei, e, i, Gef, Ge, Gf, E, I, af, fa, f;
    int GW, GX1, GZ, Gej, Gab, Gij, Gj;
    int num_j, num_i, num_e, nlinks;
    dpdfile2 F, z1, z2;
    dpdbuf4 W, WL, D, X2, Z2, Z3, lx_iajb; 
    dpdfile2 Y1, Y1new, X1, L1, mu1, GAE, GMI, lt, lx, lx_ia, lx_AB;
    dpdbuf4 L2, Z, mu2, Hx_ijab, lx_ijab;
    char lbl[32];
    double Y1_norm;
    double *X,*Y;

    sprintf(lbl, "Inhomo Y_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&Y1new, PSIF_CC_OEI, irrep, 0, 1, lbl);

    /*** Mu * L1 ***/ 

    sprintf(lbl, "%s_IA", pert);
    global_dpd_->file2_init(&mu1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->file2_axpy(&mu1, &Y1new, 2, 0);  
    global_dpd_->file2_close(&mu1);

    /*** L1 * MuBAR + L2 * MuBAR ***/

    global_dpd_->file2_init(&L1, PSIF_CC_LAMPS, 0, 0, 1, "LIA 0 -1"); 
    sprintf(lbl, "%sBAR_MI", pert);
    global_dpd_->file2_init(&mu1, PSIF_CC_OEI, irrep, 0, 0, lbl);
    global_dpd_->contract222(&mu1, &L1, &Y1new, 0, 1, -1, 1.0);
    global_dpd_->file2_close(&mu1);

    sprintf(lbl, "%sBAR_AE", pert);
    global_dpd_->file2_init(&mu1, PSIF_CC_OEI, irrep, 1, 1, lbl);
    global_dpd_->contract222(&L1, &mu1, &Y1new, 0, 1, 1.0, 1.0);
    global_dpd_->file2_close(&mu1);
    global_dpd_->file2_close(&L1);

    sprintf(lbl, "%s_IA", pert);
    global_dpd_->file2_init(&mu1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->file2_init(&lt, PSIF_CC_OEI, 0, 0, 0, "Lt_IJ");
    global_dpd_->contract222(&lt, &mu1, &Y1new, 0, 1, 2, 1.0);   //I multiplied by 2

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa"); 
    sprintf(lbl, "%sBAR_MbIj", pert, omega);
    global_dpd_->buf4_init(&mu2, PSIF_CC_LR, irrep, 10, 0, 10, 0, 0, lbl);
    global_dpd_->contract442(&mu2, &L2, &Y1new, 0, 2, -0.5, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&mu2);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    sprintf(lbl, "%sBAR_MbIj", pert);
    global_dpd_->buf4_init(&mu2, PSIF_CC_LR, irrep, 10, 0, 10, 0, 0, lbl);
    global_dpd_->contract442(&mu2, &L2, &Y1new, 0, 2, -0.5, 1.0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&mu2);

    // *** <O|[Hbar(0), X1]|0> ***
    
    sprintf(lbl, "DX_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&lx_ia, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 10, 10, 10, 10, 0, "D 2<ij|ab> - <ij|ba> (ia,jb)");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->contract422(&D, &X1, &lx_ia, 0, 0, 1, 1);
    global_dpd_->file2_axpy(&lx_ia, &Y1new, 2, 0);
    global_dpd_->file2_close(&lx_ia);
    global_dpd_->buf4_close(&D); 
    global_dpd_->file2_close(&X1);   

    // *** <O|L1(0)|[Hbar(0), X1]|0> ***

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, "Z (ij,ab)");
    global_dpd_->file2_init(&L1, PSIF_CC_LAMPS, 0, 0, 1, "LIA 0 -1"); 
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (Mn,eI)");
    global_dpd_->contract424(&W, &L1, &Z, 3, 0, 0, -1, 0);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, qpsr, 0, 5, "Z (ji,ba)");
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, "Z (ij,ab)");   //Revisit it!
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, irrep, 0, 5, 0, 5, 0, "Z (ji,ba)");  //Revisit it!
    global_dpd_->buf4_axpy(&Z2, &Z, 1);
    global_dpd_->buf4_close(&Z2);

    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);


    //MK I am OPTIMIZING IT


    //global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z3 (ij,ab)");
    //global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 5, 10, 5, 0, "WAmEf 2(mA,Ef) - (mA,fE)");
    //global_dpd_->contract424(&W, &L1, &Z2, 1, 1, 1, 1, 0);      

/*
    global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "Z3 (ij,ab)");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    //global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 10, 5, 10, 5, 0, "WAmEf 2(mA,Ef) - (mA,fE)");

        //       dpd_contract244(&LIA, &W, &Z, 1, 2, 1, 1, 0); //
        /// Out-of-core contract244 /
        GW = W.file.my_irrep;
        GZ = Z2.file.my_irrep;
        GX1 = X1.my_irrep;

        global_dpd_->file2_mat_init(&X1);
        global_dpd_->file2_mat_rd(&X1);

        for (Gej = 0; Gej < moinfo.nirreps; Gej++) {
            Gab = Gej ^ GW;
            Gij = Gab ^ GZ;

            global_dpd_->buf4_mat_irrep_init(&Z2, Gij);

            for (Ge = 0; Ge < moinfo.nirreps; Ge++) {
                Gi = Ge ^ GX1;
                Gj = GZ ^ Gab ^ Gi;

                num_j = Z2.params->qpi[Gj];
                num_i = X1.params->rowtot[Gi];
                num_e = X1.params->coltot[Ge];

                global_dpd_->buf4_mat_irrep_init_block(&W, Gej, num_j);

                for (e = 0; e < num_e; e++) {
                    E = W.params->poff[Ge] + e;
                    global_dpd_->buf4_mat_irrep_rd_block(&W, Gej, W.row_offset[Gej][E], num_j);

                    for (i = 0; i < num_i; i++) {
                        I = Z2.params->poff[Gi] + i;

                        nlinks = Z2.params->coltot[Gab] * num_j;
                        if (nlinks) {
                            C_DAXPY(nlinks, X1.matrix[Gi][i][e], &(W.matrix[Gej][0][0]), 1,
                                    &(Z2.matrix[Gij][Z2.row_offset[Gij][I]][0]), 1);
                        }
                    }
                }
                global_dpd_->buf4_mat_irrep_close_block(&W, Gej, num_j);
            }
            global_dpd_->buf4_mat_irrep_wrt(&Z2, Gij);
            global_dpd_->buf4_mat_irrep_close(&Z2, Gij);
        }

        /// End out-of-core contract244 ///
        //global_dpd_->buf4_close(&W);

*/

    global_dpd_->buf4_init(&WL, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "WL(ij,ab) + WL(ji,ba)");
    global_dpd_->buf4_axpy(&WL, &Z, 1);
    global_dpd_->buf4_close(&WL);

    global_dpd_->dot14(&X1, &Z, &Y1new, 0, 0, 1, 1);
    global_dpd_->file2_close(&X1);
    global_dpd_->file2_close(&L1); 
    global_dpd_->buf4_close(&Z);

    // *** <O|L2(0)|[Hbar(0), X1]|0> ***

    global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 0, 5, 0, 5, 0, "Z (ij,ab)");

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 0, 0, 0, 0, "WMnIj (nM,Ij)");
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract444(&W, &L2, &Z, 0, 1, 1.0, 0.0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&WL, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "WefabL2 2(ij,ba) - (ij,ab)");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->buf4_axpy(&WL, &Z, 1.0); 
    global_dpd_->dot23(&X1, &Z, &Y1new, 0, 0, 1, 1);	
    global_dpd_->buf4_close(&WL);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&Z);

    //Improve this part....
    //Build Gae*.Loovv -> Take it straight from L2???
    //############################################################
//    global_dpd_->file2_init(&GAE, PSIF_CC_LAMBDA, 0, 1, 1, "GAE");
//    global_dpd_->file2_scm(&GAE, 2); //I multiplied by 2
//    global_dpd_->file2_close(&GAE);
 
//    global_dpd_->file2_init(&GMI, PSIF_CC_LAMBDA, 0, 0, 0, "GMI");
//    global_dpd_->file2_scm(&GMI, 2); //I multiplied by 2
//    global_dpd_->file2_close(&GMI);
    //############################################################

    //tmp  =  ndot('nb,fb->nf', self.x1, self.build_Gvv(self.t2, self.l2))
    //r_y1 += ndot('inaf,nf->ia', self.Loovv, tmp) 
    //tmp  =  ndot('me,fa->mefa', self.x1, self.build_Gvv(self.t2, self.l2))
    //r_y1 += ndot('mief,mefa->ia', self.Loovv, tmp)

    global_dpd_->buf4_init(&Z, PSIF_CC_LAMBDA, 0, 0, 5, 0, 5, 0, "GAED (ij,ab)");
    global_dpd_->file2_init(&GAE, PSIF_CC_LAMBDA, 0, 1, 1, "GAE");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->contract424(&D, &GAE, &Z, 3, 1, 0, 1.0, 0.0);
    global_dpd_->file2_close(&GAE);
    global_dpd_->buf4_close(&D);    
      
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->dot24(&X1, &Z, &Y1new, 0, 0, 1, 1);
    global_dpd_->dot13(&X1, &Z, &Y1new, 0, 0, 1, 1);
    global_dpd_->buf4_close(&Z);

    //tmp  =  ndot('me,ni->meni', self.x1, self.build_Goo(self.t2, self.l2))
    //r_y1 -= ndot('meni,mnea->ia', tmp, self.Loovv)
    //tmp  =  ndot('jf,nj->fn', self.x1, self.build_Goo(self.t2, self.l2))
    //r_y1 -= ndot('inaf,fn->ia', self.Loovv, tmp)

    global_dpd_->buf4_init(&Z, PSIF_CC_LAMBDA, 0, 0, 5, 0, 5, 0, "GMID (ij,ab)");
    
    global_dpd_->file2_init(&GMI, PSIF_CC_LAMBDA, 0, 0, 0, "GMI");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->contract424(&D, &GMI, &Z, 1, 0, 1, 1, 0);
    global_dpd_->dot13(&X1, &Z, &Y1new, 0, 0, -1, 1);
    global_dpd_->dot24(&X1, &Z, &Y1new, 0, 0, -1, 1);
    global_dpd_->file2_close(&GMI);
    global_dpd_->buf4_close(&D);
    global_dpd_->file2_close(&X1);
    global_dpd_->buf4_close(&Z);    

    // Type-II L2 residual 
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 10, 10, 10, 10, 0, "LHX1Y1 (ia,jb) + (jb,ia) Residual II");
    sprintf(lbl, "X_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&X1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    global_dpd_->contract422(&L2, &X1, &Y1new, 0, 0, -1, 1);  
    global_dpd_->buf4_close(&L2);
    global_dpd_->file2_close(&X1);
    
    // *** <O|L1(0)|[Hbar(0), X2]|phi^a_i> ***

    sprintf(lbl, "Z_%s_ME", pert);
    global_dpd_->file2_init(&z1, PSIF_CC_TMP0, irrep, 0, 1, lbl);
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->file2_init(&L1, PSIF_CC_LAMPS, 0, 0, 1, "LIA 0 -1"); 
    global_dpd_->dot24(&L1, &X2, &z1, 0, 0, 2, 0); 
    
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);  
    sprintf(lbl, "Z2_%s_ME", pert);
    global_dpd_->dot23(&L1, &X2, &z1, 0, 0, -1, 1);
    global_dpd_->buf4_close(&X2);
 
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->dot24(&z1, &D, &Y1new, 0, 0, 1, 1); 
    global_dpd_->buf4_close(&D); 
    global_dpd_->file2_close(&z1);

    sprintf(lbl, "Z_%s_MN", pert);
    global_dpd_->file2_init(&z1, PSIF_CC_TMP0, irrep, 0, 0, lbl);
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->contract442(&X2, &D, &z1, 0, 0, 1.0, 0.0);
    global_dpd_->contract222(&z1, &L1, &Y1new, 1, 1, -1.0, 1.0);
    global_dpd_->file2_close(&z1);

    sprintf(lbl, "Z_%s_AE", pert);
    global_dpd_->file2_init(&z1, PSIF_CC_TMP0, irrep, 1, 1, lbl);
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->contract442(&X2, &D, &z1, 2, 2, -1.0, 0.0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&D);
    global_dpd_->contract222(&L1, &z1, &Y1new, 0, 1, 1.0, 1.0);
    global_dpd_->file2_close(&z1);
    global_dpd_->file2_close(&L1);

    // *** <O|L2(0)|[Hbar(0), X2]|0> ***
    
    // Lijab * Xijab -> Lx_IJ //
    global_dpd_->file2_init(&lx, PSIF_CC_OEI, 0, 0, 0, "Lx_IJ");
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract442(&L2, &X2, &lx, 0, 0, 1.0, 0.0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&L2);
    
    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "FME"); 
    global_dpd_->contract222(&lx, &F, &Y1new, 0, 1, -1.0, 1.0);
    global_dpd_->file2_close(&lx);
    global_dpd_->file2_close(&F);  

    // Lijab * Xijab -> Lx_AB 
    global_dpd_->file2_init(&lx, PSIF_CC_OEI, 0, 1, 1, "Lx_AB");
    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract442(&L2, &X2, &lx, 2, 2, 1.0, 0.0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&L2);

    global_dpd_->file2_init(&F, PSIF_CC_OEI, 0, 0, 1, "FME");	
    global_dpd_->contract222(&F, &lx, &Y1new, 0, 0, -1.0, 1.0);
    global_dpd_->file2_close(&lx);
    global_dpd_->file2_close(&F);

    global_dpd_->buf4_init(&lx_iajb, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LXiajb");

    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 10, 10, 10, 10, 0, "(2 LIjAb - LIjBa) (ib|ja)");
    global_dpd_->contract444(&L2, &X2, &lx_iajb, 0, 1, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf"); //Compute this part out of core
    global_dpd_->contract442(&lx_iajb, &W, &Y1new, 0, 3, -1, 1);
    global_dpd_->buf4_close(&lx_iajb);
    global_dpd_->buf4_close(&W);

//MK Optimization It is not working!
/*
    global_dpd_->buf4_init(&lx_ijab, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, "LXijab"); 
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf");
    global_dpd_->buf4_scm(&W, -1);
    // dpd_contract442(&X2, &W, &X1new, 0, 0, 1, 1); /
    // ooc code below added 7/28/05, -TDC /
    global_dpd_->file2_mat_init(&Y1new);
    global_dpd_->file2_mat_rd(&Y1new);
    for (Gam = 0; Gam < moinfo.nirreps; Gam++) {
        Gef = Gam; // W is totally symmetric //
        Gim = Gef ^ irrep;

        global_dpd_->buf4_mat_irrep_init(&lx_ijab, Gim);
        global_dpd_->buf4_mat_irrep_rd(&lx_ijab, Gim);
        global_dpd_->buf4_mat_irrep_shift13(&lx_ijab, Gim);

        for (Gi = 0; Gi < moinfo.nirreps; Gi++) {
            Ga = Gi ^ irrep;
            Gm = Ga ^ Gam;

            W.matrix[Gam] = global_dpd_->dpd_block_matrix(moinfo.occpi[Gm], W.params->coltot[Gef]);

            nrows = moinfo.occpi[Gi];
            ncols = moinfo.occpi[Gm] * W.params->coltot[Gef];

            for (A = 0; A < moinfo.virtpi[Ga]; A++) {
                a = moinfo.vir_off[Ga] + A;
                am = W.row_offset[Gam][a];

                global_dpd_->buf4_mat_irrep_rd_block(&W, Gam, am, moinfo.occpi[Gm]);

                if (nrows && ncols && moinfo.virtpi[Ga])
                    C_DGEMV('n', nrows, ncols, 1, lx_ijab.shift.matrix[Gim][Gi][0], ncols, W.matrix[Gam][0], 1, 1,
                            &(Y1new.matrix[Gi][0][A]), moinfo.virtpi[Ga]);
            }
            global_dpd_->free_dpd_block(W.matrix[Gam], moinfo.occpi[Gm], W.params->coltot[Gef]);
        }

        global_dpd_->buf4_mat_irrep_close(&lx_ijab, Gim);
    }
    global_dpd_->file2_mat_wrt(&Y1new);
    global_dpd_->file2_mat_close(&Y1new);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&lx_ijab);

*/

    // Maybe X2 x self.Hvovv first ? 
    // tmp   =  ndot('mnga,mnef->gaef',self.l2, self.x2)
    // r_y1 -=  ndot('gief,gaef->ia', self.Hvovv, tmp)

    global_dpd_->buf4_init(&lx_iajb, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LXiajb_2");

    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 10, 10, 10, 10, 0, "(2 LIjAb - LIjBa) (ia|jb)"); 
    global_dpd_->contract444(&L2, &X2, &lx_iajb, 0, 1, 1, 0);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&X2);
 
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf (am,fe)"); //Compute this part out of core
    global_dpd_->contract442(&lx_iajb, &W, &Y1new, 0, 3, -1, 1);  
    global_dpd_->buf4_close(&lx_iajb);
    global_dpd_->buf4_close(&W);

    // tmp   =  ndot('mnga,mnef->gaef',self.l2, self.x2)
    // r_y1 -=  ndot('gief,gaef->ia', self.Hvovv, tmp)

    // X2 x  Hvovv
    global_dpd_->buf4_init(&Hx_ijab, PSIF_CC_LR, irrep, 0, 11, 0, 11, 0, "Hx_ijab");

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf"); //Compute this part out of core
    global_dpd_->contract444(&X2, &W, &Hx_ijab, 0, 0, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->contract442(&Hx_ijab, &L2, &Y1new, 3, 3, -1, 1);
    global_dpd_->buf4_close(&Hx_ijab);
    global_dpd_->buf4_close(&L2);

    //tmp   =  ndot('gmae,mnef->ganf',self.Hvovv, self.x2, prefactor=2.0)
    //tmp  +=  ndot('gmea,mnef->ganf',self.Hvovv, self.x2, prefactor=-1.0)

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, irrep, 10, 10, 10, 10, 0, "(2 LIjAb - LIjBa) (ia|jb)");
    sprintf(lbl, "X_%s_IAjb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl); 
    global_dpd_->buf4_init(&lx_iajb, PSIF_CC_LR, 0, 10, 10, 10, 10, 0, "LX (ia,jb)");

    global_dpd_->contract444(&L2, &X2, &lx_iajb, 0, 0, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE) (am,fe)");	
    global_dpd_->contract442(&lx_iajb, &W, &Y1new, 0, 3, 1, 1);
    global_dpd_->buf4_close(&lx_iajb);
    global_dpd_->buf4_close(&W);


    //global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 11, 5, 11, 5, 0, "WAmEf 2(Am,Ef) - (Am,fE)"); //Make it out of core
    //global_dpd_->file2_init(&lx, PSIF_CC_OEI, irrep, 1, 1, "Lx_AB");
    //global_dpd_->dot13(&lx, &W, &Y1new, 0, 0, 1.0, 1.0); 
    //global_dpd_->buf4_close(&W);
    //global_dpd_->file2_close(&lx);    

        /* Above code replaced to remove disk-space and memory bottlenecks */
        global_dpd_->file2_init(&lx, PSIF_CC_OEI, irrep, 1, 1, "Lx_AB");
        global_dpd_->file2_mat_init(&lx);
        global_dpd_->file2_mat_rd(&lx);
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
                    C_DGEMV('t', nrows, ncols, 1, &X[W.col_offset[Gei][Gf]], ncols, lx.matrix[Ge][E], 1, 1,
                            Y1new.matrix[Gi][I], 1);
            }
            global_dpd_->buf4_mat_irrep_row_close(&W, Gei);
            free(X);
        }
        global_dpd_->buf4_close(&W);
        global_dpd_->file2_mat_wrt(&Y1new);
        global_dpd_->file2_mat_close(&Y1new);
        global_dpd_->file2_mat_close(&lx);
        global_dpd_->file2_close(&lx);

    sprintf(lbl, "X_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl); 
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->buf4_init(&lx_ijab, PSIF_CC_LR, irrep, 0, 0, 0, 0, 0, "Lx_ijkl");
    
    global_dpd_->contract444(&L2, &X2, &lx_ijab, 0, 0, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 10, 0, 10, 0, "WMnIe");
    global_dpd_->contract442(&lx_ijab, &W, &Y1new, 1, 3, 1, 1);
    global_dpd_->buf4_close(&lx_ijab);
    global_dpd_->buf4_close(&W); 

    global_dpd_->buf4_init(&lx_iajb, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LXiajb_2"); //I am reusing here
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe (Mn,eI)");

    global_dpd_->contract442(&W, &lx_iajb, &Y1new, 0, 1, 1, 1);
    global_dpd_->buf4_close(&lx_iajb);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&lx_iajb, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LXiajb_3");  

    sprintf(lbl, "X_%s_IbjA (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 10, 10, 10, 10, 0, "(2 LIjAb - LIjBa) (ib|ja)");  
    
    global_dpd_->contract444(&L2, &X2, &lx_iajb, 0, 0, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&L2);

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "WMnIe (nM,eI)");
    global_dpd_->contract442(&W, &lx_iajb, &Y1new, 0, 1, 1, 1);
    global_dpd_->buf4_close(&lx_iajb);
    global_dpd_->buf4_close(&W);

    global_dpd_->file2_init(&lx, PSIF_CC_OEI, irrep, 0, 0, "Lx_IJ");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (Mn,eI)");
    global_dpd_->dot14(&lx, &W, &Y1new, 1, 0, -1.0, 1.0);
    global_dpd_->file2_close(&lx);
    global_dpd_->buf4_close(&W);
   
    global_dpd_->buf4_init(&lx_iajb, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LXiajb");
    
    sprintf(lbl, "X_%s_IAjb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&X2, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, lbl);
    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, irrep, 10, 10, 10, 10, 0, "(2 LIjAb - LIjBa) (ia|jb)");
    global_dpd_->contract444(&L2, &X2, &lx_iajb, 1, 0, 1, 0);
    global_dpd_->buf4_close(&X2);
    global_dpd_->buf4_close(&L2);
    global_dpd_->buf4_close(&lx_iajb);

    global_dpd_->buf4_init(&lx_iajb, PSIF_CC_LR, irrep, 10, 10, 10, 10, 0, "LXiajb");
    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 0, 11, 0, 11, 0, "2WMnIe - WnMIe (Mn,eI)"); 
    global_dpd_->contract442(&W, &lx_iajb, &Y1new, 0, 1, -1, 1);   
    global_dpd_->buf4_close(&lx_iajb);
    global_dpd_->buf4_close(&W);

    /*
    Y1_norm = 0;
    Y1_norm = global_dpd_->file2_dot_self(&Y1new);
    Y1_norm = sqrt(Y1_norm);
    outfile->Printf("\tNorm of the Y1new part1i_Final.... %20.15f\n", Y1_norm);
    */

    global_dpd_->file2_close(&Y1new);

    return;
}

}  // namespace ccresponse
}  // namespace psi
