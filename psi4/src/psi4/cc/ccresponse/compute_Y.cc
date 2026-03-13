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
#include "psi4/libpsi4util/exception.h"
#include "psi4/psifiles.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace ccresponse {

void init_X(const char *pert, int irrep, double omega);
void sort_X(const char *pert, int irrep, double omega);
void cc2_sort_X(const char *pert, int irrep, double omega);
void X1_build(const char *pert, int irrep, double omega);
void X2_build(const char *pert, int irrep, double omega);
void cc2_X1_build(const char *pert, int irrep, double omega);
void cc2_X2_build(const char *pert, int irrep, double omega);
//double converged(const char *pert, int irrep, double omega);
void save_X(const char *pert, int irrep, double omega);
void print_X(const char *pert, int irrep, double omega);
void update_X(const char *pert, int irrep, double omega);
//void diis(int iter, const char *pert, int irrep, double omega);
double pseudopolar(const char *pert, int irrep, double omega);
void cleanup();
void exit_io();
void amp_write(const char *pert, int irrep, double omega);

void analyze(const char *pert, int irrep, double omega);
//Added
double pseudopolar_Y(const char *pert, int irrep, double omega);
void init_Y(const char *pert, int irrep, double omega);
void sort_Y(const char *pert, int irrep, double omega);
void Y1_inhomogenous_build(const char *pert, int irrep, double omega);
void Y2_inhomogenous_build(const char *pert, int irrep, double omega);
void Y1_homogenous_build(const char *pert, int irrep, double omega);
void Y2_homogenous_build(const char *pert, int irrep, double omega);
//void Y1_build(const char *pert, int irrep, double omega);
void save_Y(const char *pert, int irrep, double omega);
void update_Y(const char *pert, int irrep, double omega);
void G_build(const char *pert, int irrep, double omega);
double converged_Y(const char *pert, int irrep, double omega);
void diis_Y(int iter, const char *pert, int irrep, double omega);

void lambda_residuals();

void compute_Y(const char *pert, int irrep, double omega) {
    int i, iter = 0, done = 0;
    double rms, polar, X2_norm;
    char lbl[32];
    //----------------
    dpdfile2 X1, Y1, L1, GMI, GAE;
    dpdbuf4 z2, X2, Y2, L2, WL;
    double Y1_norm,Y2_norm;
    dpdfile2 Y1new;  //Just for a test
    

    timer_on("compute_Y");

    outfile->Printf("\n\tComputing Y amplitudes %s-Perturbed Wave Function (%5.3f E_h).\n", pert, omega);
    outfile->Printf("\tIter   Pseudopolarizability       RMS \n");
    outfile->Printf("\t----   --------------------   -----------\n");


//#####################################################################
/*
    global_dpd_->file2_init(&L1, PSIF_CC_LAMPS, 0, 0, 1, "LIA 0 -1");

    global_dpd_->file2_scm(&L1, 2);

    Y1_norm = 0;
    Y1_norm = global_dpd_->file2_dot_self(&L1);
    Y1_norm = sqrt(Y1_norm);
    outfile->Printf("\tNorm of L1.... %20.15f\n", Y1_norm);


    global_dpd_->file2_close(&L1);

    global_dpd_->buf4_init(&L2, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "2 LIjAb - LIjBa");
    global_dpd_->buf4_scm(&L2, 2);

    Y1_norm = 0;
    Y1_norm = global_dpd_->buf4_dot_self(&L2);
    Y1_norm = sqrt(Y1_norm);
    outfile->Printf("\tNorm of L2.... %20.15f\n", Y1_norm);

    global_dpd_->buf4_close(&L2);
*/

/*
    global_dpd_->buf4_init(&WL, PSIF_CC_LAMPS, 0, 0, 5, 0, 5, 0, "WefabL2");
    global_dpd_->buf4_scm(&WL, 2); //I multiplied by 2
    global_dpd_->buf4_close(&WL);

    global_dpd_->file2_init(&GAE, PSIF_CC_LAMBDA, 0, 1, 1, "GAE");
    global_dpd_->file2_scm(&GAE, 2); //I multiplied by 2
    global_dpd_->file2_close(&GAE);

    global_dpd_->file2_init(&GMI, PSIF_CC_LAMBDA, 0, 0, 0, "GMI");
    global_dpd_->file2_scm(&GMI, 2); //I multiplied by 2
    global_dpd_->file2_close(&GMI);

*/
//###################################################################


    //if (params.wfn == "CC2")
    //    cc2_sort_X(pert, irrep, omega);
    //else
    //sort_X(pert, irrep, omega);

    //Init Y1 and Y2
    init_Y(pert, irrep, omega);

    sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    Y2_norm = global_dpd_->buf4_dot_self(&Y2);
    global_dpd_->buf4_close(&Y2);
    Y2_norm = sqrt(Y2_norm);
    outfile->Printf("\tBefore sort ..Norm of the guessed Y2 amplitudes %20.15f\n", Y2_norm);




    sort_Y(pert, irrep, omega);

    sprintf(lbl, "Y_%s_IA (%5.3f)", pert, omega);
    global_dpd_->file2_init(&Y1, PSIF_CC_OEI, irrep, 0, 1, lbl);
    Y1_norm = global_dpd_->file2_dot_self(&Y1);
    global_dpd_->file2_close(&Y1);
    Y1_norm = sqrt(Y1_norm);
    outfile->Printf("\tNorm of the guessed Y1 amplitudes %20.15f\n", Y1_norm);

    sprintf(lbl, "Y_%s_(2IjAb-IjbA) (%5.3f)", pert, omega);
    global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
    Y2_norm = global_dpd_->buf4_dot_self(&Y2);
    global_dpd_->buf4_close(&Y2);
    Y2_norm = sqrt(Y2_norm);
    outfile->Printf("\tNorm of the guessed Y2 amplitudes %20.15f\n", Y2_norm); 

    lambda_residuals();

    polar = -2.0 * pseudopolar_Y(pert, irrep, omega);
    outfile->Printf("\t%4d   %20.12f\n", iter, polar);


    timer_on("Y1 inhomo");
    Y1_inhomogenous_build(pert, irrep, omega);
    timer_off("Y1 inhomo");
    timer_on("Y2 inhomo"); 
    Y2_inhomogenous_build(pert, irrep, omega);
    timer_off("Y2 inhomo");


    for (iter = 1; iter <= params.maxiter; iter++) {
        if (params.wfn == "CC2") {
            cc2_sort_X(pert, irrep, omega);
            cc2_X1_build(pert, irrep, omega);
            cc2_X2_build(pert, irrep, omega);
        } else {
            sort_Y(pert, irrep, omega);
 	    G_build(pert, irrep, omega);	
            Y1_homogenous_build(pert, irrep, omega);
            Y2_homogenous_build(pert, irrep, omega);
            //Y1_inhomogenous_build(pert, irrep, omega);
            //Y1_build(pert, irrep, omega);
            //X1_build(pert, irrep, omega);
            //X2_build(pert, irrep, omega);

           //Begin test
           Y1_norm = 0;
           sprintf(lbl, "New Y_%s_IA (%5.3f)", pert, omega);
           global_dpd_->file2_init(&Y1new, PSIF_CC_OEI, irrep, 0, 1, lbl);
           Y1_norm = global_dpd_->file2_dot_self(&Y1new);
           global_dpd_->file2_close(&Y1new);
           Y1_norm = sqrt(Y1_norm);
           outfile->Printf("\n\n\tNorm of the final Y1 amplitudes...... %20.15f\n", Y1_norm);

           sprintf(lbl, "New Y_%s_IjAb (%5.3f)", pert, omega);
           global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
           Y2_norm = global_dpd_->buf4_dot_self(&Y2);
           global_dpd_->buf4_close(&Y2);
           Y2_norm = sqrt(Y2_norm);
           outfile->Printf("\tNorm of the final Y2 amplitudes...... %20.15f\n", Y2_norm);
           //END TEST

        }
        update_Y(pert, irrep, omega);

           //Begin test
           Y1_norm = 0;
           sprintf(lbl, "New Y_%s_IA (%5.3f)", pert, omega);
           global_dpd_->file2_init(&Y1new, PSIF_CC_OEI, irrep, 0, 1, lbl);
           Y1_norm = global_dpd_->file2_dot_self(&Y1new);
           global_dpd_->file2_close(&Y1new);
           Y1_norm = sqrt(Y1_norm);
           outfile->Printf("\tNorm of the updated Y1 amplitudes %20.15f\n", Y1_norm);

           sprintf(lbl, "New Y_%s_IjAb (%5.3f)", pert, omega);
           global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
           Y2_norm = global_dpd_->buf4_dot_self(&Y2);
           global_dpd_->buf4_close(&Y2);
           Y2_norm = sqrt(Y2_norm);
           outfile->Printf("\tNorm of the updated Y2 amplitudes...... %20.15f\n", Y2_norm);

           //END TEST

        rms = converged_Y(pert, irrep, omega);
        if (rms <= params.convergence) {
            done = 1;
            save_Y(pert, irrep, omega);
            if (params.wfn == "CC2")
                cc2_sort_X(pert, irrep, omega);
            else
                sort_Y(pert, irrep, omega);
            outfile->Printf("\t-----------------------------------------\n");
            outfile->Printf("\t Left hand perturbed amplitudes %s-Perturbed Wfn converged to %4.3e\n", pert, rms);
            if (params.print & 2) {
                sprintf(lbl, "Y_%s_IjAb (%5.3f)", pert, omega);
                global_dpd_->buf4_init(&Y2, PSIF_CC_LR, irrep, 0, 5, 0, 5, 0, lbl);
                Y2_norm = global_dpd_->buf4_dot_self(&Y2);
                global_dpd_->buf4_close(&Y2);
                Y2_norm = sqrt(Y2_norm);
                outfile->Printf("\tNorm of the converged Y2 amplitudes %20.15f\n", Y2_norm);
                amp_write(pert, irrep, omega);
            }

            break;
        }
        if (params.diis) diis_Y(iter, pert, irrep, omega);
        save_Y(pert, irrep, omega);
        if (params.wfn == "CC2")
            cc2_sort_X(pert, irrep, omega);
        else
            sort_Y(pert, irrep, omega);

        polar = -2.0 * pseudopolar_Y(pert, irrep, omega);
        outfile->Printf("\t%4d   %20.12f    %4.3e\n", iter, polar, rms);
    }
    if (!done) {
        outfile->Printf("\tNOT CONVERGED!!!!!!!!!!!!!!!!!!!!!!!!!!\n");
/*
        dpd_close(0);
        cleanup();
        exit_io();
        throw PsiException("Failed to converge perturbed wavefunction", __FILE__, __LINE__);
*/
    }

    /* Clean up disk space */
    psio_close(PSIF_CC_DIIS_AMP, 0);
    psio_close(PSIF_CC_DIIS_ERR, 0);

    psio_open(PSIF_CC_DIIS_AMP, 0);
    psio_open(PSIF_CC_DIIS_ERR, 0);

    for (i = PSIF_CC_TMP; i <= PSIF_CC_TMP11; i++) {
        psio_close(i, 0);
        psio_open(i, 0);
    }

    if (params.analyze) analyze(pert, irrep, omega);

    /*  print_X(pert, irrep, omega); */

    timer_off("compute_Y");
}

}  // namespace ccresponse
}  // namespace psi
