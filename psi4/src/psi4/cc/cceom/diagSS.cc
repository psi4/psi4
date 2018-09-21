/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2018 The Psi4 Developers.
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
    \ingroup CCEOM
    \brief Enter brief description of file here
*/
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libqt/qt.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace cceom {
#include "psi4/physconst.h"

void sigmaSS(int index, int C_irr);
void cc2_sigmaSS(int index, int C_irr);
void precondition_SS(dpdfile2 *RIA, dpdfile2 *Ria, double eval);
void schmidt_add_SS(dpdfile2 *RIA, dpdfile2 *Ria, int C_irr, int *numCs);
void precondition_SS_RHF(dpdfile2 *RIA, double eval);
void schmidt_add_SS_RHF(dpdfile2 *RIA, int C_irr, int *numCs);
void restart_SS(double **alpha, int L, int num, int C_irr);
void dgeev_eom(int L, double **G, double *evals, double **alpha);
double norm_C1(dpdfile2 *C1A, dpdfile2 *C1B);
double norm_C1_rhf(dpdfile2 *C1A);
void scm_C1(dpdfile2 *C1A, dpdfile2 *C1B, double a);

void diagSS(int C_irr) {
    dpdfile2 Fmi, FMI, Fae, FAE, Fme, FME;
    dpdfile2 CME, Cme, C, SIA, Sia, RIA, Ria;
    dpdbuf4 CMNEF, Cmnef, CMnEf, W;
    char lbl[32], lbl2[32];
    int lwork, info, get_right_ev = 1, get_left_ev = 0;
    int L, h, i, j, k, a, C_index, errcod, keep_going = 1, numCs, iter = 0;
    double norm, tval, *lambda, *lambda_old, zero = 0.0;
    double **G, *work, *evals_complex, **alpha, **evectors_left;
    int nirreps, *openpi, *occpi, *virtpi, range;
    int *aoccpi, *avirtpi, *boccpi, *bvirtpi;
    int *converged, num_converged, num_roots;
    int begin_occ, begin_virt, end_virt, dim_SS = 0;
    int pf, cnt, irr_occ, irr_virt;

    nirreps = moinfo.nirreps;
    openpi = moinfo.openpi;
    occpi = moinfo.occpi;
    virtpi = moinfo.virtpi;

    if (params.eom_ref == 2) { /* UHF */
        aoccpi = moinfo.aoccpi;
        boccpi = moinfo.boccpi;
        avirtpi = moinfo.avirtpi;
        bvirtpi = moinfo.bvirtpi;
    }

    range = eom_params.excitation_range;
    pf = eom_params.print_singles;
    if (pf) outfile->Printf("\n\n");

    /* a bunch of tedious code to setup reasonable HOMO-LUMO guess vectors */
    C_index = 0;
    if (2 * range * range * nirreps < eom_params.cs_per_irrep[C_irr]) range = eom_params.cs_per_irrep[C_irr] / 2;

    if (!params.local) { /* guesses should already be built for local mode */

        for (cnt = 0; cnt < nirreps; ++cnt) { /* cnt loops over dp's to get C_irr */
            irr_occ = dpd_dp[C_irr][cnt][0];
            irr_virt = dpd_dp[C_irr][cnt][1];
            /* C_irr = irr_occ * irr_virt */

            if (params.eom_ref == 0) { /* ref = RHF; eom_ref = RHF */
                begin_occ = MAX(occpi[irr_occ] - range, 0);
                end_virt = MIN(range, virtpi[irr_virt]);
                for (i = begin_occ; i < occpi[irr_occ]; ++i)
                    for (a = 0; a < end_virt; ++a, ++C_index) {
                        sprintf(lbl, "%s %d", "CME", C_index);
                        global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, lbl);
                        global_dpd_->file2_mat_init(&CME);
                        CME.matrix[cnt][i][a] = 1.0 / sqrt(2.0);
                        global_dpd_->file2_mat_wrt(&CME);
                        global_dpd_->file2_mat_close(&CME);
                        global_dpd_->file2_close(&CME);
                    }
            }
            /* eom_ref = ROHF, closed shell */
            else if ((params.eom_ref < 2) && (moinfo.iopen == 0)) {
                begin_occ = MAX(occpi[irr_occ] - range, 0);
                end_virt = MIN(range, virtpi[irr_virt]);
                for (i = begin_occ; i < occpi[irr_occ]; ++i)
                    for (a = 0; a < end_virt; ++a, ++C_index) {
                        sprintf(lbl, "%s %d", "CME", C_index);
                        global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, lbl);
                        global_dpd_->file2_mat_init(&CME);
                        CME.matrix[cnt][i][a] = 1.0 / sqrt(2.0);
                        global_dpd_->file2_mat_wrt(&CME);
                        global_dpd_->file2_mat_close(&CME);
                        global_dpd_->file2_close(&CME);
                        sprintf(lbl, "%s %d", "Cme", C_index);
                        global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, lbl);
                        global_dpd_->file2_mat_init(&Cme);
                        Cme.matrix[cnt][i][a] = 1.0 / sqrt(2.0);
                        global_dpd_->file2_mat_wrt(&Cme);
                        global_dpd_->file2_mat_close(&Cme);
                        global_dpd_->file2_close(&Cme);
                        ++C_index;
                        sprintf(lbl, "%s %d", "CME", C_index);
                        global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, lbl);
                        global_dpd_->file2_mat_init(&CME);
                        CME.matrix[cnt][i][a] = 1.0 / sqrt(2.0);
                        global_dpd_->file2_mat_wrt(&CME);
                        global_dpd_->file2_mat_close(&CME);
                        global_dpd_->file2_close(&CME);
                        sprintf(lbl, "%s %d", "Cme", C_index);
                        global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, lbl);
                        global_dpd_->file2_mat_init(&Cme);
                        Cme.matrix[cnt][i][a] = -1.0 / sqrt(2.0);
                        global_dpd_->file2_mat_wrt(&Cme);
                        global_dpd_->file2_mat_close(&Cme);
                        global_dpd_->file2_close(&Cme);
                    }
            } else if (params.eom_ref == 1) { /* open-shell ROHF */
                /* alpha excitations */
                begin_occ = MAX(occpi[irr_occ] - range, 0);
                end_virt = MIN(virtpi[irr_virt] - openpi[irr_virt], range);
                for (i = begin_occ; i < occpi[irr_occ]; ++i)
                    for (a = 0; a < end_virt; ++a, ++C_index) {
                        sprintf(lbl, "%s %d", "CME", C_index);
                        global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, lbl);
                        global_dpd_->file2_mat_init(&CME);
                        CME.matrix[cnt][i][a] = 1.0;
                        global_dpd_->file2_mat_wrt(&CME);
                        global_dpd_->file2_close(&CME);
                        sprintf(lbl, "%s %d", "Cme", C_index);
                        global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, lbl);
                        global_dpd_->file2_mat_init(&Cme);
                        global_dpd_->file2_mat_wrt(&Cme);
                        global_dpd_->file2_close(&Cme);
                    }
                /* beta excitations into open shells */
                begin_occ = MAX(occpi[irr_occ] - openpi[irr_occ] - range, 0);
                begin_virt = virtpi[irr_virt] - openpi[irr_virt];
                for (i = begin_occ; i < occpi[irr_occ] - openpi[irr_occ]; ++i)
                    for (a = begin_virt; a < virtpi[irr_virt]; ++a, ++C_index) {
                        sprintf(lbl, "%s %d", "Cme", C_index);
                        global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, lbl);
                        global_dpd_->file2_mat_init(&Cme);
                        Cme.matrix[cnt][i][a] = 1.0;
                        global_dpd_->file2_mat_wrt(&Cme);
                        global_dpd_->file2_close(&Cme);
                        sprintf(lbl, "%s %d", "CME", C_index);
                        global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, lbl);
                        global_dpd_->file2_mat_init(&CME);
                        global_dpd_->file2_mat_wrt(&CME);
                        global_dpd_->file2_close(&CME);
                    }
                /* beta excitations into unoccupied orbitals */
                begin_occ = MAX(occpi[irr_occ] - openpi[irr_occ] - range, 0);
                end_virt = MIN(range - openpi[irr_virt], 0);
                for (i = begin_occ; i < occpi[irr_occ] - openpi[irr_occ]; ++i)
                    for (a = 0; a < end_virt; ++a, ++C_index) {
                        sprintf(lbl, "%s %d", "Cme", C_index);
                        global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, lbl);
                        global_dpd_->file2_mat_init(&Cme);
                        Cme.matrix[cnt][i][a] = 1.0;
                        global_dpd_->file2_mat_wrt(&Cme);
                        global_dpd_->file2_close(&Cme);
                        sprintf(lbl, "%s %d", "CME", C_index);
                        global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, lbl);
                        global_dpd_->file2_mat_init(&CME);
                        global_dpd_->file2_mat_wrt(&CME);
                        global_dpd_->file2_close(&CME);
                    }
            } else { /* UHF */
                /* alpha excitations */
                begin_occ = MAX(aoccpi[irr_occ] - range, 0);
                end_virt = MIN(avirtpi[irr_virt], range);
                for (i = begin_occ; i < aoccpi[irr_occ]; ++i)
                    for (a = 0; a < end_virt; ++a, ++C_index) {
                        sprintf(lbl, "%s %d", "CME", C_index);
                        global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, lbl);
                        global_dpd_->file2_mat_init(&CME);
                        CME.matrix[cnt][i][a] = 1.0;
                        global_dpd_->file2_mat_wrt(&CME);
                        global_dpd_->file2_close(&CME);
                        sprintf(lbl, "%s %d", "Cme", C_index);
                        global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 2, 3, lbl);
                        global_dpd_->file2_mat_init(&Cme);
                        global_dpd_->file2_mat_wrt(&Cme);
                        global_dpd_->file2_close(&Cme);
                    }
                begin_occ = MAX(boccpi[irr_occ] - range, 0);
                end_virt = MIN(bvirtpi[irr_virt], range);
                for (i = begin_occ; i < boccpi[irr_occ]; ++i)
                    for (a = 0; a < end_virt; ++a, ++C_index) {
                        sprintf(lbl, "%s %d", "CME", C_index);
                        global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, lbl);
                        global_dpd_->file2_mat_init(&CME);
                        global_dpd_->file2_mat_wrt(&CME);
                        global_dpd_->file2_close(&CME);
                        sprintf(lbl, "%s %d", "Cme", C_index);
                        global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 2, 3, lbl);
                        global_dpd_->file2_mat_init(&Cme);
                        Cme.matrix[cnt][i][a] = 1.0;
                        global_dpd_->file2_mat_wrt(&Cme);
                        global_dpd_->file2_close(&Cme);
                    }
            }
        }

        if (pf) outfile->Printf("%d initial single excitation guesses\n", C_index);
        if (C_index == 0) {
            outfile->Printf("No initial guesses obtained for %s state \n", moinfo.irr_labs[moinfo.sym ^ C_irr].c_str());
            exit(1);
        }
    } else {
        C_index = eom_params.cs_per_irrep[0];
    } /* end if(!params.local) */

    /* show initial guesses */
    /*
    for (i=0; i<C_index ;++i) {
      sprintf(lbl, "%s %d", "CME", i);
      dpd_file2_init(&CME, EOM_CME, C_irr, 0, 1, lbl);
      dpd_file2_print(&CME, outfile);
      dpd_file2_close(&CME);
      if (params.eom_ref > 0) {
        sprintf(lbl, "%s %d", "Cme", i);
        if (params.eom_ref == 1) dpd_file2_init(&Cme, EOM_Cme, C_irr, 0, 1, lbl);
        else if (params.eom_ref == 2) dpd_file2_init(&Cme, EOM_Cme, C_irr, 2, 3, lbl);
        dpd_file2_print(&Cme, outfile);
        dpd_file2_close(&Cme);
      }
    }
    */
    if (eom_params.skip_diagSS) {
        outfile->Printf("\nSkipping diagonalization of Hbar SS block.\n");
        return;
    }

    /* Setup residual vector file */
    global_dpd_->file2_init(&RIA, PSIF_EOM_R, C_irr, 0, 1, "RIA");
    global_dpd_->file2_mat_init(&RIA);
    global_dpd_->file2_mat_wrt(&RIA);
    global_dpd_->file2_close(&RIA);
    if (params.eom_ref > 0) {
        if (params.eom_ref == 1)
            global_dpd_->file2_init(&Ria, PSIF_EOM_R, C_irr, 0, 1, "Ria");
        else if (params.eom_ref == 2)
            global_dpd_->file2_init(&Ria, PSIF_EOM_R, C_irr, 2, 3, "Ria");
        global_dpd_->file2_mat_init(&Ria);
        global_dpd_->file2_mat_wrt(&Ria);
        global_dpd_->file2_close(&Ria);
    }

    /* arrays must be dimensioned with at least the final number of roots - even though
       num_roots may be limited until the first collapse by the number of good
       initial guesses obtained above. */
    L = num_roots = C_index;
    i = MAX(eom_params.cs_per_irrep[C_irr], C_index);
    converged = init_int_array(i);
    lambda_old = init_array(i);

    if (pf) outfile->Printf("\n");
    while ((keep_going == 1) && (iter < eom_params.max_iter_SS)) {
        if (pf) outfile->Printf("Iter=%-4d L=%-4d", iter + 1, L);
        keep_going = 0;
        numCs = L;
        num_converged = 0;

        /* zero S vectors */
        for (i = 0; i < L; ++i) {
            if (params.eom_ref == 0) {
                sprintf(lbl, "%s %d", "SIA", i);
                global_dpd_->file2_init(&SIA, PSIF_EOM_SIA, C_irr, 0, 1, lbl);
                global_dpd_->file2_scm(&SIA, 0.0);
                global_dpd_->file2_close(&SIA);
                if (params.full_matrix) {
                    sprintf(lbl, "%s %d", "S0", i);
                    psio_write_entry(PSIF_EOM_SIA, lbl, (char *)&zero, sizeof(double));
                }
            }
            if (params.eom_ref > 0) {
                sprintf(lbl, "%s %d", "SIA", i);
                global_dpd_->file2_init(&SIA, PSIF_EOM_SIA, C_irr, 0, 1, lbl);
                sprintf(lbl, "%s %d", "Sia", i);
                if (params.eom_ref == 1)
                    global_dpd_->file2_init(&Sia, PSIF_EOM_Sia, C_irr, 0, 1, lbl);
                else if (params.eom_ref == 2)
                    global_dpd_->file2_init(&Sia, PSIF_EOM_Sia, C_irr, 2, 3, lbl);
                scm_C1(&SIA, &Sia, 0.0);
                global_dpd_->file2_close(&SIA);
                global_dpd_->file2_close(&Sia);
            }
        }

        if (params.wfn == "EOM_CC2") {
            for (i = 0; i < L; ++i) cc2_sigmaSS(i, C_irr);
        } else {
            for (i = 0; i < L; ++i) sigmaSS(i, C_irr);
        }

        G = block_matrix(L, L);

        for (i = 0; i < L; ++i) {
            sprintf(lbl, "%s %d", "CME", i);
            global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, lbl);
            if (params.eom_ref > 0) {
                sprintf(lbl, "%s %d", "Cme", i);
                if (params.eom_ref == 1)
                    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, lbl);
                else if (params.eom_ref == 2)
                    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 2, 3, lbl);
            }
            for (j = 0; j < L; ++j) {
                if (params.eom_ref == 0) {
                    sprintf(lbl, "%s %d", "SIA", j);
                    global_dpd_->file2_init(&SIA, PSIF_EOM_SIA, C_irr, 0, 1, lbl);
                    tval = 2.0 * global_dpd_->file2_dot(&CME, &SIA);
                    global_dpd_->file2_close(&SIA);
                } else if (params.eom_ref > 0) {
                    sprintf(lbl, "%s %d", "SIA", j);
                    global_dpd_->file2_init(&SIA, PSIF_EOM_SIA, C_irr, 0, 1, lbl);
                    sprintf(lbl, "%s %d", "Sia", j);
                    if (params.eom_ref == 1)
                        global_dpd_->file2_init(&Sia, PSIF_EOM_Sia, C_irr, 0, 1, lbl);
                    else if (params.eom_ref == 2)
                        global_dpd_->file2_init(&Sia, PSIF_EOM_Sia, C_irr, 2, 3, lbl);
                    tval = global_dpd_->file2_dot(&CME, &SIA);
                    tval += global_dpd_->file2_dot(&Cme, &Sia);
                    global_dpd_->file2_close(&SIA);
                    global_dpd_->file2_close(&Sia);
                }
                G[i][j] = tval;
            }
            global_dpd_->file2_close(&CME);
            if (params.eom_ref > 0) global_dpd_->file2_close(&Cme);
        }

        /* print_mat(G, L, L, outfile); */

        lambda = init_array(L);     /* holds real part of eigenvalues of G */
        alpha = block_matrix(L, L); /* will hold eigenvectors of G */
        dgeev_eom(L, G, lambda, alpha);
        eigsort(lambda, alpha, L);
        /*  eivout(alpha, lambda, L, L, outfile); */
        free_block(G);

        if (pf) outfile->Printf("  Root    EOM Energy    Delta E     Res. Norm    Conv?\n");

        global_dpd_->file2_init(&RIA, PSIF_EOM_R, C_irr, 0, 1, "RIA");
        if (params.eom_ref == 1)
            global_dpd_->file2_init(&Ria, PSIF_EOM_R, C_irr, 0, 1, "Ria");
        else if (params.eom_ref == 2)
            global_dpd_->file2_init(&Ria, PSIF_EOM_R, C_irr, 2, 3, "Ria");

        for (k = 0; k < num_roots; ++k) {
            global_dpd_->file2_scm(&RIA, 0.0);
            if (params.eom_ref > 0) global_dpd_->file2_scm(&Ria, 0.0);
            converged[k] = 0;
            for (i = 0; i < L; ++i) {
                sprintf(lbl, "%s %d", "SIA", i);
                global_dpd_->file2_init(&SIA, PSIF_EOM_SIA, C_irr, 0, 1, lbl);
                sprintf(lbl, "%s %d", "CME", i);
                global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, lbl);
                global_dpd_->file2_axpbycz(&CME, &SIA, &RIA, -1.0 * lambda[k] * alpha[i][k], alpha[i][k], 1.0);
                global_dpd_->file2_close(&CME);
                global_dpd_->file2_close(&SIA);
                if (params.eom_ref > 0) {
                    sprintf(lbl, "%s %d", "Sia", i);
                    if (params.eom_ref == 1)
                        global_dpd_->file2_init(&Sia, PSIF_EOM_Sia, C_irr, 0, 1, lbl);
                    else if (params.eom_ref == 2)
                        global_dpd_->file2_init(&Sia, PSIF_EOM_Sia, C_irr, 2, 3, lbl);
                    sprintf(lbl, "%s %d", "Cme", i);
                    if (params.eom_ref == 1)
                        global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, lbl);
                    else if (params.eom_ref == 2)
                        global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 2, 3, lbl);
                    global_dpd_->file2_axpbycz(&Cme, &Sia, &Ria, -1.0 * lambda[k] * alpha[i][k], alpha[i][k], 1.0);
                    global_dpd_->file2_close(&Cme);
                    global_dpd_->file2_close(&Sia);
                }
            }

            if (params.eom_ref == 0)
                precondition_SS_RHF(&RIA, lambda[k]);
            else
                precondition_SS(&RIA, &Ria, lambda[k]);

            if (params.eom_ref == 0)
                norm = norm_C1_rhf(&RIA);
            else
                norm = norm_C1(&RIA, &Ria);

            if (pf) outfile->Printf("%6d%15.10lf%11.2e%12.2e", k + 1, lambda[k], lambda[k] - lambda_old[k], norm);

            if ((norm > eom_params.residual_tol_SS) ||
                (std::fabs(lambda[k] - lambda_old[k]) > eom_params.eval_tol_SS)) {
                if (pf) outfile->Printf("%7s\n", "N");
                /*
                  if (params.eom_ref == 0) precondition_SS_RHF(&RIA, lambda[k]);
                  else precondition_SS(&RIA, &Ria, lambda[k]);
                */

                /* normalize preconditioned residual */
                if (params.eom_ref == 0) {
                    norm = norm_C1_rhf(&RIA);
                    global_dpd_->file2_scm(&RIA, 1.0 / norm);
                } else {
                    norm = norm_C1(&RIA, &Ria);
                    global_dpd_->file2_scm(&RIA, 1.0 / norm);
                    global_dpd_->file2_scm(&Ria, 1.0 / norm);
                }

                if (params.eom_ref == 0)
                    schmidt_add_SS_RHF(&RIA, C_irr, &numCs);
                else
                    schmidt_add_SS(&RIA, &Ria, C_irr, &numCs);
            } else {
                if (pf) outfile->Printf("%7s\n", "Y");
                ++num_converged;
                converged[k] = 1;
            }
        }

        global_dpd_->file2_close(&RIA);
        if (params.eom_ref > 0) global_dpd_->file2_close(&Ria);

        for (i = 0; i < num_roots; ++i) lambda_old[i] = lambda[i];
        free(lambda);

        if ((iter == 2) || (L > eom_params.vectors_per_root_SS * eom_params.cs_per_irrep[C_irr])) {
            restart_SS(alpha, L, eom_params.cs_per_irrep[C_irr], C_irr);
            L = num_roots = eom_params.cs_per_irrep[C_irr];
            keep_going = 1;
            free_block(alpha);
        } else {
            if (numCs > L) { /* successfully added vectors */
                keep_going = 1;
                free_block(alpha);
                L = numCs;
            }
        }
        ++iter;
    }

    if (pf) outfile->Printf("\nLowest eigenvalues of HBar Singles-Singles Block\n");
    if (pf) outfile->Printf("Root      Excitation Energy         Total Energy\n");
    if (pf) outfile->Printf("           (eV)     (cm^-1)            (au)\n");
    if (pf)
        for (i = 0; i < num_roots; ++i)
            if (pf)
                outfile->Printf("%4d%12.3lf%12.2lf%20.10lf\n", i + 1, lambda_old[i] * pc_hartree2ev,
                                lambda_old[i] * pc_hartree2wavenumbers, lambda_old[i] + moinfo.eref + moinfo.ecc);

    free(lambda_old);
    free(converged);
    /* collapse solutions to one vector each */
    restart_SS(alpha, L, eom_params.cs_per_irrep[C_irr], C_irr);
    free_block(alpha);
    if (pf) outfile->Printf("\n");

    return;
}

void precondition_SS(dpdfile2 *RIA, dpdfile2 *Ria, double eval) {
    dpdfile2 DIA, Dia;
    int h, nirreps, i, j, a, b, ij, ab, C_irr;
    double tval;

    C_irr = RIA->my_irrep;
    nirreps = RIA->params->nirreps;

    global_dpd_->file2_mat_init(RIA);
    global_dpd_->file2_mat_rd(RIA);
    global_dpd_->file2_init(&DIA, PSIF_EOM_D, C_irr, 0, 1, "DIA");
    global_dpd_->file2_mat_init(&DIA);
    global_dpd_->file2_mat_rd(&DIA);
    for (h = 0; h < nirreps; h++)
        for (i = 0; i < RIA->params->rowtot[h]; i++)
            for (a = 0; a < RIA->params->coltot[h ^ C_irr]; a++) {
                tval = eval - DIA.matrix[h][i][a];
                if (std::fabs(tval) > 0.0001) RIA->matrix[h][i][a] /= tval;
            }
    global_dpd_->file2_mat_wrt(RIA);
    global_dpd_->file2_mat_close(RIA);
    global_dpd_->file2_close(&DIA);

    global_dpd_->file2_mat_init(Ria);
    global_dpd_->file2_mat_rd(Ria);

    if (params.eom_ref == 1)
        global_dpd_->file2_init(&Dia, PSIF_EOM_D, C_irr, 0, 1, "Dia");
    else if (params.eom_ref == 2)
        global_dpd_->file2_init(&Dia, PSIF_EOM_D, C_irr, 2, 3, "Dia");
    global_dpd_->file2_mat_init(&Dia);
    global_dpd_->file2_mat_rd(&Dia);
    for (h = 0; h < nirreps; h++)
        for (i = 0; i < Ria->params->rowtot[h]; i++)
            for (a = 0; a < Ria->params->coltot[h ^ C_irr]; a++) {
                tval = eval - Dia.matrix[h][i][a];
                if (std::fabs(tval) > 0.0001) Ria->matrix[h][i][a] /= tval;
            }
    global_dpd_->file2_mat_wrt(Ria);
    global_dpd_->file2_mat_close(Ria);
    global_dpd_->file2_close(&Dia);

    return;
}

void precondition_SS_RHF(dpdfile2 *RIA, double eval) {
    dpdfile2 DIA;
    int C_irr, h, nirreps, i, j, a, b, ij, ii, ab;
    double tval;
    psio_address next;

    /* Local correlation decs */
    int nso, nocc, nvir;
    double *T1tilde, *T1bar;

    nirreps = RIA->params->nirreps;
    C_irr = RIA->my_irrep;

    if (params.local && local.filter_singles) {
        nso = local.nso;
        nocc = local.nocc;
        nvir = local.nvir;

        local.pairdom_len = init_int_array(nocc * nocc);
        local.pairdom_nrlen = init_int_array(nocc * nocc);
        local.eps_occ = init_array(nocc);
        psio_read_entry(PSIF_CC_INFO, "Local Pair Domain Length", (char *)local.pairdom_len, nocc * nocc * sizeof(int));
        psio_read_entry(PSIF_CC_INFO, "Local Pair Domain Length (Non-redundant basis)", (char *)local.pairdom_nrlen,
                        nocc * nocc * sizeof(int));
        psio_read_entry(PSIF_CC_INFO, "Local Occupied Orbital Energies", (char *)local.eps_occ, nocc * sizeof(double));
        local.W = (double ***)malloc(nocc * nocc * sizeof(double **));
        local.V = (double ***)malloc(nocc * nocc * sizeof(double **));
        local.eps_vir = (double **)malloc(nocc * nocc * sizeof(double *));
        next = PSIO_ZERO;
        for (ij = 0; ij < nocc * nocc; ij++) {
            local.eps_vir[ij] = init_array(local.pairdom_nrlen[ij]);
            psio_read(PSIF_CC_INFO, "Local Virtual Orbital Energies", (char *)local.eps_vir[ij],
                      local.pairdom_nrlen[ij] * sizeof(double), next, &next);
        }
        next = PSIO_ZERO;
        for (ij = 0; ij < nocc * nocc; ij++) {
            local.V[ij] = block_matrix(nvir, local.pairdom_len[ij]);
            psio_read(PSIF_CC_INFO, "Local Residual Vector (V)", (char *)local.V[ij][0],
                      nvir * local.pairdom_len[ij] * sizeof(double), next, &next);
        }
        next = PSIO_ZERO;
        for (ij = 0; ij < nocc * nocc; ij++) {
            local.W[ij] = block_matrix(local.pairdom_len[ij], local.pairdom_nrlen[ij]);
            psio_read(PSIF_CC_INFO, "Local Transformation Matrix (W)", (char *)local.W[ij][0],
                      local.pairdom_len[ij] * local.pairdom_nrlen[ij] * sizeof(double), next, &next);
        }

        global_dpd_->file2_mat_init(RIA);
        global_dpd_->file2_mat_rd(RIA);

        for (i = 0; i < nocc; i++) {
            ii = i * nocc + i;

            if (!local.pairdom_len[ii]) {
                outfile->Printf("\n\tlocal_filter_T1: Pair ii = [%d] is zero-length, which makes no sense.\n", ii);
                exit(2);
            }

            T1tilde = init_array(local.pairdom_len[ii]);
            T1bar = init_array(local.pairdom_nrlen[ii]);

            /* Transform the virtuals to the redundant projected virtual basis */
            C_DGEMV('t', nvir, local.pairdom_len[ii], 1.0, &(local.V[ii][0][0]), local.pairdom_len[ii],
                    &(RIA->matrix[0][i][0]), 1, 0.0, &(T1tilde[0]), 1);

            /* Transform the virtuals to the non-redundant virtual basis */
            C_DGEMV('t', local.pairdom_len[ii], local.pairdom_nrlen[ii], 1.0, &(local.W[ii][0][0]),
                    local.pairdom_nrlen[ii], &(T1tilde[0]), 1, 0.0, &(T1bar[0]), 1);

            for (a = 0; a < local.pairdom_nrlen[ii]; a++) {
                tval = eval + local.eps_occ[i] - local.eps_vir[ii][a];
                if (std::fabs(tval) > 0.0001) T1bar[a] /= tval;
                /* else T1bar[a] = 0.0; */
            }

            /* Transform the new T1's to the redundant projected virtual basis */
            C_DGEMV('n', local.pairdom_len[ii], local.pairdom_nrlen[ii], 1.0, &(local.W[ii][0][0]),
                    local.pairdom_nrlen[ii], &(T1bar[0]), 1, 0.0, &(T1tilde[0]), 1);

            /* Transform the new T1's to the MO basis */
            C_DGEMV('n', nvir, local.pairdom_len[ii], 1.0, &(local.V[ii][0][0]), local.pairdom_len[ii], &(T1tilde[0]),
                    1, 0.0, &(RIA->matrix[0][i][0]), 1);

            free(T1bar);
            free(T1tilde);
        }

        global_dpd_->file2_mat_wrt(RIA);
        global_dpd_->file2_mat_close(RIA);

        /* Free Local Memory */
        for (i = 0; i < nocc * nocc; i++) {
            free_block(local.W[i]);
            free_block(local.V[i]);
            free(local.eps_vir[i]);
        }
        free(local.W);
        free(local.V);
        free(local.eps_vir);

        free(local.eps_occ);
        free(local.pairdom_len);
        free(local.pairdom_nrlen);
    } else {
        global_dpd_->file2_mat_init(RIA);
        global_dpd_->file2_mat_rd(RIA);
        global_dpd_->file2_init(&DIA, PSIF_EOM_D, C_irr, 0, 1, "DIA");
        global_dpd_->file2_mat_init(&DIA);
        global_dpd_->file2_mat_rd(&DIA);
        for (h = 0; h < nirreps; h++)
            for (i = 0; i < RIA->params->rowtot[h]; i++)
                for (a = 0; a < RIA->params->coltot[h ^ C_irr]; a++) {
                    tval = eval - DIA.matrix[h][i][a];
                    if (std::fabs(tval) > 0.0001) RIA->matrix[h][i][a] /= tval;
                }
        global_dpd_->file2_mat_wrt(RIA);
        global_dpd_->file2_mat_close(RIA);
        global_dpd_->file2_close(&DIA);
    }

    return;
}

void schmidt_add_SS(dpdfile2 *RIA, dpdfile2 *Ria, int C_irr, int *numCs) {
    double dotval, norm;
    int i, I;
    dpdfile2 Cme, CME;
    char CME_lbl[32], Cme_lbl[32], CMNEF_lbl[32], Cmnef_lbl[32], CMnEf_lbl[32];

    for (i = 0; i < *numCs; i++) {
        sprintf(CME_lbl, "%s %d", "CME", i);
        global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
        sprintf(Cme_lbl, "%s %d", "Cme", i);
        if (params.eom_ref == 1)
            global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, Cme_lbl);
        else if (params.eom_ref == 2)
            global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 2, 3, Cme_lbl);
        dotval = global_dpd_->file2_dot(RIA, &CME);
        dotval += global_dpd_->file2_dot(Ria, &Cme);
        global_dpd_->file2_axpy(&CME, RIA, -1.0 * dotval, 0);
        global_dpd_->file2_axpy(&Cme, Ria, -1.0 * dotval, 0);
        global_dpd_->file2_close(&CME);
        global_dpd_->file2_close(&Cme);
    }

    norm = norm_C1(RIA, Ria);

    if (norm < eom_params.schmidt_add_residual_tol) {
        return;
    } else {
        scm_C1(RIA, Ria, 1.0 / norm);
        sprintf(CME_lbl, "%s %d", "CME", *numCs);
        sprintf(Cme_lbl, "%s %d", "Cme", *numCs);
        global_dpd_->file2_copy(RIA, PSIF_EOM_CME, CME_lbl);
        global_dpd_->file2_copy(Ria, PSIF_EOM_Cme, Cme_lbl);
        ++(*numCs);
    }
    return;
}

void schmidt_add_SS_RHF(dpdfile2 *RIA, int C_irr, int *numCs) {
    double dotval, norm;
    int i, I;
    dpdfile2 CME;
    char CME_lbl[32], Cme_lbl[32];

    for (i = 0; i < *numCs; i++) {
        sprintf(CME_lbl, "%s %d", "CME", i);
        global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, CME_lbl);
        dotval = 2.0 * global_dpd_->file2_dot(RIA, &CME);
        global_dpd_->file2_axpy(&CME, RIA, -1.0 * dotval, 0);
        global_dpd_->file2_close(&CME);
    }

    norm = norm_C1_rhf(RIA);

    if (norm < eom_params.schmidt_add_residual_tol) {
        return;
    } else {
        global_dpd_->file2_scm(RIA, 1.0 / norm);
        sprintf(CME_lbl, "%s %d", "CME", *numCs);
        global_dpd_->file2_copy(RIA, PSIF_EOM_CME, CME_lbl);
        ++(*numCs);
    }
    return;
}

void restart_SS(double **alpha, int L, int num, int C_irr) {
    int i, j, I;
    char lbl[20];
    dpdfile2 C1, CME, Cme, CME2, Cme2;
    dpdbuf4 C2, CMNEF, Cmnef, CMnEf;
    double norm, dotval;

    for (I = 1; I < num; ++I) {
        for (i = 0; i < I; i++) {
            dotval = 0.0;
            for (j = 0; j < L; ++j) {
                dotval += alpha[j][i] * alpha[j][I];
            }
            for (j = 0; j < L; j++) alpha[j][I] -= dotval * alpha[j][i];
        }
        dotval = 0.0;
        for (j = 0; j < L; ++j) dotval += alpha[j][I] * alpha[j][I];
        norm = sqrt(dotval);
        for (j = 0; j < L; ++j) alpha[j][I] = alpha[j][I] / norm;
    }

    for (i = 0; i < num; ++i) {
        sprintf(lbl, "%s %d", "CME", L + i);
        global_dpd_->file2_init(&C1, PSIF_EOM_CME, C_irr, 0, 1, lbl);
        global_dpd_->file2_scm(&C1, 0.0);
        for (j = 0; j < L; ++j) {
            sprintf(lbl, "%s %d", "CME", j);
            global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, lbl);
            global_dpd_->file2_axpy(&CME, &C1, alpha[j][i], 0);
            global_dpd_->file2_close(&CME);
        }
        global_dpd_->file2_close(&C1);

        if (params.eom_ref > 0) {
            sprintf(lbl, "%s %d", "Cme", L + i);
            if (params.eom_ref == 1)
                global_dpd_->file2_init(&C1, PSIF_EOM_Cme, C_irr, 0, 1, lbl);
            else if (params.eom_ref == 2)
                global_dpd_->file2_init(&C1, PSIF_EOM_Cme, C_irr, 2, 3, lbl);
            global_dpd_->file2_scm(&C1, 0.0);
            for (j = 0; j < L; ++j) {
                sprintf(lbl, "%s %d", "Cme", j);
                if (params.eom_ref == 1)
                    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, lbl);
                else if (params.eom_ref == 2)
                    global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 2, 3, lbl);
                global_dpd_->file2_axpy(&Cme, &C1, alpha[j][i], 0);
                global_dpd_->file2_close(&Cme);
            }
            global_dpd_->file2_close(&C1);
        }
    }

    for (i = 0; i < num; ++i) {
        sprintf(lbl, "%s %d", "CME", L + i);
        global_dpd_->file2_init(&CME, PSIF_EOM_CME, C_irr, 0, 1, lbl);
        sprintf(lbl, "%s %d", "CME", i);
        global_dpd_->file2_copy(&CME, PSIF_EOM_CME, lbl);
        global_dpd_->file2_close(&CME);
        if (params.eom_ref > 0) {
            sprintf(lbl, "%s %d", "Cme", L + i);
            if (params.eom_ref == 1)
                global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 0, 1, lbl);
            else if (params.eom_ref == 2)
                global_dpd_->file2_init(&Cme, PSIF_EOM_Cme, C_irr, 2, 3, lbl);
            sprintf(lbl, "%s %d", "Cme", i);
            global_dpd_->file2_copy(&Cme, PSIF_EOM_Cme, lbl);
            global_dpd_->file2_close(&Cme);
        }
    }
    return;
}
}
}  // namespace psi
