/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2021 The Psi4 Developers.
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
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/

#include "Local.h"
#include "Params.h"
#include "MOInfo.h"
#include "psi4/cc/ccwave.h"

#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libiwl/iwl.h"
#include "psi4/libqt/qt.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/psifiles.h"
#include "psi4/libpsi4util/process.h"
#include "psi4/libmints/mintshelper.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

namespace psi {

//Constructor attempt
Local_cc::Local_cc() {
};

void Local_cc::local_init() {
    /*if (method == "PNO")
        init_pno();
    if (method == "PNO++")
        init_pnopp(basis);
    outfile->Printf("    Localization parameters ready.\n\n");*/

}

void Local_cc::init_pno() {
    outfile->Printf(" Local correlation using Pair Natural Orbitals\nCutoff value: %e \n", cutoff);
    // TODO weak pairs here
    npairs = nocc*nocc;

    dpdbuf4 T2;
    dpdbuf4 T2tilde;
    std::vector<SharedMatrix> Tij;
    std::vector<SharedMatrix> Ttij;
    SharedMatrix temp(new Matrix(nvir, nvir));
    
    std::vector<SharedMatrix> Q;

    // Check if occupied block of Fock matrix is non-diagonal (localized)
    // On the way, store the diagonal Fock matrix elements
    dpdfile2 Fij;
    global_dpd_->file2_init(&Fij, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    global_dpd_->file2_mat_init(&Fij);
    global_dpd_->file2_mat_rd(&Fij);

    double *temp_occ_array = new double[nocc];
    for(int i=0; i < nocc; i++) {
        temp_occ_array[i] = Fij.matrix[0][i][i];
    }
    global_dpd_->file2_close(&Fij);
    psio_write_entry(PSIF_CC_INFO, "Local Occupied Orbital Energies", (char *) temp_occ_array,
            nocc * sizeof(double));

    // "Vectorize" T2
    psio_tocprint(PSIF_CC_TAMPS);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS,  0, 0, 5, 0, 5, 0, "tIjAb");
    // Create T2~
    get_matvec(&T2, &Tij);
    global_dpd_->buf4_scmcopy(&T2, PSIF_CC_TMP0, "tIjAb ~", 2);
    global_dpd_->buf4_sort_axpy(&T2, PSIF_CC_TMP0, pqsr, 0, 5, "tIjAb ~", -1);
    global_dpd_->buf4_close(&T2);

    // "Vectorize" T2~
    global_dpd_->buf4_init(&T2tilde, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "tIjAb ~");
    get_matvec(&T2tilde, &Ttij);
    global_dpd_->buf4_close(&T2tilde);

    // Print check
    /* outfile->Printf("*** IJ Matrices ***\n");
    for(int ij=0; ij<Ttij.size(); ++ij)
        Ttij.at(ij)->print(); */

    // Create Density
    std::vector<SharedMatrix> Dij;
    for(int ij=0; ij < npairs; ++ij) {
        temp->zero();
        int i = ij/nocc;
        int j = ij%nocc;
        temp->gemm(0, 1, 1, Tij[ij], Ttij[ij], 0);
        temp->gemm(1, 0, 1, Tij[ij], Ttij[ij], 1);
        temp->scale(2.0 * 1/(1+(i==j)));
        temp->axpy(1.0, temp->transpose());
        temp->scale(0.5);
        Dij.push_back(temp->clone());
    }

    // Print checking
    /*outfile->Printf("\t***** Dij Matrix *****\n");
    for(int ij=0; ij < nocc*nocc; ++ij) {
        outfile->Printf("ij = %2d\n", ij);
        Dij[ij]->print();
    }*/

    // Diagonalize density
    Q = build_PNO_lists(cutoff, Dij);

    // Get semicanonical transforms
    get_semicanonical_transforms(Q);

    // Check if they can be read back in
    /*eps_pno.clear();
    next = PSIO_ZERO;
    for (int ij = 0; ij < npairs; ++ij) {
        int npno = survivor_list[ij];
        //auto eps_pno_temp = std::make_shared<Vector>(npno);
        SharedVector eps_pno_temp(new Vector(npno));
        psio_read(PSIF_CC_INFO, "Local Virtual Orbital Energies", (char *) eps_pno_temp->pointer(),
                  npno * sizeof(double), next, &next);
        eps_pno.push_back(eps_pno_temp);
    }
    outfile->Printf("Testing read-in of virtual orbital energies\n");
    for(int ij=0; ij < npairs; ++ij) {
        outfile->Printf("\nPair %d ", ij);
        int npno = survivor_list[ij];
        for(int a=0; a < npno; ++a) {
            outfile->Printf("%10.10f\t", eps_pno[ij]->get(a));
        }
    }
    next = PSIO_ZERO;
    outfile->Printf("**** Read-in Truncated Q ****\n");
    for(int ij=0; ij < npairs; ++ij) {
        int npno = survivor_list[ij];
        auto Ltest = std::make_shared<Matrix>(nvir,npno);
        psio_read(PSIF_CC_INFO, "Local Transformation Matrix Q", (char *) Ltest->pointer()[0],
        nvir * npno * sizeof(double), next, &next);
        Ltest->print();
        Ltest->zero();
    }*/

    // Assuming this memory needs to be freed
    delete[] temp_occ_array;
}

void Local_cc::init_pnopp(const double omega) {
    outfile->Printf(" Local correlation using Perturbed Pair Natural Orbitals\nCutoff value: %e\n", cutoff);
    // TODO weak pairs here
    npairs = nocc*nocc;

    std::vector<SharedMatrix> Xij;
    std::vector<SharedMatrix> Xtij;

    std::vector<SharedMatrix> Q;

    // Check if occupied block of Fock matrix is non-diagonal (localized)
    // On the way, store the diagonal Fock matrix elements
    dpdfile2 Fij;
    global_dpd_->file2_init(&Fij, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    global_dpd_->file2_mat_init(&Fij);
    global_dpd_->file2_mat_rd(&Fij);

    double *temp_occ_array = new double[nocc];
    for(int i=0; i < nocc; i++) {
        temp_occ_array[i] = Fij.matrix[0][i][i];
    }
    psio_write_entry(PSIF_CC_INFO, "Local Occupied Orbital Energies", (char *) temp_occ_array,
            nocc * sizeof(double));

    // Build the denominator of Hbar elements
    // d_ia = H_ii - H_aa
    //d_ijab = H_ii + H_jj - H_aa - H_bb

    outfile->Printf("Building denoms...\n");
    dpdfile2 FMI;
    dpdbuf4 D;
    dpdbuf4 T2;
    double *h_oo_array = new double[nocc];
    global_dpd_->file2_copy(&Fij, PSIF_CC_OEI, "FMI");
    global_dpd_->file2_close(&Fij);

    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI");
    global_dpd_->contract442(&D, &T2, &FMI, 0, 0, 1.0, 1.0);
    global_dpd_->buf4_close(&T2);
    global_dpd_->buf4_close(&D);
    global_dpd_->file2_close(&FMI);

    global_dpd_->file2_init(&FMI, PSIF_CC_OEI, 0, 0, 0, "FMI");
    global_dpd_->file2_mat_init(&FMI);
    global_dpd_->file2_mat_rd(&FMI);
    for(int i=0; i < nocc; i++) {
        h_oo_array[i] = FMI.matrix[0][i][i];
    }
    global_dpd_->file2_close(&FMI);
    
    double *h_vv_array = new double[nocc];
    dpdfile2 FAE;
    dpdfile2 Fab;
    global_dpd_->file2_init(&Fab, PSIF_CC_OEI, 0, 1, 1, "fAB");
    global_dpd_->file2_copy(&Fab, PSIF_CC_OEI, "FAE");
    global_dpd_->file2_close(&Fab);

    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
    global_dpd_->contract442(&T2, &D, &FAE, 3, 3, -1, 1);
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_close(&T2);
    global_dpd_->file2_close(&FAE);

    global_dpd_->file2_init(&FAE, PSIF_CC_OEI, 0, 1, 1, "FAE");
    global_dpd_->file2_mat_init(&FAE);
    global_dpd_->file2_mat_rd(&FAE);
    /*global_dpd_->file2_mat_print(&FAE, "outfile");*/
    for(int a=0; a < nvir; a++) {
        h_vv_array[a] = FAE.matrix[0][a][a];
    }
    global_dpd_->file2_close(&FAE);

    outfile->Printf("...done\n");

    // Build the similarity-transformed perturbation
    outfile->Printf("Building Abar...\n");

    std::vector<SharedMatrix> A;
    A.resize(3);

    // Do the AO-MO transformation of the perturbation
    char lbl[32];
    std::vector<char*> cart_list = {strdup("x"),strdup("y"),strdup("z")};

    // Do the similarity transformation with MP2 T2s
    dpdbuf4 fbar;
    for (int n=0; n < 3; n++) {
        dpdfile2 AAE;
        dpdfile2 AMI;
        sprintf(lbl, "Pertbar_%s", cart_list[n]);
        global_dpd_->buf4_init(&fbar, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, lbl);
        sprintf(lbl, "PertAE_%s", cart_list[n]);
        global_dpd_->file2_init(&AAE, PSIF_CC_OEI, 0, 1, 1, lbl);
        global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
        global_dpd_->contract424(&T2, &AAE, &fbar, 3, 1, 0, 1, 0); 
        global_dpd_->contract244(&AAE, &T2, &fbar, 1, 2, 1, 1, 1); 
        global_dpd_->buf4_close(&T2);
        global_dpd_->file2_close(&AAE);

        sprintf(lbl, "PertMI_%s", cart_list[n]);
        global_dpd_->file2_init(&AMI, PSIF_CC_OEI, 0, 0, 0, lbl);
        global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
        global_dpd_->contract424(&T2, &AMI, &fbar, 1, 0, 1, -1, 1); 
        global_dpd_->contract244(&AMI, &T2, &fbar, 0, 0, 0, -1, 1); 
        global_dpd_->buf4_close(&T2);
        global_dpd_->file2_close(&AMI);

        global_dpd_->buf4_close(&fbar);

    }
    outfile->Printf("...done\n");

    // Build X_ij as Abar/denom
    outfile->Printf("Building X_ij...\n");
    dpdbuf4 Xijab;
    std::vector<SharedMatrix> Dij;
    SharedMatrix temp(new Matrix(nvir, nvir));
    for (int ij=0; ij < npairs; ij++) {
        auto mat = std::make_shared<Matrix>(nvir, nvir);
        Dij.push_back(mat->clone());
    }
    outfile->Printf("Initialized density\n");
    for (int n=0; n < 3; n++) {
        sprintf(lbl, "Pertbar_%s", cart_list[n]);
        global_dpd_->buf4_init(&fbar, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, lbl);
        sprintf(lbl, "Xijab_%s", cart_list[n]);
        global_dpd_->buf4_init(&Xijab, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, lbl);
        global_dpd_->buf4_mat_irrep_init(&Xijab, 0);
        global_dpd_->buf4_mat_irrep_init(&fbar, 0);
        global_dpd_->buf4_mat_irrep_rd(&fbar, 0);
        for(int ij=0; ij < nocc*nocc; ++ij) {
            int i =  ij/(nocc);
            int j = ij%(nocc);
            for(int ab=0; ab < nvir*nvir; ++ab) {
                int a =  ab/(nvir);
                int b = ab%(nvir);
                outfile->Printf("Fbar before dividing: %10.10f \n Omega: %10.10f \n", fbar.matrix[0][ij][ab], omega);
                Xijab.matrix[0][ij][ab] = fbar.matrix[0][ij][ab] / (h_oo_array[i] + h_oo_array[j] - h_vv_array[a] - h_vv_array[b] + omega);
                outfile->Printf("X after dividing: %10.10f \n", Xijab.matrix[0][ij][ab]);
            }
        }
        global_dpd_->buf4_close(&fbar);
        global_dpd_->buf4_mat_irrep_wrt(&Xijab, 0);
        global_dpd_->buf4_close(&Xijab);
        outfile->Printf("Divided by denom\n");

        sprintf(lbl, "Xijab_%s", cart_list[n]);
        global_dpd_->buf4_init(&Xijab, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, lbl);
        get_matvec(&Xijab, &Xij);

        /*outfile->Printf("*** IJ Matrices ***\n");
        for(int ij=0; ij<Xij.size(); ++ij) {
            Xij.at(ij)->print();   }*/

        sprintf(lbl, "Xtijab_%s", cart_list[n]);
        global_dpd_->buf4_scmcopy(&Xijab, PSIF_CC_TMP0, lbl, 2);
        global_dpd_->buf4_sort_axpy(&Xijab, PSIF_CC_TMP0, pqsr, 0, 5, lbl, -1);
        global_dpd_->buf4_close(&Xijab);

        global_dpd_->buf4_init(&Xijab, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, lbl);
        get_matvec(&Xijab, &Xtij);
        global_dpd_->buf4_close(&Xijab);
        outfile->Printf("Got X_ij\n");

        // Create density
        std::vector<SharedMatrix> Dtemp;
        for(int ij=0; ij < npairs; ++ij) {
            int i = ij/nocc;
            int j = ij%nocc;
            temp->gemm(0, 1, 1, Xij[ij], Xtij[ij], 0);
            temp->gemm(1, 0, 1, Xij[ij], Xtij[ij], 1);
            temp->scale(2.0 * 1/(1+(i==j)));
            temp->axpy(1.0, temp->transpose());
            temp->scale(0.5);
            Dtemp.push_back(temp->clone());
        }

        for(int ij=0; ij < npairs; ++ij) {
            Dij[ij]->add(Dtemp[ij]);
        }
    }
    outfile->Printf("...done\n");
    for(int ij=0; ij < npairs; ++ij) {
        Dij[ij]->scale(1.0/3.0);
    }
    // Diagonalize density
    Q = build_PNO_lists(cutoff, Dij);
    outfile->Printf("Density diagonalized successfully.\n");

    // Get semicanonical transforms
    get_semicanonical_transforms(Q);

    // Assuming this memory needs to be freed
    delete[] temp_occ_array;
}
void Local_cc::get_matvec(dpdbuf4 *buf_obj, std::vector<SharedMatrix> *matvec) {
    
    auto mat = std::make_shared<Matrix>(nvir, nvir);

    global_dpd_->buf4_mat_irrep_init(buf_obj,0);
    global_dpd_->buf4_mat_irrep_rd(buf_obj, 0);

    for(int ij=0; ij < nocc*nocc; ++ij) {
        for(int ab=0; ab < nvir*nvir; ++ab) {
            int a =  ab/(nvir);
            int b = ab%(nvir);
            // outfile->Printf("ij: %d, ab: %d a: %d b: %d \nT2 val(ij,ab): %f \n", ij, ab, a, b, buf_obj->matrix[0][ij][ab]);
            mat->set(a, b, buf_obj->matrix[0][ij][ab]);
        }
        matvec->push_back(mat->clone());
    }

    global_dpd_->buf4_mat_irrep_close(buf_obj, 0);
}

void Local_cc::get_semicanonical_transforms(std::vector<SharedMatrix> Q) {
    
    std::vector<SharedMatrix> L;
    std::vector<SharedVector> eps_pno;

    // Read in virtual block of Fock matrix
    dpdfile2 Fab;
    global_dpd_->file2_init(&Fab, PSIF_CC_OEI, 0, 1, 1, "fAB");
    global_dpd_->file2_mat_init(&Fab);
    global_dpd_->file2_mat_rd(&Fab);
    auto Fvir = std::make_shared<Matrix>(nvir, nvir);
    for(int a=0; a < nvir; ++a) {
        for(int b=0; b < nvir; ++b) {
            Fvir->set(a, b, Fab.matrix[0][a][b]);
        }
    }
    global_dpd_->file2_close(&Fab);

    // Transform F_vir to PNO basis
    for(int ij=0; ij < npairs; ++ij) {
        int npno = Q[ij]->colspi(0);
        auto atemp = std::make_shared<Matrix>(nvir, npno);
        auto Fpno = std::make_shared<Matrix>(npno, npno);
        atemp->gemm(0, 0, 1, Fvir, Q[ij], 0);
        Fpno->gemm(1, 0, 1, Q[ij], atemp, 1);
    // Diagonalize to obtain L, eps_vir
        auto eps = std::make_shared<Vector>(npno);
        auto evecs = std::make_shared<Matrix>(npno, npno);
        Fpno->diagonalize(evecs, eps, ascending);
        L.push_back(evecs->clone());
        eps_pno.push_back(eps);
    }

    // Write L and eps_pno to file
    // for use during CC iterations
    psio_address next;
    next = PSIO_ZERO;
    for(int ij=0; ij < npairs; ++ij) {
        int npno = Q[ij]->colspi(0);
        psio_write(PSIF_CC_INFO, "Semicanonical Transformation Matrix L", (char *) L[ij]->pointer()[0],
                npno * npno * sizeof(double), next, &next);
    }
    next = PSIO_ZERO;
    for(int ij=0; ij < npairs; ++ij) {
        int npno = Q[ij]->colspi(0);
        psio_write(PSIF_CC_INFO, "Local Virtual Orbital Energies", (char *) eps_pno[ij]->pointer(),
                npno * sizeof(double), next, &next);
    }

}

std::vector<SharedMatrix> Local_cc::build_PNO_lists(double cutoff, std::vector<SharedMatrix> Dij) {
    std::vector<SharedMatrix> Q_full(npairs);
    std::vector<SharedVector> occ_num(npairs);
    std::vector<SharedMatrix> Q;
    std::vector<int> survivor_list;

    for(int ij=0; ij < npairs; ++ij) {
        occ_num[ij] = std::make_shared<Vector>(nvir);
        Q_full[ij]   = std::make_shared<Matrix>(nvir, nvir);
    }

    // Diagonalize D
    for(int ij=0; ij < npairs; ++ij) {
        Dij[ij]->diagonalize(Q_full[ij], occ_num[ij], descending);
    }

    // Identify survivors
    double abs_occ;
    int survivors;
    for(int ij=0; ij < npairs; ++ij) {
        survivors = 0;
        for(int a=0; a < nvir; ++a) {
            abs_occ = fabs(occ_num[ij]->get(a));
            if( abs_occ >= cutoff) {
                survivors += 1;
            }
        }
        //outfile->Printf("Survivors: %i \t", survivors);
        survivor_list.push_back(survivors);
    }

    // Compute stats
    int total_pno = 0;
    double t2_ratio = 0.0;
    outfile->Printf("\nSurvivor list: "); 
    for(int ij=0; ij < npairs; ++ij) {
        outfile->Printf("%d\t", survivor_list[ij]);
        total_pno += survivor_list[ij];
        t2_ratio += pow(survivor_list[ij],2);
    }
    outfile->Printf("\nT2 ratio: %10.10lf \n", t2_ratio);
    double avg_pno = static_cast<float>(total_pno) / static_cast<float>(npairs);
    t2_ratio /= (nocc*nocc*nvir*nvir);

    // Print stats
    outfile->Printf("Total number of PNOs: %i \n", total_pno);
    outfile->Printf("Average number of PNOs: %10.10lf \n", avg_pno);
    outfile->Printf("T2 ratio: %10.10lf \n", t2_ratio);

    // Truncate Q
    // If I understood Slices I wouldn't need to use
    // for loops to set individual matrix elements
    for(int ij=0; ij < npairs; ++ij) {
        int npno = survivor_list[ij];
        auto qtemp = std::make_shared<Matrix>(nvir, npno);
        for(int a=0; a < nvir; ++a) {
            for(int aij=0; aij < npno; ++aij) {
                qtemp->set(a, aij, Q_full[ij]->get(a, aij));
            }
        }
        Q.push_back(qtemp->clone());
    }

    // Print check Q
    /*outfile->Printf("**** Truncated Q ****\n");
    for (auto &qel : Q) {
        qel->print();
    }*/

    // Write PNO dimensions, Q to file
    int *survivors_list = new int[npairs];
    for (int ij=0; ij < npairs; ij++) {
        survivors_list[ij] = survivor_list[ij];
    }
    psio_write_entry(PSIF_CC_INFO, "PNO dimensions", (char *) survivors_list, npairs * sizeof(int));
    // Segfaulting for PNO++
    psio_address next;
    next = PSIO_ZERO;
    for(int ij=0; ij < npairs; ++ij) {
        int npno = survivor_list[ij];
        psio_write(PSIF_CC_INFO, "Local Transformation Matrix Q", (char *) Q[ij]->pointer()[0],
                nvir * npno * sizeof(double), next, &next);
    }
    outfile->Printf("Wrote PNOs to file\n");
    
    delete[] survivors_list;
    return Q;
}

void Local_cc::local_done() { outfile->Printf("    Local parameters free.\n"); }

void Local_cc::local_filter_T1(dpdfile2 *T1) {
    int ii;
    psio_address next;
    npairs = nocc*nocc;

    // Switching to cleaner std::vector for memory allocation
    std::vector<SharedMatrix> Q;
    std::vector<SharedMatrix> L;
    std::vector<SharedVector> eps_pno;

    //survivor_list.resize(npairs);
    int *survivors_list = new int[npairs];
    double *occ_eps = new double[nocc];
    /*   local.weak_pairs = init_int_array(nocc*nocc); */
    psio_read_entry(PSIF_CC_INFO, "Local Occupied Orbital Energies", (char *) occ_eps, nocc * sizeof(double));
    psio_read_entry(PSIF_CC_INFO, "PNO dimensions", (char *) survivors_list, npairs * sizeof(int));
    //survivor_list.insert(survivor_list.begin(), std::begin(survivors), std::end(survivors));

    /*outfile->Printf("\nChecking read_in of survivor_list, nocc=%d\n", nocc);
    for (auto sel : survivors) {
        outfile->Printf("%d\t", sel);
    }
    outfile->Printf("End\n");*/

    next = PSIO_ZERO;
    for (int ij = 0; ij < npairs; ++ij) {
        int npno = survivors_list[ij];
        auto eps_pno_temp = std::make_shared<Vector>(npno);
        psio_read(PSIF_CC_INFO, "Local Virtual Orbital Energies", (char *)eps_pno_temp->pointer(),
                  npno * sizeof(double), next, &next);
        eps_pno.push_back(eps_pno_temp);
    }
    next = PSIO_ZERO;
    for (int ij = 0; ij < npairs; ++ij) {
        int npno = survivors_list[ij];
        auto qtemp = std::make_shared<Matrix>(nvir, npno);
        psio_read(PSIF_CC_INFO, "Local Transformation Matrix Q", (char *)qtemp->pointer()[0],
                  sizeof(double) * nvir * npno, next, &next);
        Q.push_back(qtemp->clone());
    }
    /*outfile->Printf("Testing read-in of Q\n");
    for (auto &qel : Q) {
        qel->print();
    }*/
    next = PSIO_ZERO;
    for (int ij = 0; ij < nocc * nocc; ij++) {
        int npno = survivors_list[ij];
        auto ltemp = std::make_shared<Matrix>(npno, npno);
        psio_read(PSIF_CC_INFO, "Semicanonical Transformation Matrix L", (char *)ltemp->pointer()[0],
                  sizeof(double) * npno * npno, next, &next);
        L.push_back(ltemp->clone());
    }

    global_dpd_->file2_mat_init(T1);
    global_dpd_->file2_mat_rd(T1);

    double *T1tilde, *T1bar;
    double **qtemp, **ltemp;

    for (int i = 0; i < nocc; i++) {
        ii = i * nocc + i; /* diagonal element of pair matrices */
        int npno = survivors_list[ii];

        // Print checking
        /*Vector T1vec(nvir);
        for(int a=0; a < nvir; ++a) {
            T1vec.set(a, T1->matrix[0][i][a]);
        }
        outfile->Printf("**** T1[%d] vector before denom applied ****", i);
        T1vec.print();*/
        T1tilde = init_array(npno);
        T1bar = init_array(npno);
        qtemp = Q[ii]->to_block_matrix();
        ltemp = L[ii]->to_block_matrix();

        /* Transform the virtuals to the redundant projected virtual basis */
        // T1tilde->gemv(1, 1, Q[ii].get(), &T1vec, 0);
        C_DGEMV('t', nvir, npno, 1.0, qtemp[0], npno, &(T1->matrix[0][i][0]), 1.0, 0.0, &(T1tilde[0]), 1);

        /* Transform the virtuals to the non-redundant virtual basis */
        // T1bar->gemv(1, 1, L[ii].get(), T1tilde.get(), 0);
        C_DGEMV('t', npno, npno, 1.0, ltemp[0], npno, &(T1tilde[0]), 1.0, 0.0, &(T1bar[0]), 1);

        /* Apply the denominators */
        for (int a = 0; a < npno; a++) {
            //outfile->Printf("T1 Residual before denom: %10.15f \n", T1bar[a]);
            T1bar[a] /= occ_eps[i] - eps_pno[ii]->get(a);
            //outfile->Printf("T1 Residual after denom: %10.15f \n", T1bar[a]);
        }

        /* Transform the new T1's to the redundant projected virtual basis */
        // T1tilde->gemv(0, 1, L[ii].get(), T1bar.get(), 0);
        C_DGEMV('n', npno, npno, 1.0, ltemp[0], npno, &(T1bar[0]), 1.0, 0.0, &(T1tilde[0]), 1);

        /* Transform the new T1's to the MO basis */
        // T1vec.gemv(0, 1, Q[ii].get(), T1tilde.get(), 0);
        C_DGEMV('n', nvir, npno, 1.0, qtemp[0], npno, &(T1tilde[0]), 1.0, 0.0, &(T1->matrix[0][i][0]), 1);

    }

    global_dpd_->file2_mat_wrt(T1);
    global_dpd_->file2_mat_close(T1);

    delete[] survivors_list;
    delete[] occ_eps;
    free(T1tilde);
    free(T1bar);
}

void Local_cc::local_filter_T2(dpdbuf4 *T2) {
    psio_address next;
    npairs = nocc*nocc;

    // Switching to cleaner std::vector for memory allocation
    std::vector<SharedMatrix> Q;
    std::vector<SharedMatrix> L;
    std::vector<SharedVector> eps_pno;

    //survivor_list.resize(npairs, 0);
    int *survivors_list = new int[npairs];
    double *occ_eps = new double[nocc];

    psio_read_entry(PSIF_CC_INFO, "Local Occupied Orbital Energies", (char *) occ_eps, nocc * sizeof(double));
    psio_read_entry(PSIF_CC_INFO, "PNO dimensions", (char *) survivors_list, npairs * sizeof(int));
    /*outfile->Printf("Testing read-in of eps_occ\n");
    for (int i=0; i < nocc; i++) {
        outfile->Printf("%5.15f\t",occ_eps[i]);
    }
    for (auto qel : eps_occ) {
        outfile->Printf("%20.10f\t", qel);;
    }*/

    next = PSIO_ZERO;
    for (int ij = 0; ij < npairs; ++ij) {
        int npno = survivors_list[ij];
        auto eps_pno_temp = std::make_shared<Vector>(npno);
        psio_read(PSIF_CC_INFO, "Local Virtual Orbital Energies", (char *)eps_pno_temp->pointer(),
                  npno * sizeof(double), next, &next);
        eps_pno.push_back(eps_pno_temp);
    }
    /*outfile->Printf("Testing read-in of virtual orbital energies\n");
    for(int ij=0; ij < npairs; ++ij) {
        outfile->Printf("\nPair %d ", ij);
        int npno = survivors[ij];
        for(int a=0; a < npno; ++a) {
            outfile->Printf("%10.10f\t", eps_pno[ij]->get(a));
        }
    }*/
    next = PSIO_ZERO;
    for (int ij = 0; ij < npairs; ij++) {
        int npno = survivors_list[ij];
        auto qtemp = std::make_shared<Matrix>(nvir, npno);
        psio_read(PSIF_CC_INFO, "Local Transformation Matrix Q", (char *)qtemp->pointer()[0],
                  sizeof(double) * nvir * npno, next, &next);
        Q.push_back(qtemp->clone());
    }
    next = PSIO_ZERO;
    for (int ij = 0; ij < npairs; ij++) {
        int npno = survivors_list[ij];
        auto ltemp = std::make_shared<Matrix>(npno, npno);
        psio_read(PSIF_CC_INFO, "Semicanonical Transformation Matrix L", (char *)ltemp->pointer()[0],
                  sizeof(double) * npno * npno, next, &next);
        L.push_back(ltemp->clone());
    }
    /*outfile->Printf("Testing read-in of L\n");
    for (auto &qel : L) {
        qel->print();
    }*/

    /* Grab the MO-basis T2's */
    global_dpd_->buf4_mat_irrep_init(T2, 0);
    global_dpd_->buf4_mat_irrep_rd(T2, 0);

    for (int ij = 0; ij < npairs; ++ij) {
        int npno = survivors_list[ij];
        auto T2temp = std::make_shared<Matrix>(nvir, nvir);
        auto atemp = std::make_shared<Matrix>(nvir, npno);
        auto btemp = std::make_shared<Matrix>(npno, npno);
        auto T2bar = std::make_shared<Matrix>(npno, npno);

        for(int ab=0; ab < nvir*nvir; ++ab) {
            int a = ab / nvir;
            int b = ab % nvir;
            T2temp->set(a, b, T2->matrix[0][ij][ab]);
        }

        //outfile->Printf("T2 residual before change of basis\n");
        //T2temp->print();
        /* Transform the virtuals to the redundant projected virtual basis */
        atemp->gemm(0, 0, 1, T2temp, Q[ij], 0);
        T2bar->gemm(1, 0, 1, Q[ij], atemp, 0);

        /* Transform the virtuals to the non-redundant virtual basis */
        btemp->gemm(0, 0, 1, T2bar, L[ij], 0);
        T2bar->gemm(1, 0, 1, L[ij], btemp, 0);

        /* Divide the new amplitudes by the denominators */
        int i = ij / nocc;
        int j = ij % nocc;
        for (int a = 0; a < npno; a++) {
            for (int b = 0; b < npno; b++) {
                double val = T2bar->get(a,b);
                T2bar->set(a, b, val * 1.0 /
                    (occ_eps[i] + occ_eps[j] - eps_pno[ij]->get(a) - eps_pno[ij]->get(b)));
            }
        }
        /* Transform the new T2's to the redundant virtual basis */
        btemp->gemm(0, 0, 1, L[ij], T2bar, 0);
        T2bar->gemm(0, 1, 1, btemp, L[ij], 0);

        /* Transform the new T2's to the MO basis */
        atemp->gemm(0, 0, 1, Q[ij], T2bar, 0);
        T2temp->gemm(0, 1, 1, atemp, Q[ij], 0);

        for(int ab=0; ab < nvir*nvir; ++ab) {
            int a = ab / nvir;
            int b = ab % nvir;
            T2->matrix[0][ij][ab] = T2temp->get(a,b);

        }
    }

    /* Write the updated MO-basis T2's to disk */
    global_dpd_->buf4_mat_irrep_wrt(T2, 0);
    global_dpd_->buf4_mat_irrep_close(T2, 0);

    delete[] survivors_list;
    delete[] occ_eps;

}

void Local_cc::init_filter_T2() {
    dpdbuf4 T2;
    dpdbuf4 D;

    //Re-initialize T2s, but with local filter applied
    global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D <ij|ab>");
    global_dpd_->buf4_copy(&D, PSIF_CC_TAMPS, "tIjAb");
    global_dpd_->buf4_close(&D);
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tIjAb");
    local_filter_T2(&T2);
    global_dpd_->buf4_close(&T2);
}
void Local_cc::mp2_pair_energy() {
    
}

// Destructor attempt
Local_cc::~Local_cc() {
};

} // namespace psi
