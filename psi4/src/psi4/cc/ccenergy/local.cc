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

#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libiwl/iwl.h"
#include "psi4/libqt/qt.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/psifiles.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/local.h"
#include "psi4/libmints/basisset.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

namespace psi {

//Constructor attempt
Local_cc::Local_cc() {
};

void Local_cc::local_init() {

    //TODO: compute mp2 energy and store weak pairs
    weak_pair_energy = 0.0;
    weak_pairs = init_int_array(nocc * nocc);

    npairs = nocc*nocc;

    if (method == "PNO")
        init_pno();
    if (method == "PNO++")
        outfile->Printf(" Local correlation using Perturbed Pair Natural Orbitals\nCutoff value: %e", cutoff);
    outfile->Printf("    Localization parameters ready.\n\n");

}

void Local_cc::init_pno() {
    outfile->Printf(" Local correlation using Pair Natural Orbitals\nCutoff value: %e \n", cutoff);
    dpdbuf4 T2;
    dpdbuf4 T2tilde;
    dpdbuf4 D;
    std::vector<SharedMatrix> Tij;
    std::vector<SharedMatrix> Ttij;
    SharedMatrix temp(new Matrix(nvir, nvir));


    // Check if occupied block of Fock matrix is non-diagonal (localized)
    // On the way, store the occupied orbital energies
    dpdfile2 Fij;
    global_dpd_->file2_init(&Fij, PSIF_CC_OEI, 0, 1, 1, "fIJ");
    global_dpd_->file2_mat_init(&Fij);
    global_dpd_->file2_mat_rd(&Fij);
    auto Focc = std::make_shared<Matrix>(nocc, nocc);
    for(int i=0; i < nocc; ++i) {
        for(int j=0; j < nocc; ++j) {
            Focc->set(i, j, Fij.matrix[0][i][j]);
        }
    }
    global_dpd_->file2_close(&Fij);
    auto local_occ_eps = std::make_shared<Vector>(nocc);
    auto evecs = std::make_shared<Matrix>(nocc, nocc);
    Focc->diagonalize(evecs, local_occ_eps, ascending);
    outfile->Printf("Here are original local occ orbital energies\n");
    for(int i=0; i < nocc; i++) {
        outfile->Printf("%f ", local_occ_eps->get(i));
    }
    psio_write_entry(PSIF_CC_INFO, "Local Occupied Orbital Energies", (char *) local_occ_eps->pointer(),
            nocc * sizeof(double));
    /* outfile->Printf("*** Vir block of F ***");
    global_dpd_->file2_mat_print(&Fij, "outfile");
    global_dpd_->file2_close(&Fij);*/

    // "Vectorize" T2
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS,  0, 0, 5, 0, 5, 0, "tIjAb");
    /*outfile->Printf("*** T2s ***");
    global_dpd_->buf4_print(&T2, "outfile", 1);*/
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
    temp->zero();
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

    // Diagonalize Density
    std::vector<SharedMatrix> Q_full(npairs);
    std::vector<SharedVector> occ_num(npairs);
    for(int ij=0; ij < npairs; ++ij) {
        occ_num[ij] = std::make_shared<Vector>(nvir);
        Q_full[ij]   = std::make_shared<Matrix>(nvir, nvir);
    }

    for(int ij=0; ij < npairs; ++ij) {
        Dij[ij]->diagonalize(Q_full[ij], occ_num[ij], descending);
    }

    // Print checking
    //for(int ij=0; ij < npairs; ++ij) {

    /*outfile->Printf( "Pair: %d\n", 1);
    for(int a=0; a < nvir; ++a) {
        outfile->Printf("%20.12lf\n", occ_num[1]->get(a));
    }*/

    // Identify survivors
    double abs_occ;
    int survivors;
    for(int ij=0; ij < npairs; ++ij) {
        survivors = 0;
        for(int a=0; a < nvir; ++a) {
            abs_occ = abs(occ_num[ij]->get(a));
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
        qtemp->zero();
    }

    // Print check Q
    /*outfile->Printf("**** Truncated Q ****\n");
    for (auto &qel : Q) {
        qel->print();
    }*/

    // Get semicanonical transforms
    get_semicanonical_transforms();
    // Print check L
    /*outfile->Printf("**** Truncated L ****\n");
    for (auto &qel : L) {
        qel->print();
    }*/

    // Write Q, L, eps_pno to file
    psio_address next;
    psio_write_entry(PSIF_CC_INFO, "PNO dimensions", (char *) &survivor_list, npairs * sizeof(int));
    next = PSIO_ZERO;
    for(int ij=0; ij < npairs; ++ij) {
        int npno = survivor_list[ij];
        psio_write(PSIF_CC_INFO, "Local Transformation Matrix Q", (char *) Q[ij]->pointer()[0],
                nvir * npno * sizeof(double), next, &next);
    }
    next = PSIO_ZERO;
    for(int ij=0; ij < npairs; ++ij) {
        int npno = survivor_list[ij];
        psio_write(PSIF_CC_INFO, "Semicanonical Transformation Matrix L", (char *) L[ij]->pointer()[0],
                npno * npno * sizeof(double), next, &next);
    }
    next = PSIO_ZERO;
    for(int ij=0; ij < npairs; ++ij) {
        int npno = survivor_list[ij];
        psio_write(PSIF_CC_INFO, "Local Virtual Orbital Energies", (char *) eps_pno[ij]->pointer(),
                npno * sizeof(double), next, &next);
    }
    next = PSIO_ZERO;
    for(int ij=0; ij < npairs; ++ij) {
        int npno = survivor_list[ij];
        psio_write(PSIF_CC_INFO, "Local Survivor List", (char *) eps_pno[ij]->pointer(),
                npno * sizeof(double), next, &next);
    }

    // Check if they can be read back in
    /*next = PSIO_ZERO;
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
    /*free(&Q);
    free(&L);
    free(&eps_pno);*/
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

void Local_cc::get_semicanonical_transforms() {
    
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
        int npno = survivor_list[ij];
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
}

void Local_cc::local_done() { outfile->Printf("    Local parameters free.\n"); }

void Local_cc::local_filter_T1(dpdfile2 *T1) {
    int ii;
    psio_address next;
    npairs = nocc*nocc;

    /*   local.weak_pairs = init_int_array(nocc*nocc); */
    eps_occ = init_array(nocc);
    double *survivors = init_array(npairs);
    psio_read_entry(PSIF_CC_INFO, "Local Occupied Orbital Energies", (char *) eps_occ, nocc * sizeof(double));
    psio_read_entry(PSIF_CC_INFO, "PNO dimensions", (char *) &survivor_list, npairs * sizeof(int));
    //survivor_list.insert(survivor_list.begin(), std::begin(survivors), std::end(survivors));

    // Is there an issue with using local_ variables?
    std::vector<SharedVector> eps_pno;
    std::vector<SharedMatrix> Q;
    std::vector<SharedMatrix> L;

    next = PSIO_ZERO;
    for (int ij = 0; ij < npairs; ++ij) {
        int npno = survivor_list[ij];
        auto eps_pno_temp = std::make_shared<Vector>(npno);
        psio_read(PSIF_CC_INFO, "Local Virtual Orbital Energies", (char *)eps_pno_temp->pointer(),
                  npno * sizeof(double), next, &next);
        eps_pno.push_back(eps_pno_temp);
    }
    next = PSIO_ZERO;
    for (int ij = 0; ij < npairs; ++ij) {
        int npno = survivor_list[ij];
        auto qtemp = std::make_shared<Matrix>(nvir, npno);
        psio_read(PSIF_CC_INFO, "Local Transformation Matrix Q", (char *)qtemp->pointer()[0],
                  sizeof(double) * nvir * npno, next, &next);
        Q.push_back(qtemp->clone());
        qtemp->zero();
    }
    /*outfile->Printf("Testing read-in of Q\n");
    for (auto &qel : Q) {
        qel->print();
    }*/
    next = PSIO_ZERO;
    for (int ij = 0; ij < nocc * nocc; ij++) {
        int npno = survivor_list[ij];
        auto ltemp = std::make_shared<Matrix>(npno, npno);
        psio_read(PSIF_CC_INFO, "Semicanonical Transformation Matrix L", (char *)ltemp->pointer()[0],
                  sizeof(double) * npno * npno, next, &next);
        L.push_back(ltemp->clone());
        ltemp->zero();
    }

    global_dpd_->file2_mat_init(T1);
    global_dpd_->file2_mat_rd(T1);

    double *T1tilde, *T1bar;
    double **qtemp, **ltemp;

    for (int i = 0; i < nocc; i++) {
        ii = i * nocc + i; /* diagonal element of pair matrices */
        int npno = survivor_list[ii];

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
            T1bar[a] /= eps_occ[i] - eps_pno[ii]->get(a);
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

    /*for (int ij = 0; ij < npairs; ij++) {
        free_block(Q[ij]);
        free_block(L[ij]);
        free(eps_pno[ij]);
    }
    free(Q);
    free(L);
    free(eps_pno);
    free(eps_occ);
    free(weak_pairs); */
}

void Local_cc::local_filter_T2(dpdbuf4 *T2) {
    psio_address next;
    npairs = nocc*nocc;

    eps_occ = init_array(nocc);
    psio_read_entry(PSIF_CC_INFO, "Local Occupied Orbital Energies", (char *) eps_occ, nocc * sizeof(double));
    psio_read_entry(PSIF_CC_INFO, "PNO dimensions", (char *) &survivor_list, npairs * sizeof(int));

    // Is there an issue with using local_ variables?
    std::vector<SharedVector> eps_pno;
    std::vector<SharedMatrix> Q;
    std::vector<SharedMatrix> L;

    next = PSIO_ZERO;
    for (int ij = 0; ij < npairs; ++ij) {
        int npno = survivor_list[ij];
        auto eps_pno_temp = std::make_shared<Vector>(npno);
        psio_read(PSIF_CC_INFO, "Local Virtual Orbital Energies", (char *)eps_pno_temp->pointer(),
                  npno * sizeof(double), next, &next);
        eps_pno.push_back(eps_pno_temp);
    }
    /*outfile->Printf("Testing read-in of virtual orbital energies\n");
    for(int ij=0; ij < npairs; ++ij) {
        outfile->Printf("\nPair %d ", ij);
        int npno = survivor_list[ij];
        for(int a=0; a < npno; ++a) {
            outfile->Printf("%f ", local_.eps_pno[ij]->get(a));
        }
    }*/
    next = PSIO_ZERO;
    for (int ij = 0; ij < npairs; ij++) {
        int npno = survivor_list[ij];
        auto qtemp = std::make_shared<Matrix>(nvir, npno);
        psio_read(PSIF_CC_INFO, "Local Transformation Matrix Q", (char *)qtemp->pointer()[0],
                  sizeof(double) * nvir * npno, next, &next);
        Q.push_back(qtemp->clone());
        qtemp->zero();
    }
    next = PSIO_ZERO;
    for (int ij = 0; ij < nocc * nocc; ij++) {
        int npno = survivor_list[ij];
        auto ltemp = std::make_shared<Matrix>(npno, npno);
        psio_read(PSIF_CC_INFO, "Semicanonical Transformation Matrix L", (char *)ltemp->pointer()[0],
                  sizeof(double) * npno * npno, next, &next);
        L.push_back(ltemp->clone());
        ltemp->zero();
    }

    /* Grab the MO-basis T2's */
    global_dpd_->buf4_mat_irrep_init(T2, 0);
    global_dpd_->buf4_mat_irrep_rd(T2, 0);

    for (int ij = 0; ij < npairs; ++ij) {
        int npno = survivor_list[ij];
        auto atemp = std::make_shared<Matrix>(nvir, npno);
        auto T2tilde = std::make_shared<Matrix>(npno, npno);
        auto T2temp = std::make_shared<Matrix>(nvir, nvir);

        for(int ab=0; ab < nvir*nvir; ++ab) {
            int a = ab / nvir;
            int b = ab % nvir;
            T2temp->set(a, b, T2->matrix[0][ij][ab]);
        }

        /* Transform the virtuals to the redundant projected virtual basis */
        atemp->gemm(0, 0, 1, T2temp, Q[ij], 0);
        T2tilde->gemm(1, 0, 1, Q[ij], atemp, 0);

        auto btemp = std::make_shared<Matrix>(npno, npno);
        auto T2bar = std::make_shared<Matrix>(npno, npno);
        /* Transform the virtuals to the non-redundant virtual basis */
        btemp->gemm(0, 0, 1, T2tilde, L[ij], 0);
        T2bar->gemm(1, 0, 1, L[ij], btemp, 0);

        /* Divide the new amplitudes by the denominators */
        int i = ij / nocc;
        int j = ij % nocc;
        for (int a = 0; a < npno; a++) {
            for (int b = 0; b < npno; b++) {
                T2bar->set(a, b, T2bar->get(a,b) /
                    (eps_occ[i] + eps_occ[j] - eps_pno[ij]->get(a) - eps_pno[ij]->get(b)));
            }
        }
        /* Transform the new T2's to the redundant virtual basis */
        btemp->gemm(0, 0, 1, L[ij], T2bar, 0);
        T2tilde->gemm(0, 1, 1, btemp, L[ij], 0);

        /* Transform the new T2's to the MO basis */
        atemp->gemm(0, 0, 1, Q[ij], T2tilde, 0);
        T2temp->gemm(0, 1, 1, atemp, Q[ij], 0);

        for(int ab=0; ab < nvir*nvir; ++ab) {
            int a = ab / nvir;
            int b = ab % nvir;
            T2->matrix[0][ij][ab] = T2temp->get(a,b);

        }
    }

    /*free_block(T2tilde);
    free_block(T2bar);*/

    /* Write the updated MO-basis T2's to disk */
    global_dpd_->buf4_mat_irrep_wrt(T2, 0);
    global_dpd_->buf4_mat_irrep_close(T2, 0);

}

// Destructor attempt
Local_cc::~Local_cc() {
};

} // namespace psi
