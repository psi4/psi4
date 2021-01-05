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
    \ingroup CCLAMBDA
    \brief Enter brief description of file here
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libiwl/iwl.h"
#include "psi4/libqt/qt.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libpsi4util/exception.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libmints/vector.h"
#include "psi4/libmints/local.h"
#include "psi4/libmints/basisset.h"
#include "psi4/psifiles.h"
#include "MOInfo.h"
#include "Params.h"
#include "Local.h"
#include "cclambda.h"
#define EXTERN
#include "globals.h"

namespace psi {
namespace cclambda {

/*!
** local_init(): Set up parameters of local excitation domains.
**
** The orbital domains constructed here are based on those described
** in Broughton and Pulay, J. Comp. Chem. 14, 736-740 (1993).  The
** localization of the occupied orbitals is done elsewhere (see the
** program "local").  Pair domains are defined as the union of pairs
** of single occupied orbital domains.  "Weak pairs", which are
** defined as pair domains whose individual occupied orbital domains
** have no atoms in common, are identified (cf. int *weak_pairs).
**
** TDC, Jan-June 2002
*/

void CCLambdaWavefunction::local_init() {
    local.nso = moinfo.nso;
    local.nocc = moinfo.occpi[0];  /* active doubly occupied orbitals */
    local.nvir = moinfo.virtpi[0]; /* active virtual orbitals */

    outfile->Printf("\tLocalization parameters ready.\n\n");
}

void CCLambdaWavefunction::local_done() { outfile->Printf("\tLocal parameters free.\n"); }

void CCLambdaWavefunction::local_filter_T1(dpdfile2 *T1) {
    int i, a, ij, ii;
    int nocc, nvir;
    double *T1tilde, *T1bar;
    psio_address next;

    nocc = local.nocc;
    nvir = local.nvir;

    if (local.method == "PNO") {
        pno_filter_T1(T1);
    }
    else {

        local.pairdom_len = init_int_array(nocc * nocc);
        local.pairdom_nrlen = init_int_array(nocc * nocc);
        local.eps_occ = init_array(nocc);
        psio_read_entry(PSIF_CC_INFO, "Local Pair Domain Length", (char *)local.pairdom_len, sizeof(int) * nocc * nocc);
        psio_read_entry(PSIF_CC_INFO, "Local Pair Domain NR Length", (char *)local.pairdom_nrlen,
                        sizeof(int) * nocc * nocc);
        psio_read_entry(PSIF_CC_INFO, "Local Occupied Orbital Energies", (char *)local.eps_occ, nocc * sizeof(double));

        local.W = (double ***)malloc(sizeof(double **) * nocc * nocc);
        local.V = (double ***)malloc(sizeof(double **) * nocc * nocc);
        local.eps_vir = (double **)malloc(sizeof(double *) * nocc * nocc);
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
                      sizeof(double) * nvir * local.pairdom_len[ij], next, &next);
        }
        next = PSIO_ZERO;
        for (ij = 0; ij < nocc * nocc; ij++) {
            local.W[ij] = block_matrix(local.pairdom_len[ij], local.pairdom_nrlen[ij]);
            psio_read(PSIF_CC_INFO, "Local Transformation Matrix (W)", (char *)local.W[ij][0],
                      sizeof(double) * local.pairdom_len[ij] * local.pairdom_nrlen[ij], next, &next);
        }

        global_dpd_->file2_mat_init(T1);
        global_dpd_->file2_mat_rd(T1);

        for (i = 0; i < nocc; i++) {
            ii = i * nocc + i; /* diagonal element of pair matrices */

            if (!local.pairdom_len[ii]) {
                outfile->Printf("\n\tlocalfilter_T1: Pair ii = [%d] is zero-length, which makes no sense.\n", ii);
                throw PsiException("cclambda: error", __FILE__, __LINE__);
            }

            T1tilde = init_array(local.pairdom_len[ii]);
            T1bar = init_array(local.pairdom_nrlen[ii]);

            /* Transform the virtuals to the redundant projected virtual basis */
            C_DGEMV('t', nvir, local.pairdom_len[ii], 1.0, &(local.V[ii][0][0]), local.pairdom_len[ii],
                    &(T1->matrix[0][i][0]), 1, 0.0, &(T1tilde[0]), 1);

            /* Transform the virtuals to the non-redundant virtual basis */
            C_DGEMV('t', local.pairdom_len[ii], local.pairdom_nrlen[ii], 1.0, &(local.W[ii][0][0]), local.pairdom_nrlen[ii],
                    &(T1tilde[0]), 1, 0.0, &(T1bar[0]), 1);

            /* Apply the denominators */
            for (a = 0; a < local.pairdom_nrlen[ii]; a++) T1bar[a] /= (local.eps_occ[i] - local.eps_vir[ii][a]);

            /* Transform the new T1's to the redundant projected virtual basis */
            C_DGEMV('n', local.pairdom_len[ii], local.pairdom_nrlen[ii], 1.0, &(local.W[ii][0][0]), local.pairdom_nrlen[ii],
                    &(T1bar[0]), 1, 0.0, &(T1tilde[0]), 1);

            /* Transform the new T1's to the MO basis */
            C_DGEMV('n', nvir, local.pairdom_len[ii], 1.0, &(local.V[ii][0][0]), local.pairdom_len[ii], &(T1tilde[0]), 1,
                    0.0, &(T1->matrix[0][i][0]), 1);

            free(T1bar);
            free(T1tilde);
        }

        global_dpd_->file2_mat_wrt(T1);
        global_dpd_->file2_mat_close(T1);

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
    }
}

void CCLambdaWavefunction::pno_filter_T1(dpdfile2 *T1) {
    int ii;
    psio_address next;

    auto nocc = local.nocc;
    auto nvir = local.nvir;
    auto npairs = local.npairs;
    auto survivor_list = local.survivor_list;

    /*   local.weak_pairs = init_int_array(nocc*nocc); */
    local.eps_occ = init_array(nocc);
    psio_read_entry(PSIF_CC_INFO, "Local Occupied Orbital Energies", (char *) local.eps_occ, nocc * sizeof(double));

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
            T1bar[a] /= local.eps_occ[i] - eps_pno[ii]->get(a);
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
        free_block(local.Q[ij]);
        free_block(local.L[ij]);
        free(local.eps_pno[ij]);
    }
    free(local.Q);
    free(local.L);
    free(local.eps_pno);
    free(local.eps_occ);
    free(local.weak_pairs); */
}

void CCLambdaWavefunction::local_filter_T2(dpdbuf4 *T2) {
    int ij, i, j, a, b;
    int nso, nocc, nvir;
    double **X1, **X2, **T2tilde, **T2bar;
    psio_address next;

    nso = local.nso;
    nocc = local.nocc;
    nvir = local.nvir;

    if (local.method == "PNO") {
        pno_filter_T2(T2);
    }
    else {

        local.pairdom_len = init_int_array(nocc * nocc);
        local.pairdom_nrlen = init_int_array(nocc * nocc);
        local.weak_pairs = init_int_array(nocc * nocc);
        local.eps_occ = init_array(nocc);
        psio_read_entry(PSIF_CC_INFO, "Local Pair Domain Length", (char *)local.pairdom_len, sizeof(int) * nocc * nocc);
        psio_read_entry(PSIF_CC_INFO, "Local Pair Domain NR Length", (char *)local.pairdom_nrlen,
                        sizeof(int) * nocc * nocc);
        psio_read_entry(PSIF_CC_INFO, "Local Occupied Orbital Energies", (char *)local.eps_occ, nocc * sizeof(double));
        psio_read_entry(PSIF_CC_INFO, "Local Weak Pairs", (char *)local.weak_pairs, sizeof(int) * nocc * nocc);
        local.W = (double ***)malloc(sizeof(double **) * nocc * nocc);
        local.V = (double ***)malloc(sizeof(double **) * nocc * nocc);
        local.eps_vir = (double **)malloc(sizeof(double *) * nocc * nocc);
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
                      sizeof(double) * nvir * local.pairdom_len[ij], next, &next);
        }
        next = PSIO_ZERO;
        for (ij = 0; ij < nocc * nocc; ij++) {
            local.W[ij] = block_matrix(local.pairdom_len[ij], local.pairdom_nrlen[ij]);
            psio_read(PSIF_CC_INFO, "Local Transformation Matrix (W)", (char *)local.W[ij][0],
                      sizeof(double) * local.pairdom_len[ij] * local.pairdom_nrlen[ij], next, &next);
        }

        /* Grab the MO-basis T2's */
        global_dpd_->buf4_mat_irrep_init(T2, 0);
        global_dpd_->buf4_mat_irrep_rd(T2, 0);

        X1 = block_matrix(nso, nvir);
        X2 = block_matrix(nvir, nso);
        T2tilde = block_matrix(nso, nso);
        T2bar = block_matrix(nvir, nvir);

        for (i = 0, ij = 0; i < nocc; i++) {
            for (j = 0; j < nocc; j++, ij++) {
                if (!local.weak_pairs[ij]) {
                    /* Transform the virtuals to the redundant projected virtual basis */
                    C_DGEMM('t', 'n', local.pairdom_len[ij], nvir, nvir, 1.0, &(local.V[ij][0][0]), local.pairdom_len[ij],
                            &(T2->matrix[0][ij][0]), nvir, 0.0, &(X1[0][0]), nvir);
                    C_DGEMM('n', 'n', local.pairdom_len[ij], local.pairdom_len[ij], nvir, 1.0, &(X1[0][0]), nvir,
                            &(local.V[ij][0][0]), local.pairdom_len[ij], 0.0, &(T2tilde[0][0]), nso);

                    /* Transform the virtuals to the non-redundant virtual basis */
                    C_DGEMM('t', 'n', local.pairdom_nrlen[ij], local.pairdom_len[ij], local.pairdom_len[ij], 1.0,
                            &(local.W[ij][0][0]), local.pairdom_nrlen[ij], &(T2tilde[0][0]), nso, 0.0, &(X2[0][0]), nso);
                    C_DGEMM('n', 'n', local.pairdom_nrlen[ij], local.pairdom_nrlen[ij], local.pairdom_len[ij], 1.0,
                            &(X2[0][0]), nso, &(local.W[ij][0][0]), local.pairdom_nrlen[ij], 0.0, &(T2bar[0][0]), nvir);

                    /* Divide the new amplitudes by the denominators */
                    for (a = 0; a < local.pairdom_nrlen[ij]; a++) {
                        for (b = 0; b < local.pairdom_nrlen[ij]; b++) {
                            T2bar[a][b] /=
                                (local.eps_occ[i] + local.eps_occ[j] - local.eps_vir[ij][a] - local.eps_vir[ij][b]);
                        }
                    }

                    /* Transform the new T2's to the redundant virtual basis */
                    C_DGEMM('n', 'n', local.pairdom_len[ij], local.pairdom_nrlen[ij], local.pairdom_nrlen[ij], 1.0,
                            &(local.W[ij][0][0]), local.pairdom_nrlen[ij], &(T2bar[0][0]), nvir, 0.0, &(X1[0][0]), nvir);
                    C_DGEMM('n', 't', local.pairdom_len[ij], local.pairdom_len[ij], local.pairdom_nrlen[ij], 1.0,
                            &(X1[0][0]), nvir, &(local.W[ij][0][0]), local.pairdom_nrlen[ij], 0.0, &(T2tilde[0][0]), nso);

                    /* Transform the new T2's to the MO basis */
                    C_DGEMM('n', 'n', nvir, local.pairdom_len[ij], local.pairdom_len[ij], 1.0, &(local.V[ij][0][0]),
                            local.pairdom_len[ij], &(T2tilde[0][0]), nso, 0.0, &(X2[0][0]), nso);
                    C_DGEMM('n', 't', nvir, nvir, local.pairdom_len[ij], 1.0, &(X2[0][0]), nso, &(local.V[ij][0][0]),
                            local.pairdom_len[ij], 0.0, &(T2->matrix[0][ij][0]), nvir);
                } else /* This must be a neglected weak pair; force it to zero */
                    memset((void *)T2->matrix[0][ij], 0, sizeof(double) * nvir * nvir);
            }
        }
        free_block(X1);
        free_block(X2);
        free_block(T2tilde);
        free_block(T2bar);

        /* Write the updated MO-basis T2's to disk */
        global_dpd_->buf4_mat_irrep_wrt(T2, 0);
        global_dpd_->buf4_mat_irrep_close(T2, 0);

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
        free(local.weak_pairs);
    }
}

void CCLambdaWavefunction::pno_filter_T2(dpdbuf4 *T2) {
    psio_address next;

    auto nso = local.nso;
    auto nocc = local.nocc;
    auto nvir = local.nvir;
    auto npairs = local.npairs;

    local.eps_occ = init_array(nocc);
    psio_read_entry(PSIF_CC_INFO, "Local Occupied Orbital Energies", (char *) local.eps_occ, nocc * sizeof(double));

    // Is there an issue with using local variables?
    std::vector<SharedVector> eps_pno;
    std::vector<SharedMatrix> Q;
    std::vector<SharedMatrix> L;

    next = PSIO_ZERO;
    for (int ij = 0; ij < npairs; ++ij) {
        int npno = local.survivor_list[ij];
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
            outfile->Printf("%f ", local.eps_pno[ij]->get(a));
        }
    }*/
    next = PSIO_ZERO;
    for (int ij = 0; ij < npairs; ij++) {
        int npno = local.survivor_list[ij];
        auto qtemp = std::make_shared<Matrix>(nvir, npno);
        psio_read(PSIF_CC_INFO, "Local Transformation Matrix Q", (char *)qtemp->pointer()[0],
                  sizeof(double) * nvir * npno, next, &next);
        Q.push_back(qtemp->clone());
        qtemp->zero();
    }
    next = PSIO_ZERO;
    for (int ij = 0; ij < nocc * nocc; ij++) {
        int npno = local.survivor_list[ij];
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
        int npno = local.survivor_list[ij];
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
                    (local.eps_occ[i] + local.eps_occ[j] - eps_pno[ij]->get(a) - eps_pno[ij]->get(b)));
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

}  // namespace cclambda
}  // namespace psi
