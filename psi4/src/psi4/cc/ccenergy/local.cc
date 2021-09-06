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
namespace ccenergy {

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

void CCEnergyWavefunction::local_init() {
    local_.nso = moinfo_.nso;
    local_.nocc = moinfo_.occpi[0];  /* active doubly occupied orbitals */
    local_.nvir = moinfo_.virtpi[0]; /* active virtual orbitals */
    auto nocc = local_.nocc;

    local_.weak_pair_energy = 0.0;

    local_.weak_pairs = init_int_array(nocc * nocc);
    // The following used to read "Local Weak Pairs", not generated
    // in the current version of Psi4
    // psio_read_entry(PSIF_CC_INFO, "Local Weak Pairs", (char *)local_.weak_pairs,
    //                sizeof(int) * local_.nocc * local_.nocc);
    
    if (local_.method == "PNO")
        init_pno();
    if (local_.method == "PNO++")
        outfile->Printf(" Local correlation using Perturbed Pair Natural Orbitals\nCutoff value: %e", local_.cutoff);
    outfile->Printf("    Localization parameters ready.\n\n");
}

void CCEnergyWavefunction::init_pno() {
    outfile->Printf(" Local correlation using Pair Natural Orbitals\nCutoff value: %e \n", local_.cutoff);
    auto nocc = local_.nocc;
    auto nvir = local_.nvir;
    auto pno_cut = local_.cutoff;
    auto npairs = nocc*nocc;
    dpdbuf4 T2;
    dpdbuf4 T2tilde;
    dpdbuf4 D;
    std::vector<SharedMatrix> Tij;
    std::vector<SharedMatrix> Ttij;
    SharedMatrix temp(new Matrix(nvir, nvir));

    // Create T2_tilde and "vectorize" it
    global_dpd_->buf4_init(&T2tilde, PSIF_CC_DINTS, 0, 0, 5, 0, 5, 0, "D 2<ij|ab> - <ij|ba>");
    global_dpd_->buf4_init(&D, PSIF_CC_DENOM, 0, 0, 5, 0, 5, 0, "dIjAb");
    global_dpd_->buf4_dirprd(&D, &T2tilde);
    get_matvec(&T2tilde, &Ttij);
    global_dpd_->buf4_close(&T2tilde);
    global_dpd_->buf4_close(&D);

    //Print check
    // This doesn't work; why?
    amp_write();

    // "Vectorize" T2
    global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS,  0, 0, 5, 0, 5, 0, "tIjAb");
    get_matvec(&T2, &Tij);
    global_dpd_->buf4_close(&T2);

    // Create Density
    std::vector<SharedMatrix> Dij;
    temp->zero();
    for(int ij=0; ij < npairs; ++ij) {
        temp->zero();
        int i = ij/nocc;
        int j = ij/nocc;
        temp->gemm(0, 1, 1, Tij[ij], Ttij[ij], 0);
        temp->gemm(1, 0, 1, Tij[ij], Ttij[ij], 1);
        temp->scale(2.0 * 1/(1+(i==j)));
        Dij.push_back(temp->clone());
    }

    // Print checking

    /* outfile->Printf("\t***** Dij Matrix *****\n");
    for(int ij=0; ij < nocc*nocc; ++ij) {
        int i = ij/nocc;
        int j = ij%nocc;
        outfile->Printf("ij = %2d, i = %2d, j = %2d", ij, i, j);
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

    /*for(int ij=0; ij < npairs; ++ij) {
        outfile->Printf( "Pair: %d\n", ij);
        for(int a=0; a < nvir; ++a) {
            outfile->Printf("%20.12lf\n", occ_num[ij]->get(a));
        }
    }*/

    // Identify survivors
    double abs_occ;
    int survivors;
    auto survivor_list = std::make_shared<Vector>(npairs);
    for(int ij=0; ij < npairs; ++ij) {
        survivors = 0;
        for(int a=0; a < nvir; ++a) {
            abs_occ = fabs(occ_num[ij]->get(a));
            if( abs_occ >= pno_cut) {
                ++survivors;
            }
        }
        survivor_list->set(ij, survivors);
    }

    // Compute stats
    int total_pno = 0;
    double t2_ratio = 0.0;
    for(int ij=0; ij < npairs; ++ij) {
        total_pno += survivor_list->get(ij);
        t2_ratio += survivor_list->get(ij)*survivor_list->get(ij);
    }
    double avg_pno = total_pno / npairs;
    t2_ratio /= (nocc*nocc*nvir*nvir);

    // Print stats
    outfile->Printf("\nTotal number of PNOs: %i \n", total_pno);
    outfile->Printf("Average number of PNOs: %d \n", avg_pno);
    outfile->Printf("T2 ratio: %d \n", t2_ratio);

}


void CCEnergyWavefunction::get_matvec(dpdbuf4 *buf_obj, std::vector<SharedMatrix> *matvec) {
    
    auto nocc = local_.nocc;
    auto nvir = local_.nvir;
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
    // Print check
    outfile->Printf("*** IJ Matrices ***\n");
    for(int ij=0; ij<matvec->size(); ++ij)
        matvec->at(ij)->print();
}

void CCEnergyWavefunction::local_done() { outfile->Printf("    Local parameters free.\n"); }

void CCEnergyWavefunction::local_filter_T1(dpdfile2 *T1) {
    int ii;
    double *T1tilde, *T1bar;
    psio_address next;

    auto nocc = local_.nocc;
    auto nvir = local_.nvir;

    /*   local.weak_pairs = init_int_array(nocc*nocc); */
    local_.pairdom_len = init_int_array(nocc * nocc);
    local_.pairdom_nrlen = init_int_array(nocc * nocc);
    local_.eps_occ = init_array(nocc);
    /*   psio_read_entry(CC_INFO, "Local Weak Pairs", (char *) local.weak_pairs, */
    /* 		  local.nocc*local.nocc*sizeof(int)); */
    psio_read_entry(PSIF_CC_INFO, "Local Pair Domain Length", (char *)local_.pairdom_len, sizeof(int) * nocc * nocc);
    psio_read_entry(PSIF_CC_INFO, "Local Pair Domain NR Length", (char *)local_.pairdom_nrlen,
                    sizeof(int) * nocc * nocc);
    psio_read_entry(PSIF_CC_INFO, "Local Occupied Orbital Energies", (char *)local_.eps_occ, nocc * sizeof(double));

    local_.W = (double ***)malloc(sizeof(double **) * nocc * nocc);
    local_.V = (double ***)malloc(sizeof(double **) * nocc * nocc);
    local_.eps_vir = (double **)malloc(sizeof(double *) * nocc * nocc);
    next = PSIO_ZERO;
    for (int ij = 0; ij < nocc * nocc; ij++) {
        local_.eps_vir[ij] = init_array(local_.pairdom_nrlen[ij]);
        psio_read(PSIF_CC_INFO, "Local Virtual Orbital Energies", (char *)local_.eps_vir[ij],
                  local_.pairdom_nrlen[ij] * sizeof(double), next, &next);
    }
    next = PSIO_ZERO;
    for (int ij = 0; ij < nocc * nocc; ij++) {
        local_.V[ij] = block_matrix(nvir, local_.pairdom_len[ij]);
        psio_read(PSIF_CC_INFO, "Local Residual Vector (V)", (char *)local_.V[ij][0],
                  sizeof(double) * nvir * local_.pairdom_len[ij], next, &next);
    }
    next = PSIO_ZERO;
    for (int ij = 0; ij < nocc * nocc; ij++) {
        local_.W[ij] = block_matrix(local_.pairdom_len[ij], local_.pairdom_nrlen[ij]);
        psio_read(PSIF_CC_INFO, "Local Transformation Matrix (W)", (char *)local_.W[ij][0],
                  sizeof(double) * local_.pairdom_len[ij] * local_.pairdom_nrlen[ij], next, &next);
    }

    global_dpd_->file2_mat_init(T1);
    global_dpd_->file2_mat_rd(T1);

    for (int i = 0; i < nocc; i++) {
        ii = i * nocc + i; /* diagonal element of pair matrices */

        if (!local_.pairdom_len[ii]) {
            outfile->Printf("\n    local_filter_T1: Pair ii = [%d] is zero-length, which makes no sense.\n", ii);
            throw PsiException("local_filter_T1: Pair ii is zero-length, which makes no sense.", __FILE__, __LINE__);
        }

        T1tilde = init_array(local_.pairdom_len[ii]);
        T1bar = init_array(local_.pairdom_nrlen[ii]);

        /* Transform the virtuals to the redundant projected virtual basis */
        C_DGEMV('t', nvir, local_.pairdom_len[ii], 1.0, &(local_.V[ii][0][0]), local_.pairdom_len[ii],
                &(T1->matrix[0][i][0]), 1, 0.0, &(T1tilde[0]), 1);

        /* Transform the virtuals to the non-redundant virtual basis */
        C_DGEMV('t', local_.pairdom_len[ii], local_.pairdom_nrlen[ii], 1.0, &(local_.W[ii][0][0]),
                local_.pairdom_nrlen[ii], &(T1tilde[0]), 1, 0.0, &(T1bar[0]), 1);

        /* Apply the denominators */
        for (int a = 0; a < local_.pairdom_nrlen[ii]; a++) T1bar[a] /= (local_.eps_occ[i] - local_.eps_vir[ii][a]);

        /* Transform the new T1's to the redundant projected virtual basis */
        C_DGEMV('n', local_.pairdom_len[ii], local_.pairdom_nrlen[ii], 1.0, &(local_.W[ii][0][0]),
                local_.pairdom_nrlen[ii], &(T1bar[0]), 1, 0.0, &(T1tilde[0]), 1);

        /* Transform the new T1's to the MO basis */
        C_DGEMV('n', nvir, local_.pairdom_len[ii], 1.0, &(local_.V[ii][0][0]), local_.pairdom_len[ii], &(T1tilde[0]), 1,
                0.0, &(T1->matrix[0][i][0]), 1);

        free(T1bar);
        free(T1tilde);
    }

    global_dpd_->file2_mat_wrt(T1);
    global_dpd_->file2_mat_close(T1);

    for (int ij = 0; ij < nocc * nocc; ij++) {
        free_block(local_.W[ij]);
        free_block(local_.V[ij]);
        free(local_.eps_vir[ij]);
    }
    free(local_.W);
    free(local_.V);
    free(local_.eps_vir);

    free(local_.eps_occ);
    free(local_.pairdom_len);
    free(local_.pairdom_nrlen);
    /*   free(local.weak_pairs); */
}

void CCEnergyWavefunction::local_filter_T2(dpdbuf4 *T2) {
    psio_address next;

    auto nso = local_.nso;
    auto nocc = local_.nocc;
    auto nvir = local_.nvir;

    /*   local.weak_pairs = init_int_array(nocc*nocc); */
    local_.pairdom_len = init_int_array(nocc * nocc);
    local_.pairdom_nrlen = init_int_array(nocc * nocc);
    local_.eps_occ = init_array(nocc);
    /*   psio_read_entry(CC_INFO, "Local Weak Pairs", (char *) local.weak_pairs, */
    /* 		  local.nocc*local.nocc*sizeof(int)); */
    psio_read_entry(PSIF_CC_INFO, "Local Pair Domain Length", (char *)local_.pairdom_len, sizeof(int) * nocc * nocc);
    psio_read_entry(PSIF_CC_INFO, "Local Pair Domain NR Length", (char *)local_.pairdom_nrlen,
                    sizeof(int) * nocc * nocc);
    psio_read_entry(PSIF_CC_INFO, "Local Occupied Orbital Energies", (char *)local_.eps_occ, nocc * sizeof(double));

    local_.W = (double ***)malloc(sizeof(double **) * nocc * nocc);
    local_.V = (double ***)malloc(sizeof(double **) * nocc * nocc);
    local_.eps_vir = (double **)malloc(sizeof(double *) * nocc * nocc);
    next = PSIO_ZERO;
    for (int ij = 0; ij < nocc * nocc; ij++) {
        local_.eps_vir[ij] = init_array(local_.pairdom_nrlen[ij]);
        psio_read(PSIF_CC_INFO, "Local Virtual Orbital Energies", (char *)local_.eps_vir[ij],
                  local_.pairdom_nrlen[ij] * sizeof(double), next, &next);
    }
    next = PSIO_ZERO;
    for (int ij = 0; ij < nocc * nocc; ij++) {
        local_.V[ij] = block_matrix(nvir, local_.pairdom_len[ij]);
        psio_read(PSIF_CC_INFO, "Local Residual Vector (V)", (char *)local_.V[ij][0],
                  sizeof(double) * nvir * local_.pairdom_len[ij], next, &next);
    }
    next = PSIO_ZERO;
    for (int ij = 0; ij < nocc * nocc; ij++) {
        local_.W[ij] = block_matrix(local_.pairdom_len[ij], local_.pairdom_nrlen[ij]);
        psio_read(PSIF_CC_INFO, "Local Transformation Matrix (W)", (char *)local_.W[ij][0],
                  sizeof(double) * local_.pairdom_len[ij] * local_.pairdom_nrlen[ij], next, &next);
    }

    /* Grab the MO-basis T2's */
    global_dpd_->buf4_mat_irrep_init(T2, 0);
    global_dpd_->buf4_mat_irrep_rd(T2, 0);

    auto X1 = block_matrix(nso, nvir);
    auto X2 = block_matrix(nvir, nso);
    auto T2tilde = block_matrix(nso, nso);
    auto T2bar = block_matrix(nvir, nvir);

    for (int i = 0, ij = 0; i < nocc; i++) {
        for (int j = 0; j < nocc; j++, ij++) {
            if (!local_.weak_pairs[ij]) {
                /* Transform the virtuals to the redundant projected virtual basis */
                C_DGEMM('t', 'n', local_.pairdom_len[ij], nvir, nvir, 1.0, &(local_.V[ij][0][0]),
                        local_.pairdom_len[ij], &(T2->matrix[0][ij][0]), nvir, 0.0, &(X1[0][0]), nvir);
                C_DGEMM('n', 'n', local_.pairdom_len[ij], local_.pairdom_len[ij], nvir, 1.0, &(X1[0][0]), nvir,
                        &(local_.V[ij][0][0]), local_.pairdom_len[ij], 0.0, &(T2tilde[0][0]), nso);

                /* Transform the virtuals to the non-redundant virtual basis */
                C_DGEMM('t', 'n', local_.pairdom_nrlen[ij], local_.pairdom_len[ij], local_.pairdom_len[ij], 1.0,
                        &(local_.W[ij][0][0]), local_.pairdom_nrlen[ij], &(T2tilde[0][0]), nso, 0.0, &(X2[0][0]), nso);
                C_DGEMM('n', 'n', local_.pairdom_nrlen[ij], local_.pairdom_nrlen[ij], local_.pairdom_len[ij], 1.0,
                        &(X2[0][0]), nso, &(local_.W[ij][0][0]), local_.pairdom_nrlen[ij], 0.0, &(T2bar[0][0]), nvir);

                /* Divide the new amplitudes by the denominators */
                for (int a = 0; a < local_.pairdom_nrlen[ij]; a++)
                    for (int b = 0; b < local_.pairdom_nrlen[ij]; b++)
                        T2bar[a][b] /=
                            (local_.eps_occ[i] + local_.eps_occ[j] - local_.eps_vir[ij][a] - local_.eps_vir[ij][b]);

                /* Transform the new T2's to the redundant virtual basis */
                C_DGEMM('n', 'n', local_.pairdom_len[ij], local_.pairdom_nrlen[ij], local_.pairdom_nrlen[ij], 1.0,
                        &(local_.W[ij][0][0]), local_.pairdom_nrlen[ij], &(T2bar[0][0]), nvir, 0.0, &(X1[0][0]), nvir);
                C_DGEMM('n', 't', local_.pairdom_len[ij], local_.pairdom_len[ij], local_.pairdom_nrlen[ij], 1.0,
                        &(X1[0][0]), nvir, &(local_.W[ij][0][0]), local_.pairdom_nrlen[ij], 0.0, &(T2tilde[0][0]), nso);

                /* Transform the new T2's to the MO basis */
                C_DGEMM('n', 'n', nvir, local_.pairdom_len[ij], local_.pairdom_len[ij], 1.0, &(local_.V[ij][0][0]),
                        local_.pairdom_len[ij], &(T2tilde[0][0]), nso, 0.0, &(X2[0][0]), nso);
                C_DGEMM('n', 't', nvir, nvir, local_.pairdom_len[ij], 1.0, &(X2[0][0]), nso, &(local_.V[ij][0][0]),
                        local_.pairdom_len[ij], 0.0, &(T2->matrix[0][ij][0]), nvir);
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

    for (int i = 0; i < nocc * nocc; i++) {
        free_block(local_.W[i]);
        free_block(local_.V[i]);
        free(local_.eps_vir[i]);
    }
    free(local_.W);
    free(local_.V);
    free(local_.eps_vir);

    free(local_.eps_occ);
    free(local_.pairdom_len);
    free(local_.pairdom_nrlen);
    /*   free(local.weak_pairs); */
}
}  // namespace ccenergy
}  // namespace psi
