/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2022 The Psi4 Developers.
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

#include "Params.h"
#include "MOInfo.h"
#include "psi4/cc/ccwave.h"

#include "psi4/libpsi4util/PsiOutStream.h"
#include "psi4/libciomr/libciomr.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libiwl/iwl.h"
#include "psi4/libqt/qt.h"
#include "psi4/libmints/wavefunction.h"
#include "psi4/psifiles.h"

#include "psi4/libmints/wavefunction.h"
#include "psi4/libtrans/mospace.h"
#include "psi4/libmints/matrix.h"

#include <cstdio>
#include <cstdlib>
#include <cmath>

namespace psi {
namespace ccenergy {

/* rotate(): Rotate the orbitals using a linear transformation matrix
** built from converged T1 amplitudes.  I still need to add spin-restricted
** Brueckner rotations from my 1997 paper.
**
** TDC, 5/03
*/

bool CCEnergyWavefunction::rotate() {
    dpdfile2 T1;
    double **U, **S, **X;
    double *evals, *work, **MO_S;
    double **scf, **scf_new, **scf_a, **scf_b;
    double **scf_orig, **scf_a_orig, **scf_b_orig;
    double **D, **D_a, **D_b;          /* SCF densities */
    double **fock, **fock_a, **fock_b; /* Fock matrices (SO or MO basis) */
    double ***Foo, ***Fvv;             /* occ-occ and vir-vir block of Fock matrix */
    int phase_ok = 1, max_col;

    auto nso = moinfo_.nso;
    auto nmo = moinfo_.nmo;
    auto offset = init_int_array(nirrep_);
    for (int h = 1; h < nirrep_; h++) offset[h] = offset[h - 1] + moinfo_.orbspi[h - 1];

    // => Check if orbitals are converged <=
    // ==> Compute max T amplitude (infinity norm) <==
    auto max = 0.0;
    if (params_.ref == 0) { /** RHF **/
        global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
        Matrix t1(&T1);
        max = t1.absmax();
        global_dpd_->file2_close(&T1);
    } else if (params_.ref == 2) { /** UHF **/

        global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
        Matrix t1_a(&T1);
        max = t1_a.absmax();
        global_dpd_->file2_close(&T1);

        global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
        Matrix t1_b(&T1);
        max = std::max(max, t1_b.absmax());
        global_dpd_->file2_close(&T1);
    }

    // ==> Check max T amplitude against threshhold <==
    if (std::fabs(max) <= params_.bconv) {
        outfile->Printf("    Brueckner orbitals converged.  Maximum T1 = %15.12f\n", std::fabs(max));
        return true;
    } else
        outfile->Printf("    Rotating orbitals.  Maximum T1 = %15.12f\n", std::fabs(max));

    /* grab the SO-basis overlap integrals for later use */
    auto SO_S = S_->to_block_matrix();

    int stat = 0;
    if (params_.ref == 0) { /* RHF */

        U = block_matrix(nmo, nmo);
        for (int i = 0; i < nmo; i++) U[i][i] = 1.0;

        global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
        global_dpd_->file2_mat_init(&T1);
        global_dpd_->file2_mat_rd(&T1);
        for (int h = 0; h < nirrep_; h++) {
            for (int i = 0; i < moinfo_.occpi[h]; i++) {
                auto ii = moinfo_.qt2pitzer[moinfo_.qt_occ[i] + moinfo_.occ_off[h]];
                for (int a = 0; a < moinfo_.virtpi[h]; a++) {
                    auto aa = moinfo_.qt2pitzer[moinfo_.qt_vir[a] + moinfo_.vir_off[h]];

                    U[ii][aa] = T1.matrix[h][i][a];
                    U[aa][ii] = -T1.matrix[h][i][a];
                }
            }
        }
        global_dpd_->file2_mat_close(&T1);
        global_dpd_->file2_close(&T1);

        scf = Ca_->to_block_matrix();
        scf_orig = Ca_->to_block_matrix();
        scf_new = block_matrix(nso, nmo);
        C_DGEMM('n', 't', nso, nmo, nmo, 1, &(scf[0][0]), nmo, &(U[0][0]), nmo, 0, &(scf_new[0][0]), nmo);
        free_block(U);
        free_block(scf);

        /* transform the overlap into the new MO basis */
        MO_S = block_matrix(nmo, nmo);
        X = block_matrix(nso, nso);
        C_DGEMM('t', 'n', nmo, nso, nso, 1, &(scf_new[0][0]), nmo, &(SO_S[0][0]), nso, 0, &(X[0][0]), nso);
        C_DGEMM('n', 'n', nmo, nmo, nso, 1, &(X[0][0]), nso, &(scf_new[0][0]), nmo, 0, &(MO_S[0][0]), nmo);
        free_block(X);

        /* build S^-1/2 for this basis */
        evals = init_array(nmo);
        work = init_array(nmo * 3);
        if ((stat = C_DSYEV('v', 'u', nmo, &(MO_S[0][0]), nmo, evals, work, nmo * 3))) {
            outfile->Printf("rotate(): Error in overlap diagonalization. stat = %d\n", stat);
            throw PsiException("rotate(): Error in overlap diagonalization.", __FILE__, __LINE__);
        }
        S = block_matrix(nmo, nmo);
        for (int i = 0; i < nmo; i++) {
            if (std::fabs(evals[i]) > 1e-8)
                S[i][i] = 1 / sqrt(evals[i]);
            else
                S[i][i] = 0.0;
        }
        free(evals);
        free(work);
        X = block_matrix(nmo, nmo);
        C_DGEMM('t', 'n', nmo, nmo, nmo, 1, &(MO_S[0][0]), nso, &(S[0][0]), nmo, 0, &(X[0][0]), nmo);
        C_DGEMM('n', 'n', nmo, nmo, nmo, 1, &(X[0][0]), nmo, &(MO_S[0][0]), nso, 0, &(S[0][0]), nmo);
        free_block(X);

        /* orthogonalize the new MO basis */
        scf = block_matrix(nso, nmo);
        C_DGEMM('n', 'n', nmo, nmo, nmo, 1, &(scf_new[0][0]), nmo, &(S[0][0]), nmo, 0, &(scf[0][0]), nmo);
        free_block(S);
        free_block(MO_S);
        free_block(scf_new);

        /* build the SO-basis density for the new MOs */
        D = block_matrix(nso, nso);
        for (int h = 0; h < nirrep_; h++)
            for (int p = offset[h]; p < offset[h] + moinfo_.orbspi[h]; p++)
                for (int q = offset[h]; q < offset[h] + moinfo_.orbspi[h]; q++)
                    for (int i = offset[h]; i < offset[h] + moinfo_.frdocc[h] + moinfo_.occpi[h]; i++)
                        D[p][q] += scf[p][i] * scf[q][i];

        /* build the SO-basis Fock matrix */
        fock = block_matrix(nso, nso);
        rhf_fock_build(fock, D);
        free_block(D);

        Fa_->set(fock);
        Fb_->set(fock);

        /*
    outfile->Printf( "\n    SO-basis Fock matrix:\n");
    mat_print(fock, nso, nso, outfile);
    */

        /* transform the fock matrix to the new MO basis */
        X = block_matrix(nso, nso);
        C_DGEMM('n', 'n', nso, nmo, nso, 1.0, &(fock[0][0]), nso, &(scf[0][0]), nmo, 0, &(X[0][0]), nso);
        C_DGEMM('t', 'n', nmo, nmo, nso, 1.0, &(scf[0][0]), nmo, &(X[0][0]), nso, 0, &(fock[0][0]), nso);
        free_block(X);

        /*
    outfile->Printf( "\n    MO-basis Fock matrix:\n");
    mat_print(fock, nmo, nmo, outfile);
    */

        /* extract the occ-occ and vir-vir block of the Fock matrix */
        Foo = (double ***)malloc(nirrep_ * sizeof(double **));
        Fvv = (double ***)malloc(nirrep_ * sizeof(double **));
        X = block_matrix(nmo, nmo);
        for (int h = 0; h < nirrep_; h++) {
            /* leave the frozen orbitals alone */
            for (int i = offset[h]; i < offset[h] + moinfo_.frdocc[h]; i++) X[i][i] = 1.0;
            for (int end = offset[h] + moinfo_.orbspi[h], start = end - moinfo_.fruocc[h],
                 i = start; i < end; i++)
                X[i][i] = 1.0;

            Foo[h] = block_matrix(moinfo_.occpi[h], moinfo_.occpi[h]);
            Fvv[h] = block_matrix(moinfo_.virtpi[h], moinfo_.virtpi[h]);

            for (int i = offset[h] + moinfo_.frdocc[h], I = 0; i < offset[h] + moinfo_.frdocc[h] + moinfo_.occpi[h];
                 i++, I++)
                for (int j = offset[h] + moinfo_.frdocc[h], J = 0; j < offset[h] + moinfo_.frdocc[h] + moinfo_.occpi[h];
                     j++, J++)
                    Foo[h][I][J] = fock[i][j];

            for (int start = offset[h] + moinfo_.frdocc[h] + moinfo_.occpi[h],
                 end = start + moinfo_.virtpi[h], a = start, A = 0; a < end; a++, A++)
                for (int b = start, B = 0; b < end; b++, B++)
                    Fvv[h][A][B] = fock[a][b];

            /*
      outfile->Printf( "\n    Occ-occ Fock matrix for irrep %d:\n", h);
      mat_print(Foo[h], moinfo.occpi[h], moinfo.occpi[h], outfile);

      outfile->Printf( "\n    Vir-vir Fock matrix for irrep %d:\n", h);
      mat_print(Fvv[h], moinfo.virtpi[h], moinfo.virtpi[h], outfile);
      */

            if (moinfo_.occpi[h]) {
                evals = init_array(moinfo_.occpi[h]);
                work = init_array(3 * moinfo_.occpi[h]);
                if ((stat = C_DSYEV('v', 'u', moinfo_.occpi[h], &(Foo[h][0][0]), moinfo_.occpi[h], evals, work,
                                    moinfo_.occpi[h] * 3))) {
                    outfile->Printf("rotate(): Error in Foo[%1d] diagonalization. stat = %d\n", h, stat);
                    throw PsiException("rotate(): Error in Foo diagonalization.", __FILE__, __LINE__);
                }
                free(evals);
                free(work);

                /*
    outfile->Printf( "\n    Eigenfunctions of Occ-occ Fock matrix for irrep %1d:\n", h);
    mat_print(Foo[h], moinfo.occpi[h], moinfo.occpi[h], outfile);
    */

                for (int i = offset[h] + moinfo_.frdocc[h], I = 0; i < offset[h] + moinfo_.frdocc[h] + moinfo_.occpi[h];
                     i++, I++)
                    for (int j = offset[h] + moinfo_.frdocc[h], J = 0;
                         j < offset[h] + moinfo_.frdocc[h] + moinfo_.occpi[h]; j++, J++)
                        X[i][j] = Foo[h][J][I];
            }

            if (moinfo_.virtpi[h]) {
                evals = init_array(moinfo_.virtpi[h]);
                work = init_array(3 * moinfo_.virtpi[h]);
                if ((stat = C_DSYEV('v', 'u', moinfo_.virtpi[h], &(Fvv[h][0][0]), moinfo_.virtpi[h], evals, work,
                                    moinfo_.virtpi[h] * 3))) {
                    outfile->Printf("rotate(): Error in Fvv[%1d] diagonalization. stat = %d\n", h, stat);
                    throw PsiException("rotate(): Error in Foo diagonalization.", __FILE__, __LINE__);
                }
                free(evals);
                free(work);

                /*
    outfile->Printf( "\n    Eigenfunctions of Vir-vir Fock matrix for irrep %1d:\n", h);
    mat_print(Fvv[h], moinfo.virtpi[h], moinfo.virtpi[h], outfile);
    */

                for (int start = offset[h] + moinfo_.frdocc[h] + moinfo_.occpi[h],
                     end = start + moinfo_.virtpi[h],
                     a = start, A = 0; a < end; a++, A++)
                    for (int b = start, B = 0; b < end; b++, B++)
                        X[a][b] = Fvv[h][B][A];
            }

            free_block(Foo[h]);
            free_block(Fvv[h]);
        }
        free(Foo);
        free(Fvv);
        free_block(fock);

        /* semicanonicalization of the basis */
        /*
    outfile->Printf( "\n    Semicanonical transformation matrix:\n");
    mat_print(X, nmo, nmo, outfile);
    */

        scf_new = block_matrix(nso, nmo);
        C_DGEMM('n', 'n', nso, nmo, nmo, 1, &(scf[0][0]), nmo, &(X[0][0]), nmo, 0, &(scf_new[0][0]), nmo);
        free_block(X);
        free_block(scf);

        /* Reorder new MO's to Pitzer and write to wfn */
        /*
    outfile->Printf( "\n    Semicanonical Brueckner orbitals (Pitzer order):\n");
    mat_print(scf_new, nso, nmo, outfile);
    */

        /* correct orbital phases for amplitude restarts */
        MO_S = block_matrix(nmo, nmo);
        X = block_matrix(nso, nmo);
        C_DGEMM('n', 'n', nso, nmo, nso, 1, &(SO_S[0][0]), nso, &(scf_new[0][0]), nmo, 0, &(X[0][0]), nmo);
        C_DGEMM('t', 'n', nmo, nmo, nso, 1, &(scf_orig[0][0]), nmo, &(X[0][0]), nmo, 0, &(MO_S[0][0]), nmo);
        free_block(X);

        for (int p = 0; p < nmo; p++) {
            max = 0.0;
            for (int q = 0; q < nmo; q++) {
                if (std::fabs(MO_S[p][q]) > max) {
                    max = std::fabs(MO_S[p][q]);
                    max_col = q;
                }
            }
            if (max_col != p) phase_ok = 0;
        }

        if (phase_ok) {
            for (int p = 0; p < nmo; p++) {
                if (MO_S[p][p] < 0.0) {
                    for (int q = 0; q < nso; q++) scf_new[q][p] *= -1.0;
                }
            }
        }

        free_block(MO_S);

        /*
    outfile->Printf( "\n    Original SCF MOs:\n");
    mat_print(scf_orig, nso, nmo, outfile);
    outfile->Printf( "\n    New SCF MOs:\n");
    mat_print(scf_new, nso, nmo, outfile);
    */

        Ca_->set(scf_new);
        Cb_->set(scf_new);

        free_block(scf_new);
        free_block(scf_orig);

    } else if (params_.ref == 2) { /* UHF */

        /* AA block */
        U = block_matrix(nmo, nmo);
        for (int i = 0; i < nmo; i++) U[i][i] = 1.0;

        global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
        global_dpd_->file2_mat_init(&T1);
        global_dpd_->file2_mat_rd(&T1);
        for (int h = 0; h < nirrep_; h++) {
            for (int i = 0; i < moinfo_.aoccpi[h]; i++) {
                auto ii = moinfo_.qt2pitzer_a[moinfo_.qt_aocc[i] + moinfo_.aocc_off[h]];
                for (int a = 0; a < moinfo_.avirtpi[h]; a++) {
                    auto aa = moinfo_.qt2pitzer_a[moinfo_.qt_avir[a] + moinfo_.avir_off[h]];

                    U[ii][aa] = T1.matrix[h][i][a];
                    U[aa][ii] = -T1.matrix[h][i][a];
                }
            }
        }
        global_dpd_->file2_mat_close(&T1);
        global_dpd_->file2_close(&T1);

        scf = Ca_->to_block_matrix();
        scf_a_orig = Ca_->to_block_matrix();

        scf_new = block_matrix(nso, nmo);
        C_DGEMM('n', 't', nso, nmo, nmo, 1, &(scf[0][0]), nmo, &(U[0][0]), nmo, 0, &(scf_new[0][0]), nmo);
        free_block(U);
        free_block(scf);

        MO_S = block_matrix(nmo, nmo);

        /* transform the overlap into the new alpha MO basis */
        X = block_matrix(nso, nso);
        C_DGEMM('t', 'n', nmo, nso, nso, 1, &(scf_new[0][0]), nmo, &(SO_S[0][0]), nso, 0, &(X[0][0]), nso);
        C_DGEMM('n', 'n', nmo, nmo, nso, 1, &(X[0][0]), nso, &(scf_new[0][0]), nmo, 0, &(MO_S[0][0]), nmo);
        free_block(X);

        evals = init_array(nmo);
        work = init_array(nmo * 3);
        if ((stat = C_DSYEV('v', 'u', nmo, &(MO_S[0][0]), nmo, evals, work, nmo * 3))) {
            outfile->Printf("rotate(): Error in overlap diagonalization. stat = %d\n", stat);
            throw PsiException("rotate(): Error in overlap diagonalization.", __FILE__, __LINE__);
        }

        /* build S^-1/2 for this basis */
        S = block_matrix(nmo, nmo);
        for (int i = 0; i < nmo; i++) {
            if (std::fabs(evals[i]) > 1e-8)
                S[i][i] = 1 / sqrt(evals[i]);
            else
                S[i][i] = 0.0;
        }
        free(evals);
        free(work);
        X = block_matrix(nmo, nmo);
        C_DGEMM('t', 'n', nmo, nmo, nmo, 1, &(MO_S[0][0]), nso, &(S[0][0]), nmo, 0, &(X[0][0]), nmo);
        C_DGEMM('n', 'n', nmo, nmo, nmo, 1, &(X[0][0]), nmo, &(MO_S[0][0]), nso, 0, &(S[0][0]), nmo);
        free_block(X);

        /* orthogonalize the basis */
        scf_a = block_matrix(nso, nmo);
        C_DGEMM('n', 'n', nmo, nmo, nmo, 1, &(scf_new[0][0]), nmo, &(S[0][0]), nmo, 0, &(scf_a[0][0]), nmo);
        free_block(S);
        free_block(MO_S);
        free_block(scf_new);

        /* BB block */
        U = block_matrix(nmo, nmo);
        for (int i = 0; i < nmo; i++) U[i][i] = 1.0;

        global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
        global_dpd_->file2_mat_init(&T1);
        global_dpd_->file2_mat_rd(&T1);
        for (int h = 0; h < nirrep_; h++) {
            for (int i = 0; i < moinfo_.boccpi[h]; i++) {
                auto ii = moinfo_.qt2pitzer_b[moinfo_.qt_bocc[i] + moinfo_.bocc_off[h]];
                for (int a = 0; a < moinfo_.bvirtpi[h]; a++) {
                    auto aa = moinfo_.qt2pitzer_b[moinfo_.qt_bvir[a] + moinfo_.bvir_off[h]];

                    U[ii][aa] = T1.matrix[h][i][a];
                    U[aa][ii] = -T1.matrix[h][i][a];
                }
            }
        }
        global_dpd_->file2_mat_close(&T1);
        global_dpd_->file2_close(&T1);

        scf = Cb_->to_block_matrix();
        scf_b_orig = Cb_->to_block_matrix();

        scf_new = block_matrix(nso, nmo);
        C_DGEMM('n', 't', nso, nmo, nmo, 1, &(scf[0][0]), nmo, &(U[0][0]), nmo, 0, &(scf_new[0][0]), nmo);
        free_block(U);
        free_block(scf);

        MO_S = block_matrix(nmo, nmo);

        /* transform the overlap into the new beta MO basis */
        X = block_matrix(nso, nso);
        C_DGEMM('t', 'n', nmo, nso, nso, 1, &(scf_new[0][0]), nmo, &(SO_S[0][0]), nso, 0, &(X[0][0]), nso);
        C_DGEMM('n', 'n', nmo, nmo, nso, 1, &(X[0][0]), nso, &(scf_new[0][0]), nmo, 0, &(MO_S[0][0]), nmo);
        free_block(X);

        evals = init_array(nmo);
        work = init_array(nmo * 3);
        if ((stat = C_DSYEV('v', 'u', nmo, &(MO_S[0][0]), nmo, evals, work, nmo * 3))) {
            outfile->Printf("rotate(): Error in overlap diagonalization. stat = %d\n", stat);
            throw PsiException("rotate(): Error in Foo diagonalization.", __FILE__, __LINE__);
        }

        /* build S^-1/2 for this basis */
        S = block_matrix(nmo, nmo);
        for (int i = 0; i < nmo; i++) {
            if (std::fabs(evals[i]) > 1e-8)
                S[i][i] = 1 / sqrt(evals[i]);
            else
                S[i][i] = 0.0;
        }
        free(evals);
        free(work);
        X = block_matrix(nmo, nmo);
        C_DGEMM('t', 'n', nmo, nmo, nmo, 1, &(MO_S[0][0]), nso, &(S[0][0]), nmo, 0, &(X[0][0]), nmo);
        C_DGEMM('n', 'n', nmo, nmo, nmo, 1, &(X[0][0]), nmo, &(MO_S[0][0]), nso, 0, &(S[0][0]), nmo);
        free_block(X);

        /* orthogonalize the basis */
        scf_b = block_matrix(nso, nmo);
        C_DGEMM('n', 'n', nmo, nmo, nmo, 1, &(scf_new[0][0]), nmo, &(S[0][0]), nmo, 0, &(scf_b[0][0]), nmo);
        free_block(S);
        free_block(MO_S);
        free_block(scf_new);

        /* build the SO-basis alpha and beta densities for the new MOs */
        D_a = block_matrix(nso, nso);
        D_b = block_matrix(nso, nso);
        for (int h = 0; h < nirrep_; h++)
            for (int p = offset[h]; p < offset[h] + moinfo_.orbspi[h]; p++)
                for (int q = offset[h]; q < offset[h] + moinfo_.orbspi[h]; q++) {
                    for (int i = offset[h]; i < offset[h] + moinfo_.frdocc[h] + moinfo_.aoccpi[h]; i++)
                        D_a[p][q] += scf_a[p][i] * scf_a[q][i];
                    for (int i = offset[h]; i < offset[h] + moinfo_.frdocc[h] + moinfo_.boccpi[h]; i++)
                        D_b[p][q] += scf_b[p][i] * scf_b[q][i];
                }

        /* build the alpha and beta SO-basis Fock matrices */
        fock_a = block_matrix(nso, nso);
        fock_b = block_matrix(nso, nso);
        uhf_fock_build(fock_a, fock_b, D_a, D_b);
        free_block(D_a);
        free_block(D_b);

        Fa_->set(fock_a);
        Fb_->set(fock_b);

        /* transform the fock matrices to the new alpha and beta MO bases */
        X = block_matrix(nso, nso);
        C_DGEMM('n', 'n', nso, nmo, nso, 1.0, &(fock_a[0][0]), nso, &(scf_a[0][0]), nmo, 0, &(X[0][0]), nso);
        C_DGEMM('t', 'n', nmo, nmo, nso, 1.0, &(scf_a[0][0]), nmo, &(X[0][0]), nso, 0, &(fock_a[0][0]), nso);

        C_DGEMM('n', 'n', nso, nmo, nso, 1.0, &(fock_b[0][0]), nso, &(scf_b[0][0]), nmo, 0, &(X[0][0]), nso);
        C_DGEMM('t', 'n', nmo, nmo, nso, 1.0, &(scf_b[0][0]), nmo, &(X[0][0]), nso, 0, &(fock_b[0][0]), nso);
        free_block(X);

        /** alpha Fock semicanonicalization **/

        Foo = (double ***)malloc(nirrep_ * sizeof(double **));
        Fvv = (double ***)malloc(nirrep_ * sizeof(double **));
        X = block_matrix(nmo, nmo);
        for (int h = 0; h < nirrep_; h++) {
            /* leave the frozen orbitals alone */
            for (int i = offset[h]; i < offset[h] + moinfo_.frdocc[h]; i++) X[i][i] = 1.0;
            for (int end = offset[h] + moinfo_.orbspi[h], start = end - moinfo_.fruocc[h],
                 i = start; i < end; i++)
                X[i][i] = 1.0;

            Foo[h] = block_matrix(moinfo_.aoccpi[h], moinfo_.aoccpi[h]);
            Fvv[h] = block_matrix(moinfo_.avirtpi[h], moinfo_.avirtpi[h]);

            for (int i = offset[h] + moinfo_.frdocc[h], I = 0; i < offset[h] + moinfo_.frdocc[h] + moinfo_.aoccpi[h];
                 i++, I++)
                for (int j = offset[h] + moinfo_.frdocc[h], J = 0;
                     j < offset[h] + moinfo_.frdocc[h] + moinfo_.aoccpi[h]; j++, J++)
                    Foo[h][I][J] = fock_a[i][j];

            for (int start = offset[h] + moinfo_.frdocc[h] + moinfo_.aoccpi[h],
                 end = start + moinfo_.avirtpi[h], a = start, A = 0; a < end; a++, A++)
                for (int b = start, B = 0; b < end; b++, B++)
                    Fvv[h][A][B] = fock_a[a][b];

            if (moinfo_.aoccpi[h]) {
                evals = init_array(moinfo_.aoccpi[h]);
                work = init_array(3 * moinfo_.aoccpi[h]);
                if ((stat = C_DSYEV('v', 'u', moinfo_.aoccpi[h], &(Foo[h][0][0]), moinfo_.aoccpi[h], evals, work,
                                    moinfo_.aoccpi[h] * 3))) {
                    outfile->Printf("rotate(): Error in alpha Foo[%1d] diagonalization. stat = %d\n", h, stat);
                    throw PsiException("rotate(): Error in Foo diagonalization.", __FILE__, __LINE__);
                }
                free(evals);
                free(work);

                for (int i = offset[h] + moinfo_.frdocc[h], I = 0;
                     i < offset[h] + moinfo_.frdocc[h] + moinfo_.aoccpi[h]; i++, I++)
                    for (int j = offset[h] + moinfo_.frdocc[h], J = 0;
                         j < offset[h] + moinfo_.frdocc[h] + moinfo_.aoccpi[h]; j++, J++)
                        X[i][j] = Foo[h][J][I];
            }

            if (moinfo_.avirtpi[h]) {
                evals = init_array(moinfo_.avirtpi[h]);
                work = init_array(3 * moinfo_.avirtpi[h]);
                if ((stat = C_DSYEV('v', 'u', moinfo_.avirtpi[h], &(Fvv[h][0][0]), moinfo_.avirtpi[h], evals, work,
                                    moinfo_.avirtpi[h] * 3))) {
                    outfile->Printf("rotate(): Error in alpha Fvv[%1d] diagonalization. stat = %d\n", h, stat);
                    throw PsiException("rotate(): Error in alpha Fvv diagonalization.", __FILE__, __LINE__);
                }
                free(evals);
                free(work);

                for (int start = offset[h] + moinfo_.frdocc[h] + moinfo_.aoccpi[h],
                     end = start + moinfo_.avirtpi[h],
                     a = start, A = 0; a < end; a++, A++)
                    for (int b = start, B = 0; b < end; b++, B++)
                        X[a][b] = Fvv[h][B][A];
            }

            free_block(Foo[h]);
            free_block(Fvv[h]);
        }
        free(Foo);
        free(Fvv);
        free_block(fock_a);

        scf_new = block_matrix(nso, nmo);
        C_DGEMM('n', 'n', nso, nmo, nmo, 1, &(scf_a[0][0]), nmo, &(X[0][0]), nmo, 0, &(scf_new[0][0]), nmo);
        free_block(X);
        free_block(scf_a);

        /* correct orbital phases for amplitude restarts */
        MO_S = block_matrix(nmo, nmo);
        X = block_matrix(nso, nmo);
        C_DGEMM('n', 'n', nso, nmo, nso, 1, &(SO_S[0][0]), nso, &(scf_new[0][0]), nmo, 0, &(X[0][0]), nmo);
        C_DGEMM('t', 'n', nmo, nmo, nso, 1, &(scf_a_orig[0][0]), nmo, &(X[0][0]), nmo, 0, &(MO_S[0][0]), nmo);
        free_block(X);

        for (int p = 0; p < nmo; p++) {
            max = 0.0;
            for (int q = 0; q < nmo; q++) {
                if (std::fabs(MO_S[p][q]) > max) {
                    max = std::fabs(MO_S[p][q]);
                    max_col = q;
                }
            }
            if (max_col != p) phase_ok = 0;
        }

        if (phase_ok) {
            for (int p = 0; p < nmo; p++) {
                if (MO_S[p][p] < 0.0) {
                    for (int q = 0; q < nso; q++) scf_new[q][p] *= -1.0;
                }
            }
        }

        free_block(MO_S);

        Ca_->set(scf_new);

        free_block(scf_new);
        free_block(scf_a_orig);

        /** beta Fock semicanonicalization **/

        Foo = (double ***)malloc(nirrep_ * sizeof(double **));
        Fvv = (double ***)malloc(nirrep_ * sizeof(double **));
        X = block_matrix(nmo, nmo);
        for (int h = 0; h < nirrep_; h++) {
            /* leave the frozen orbitals alone */
            for (int i = offset[h]; i < offset[h] + moinfo_.frdocc[h]; i++) X[i][i] = 1.0;
            for (int end = offset[h] + moinfo_.orbspi[h], start = end - moinfo_.fruocc[h],
                 i = start; i < end; i++)
                X[i][i] = 1.0;

            Foo[h] = block_matrix(moinfo_.boccpi[h], moinfo_.boccpi[h]);
            Fvv[h] = block_matrix(moinfo_.bvirtpi[h], moinfo_.bvirtpi[h]);

            for (int i = offset[h] + moinfo_.frdocc[h], I = 0; i < offset[h] + moinfo_.frdocc[h] + moinfo_.boccpi[h];
                 i++, I++)
                for (int j = offset[h] + moinfo_.frdocc[h], J = 0;
                     j < offset[h] + moinfo_.frdocc[h] + moinfo_.boccpi[h]; j++, J++)
                    Foo[h][I][J] = fock_b[i][j];

            for (int start = offset[h] + moinfo_.frdocc[h] + moinfo_.boccpi[h],
                 end = start + moinfo_.bvirtpi[h],
                 a = start, A = 0; a < end; a++, A++)
                for (int b = start, B = 0; b < end; b++, B++)
                    Fvv[h][A][B] = fock_b[a][b];

            if (moinfo_.boccpi[h]) {
                evals = init_array(moinfo_.boccpi[h]);
                work = init_array(3 * moinfo_.boccpi[h]);
                if ((stat = C_DSYEV('v', 'u', moinfo_.boccpi[h], &(Foo[h][0][0]), moinfo_.boccpi[h], evals, work,
                                    moinfo_.boccpi[h] * 3))) {
                    outfile->Printf("rotate(): Error in alpha Foo[%1d] diagonalization. stat = %d\n", h, stat);
                    throw PsiException("rotate(): Error in alpha Foo diagonalization.", __FILE__, __LINE__);
                }
                free(evals);
                free(work);

                for (int i = offset[h] + moinfo_.frdocc[h], I = 0;
                     i < offset[h] + moinfo_.frdocc[h] + moinfo_.boccpi[h]; i++, I++)
                    for (int j = offset[h] + moinfo_.frdocc[h], J = 0;
                         j < offset[h] + moinfo_.frdocc[h] + moinfo_.boccpi[h]; j++, J++)
                        X[i][j] = Foo[h][J][I];
            }

            if (moinfo_.bvirtpi[h]) {
                evals = init_array(moinfo_.bvirtpi[h]);
                work = init_array(3 * moinfo_.bvirtpi[h]);
                if ((stat = C_DSYEV('v', 'u', moinfo_.bvirtpi[h], &(Fvv[h][0][0]), moinfo_.bvirtpi[h], evals, work,
                                    moinfo_.bvirtpi[h] * 3))) {
                    outfile->Printf("rotate(): Error in alpha Fvv[%1d] diagonalization. stat = %d\n", h, stat);
                    throw PsiException("rotate(): Error in alpha Fvv diagonalization.", __FILE__, __LINE__);
                }
                free(evals);
                free(work);

                for (int start = offset[h] + moinfo_.frdocc[h] + moinfo_.boccpi[h],
                     end = start + moinfo_.bvirtpi[h],
                     a = start, A = 0; a < end; a++, A++)
                    for (int b = start, B = 0; b < end; b++, B++)
                        X[a][b] = Fvv[h][B][A];
            }

            free_block(Foo[h]);
            free_block(Fvv[h]);
        }
        free(Foo);
        free(Fvv);
        free_block(fock_b);

        scf_new = block_matrix(nso, nmo);
        C_DGEMM('n', 'n', nso, nmo, nmo, 1, &(scf_b[0][0]), nmo, &(X[0][0]), nmo, 0, &(scf_new[0][0]), nmo);
        free_block(X);
        free_block(scf_b);

        /* correct orbital phases for amplitude restarts */
        MO_S = block_matrix(nmo, nmo);
        X = block_matrix(nso, nmo);
        C_DGEMM('n', 'n', nso, nmo, nso, 1, &(SO_S[0][0]), nso, &(scf_new[0][0]), nmo, 0, &(X[0][0]), nmo);
        C_DGEMM('t', 'n', nmo, nmo, nso, 1, &(scf_b_orig[0][0]), nmo, &(X[0][0]), nmo, 0, &(MO_S[0][0]), nmo);
        free_block(X);

        for (int p = 0; p < nmo; p++) {
            max = 0.0;
            for (int q = 0; q < nmo; q++) {
                if (std::fabs(MO_S[p][q]) > max) {
                    max = std::fabs(MO_S[p][q]);
                    max_col = q;
                }
            }
            if (max_col != p) phase_ok = 0;
        }

        if (phase_ok) {
            for (int p = 0; p < nmo; p++) {
                if (MO_S[p][p] < 0.0) {
                    for (int q = 0; q < nso; q++) scf_new[q][p] *= -1.0;
                }
            }
        }

        free_block(MO_S);

        Cb_->set(scf_new);

        free_block(scf_new);
        free_block(scf_b_orig);
    }

    free_block(SO_S);
    free(offset);

    return false;
}
}  // namespace ccenergy
}  // namespace psi
