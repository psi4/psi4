/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2017 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

/*! \file
    \ingroup DPD
    \brief Enter brief description of file here
*/

/* T3_RHF_AAA(): Computes all T3(IJK,ABC) amplitudes for a given I, J,
** K combination for input C2, F, and E intermediates.  This function
** will work for AAA or BBB spin cases, with either RHF/ROHF or UHF
** orbitals.
**
** Arguments:
**
**   double ***W1: The target triples amplitudes in an nirreps x AB x
**   C array.  The memory for this must be allocated externally.
**
**   int nirreps: Number of irreps.
**
**   int I: Absolute index of orbital I.
**
**   int Gi: Irrep of I.
**
**   int J: Absolute index of orbital J.
**
**   int Gj: Irrep of J.
**
**   int K: Absolute index of orbital K.
**
**   int Gk: Irrep of K.
**
**   dpdbuf4 *C2: Pointer to dpd buffer for double excitation amps,
**   ordered (IJ,AB).
**
**   dpdbuf4 *F: Pointer to dpd buffer for three-virtual-index
**   intermediate, ordered (IA,BC).
**
**   dpdbuf4 *E: Pointer to dpd buffer for three-occupied-index
**   intermediate, ordered (IJ,KA).
**
**   dpdfile2 *fIJ: Pointer to the dpd file2 for the occ-occ block of
**   the Fock matrix (or other appropriate one-electron operator).
**
**   dpdfile2 *fAB: Pointer to the dpd file2 for the vir-vir block of
**   the Fock matrix (or other appropriate one-electron operator).
**
**   int *occpi: Number of occupied orbitals per irrep lookup array.
**
**   int *occ_off: Offset lookup for translating between absolute and
**   relative orbital indices for occupied space.
**
**   int *virtpi: Number of virtual orbitals per irrep lookup array.
**
**   int *vir_off: Offset lookup for translating between absolute and
**   relative orbital indices for virtual space.
**
**   double omega: a constant to add to the final denominators -
**   needed for CC3 EOM
**
** TDC, July 2004
*/

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "psi4/libqt/qt.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/psifiles.h"

namespace psi {

void DPD::T3_AAA(double ***W1, int nirreps, int I, int Gi, int J, int Gj, int K, int Gk,
                 dpdbuf4 *C2, dpdbuf4 *F, dpdbuf4 *E, dpdfile2 *fIJ, dpdfile2 *fAB,
                 int *occpi, int *occ_off, int *virtpi, int *vir_off, double omega)
{
    int h;
    int i, j, k;
    int ij, ji, ik, ki, jk, kj;
    int Gij, Gji, Gik, Gki, Gjk, Gkj, Gijk;
    int Ga, Gb, Gc;
    int Gd, Gl;
    int Gid, Gjd, Gkd;
    int Gab, Gcb, Gca;
    int Gla, Glb, Glc;
    int Gil, Gjl, Gkl;
    int a, b, c, A, B, C;
    int ab;
    int cd, bd, ad;
    int id, jd, kd;
    int la, lb, lc;
    int il, jl, kl;
    int nrows, ncols, nlinks;
    double dijk, denom;
    double ***W2;
    int GE, GF, GC, GX3;

    GC = C2->file.my_irrep;
    /* F and E are assumed to have same irrep */
    GF = GE =  F->file.my_irrep;
    GX3 = GC^GF;

    file2_mat_init(fIJ);
    file2_mat_init(fAB);
    file2_mat_rd(fIJ);
    file2_mat_rd(fAB);

    for(h=0; h < nirreps; h++) {
        buf4_mat_irrep_init(C2, h);
        buf4_mat_irrep_rd(C2, h);

        buf4_mat_irrep_init(E, h);
        buf4_mat_irrep_rd(E, h);
    }

    i = I - occ_off[Gi];
    j = J - occ_off[Gj];
    k = K - occ_off[Gk];

    Gij = Gji = Gi ^ Gj;
    Gik = Gki = Gi ^ Gk;
    Gjk = Gkj = Gj ^ Gk;
    Gijk = Gi ^ Gj ^ Gk;

    ij = C2->params->rowidx[I][J];
    ji = C2->params->rowidx[J][I];
    jk = C2->params->rowidx[J][K];
    kj = C2->params->rowidx[K][J];
    ik = C2->params->rowidx[I][K];
    ki = C2->params->rowidx[K][I];

    dijk = 0.0;
    if(fIJ->params->rowtot[Gi]) dijk += fIJ->matrix[Gi][i][i];
    if(fIJ->params->rowtot[Gj]) dijk += fIJ->matrix[Gj][j][j];
    if(fIJ->params->rowtot[Gk]) dijk += fIJ->matrix[Gk][k][k];

    W2 = (double ***) malloc(nirreps * sizeof(double **));

    for(Gab=0; Gab < nirreps; Gab++) {
        Gc = Gab ^ Gijk ^ GX3; /* changed */

        W2[Gab] = dpd_block_matrix(F->params->coltot[Gab], virtpi[Gc]);

        if(F->params->coltot[Gab] && virtpi[Gc]) {
            ::memset(W1[Gab][0], 0, F->params->coltot[Gab]*virtpi[Gc]*sizeof(double));
        }
    }

    for(Gd=0; Gd < nirreps; Gd++) {

        /* +t_kjcd * F_idab */
        Gid = Gi ^ Gd;
        Gab = Gid ^ GF; /* changed */

        Gc = Gjk ^ Gd ^ GC; /* changed */

        cd = C2->col_offset[Gjk][Gc];
        id = F->row_offset[Gid][I];

        F->matrix[Gid] = dpd_block_matrix(virtpi[Gd], F->params->coltot[Gid^GF]);
        buf4_mat_irrep_rd_block(F, Gid, id, virtpi[Gd]);

        nrows = F->params->coltot[Gid^GF];
        ncols = virtpi[Gc];
        nlinks = virtpi[Gd];

        if(nrows && ncols && nlinks)
            C_DGEMM('t','t',nrows, ncols, nlinks, 1.0, F->matrix[Gid][0], nrows,
                    &(C2->matrix[Gjk][kj][cd]), nlinks, 1.0, W1[Gab][0], ncols);

        free_dpd_block(F->matrix[Gid], virtpi[Gd], F->params->coltot[Gid^GF]);

        /* +t_ikcd * F_jdab */
        Gjd = Gj ^ Gd;
        Gab = Gjd ^ GF; /* changed */
        Gc = Gik ^ Gd ^ GC;  /* changed */

        cd = C2->col_offset[Gik][Gc];
        jd = F->row_offset[Gjd][J];

        F->matrix[Gjd] = dpd_block_matrix(virtpi[Gd], F->params->coltot[Gjd^GF]);
        buf4_mat_irrep_rd_block(F, Gjd, jd, virtpi[Gd]);

        nrows = F->params->coltot[Gjd^GF];
        ncols = virtpi[Gc];
        nlinks = virtpi[Gd];

        if(nrows && ncols && nlinks)
            C_DGEMM('t','t',nrows, ncols, nlinks, 1.0, F->matrix[Gjd][0], nrows,
                    &(C2->matrix[Gik][ik][cd]), nlinks, 1.0, W1[Gab][0], ncols);

        free_dpd_block(F->matrix[Gjd], virtpi[Gd], F->params->coltot[Gjd^GF]);

        /* -t_ijcd * F_kdab */
        Gkd = Gk ^ Gd; /*changed */
        Gab = Gkd ^ GF;
        Gc = Gij ^ Gd ^ GC;

        cd = C2->col_offset[Gij][Gc];
        kd = F->row_offset[Gkd][K];

        F->matrix[Gkd] = dpd_block_matrix(virtpi[Gd], F->params->coltot[Gkd^GF]);
        buf4_mat_irrep_rd_block(F, Gkd, kd, virtpi[Gd]);

        nrows = F->params->coltot[Gkd^GF];
        ncols = virtpi[Gc];
        nlinks = virtpi[Gd];

        if(nrows && ncols && nlinks)
            C_DGEMM('t', 't', nrows, ncols, nlinks, -1.0, F->matrix[Gkd][0], nrows,
                    &(C2->matrix[Gij][ij][cd]), nlinks, 1.0, W1[Gab][0], ncols);

        free_dpd_block(F->matrix[Gkd], virtpi[Gd], F->params->coltot[Gkd^GF]);

    }

    for(Gl=0; Gl < nirreps; Gl++) {

        /* -t_ilab * E_jklc */
        Gil = Gi ^ Gl; /* changed */
        Gab = Gil ^ GC;
        Gc = Gjk ^ Gl ^ GE;

        lc = E->col_offset[Gjk][Gl];
        il = C2->row_offset[Gil][I];

        nrows = C2->params->coltot[Gil^GC];
        ncols = virtpi[Gc];
        nlinks = occpi[Gl];

        if(nrows && ncols && nlinks)
            C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0, C2->matrix[Gil][il], nrows,
                    &(E->matrix[Gjk][jk][lc]), ncols, 1.0, W1[Gab][0], ncols);

        /* +t_jlab * E_iklc */
        Gjl = Gj ^ Gl; /* changed */
        Gab = Gjl ^ GC;
        Gc = Gik ^ Gl ^ GE;

        lc = E->col_offset[Gik][Gl];
        jl = C2->row_offset[Gjl][J];

        nrows = C2->params->coltot[Gjl^GC];
        ncols = virtpi[Gc];
        nlinks = occpi[Gl];

        if(nrows && ncols && nlinks)
            C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, C2->matrix[Gjl][jl], nrows,
                    &(E->matrix[Gik][ik][lc]), ncols, 1.0, W1[Gab][0], ncols);

        /* +t_klab * E_jilc */
        Gkl = Gk ^ Gl; /* changed! */
        Gab = Gkl ^ GC;
        Gc = Gji ^ Gl ^ GE;

        lc = E->col_offset[Gji][Gl];
        kl = C2->row_offset[Gkl][K];

        nrows = C2->params->coltot[Gkl^GC];
        ncols = virtpi[Gc];
        nlinks = occpi[Gl];

        if(nrows && ncols && nlinks)
            C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, C2->matrix[Gkl][kl], nrows,
                    &(E->matrix[Gji][ji][lc]), ncols, 1.0, W1[Gab][0], ncols);

    }

    for(Gd=0; Gd < nirreps; Gd++) {
        /* +t_kjbd * F_idca */
        Gid = Gi ^ Gd; /* changed */
        Gca = Gid ^ GF;
        Gb = Gjk ^ Gd ^ GC;

        bd = C2->col_offset[Gjk][Gb];
        id = F->row_offset[Gid][I];

        F->matrix[Gid] = dpd_block_matrix(virtpi[Gd], F->params->coltot[Gid^GF]);
        buf4_mat_irrep_rd_block(F, Gid, id, virtpi[Gd]);

        nrows = F->params->coltot[Gid^GF];
        ncols = virtpi[Gb];
        nlinks = virtpi[Gd];

        if(nrows && ncols && nlinks)
            C_DGEMM('t','t',nrows, ncols, nlinks, 1.0, F->matrix[Gid][0], nrows,
                    &(C2->matrix[Gjk][kj][bd]), nlinks, 1.0, W2[Gca][0], ncols);

        free_dpd_block(F->matrix[Gid], virtpi[Gd], F->params->coltot[Gid^GF]);

        /* +t_ikbd * F_jdca */
        Gjd = Gj ^ Gd;
        Gca = Gjd ^ GF ;
        Gb = Gik ^ Gd ^ GC;

        bd = C2->col_offset[Gik][Gb];
        jd = F->row_offset[Gjd][J];

        F->matrix[Gjd] = dpd_block_matrix(virtpi[Gd], F->params->coltot[Gjd^GF]);
        buf4_mat_irrep_rd_block(F, Gjd, jd, virtpi[Gd]);

        nrows = F->params->coltot[Gjd^GF];
        ncols = virtpi[Gb];
        nlinks = virtpi[Gd];

        if(nrows && ncols && nlinks)
            C_DGEMM('t','t',nrows, ncols, nlinks, 1.0, F->matrix[Gjd][0], nrows,
                    &(C2->matrix[Gik][ik][bd]), nlinks, 1.0, W2[Gca][0], ncols);

        free_dpd_block(F->matrix[Gjd], virtpi[Gd], F->params->coltot[Gjd^GF]);

        /* -t_ijbd * F_kdca */
        Gkd = Gk ^ Gd;
        Gca = Gkd ^ GF;
        Gb = Gij ^ Gd ^ GC;

        bd = C2->col_offset[Gij][Gb];
        kd = F->row_offset[Gkd][K];

        F->matrix[Gkd] = dpd_block_matrix(virtpi[Gd], F->params->coltot[Gkd^GF]);
        buf4_mat_irrep_rd_block(F, Gkd, kd, virtpi[Gd]);

        nrows = F->params->coltot[Gkd^GF];
        ncols = virtpi[Gb];
        nlinks = virtpi[Gd];

        if(nrows && ncols && nlinks)
            C_DGEMM('t','t',nrows, ncols, nlinks, -1.0, F->matrix[Gkd][0], nrows,
                    &(C2->matrix[Gij][ij][bd]), nlinks, 1.0, W2[Gca][0], ncols);

        free_dpd_block(F->matrix[Gkd], virtpi[Gd], F->params->coltot[Gkd^GF]);
    }

    for(Gl=0; Gl < nirreps; Gl++) {
        /* -t_ilca * E_jklb */
        Gil = Gi ^ Gl;
        Gca = Gil ^ GC;
        Gb = Gjk ^ Gl ^ GE;

        lb = E->col_offset[Gjk][Gl];
        il = C2->row_offset[Gil][I];

        nrows = C2->params->coltot[Gil^GC];
        ncols = virtpi[Gb];
        nlinks = occpi[Gl];

        if(nrows && ncols && nlinks)
            C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0, C2->matrix[Gil][il], nrows,
                    &(E->matrix[Gjk][jk][lb]), ncols, 1.0, W2[Gca][0], ncols);

        /* +t_jlca * E_iklb */
        Gjl = Gj ^ Gl;
        Gca = Gjl ^ GC;
        Gb = Gik ^ Gl^ GE;

        lb = E->col_offset[Gik][Gl];
        jl = C2->row_offset[Gjl][J];

        nrows = C2->params->coltot[Gjl^GC];
        ncols = virtpi[Gb];
        nlinks = occpi[Gl];

        if(nrows && ncols && nlinks)
            C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, C2->matrix[Gjl][jl], nrows,
                    &(E->matrix[Gik][ik][lb]), ncols, 1.0, W2[Gca][0], ncols);

        /* +t_klca * E_jilb */
        Gkl = Gk ^ Gl;
        Gca = Gkl ^ GC;
        Gb = Gji ^ Gl ^ GE;

        lb = E->col_offset[Gji][Gl];
        kl = C2->row_offset[Gkl][K];

        nrows = C2->params->coltot[Gkl^GC];
        ncols = virtpi[Gb];
        nlinks = occpi[Gl];

        if(nrows && ncols && nlinks)
            C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, C2->matrix[Gkl][kl], nrows,
                    &(E->matrix[Gji][ji][lb]), ncols, 1.0, W2[Gca][0], ncols);
    }

    sort_3d(W2, W1, nirreps, Gijk^GX3, F->params->coltot, F->params->colidx,
            F->params->colorb, F->params->rsym, F->params->ssym, vir_off,
            vir_off, virtpi, vir_off, F->params->colidx, bca, 1);

    for(Gab=0; Gab < nirreps; Gab++) {
        Gc = Gab ^ Gijk ^ GX3; /* changed */
        if(F->params->coltot[Gab] && virtpi[Gc]) {
            ::memset(W2[Gab][0], 0, F->params->coltot[Gab]*virtpi[Gc]*sizeof(double));
        }
    }

    for(Gd=0; Gd < nirreps; Gd++) {
        /* -t_kjad * F_idcb */
        Gid = Gi ^ Gd;
        Gcb = Gid ^ GF;
        Ga = Gkj ^ Gd ^ GC;

        ad = C2->col_offset[Gkj][Ga];
        id = F->row_offset[Gid][I];

        F->matrix[Gid] = dpd_block_matrix(virtpi[Gd], F->params->coltot[Gid^GF]);
        buf4_mat_irrep_rd_block(F, Gid, id, virtpi[Gd]);

        nrows = F->params->coltot[Gid^GF];
        ncols = virtpi[Ga];
        nlinks = virtpi[Gd];

        if(nrows && ncols && nlinks)
            C_DGEMM('t','t',nrows, ncols, nlinks, -1.0, F->matrix[Gid][0], nrows,
                    &(C2->matrix[Gkj][kj][ad]), nlinks, 1.0, W2[Gcb][0], ncols);

        free_dpd_block(F->matrix[Gid], virtpi[Gd], F->params->coltot[Gid^GF]);

        /* -t_ikad * F_jdcb */
        Gjd = Gj ^ Gd;
        Gcb = Gjd ^ GF;
        Ga = Gik ^ Gd ^ GC;

        ad = C2->col_offset[Gik][Ga];
        jd = F->row_offset[Gjd][J];

        F->matrix[Gjd] = dpd_block_matrix(virtpi[Gd], F->params->coltot[Gjd^GF]);
        buf4_mat_irrep_rd_block(F, Gjd, jd, virtpi[Gd]);

        nrows = F->params->coltot[Gjd^GF];
        ncols = virtpi[Ga];
        nlinks = virtpi[Gd];

        if(nrows && ncols && nlinks)
            C_DGEMM('t','t',nrows, ncols, nlinks, -1.0, F->matrix[Gjd][0], nrows,
                    &(C2->matrix[Gik][ik][ad]), nlinks, 1.0, W2[Gcb][0], ncols);

        free_dpd_block(F->matrix[Gjd], virtpi[Gd], F->params->coltot[Gjd^GF]);

        /* +t_ijad * F_kdcb */
        Gkd = Gk ^ Gd;
        Gcb = Gkd ^ GF;
        Ga = Gij ^ Gd ^ GC;

        ad = C2->col_offset[Gij][Ga];
        kd = F->row_offset[Gkd][K];

        F->matrix[Gkd] = dpd_block_matrix(virtpi[Gd], F->params->coltot[Gkd^GF]);
        buf4_mat_irrep_rd_block(F, Gkd, kd, virtpi[Gd]);

        nrows = F->params->coltot[Gkd^GF];
        ncols = virtpi[Ga];
        nlinks = virtpi[Gd];

        if(nrows && ncols && nlinks)
            C_DGEMM('t','t',nrows, ncols, nlinks, 1.0, F->matrix[Gkd][0], nrows,
                    &(C2->matrix[Gij][ij][ad]), nlinks, 1.0, W2[Gcb][0], ncols);

        free_dpd_block(F->matrix[Gkd], virtpi[Gd], F->params->coltot[Gkd^GF]);

    }

    for(Gl=0; Gl < nirreps; Gl++) {
        /* +t_ilcb * E_jkla */
        Gil = Gi ^ Gl;
        Gcb  = Gil ^ GC;
        Ga = Gjk ^ Gl ^ GE;

        la = E->col_offset[Gjk][Gl];
        il = C2->row_offset[Gil][I];

        nrows = C2->params->coltot[Gil^GC];
        ncols = virtpi[Ga];
        nlinks = occpi[Gl];

        if(nrows && ncols && nlinks)
            C_DGEMM('t', 'n', nrows, ncols, nlinks, 1.0, C2->matrix[Gil][il], nrows,
                    &(E->matrix[Gjk][jk][la]), ncols, 1.0, W2[Gcb][0], ncols);

        /* -t_jlcb * E_ikla */
        Gjl = Gj ^ Gl;
        Gcb = Gjl ^ GC;
        Ga = Gik ^ Gl ^ GE;

        la = E->col_offset[Gik][Gl];
        jl = C2->row_offset[Gjl][J];

        nrows = C2->params->coltot[Gjl^GC];
        ncols = virtpi[Ga];
        nlinks = occpi[Gl];

        if(nrows && ncols && nlinks)
            C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0, C2->matrix[Gjl][jl], nrows,
                    &(E->matrix[Gik][ik][la]), ncols, 1.0, W2[Gcb][0], ncols);

        /* -t_klcb * E_jila */
        Gkl = Gk ^ Gl;
        Gcb = Gkl ^ GC;
        Ga = Gji ^ Gl ^ GE;

        la = E->col_offset[Gji][Gl];
        kl = C2->row_offset[Gkl][K];

        nrows = C2->params->coltot[Gkl^GC];
        ncols = virtpi[Ga];
        nlinks = occpi[Gl];

        if(nrows && ncols && nlinks)
            C_DGEMM('t', 'n', nrows, ncols, nlinks, -1.0, C2->matrix[Gkl][kl], nrows,
                    &(E->matrix[Gji][ji][la]), ncols, 1.0, W2[Gcb][0], ncols);

    }

    sort_3d(W2, W1, nirreps, Gijk^GX3, F->params->coltot, F->params->colidx,
            F->params->colorb, F->params->rsym, F->params->ssym, vir_off,
            vir_off, virtpi, vir_off, F->params->colidx, cba, 1);

    for(Gab=0; Gab < nirreps; Gab++) {
        Gc = Gab ^ Gijk ^ GX3; /* assumes totally symmetric! */

        for(ab=0; ab < F->params->coltot[Gab]; ab++) {
            A = F->params->colorb[Gab][ab][0];
            B = F->params->colorb[Gab][ab][1];
            Ga = F->params->rsym[A];
            Gb = F->params->ssym[B];
            a = A - vir_off[Ga];
            b = B - vir_off[Gb];

            for(c=0; c < virtpi[Gc]; c++) {
                C = vir_off[Gc] + c;

                denom = dijk;
                if(fAB->params->rowtot[Ga]) denom -= fAB->matrix[Ga][a][a];
                if(fAB->params->rowtot[Gb]) denom -= fAB->matrix[Gb][b][b];
                if(fAB->params->rowtot[Gc]) denom -= fAB->matrix[Gc][c][c];

                W1[Gab][ab][c] /= (omega + denom);

            } /* c */
        } /* ab */
    } /* Gab */

    for(Gab=0; Gab < nirreps; Gab++) {
        Gc = Gab ^ Gijk ^ GX3; /* changed */
        free_dpd_block(W2[Gab], F->params->coltot[Gab], virtpi[Gc]);
    }
    free(W2);

    file2_mat_close(fIJ);
    file2_mat_close(fAB);

    for(h=0; h < nirreps; h++) {
        buf4_mat_irrep_close(C2, h);
        buf4_mat_irrep_close(E, h);
    }
}

}
