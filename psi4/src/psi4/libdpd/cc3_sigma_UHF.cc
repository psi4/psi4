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
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libqt/qt.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/psifiles.h"
#include "psi4/libparallel/ParallelPrinter.h"
namespace psi {

/*
  This function computes contributions to singles and doubles of
  matrix elements of triples:
    SIA   <-- <S|(Dints)           <T|(Wmbij,Wabei) CMNEF |0> |T> / (w-wt)
    SIjAb <-- <D|(FME,WAmEf,WMnIe) <T|(Wmbij,Wabei) CMNEF |0> |T> / (w-wt)
  These are used to make X3 quantity in T3_RHF:
    CIjAb, WAbEi, WMbIj, fIJ, fAB, omega
*/

void DPD::cc3_sigma_UHF_AAA(dpdbuf4 *CMNEF, dpdbuf4 *WABEI, dpdbuf4 *WMBIJ,
                            int do_singles, dpdbuf4 *Dints_anti, dpdfile2 *SIA, int do_doubles, dpdfile2 *FME,
                            dpdbuf4 *WMAFE, dpdbuf4 *WMNIE, dpdbuf4 *SIJAB, int *aoccpi, int *aocc_off,
                            int *avirtpi, int *avir_off, double omega, std::string out)
{
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
            std::shared_ptr<OutFile>(new OutFile(out)));
   int h, nirreps;
    int *occ_off, *occpi, *vir_off, *virtpi;
    int Gi, Gj, Gk, Gijk, Ga, Gb, Gc;
    int i, j, k, I, J, K;
    int a, b, c, A, B, C;
    int Gij, ij, Gab, ab, Gjk, jk;
    double ***W1, **Z;
    dpdfile2 fIJ, fAB;
    int nrows, ncols, nlinks;
    int Gd, d, cd, dc, Gid, id, DD;
    int Gm, m, Gmi, mi, im, mc, M;
    int GS, GH, GU, GC, GX3;

    GC = CMNEF->file.my_irrep;

    GU = WMBIJ->file.my_irrep;

    if (do_singles)
        GH = Dints_anti->file.my_irrep;
    else if (do_doubles)
        GH = WMAFE->file.my_irrep;

    GS = SIJAB->file.my_irrep;

    GX3 = GU^GC;

    nirreps = CMNEF->params->nirreps;
    occpi = aoccpi;
    occ_off = aocc_off;
    virtpi = avirtpi;
    vir_off = avir_off;

    file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");

    if (do_singles) {
        file2_mat_init(SIA);
        file2_mat_rd(SIA);

        for(h=0; h < nirreps; h++) {
            buf4_mat_irrep_init(Dints_anti, h);
            buf4_mat_irrep_rd(Dints_anti, h);
        }
    }

    if (do_doubles) {
        file2_mat_init(FME);
        file2_mat_rd(FME);

        for(h=0; h < nirreps; h++) {
            buf4_mat_irrep_init(SIJAB, h);
            buf4_mat_irrep_rd(SIJAB, h);
        }
        for(h=0; h < nirreps; h++) {
            buf4_mat_irrep_init(WMNIE, h);
            buf4_mat_irrep_rd(WMNIE, h);
        }
    }

    /* target T3 amplitudes go in here */
    W1 = (double ***) malloc(nirreps * sizeof(double **));

    for(Gi=0; Gi < nirreps; Gi++) {
        for(Gj=0; Gj < nirreps; Gj++) {
            Gij = Gi ^ Gj;
            for(Gk=0; Gk < nirreps; Gk++) {
                Gijk = Gi ^ Gj ^ Gk;
                Gjk = Gj ^ Gk;

                for(Gab=0; Gab < nirreps; Gab++) {
                    /* changed */
                    Gc = Gab ^ Gijk ^ GX3;
                    W1[Gab] = dpd_block_matrix(WABEI->params->coltot[Gab], virtpi[Gc]);
                }

                for(i=0; i < occpi[Gi]; i++) {
                    I = occ_off[Gi] + i;
                    for(j=0; j < occpi[Gj]; j++) {
                        J = occ_off[Gj] + j;
                        for(k=0; k < occpi[Gk]; k++) {
                            K = occ_off[Gk] + k;

                            T3_AAA(W1, nirreps, I, Gi, J, Gj, K, Gk, CMNEF, WABEI, WMBIJ, &fIJ, &fAB,
                                   occpi, occ_off, virtpi, vir_off, omega);

                            if (do_singles) {
                                /* t_KC <-- 1/4 t_IJKABC <IJ||AB> */

                                Gc = Gk ^ GS;  /* changed */
                                Gab = Gij ^ GH ;  /* changed */

                                ij = Dints_anti->params->rowidx[I][J];

                                nrows = Dints_anti->params->coltot[Gij^GH]; /* changed */
                                ncols = virtpi[Gc];

                                if(nrows && ncols)
                                    C_DGEMV('t', nrows, ncols, 0.25, W1[Gab][0], ncols, Dints_anti->matrix[Gij][ij], 1,
                                            1.0, SIA->matrix[Gk][k], 1);
                            }

                            if (do_doubles) {
                                /* t_IJAB <-- t_IJKABC F_KC */
                                Gc = Gk ^ GH;
                                Gab = Gij ^ GS;

                                nrows = SIJAB->params->coltot[Gij^GS];
                                ncols = virtpi[Gc];
                                ij = SIJAB->params->rowidx[I][J];

                                if(nrows && ncols)
                                    C_DGEMV('n', nrows, ncols, 1.0, W1[Gab][0], ncols, FME->matrix[Gk][k], 1,
                                            1.0, SIJAB->matrix[Gij][ij], 1);

                                /* t_JKDC <-- +1/2 t_IJKABC W_IDAB */
                                /* t_JKCD <-- -1/2 t_IJKABC W_IDAB */
                                jk = SIJAB->params->rowidx[J][K];
                                for(Gd=0; Gd < nirreps; Gd++) {
                                    Gid = Gi ^ Gd;
                                    Gab = Gid ^ GH;
                                    Gc = Gab ^ Gijk ^ GX3;

                                    id = WMAFE->row_offset[Gid][I];

                                    Z = dpd_block_matrix(virtpi[Gc],virtpi[Gd]);
                                    WMAFE->matrix[Gid] = dpd_block_matrix(virtpi[Gd], WMAFE->params->coltot[Gid^GH]);
                                    buf4_mat_irrep_rd_block(WMAFE, Gid, id, virtpi[Gd]);

                                    nrows = virtpi[Gc];
                                    ncols = virtpi[Gd];
                                    nlinks = WMAFE->params->coltot[Gid^GH];

                                    if(nrows && ncols && nlinks)
                                        C_DGEMM('t', 't', nrows, ncols, nlinks, 0.5, W1[Gab][0], nrows,
                                                WMAFE->matrix[Gid][0], nlinks, 0.0, Z[0], ncols);

                                    for(c=0; c < virtpi[Gc]; c++) {
                                        C = vir_off[Gc] + c;
                                        for(d=0; d < virtpi[Gd]; d++) {
                                            DD = vir_off[Gd] + d;
                                            cd = SIJAB->params->colidx[C][DD];
                                            dc = SIJAB->params->colidx[DD][C];
                                            SIJAB->matrix[Gjk][jk][dc] += Z[c][d];
                                            SIJAB->matrix[Gjk][jk][cd] += -Z[c][d];
                                        }
                                    }
                                    free_dpd_block(WMAFE->matrix[Gid], virtpi[Gd], WMAFE->params->coltot[Gid^GH]);
                                    free_dpd_block(Z,virtpi[Gc],virtpi[Gd]);
                                }


                                /* S_MIAB <-- +1/2 t_IJKABC W_JKMC */
                                /* S_IMAB <-- -1/2 t_IJKABC W_JKMC */
                                jk = WMNIE->params->rowidx[J][K];
                                for(Gm=0; Gm < nirreps; Gm++) {
                                    Gmi = Gm ^ Gi;
                                    Gab = Gmi ^ GS;
                                    Gc = Gab ^ Gijk ^ GX3;

                                    mc = WMNIE->col_offset[Gjk][Gm];

                                    nrows = WABEI->params->coltot[Gab];
                                    ncols = occpi[Gm];
                                    nlinks = virtpi[Gc];

                                    Z = dpd_block_matrix(nrows, ncols);

                                    if(nrows && ncols && nlinks)
                                        C_DGEMM('n', 't', nrows, ncols, nlinks, 0.5, W1[Gab][0], nlinks,
                                                &(WMNIE->matrix[Gjk][jk][mc]), nlinks, 0.0, Z[0], ncols);

                                    for(m=0; m < ncols; m++) {
                                        M = occ_off[Gm] + m;
                                        mi = SIJAB->params->rowidx[M][I];
                                        im = SIJAB->params->rowidx[I][M];
                                        for(ab=0; ab < nrows; ab++) {
                                            SIJAB->matrix[Gmi][mi][ab] += Z[ab][m];
                                            SIJAB->matrix[Gmi][im][ab] -= Z[ab][m];
                                        }
                                    }

                                    free_dpd_block(Z, nrows, ncols);
                                }
                            } /* end do_doubles */

                        } /* k */
                    } /* j */
                } /* i */

                for(Gab=0; Gab < nirreps; Gab++) {
                    /* This will need to change for non-totally-symmetric cases */
                    Gc = Gab ^ Gijk;
                    free_dpd_block(W1[Gab], WABEI->params->coltot[Gab], virtpi[Gc]);
                }

            } /* Gk */
        } /* Gj */
    } /* Gi */

    free(W1);

    file2_close(&fIJ);
    file2_close(&fAB);

    if (do_singles) {

        file2_mat_wrt(SIA);
        file2_mat_close(SIA);

        for(h=0; h < nirreps; h++)
            buf4_mat_irrep_close(Dints_anti, h);
    }

    if (do_doubles) {

        file2_mat_close(FME);

        for(h=0; h < nirreps; h++)
            buf4_mat_irrep_close(WMNIE, h);

        for(h=0; h < nirreps; h++) {
            buf4_mat_irrep_wrt(SIJAB, h);
            buf4_mat_irrep_close(SIJAB, h);
        }
    }
}

void DPD::cc3_sigma_UHF_BBB(dpdbuf4 *Cmnef, dpdbuf4 *Wabei, dpdbuf4 *Wmbij,
                            int do_singles, dpdbuf4 *Dijab_anti, dpdfile2 *Sia, int do_doubles,
                            dpdfile2 *Fme, dpdbuf4 *Wmafe, dpdbuf4 *Wmnie, dpdbuf4 *Sijab,
                            int *boccpi, int *bocc_off, int *bvirtpi, int *bvir_off,
                            double omega, std::string out)
{
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
            std::shared_ptr<OutFile>(new OutFile(out)));
   int h, nirreps;
    int *occ_off, *occpi, *vir_off, *virtpi;
    int Gi, Gj, Gk, Gijk, Ga, Gb, Gc;
    int i, j, k, I, J, K;
    int a, b, c, A, B, C;
    int Gij, ij, Gab, ab, Gjk, jk;
    double ***W1, **Z;
    dpdfile2 fIJ, fAB;
    int nrows, ncols, nlinks;
    int Gd, d, cd, dc, Gid, id, DD;
    int Gm, m, Gmi, mi, im, mc, M;
    int GS, GH, GU, GC, GX3;

    GC = Cmnef->file.my_irrep;
    GU = Wmbij->file.my_irrep;

    if (do_singles)
        GH = Dijab_anti->file.my_irrep;
    else if (do_doubles)
        GH = Wmafe->file.my_irrep;

    GS = Sijab->file.my_irrep;
    GX3 = GU^GC;

    nirreps = Cmnef->params->nirreps;
    occpi = boccpi;
    occ_off = bocc_off;
    virtpi = bvirtpi;
    vir_off = bvir_off;

    file2_init(&fIJ, PSIF_CC_OEI, 0, 2, 2, "fij");
    file2_init(&fAB, PSIF_CC_OEI, 0, 3, 3, "fab");

    if (do_singles) {

        file2_mat_init(Sia);
        file2_mat_rd(Sia);

        for(h=0; h < nirreps; h++) {
            buf4_mat_irrep_init(Dijab_anti, h);
            buf4_mat_irrep_rd(Dijab_anti, h);
        }
    }

    if (do_doubles) {

        file2_mat_init(Fme);
        file2_mat_rd(Fme);

        for(h=0; h < nirreps; h++) {
            buf4_mat_irrep_init(Sijab, h);
            buf4_mat_irrep_rd(Sijab, h);
        }

        for(h=0; h < nirreps; h++) {
            buf4_mat_irrep_init(Wmnie, h);
            buf4_mat_irrep_rd(Wmnie, h);
        }
    }

    /* target T3 amplitudes go in here */
    W1 = (double ***) malloc(nirreps * sizeof(double **));

    for(Gi=0; Gi < nirreps; Gi++) {
        for(Gj=0; Gj < nirreps; Gj++) {
            Gij = Gi ^ Gj;
            for(Gk=0; Gk < nirreps; Gk++) {
                Gijk = Gi ^ Gj ^ Gk;
                Gjk = Gj ^ Gk;

                for(Gab=0; Gab < nirreps; Gab++) {
                    Gc = Gab ^ Gijk ^ GX3;
                    W1[Gab] = dpd_block_matrix(Wabei->params->coltot[Gab], virtpi[Gc]);
                }

                for(i=0; i < occpi[Gi]; i++) {
                    I = occ_off[Gi] + i;
                    for(j=0; j < occpi[Gj]; j++) {
                        J = occ_off[Gj] + j;
                        for(k=0; k < occpi[Gk]; k++) {
                            K = occ_off[Gk] + k;

                            T3_AAA(W1, nirreps, I, Gi, J, Gj, K, Gk, Cmnef, Wabei, Wmbij, &fIJ, &fAB,
                                   occpi, occ_off, virtpi, vir_off, omega);

                            if (do_singles) {

                                /* S_kc <-- 1/4 t_ijkabc <ij||ab> */
                                Gc = Gk ^ GS;
                                Gab = Gij ^ GH;

                                ij = Dijab_anti->params->rowidx[I][J];

                                nrows = Dijab_anti->params->coltot[Gij^GH];
                                ncols = virtpi[Gc];

                                if(nrows && ncols)
                                    C_DGEMV('t', nrows, ncols, 0.25, W1[Gab][0], ncols, Dijab_anti->matrix[Gij][ij], 1,
                                            1.0, Sia->matrix[Gk][k], 1);
                            }

                            if (do_doubles) {

                                /* S_ijab <-- t_ijkabc F_kc */
                                Gc = Gk ^ GH;
                                Gab = Gij ^ GS;

                                nrows = Sijab->params->coltot[Gij^GS];
                                ncols = virtpi[Gc];
                                ij = Sijab->params->rowidx[I][J];

                                if(nrows && ncols)
                                    C_DGEMV('n', nrows, ncols, 1.0, W1[Gab][0], ncols, Fme->matrix[Gk][k], 1,
                                            1.0, Sijab->matrix[Gij][ij], 1);

                                /* S_jkdc <-- 1/2 t_ijkabc W_idab */
                                /* S_jkcd <-- -1/2 t_ijkabc W_idab */
                                jk = Sijab->params->rowidx[J][K];
                                for(Gd=0; Gd < nirreps; Gd++) {

                                    Gid = Gi ^ Gd;
                                    Gab = Gid ^ GH;
                                    Gc = Gab ^ Gijk ^ GX3;

                                    id = Wmafe->row_offset[Gid][I];

                                    Z = dpd_block_matrix(virtpi[Gc],virtpi[Gd]);
                                    Wmafe->matrix[Gid] = dpd_block_matrix(virtpi[Gd], Wmafe->params->coltot[Gid^GH]);
                                    buf4_mat_irrep_rd_block(Wmafe, Gid, id, virtpi[Gd]);

                                    nrows = virtpi[Gc];
                                    ncols = virtpi[Gd];
                                    nlinks = Wmafe->params->coltot[Gid^GH];

                                    if(nrows && ncols && nlinks)
                                        C_DGEMM('t', 't', nrows, ncols, nlinks, 0.5, W1[Gab][0], nrows,
                                                Wmafe->matrix[Gid][0], nlinks, 0.0, Z[0], ncols);

                                    for(c=0; c < virtpi[Gc]; c++) {
                                        C = vir_off[Gc] + c;
                                        for(d=0; d < virtpi[Gd]; d++) {
                                            DD = vir_off[Gd] + d;
                                            cd = Sijab->params->colidx[C][DD];
                                            dc = Sijab->params->colidx[DD][C];
                                            Sijab->matrix[Gjk][jk][dc] += Z[c][d];
                                            Sijab->matrix[Gjk][jk][cd] += -Z[c][d];
                                        }
                                    }
                                    free_dpd_block(Wmafe->matrix[Gid], virtpi[Gd], Wmafe->params->coltot[Gid^GH]);
                                    free_dpd_block(Z,virtpi[Gc],virtpi[Gd]);
                                }

                                /* S_miab <-- +1/2 t_ijkabc W_jkmc */
                                /* S_imab <-- -1/2 t_ijkabc W_jkmc */
                                jk = Wmnie->params->rowidx[J][K];
                                for(Gm=0; Gm < nirreps; Gm++) {
                                    Gmi = Gm ^ Gi;
                                    Gab = Gmi ^ GS;
                                    Gc = Gab ^ Gijk ^ GX3;

                                    mc = Wmnie->col_offset[Gjk][Gm];

                                    nrows = Wabei->params->coltot[Gab];
                                    ncols = occpi[Gm];
                                    nlinks = virtpi[Gc];

                                    Z = dpd_block_matrix(nrows, ncols);

                                    if(nrows && ncols && nlinks)
                                        C_DGEMM('n', 't', nrows, ncols, nlinks, 0.5, W1[Gab][0], nlinks,
                                                &(Wmnie->matrix[Gjk][jk][mc]), nlinks, 0.0, Z[0], ncols);

                                    for(m=0; m < ncols; m++) {
                                        M = occ_off[Gm] + m;
                                        mi = Sijab->params->rowidx[M][I];
                                        im = Sijab->params->rowidx[I][M];
                                        for(ab=0; ab < nrows; ab++) {
                                            Sijab->matrix[Gmi][mi][ab] += Z[ab][m];
                                            Sijab->matrix[Gmi][im][ab] -= Z[ab][m];
                                        }
                                    }

                                    free_dpd_block(Z, nrows, ncols);
                                }
                            } /* end do_doubles */

                        } /* k */
                    } /* j */
                } /* i */

                for(Gab=0; Gab < nirreps; Gab++) {
                    /* This will need to change for non-totally-symmetric cases */
                    Gc = Gab ^ Gijk;
                    free_dpd_block(W1[Gab], Wabei->params->coltot[Gab], virtpi[Gc]);
                }

            } /* Gk */
        } /* Gj */
    } /* Gi */

    free(W1);

    file2_close(&fIJ);
    file2_close(&fAB);

    if (do_singles) {

        file2_mat_wrt(Sia);
        file2_mat_close(Sia);

        for(h=0; h < nirreps; h++)
            buf4_mat_irrep_close(Dijab_anti, h);
    }

    if (do_doubles) {

        file2_mat_close(Fme);

        for(h=0; h < nirreps; h++)
            buf4_mat_irrep_close(Wmnie, h);

        for(h=0; h < nirreps; h++) {
            buf4_mat_irrep_wrt(Sijab, h);
            buf4_mat_irrep_close(Sijab, h);
        }
    }
}

void DPD::cc3_sigma_UHF_AAB(dpdbuf4 *C2AA, dpdbuf4 *C2AB, dpdbuf4 *C2BA,
                            dpdbuf4 *FAA, dpdbuf4 *FAB, dpdbuf4 *FBA,
                            dpdbuf4 *EAA, dpdbuf4 *EAB, dpdbuf4 *EBA,
                            int do_singles, dpdbuf4 *DAA, dpdbuf4 *DAB, dpdfile2 *SIA, dpdfile2 *Sia,
                            int do_doubles, dpdfile2 *FME, dpdfile2 *Fme,
                            dpdbuf4 *WMAFE, dpdbuf4 *WMaFe, dpdbuf4 *WmAfE,
                            dpdbuf4 *WMNIE, dpdbuf4 *WMnIe, dpdbuf4 *WmNiE,
                            dpdbuf4 *SIJAB, dpdbuf4 *SIjAb, int *aoccpi, int *aocc_off, int *boccpi,
                            int *bocc_off, int *avirtpi, int *avir_off, int *bvirtpi, int *bvir_off,
                            double omega, std::string out)
{
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
            std::shared_ptr<OutFile>(new OutFile(out)));
   int h, nirreps;
    int Gi, Gj, Gk, Gijk;
    int Ga, Gb, Gc, Gab;
    int Gij, ij, ji, Gjk, jk, Gbc, bc, Gcb, cb;
    int Gd, Gkd, kd, d, DD, ad, da, dc, Gid, id, Gac, ac, bd;
    int i, j, k, I, J, K;
    int a, b, c, A, B, C;
    int ab;
    double ***W1, ***W2, ***W3;
    dpdfile2 fIJ, fAB, fij, fab;
    int nrows, ncols, nlinks;
    int **W_offset, offset;
    double **Z;
    int Gm, m, Gmi, mi, im, mc, M;
    int Gmk, mk, ma, kj, Gim;
    int GS, GH, GU, GC, GX3;

    GC = C2AA->file.my_irrep;
    GU = EAA->file.my_irrep;

    if (do_singles)
        GH = DAA->file.my_irrep;
    else if (do_doubles)
        GH = WMAFE->file.my_irrep;

    GS = SIJAB->file.my_irrep;
    GX3 = GU^GC;

    nirreps = C2AA->params->nirreps;

    W_offset = init_int_matrix(nirreps, nirreps);
    for(Gab=0; Gab < nirreps; Gab++) {
        for(Ga=0,offset=0; Ga < nirreps; Ga++) {
            Gb = Ga ^ Gab;
            W_offset[Gab][Ga] = offset;
            offset += avirtpi[Ga] * avirtpi[Gb];
        }
    }

    if (do_singles) {
        file2_mat_init(SIA);
        file2_mat_rd(SIA);
        file2_mat_init(Sia);
        file2_mat_rd(Sia);

        for(h=0; h < nirreps; h++) {
            buf4_mat_irrep_init(DAA, h);
            buf4_mat_irrep_rd(DAA, h);
            buf4_mat_irrep_init(DAB, h);
            buf4_mat_irrep_rd(DAB, h);
        }
    }

    if (do_doubles) {
        file2_mat_init(FME);
        file2_mat_rd(FME);
        file2_mat_init(Fme);
        file2_mat_rd(Fme);

        for(h=0; h < nirreps; h++) {
            buf4_mat_irrep_init(WMnIe, h);
            buf4_mat_irrep_rd(WMnIe, h);

            buf4_mat_irrep_init(WMNIE, h);
            buf4_mat_irrep_rd(WMNIE, h);

            buf4_mat_irrep_init(WmNiE, h);
            buf4_mat_irrep_rd(WmNiE, h);
        }
        for(h=0; h < nirreps; h++) {
            buf4_mat_irrep_init(SIJAB, h);
            buf4_mat_irrep_rd(SIJAB, h);
        }
        for(h=0; h < nirreps; h++) {
            buf4_mat_irrep_init(SIjAb, h);
            buf4_mat_irrep_rd(SIjAb, h);
        }
    }

    file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
    file2_init(&fij, PSIF_CC_OEI, 0, 2, 2, "fij");
    file2_init(&fab, PSIF_CC_OEI, 0, 3, 3, "fab");

    /* target T3 amplitudes go in here */
    W1 = (double ***) malloc(nirreps * sizeof(double **));
    W2 = (double ***) malloc(nirreps * sizeof(double **));
    W3 = (double ***) malloc(nirreps * sizeof(double **));

    for(Gi=0; Gi < nirreps; Gi++) {
        for(Gj=0; Gj < nirreps; Gj++) {
            Gij = Gi ^ Gj;
            for(Gk=0; Gk < nirreps; Gk++) {
                Gijk = Gi ^ Gj ^ Gk;
                Gjk = Gj ^ Gk;

                for(Gab=0; Gab < nirreps; Gab++) {
                    Gc = Gab ^ Gijk ^ GX3;
                    W1[Gab] = dpd_block_matrix(FAA->params->coltot[Gab], bvirtpi[Gc]);
                }
                for(Ga=0; Ga < nirreps; Ga++) {
                    Gcb = Ga ^ Gijk ^ GX3;
                    W2[Ga] = dpd_block_matrix(avirtpi[Ga], WmAfE->params->coltot[Gcb]);  /* alpha-beta-alpha */
                }
                for(Gb=0; Gb < nirreps; Gb++) {
                    Gac = Gb ^ Gijk ^ GX3;
                    W3[Gb] = dpd_block_matrix(avirtpi[Gb], WMaFe->params->coltot[Gac]);  /* alpha-alpha-beta */
                }

                for(i=0; i < aoccpi[Gi]; i++) {
                    I = aocc_off[Gi] + i;
                    for(j=0; j < aoccpi[Gj]; j++) {
                        J = aocc_off[Gj] + j;
                        for(k=0; k < boccpi[Gk]; k++) {
                            K = bocc_off[Gk] + k;

                            T3_AAB(W1, nirreps, I, Gi, J, Gj, K, Gk, C2AA, C2AB, C2BA,
                                   FAA, FAB, FBA, EAA, EAB, EBA, &fIJ, &fij, &fAB, &fab,
                                   aoccpi, aocc_off, boccpi, bocc_off, avirtpi, avir_off, bvirtpi, bvir_off, omega);

                            if (do_singles) {

                                /* S_kc <-- 1/4 t_IJkABc <IJ||AB> */

                                Gc = Gk ^ GS;
                                Gab = Gij ^ GH;

                                ij = DAA->params->rowidx[I][J];

                                nrows = DAA->params->coltot[Gij^GH];
                                ncols = bvirtpi[Gc];

                                if(nrows && ncols)
                                    C_DGEMV('t', nrows, ncols, 0.25, W1[Gab][0], ncols, DAA->matrix[Gij][ij], 1,
                                            1.0, Sia->matrix[Gk][k], 1);

                                /* S_IA <-- t_IJkABc <Jk|Bc> */

                                Ga = Gi ^ GS;
                                Gbc = Gjk ^ GH;

                                jk = DAB->params->rowidx[J][K];

                                for(Gab=0; Gab < nirreps; Gab++) {
                                    Gb = Ga ^ Gab;
                                    Gc = Gb ^ Gbc;

                                    ab = W_offset[Gab][Ga];
                                    bc = DAB->col_offset[Gjk][Gb];

                                    nrows = avirtpi[Ga];
                                    ncols = avirtpi[Gb] * bvirtpi[Gc];

                                    if(nrows && ncols)
                                        C_DGEMV('n', nrows, ncols, 1.0, W1[Gab][ab], ncols, &(DAB->matrix[Gjk][jk][bc]), 1,
                                                1.0, SIA->matrix[Gi][i], 1);
                                }
                            }

                            if (do_doubles) {
                                /* S_IJAB <-- t_IJkABc F_kc */
                                Gc = Gk ^ GH;
                                Gab = Gij ^ GS;

                                nrows = SIJAB->params->coltot[Gij^GS];
                                ncols = bvirtpi[Gc];
                                ij = SIJAB->params->rowidx[I][J];

                                if(nrows && ncols)
                                    C_DGEMV('n', nrows, ncols, 1.0, W1[Gab][0], ncols, Fme->matrix[Gk][k], 1,
                                            1.0, SIJAB->matrix[Gij][ij], 1);

                                /* S_JkBc <-- t_IJkABc F_IA */
                                Ga = Gi ^ GH;
                                Gbc = Gjk ^ GS;

                                jk = C2AB->params->rowidx[J][K];

                                for(Gab=0; Gab < nirreps; Gab++) {
                                    Gb = Ga ^ Gab;
                                    Gc = Gb ^ Gbc;

                                    ab = W_offset[Gab][Ga];
                                    bc = SIjAb->col_offset[Gjk][Gb];

                                    nrows = avirtpi[Ga];
                                    ncols = avirtpi[Gb] * bvirtpi[Gc];

                                    if(nrows && ncols)
                                        C_DGEMV('t', nrows, ncols, 1.0, W1[Gab][ab], ncols, FME->matrix[Gi][i], 1,
                                                1.0, &(SIjAb->matrix[Gjk][jk][bc]), 1);
                                }

                                /* S_JIDA <-- t_IJkABc W_kDcB */
                                /* sort W1(AB,c) to W2(A,cB) */
                                for(Gab=0; Gab < nirreps; Gab++) {
                                    Gc = Gab ^ Gijk ^ GX3;
                                    for(ab=0; ab < FAA->params->coltot[Gab]; ab++) {
                                        A = FAA->params->colorb[Gab][ab][0];
                                        B = FAA->params->colorb[Gab][ab][1];
                                        Ga = FAA->params->rsym[A];
                                        a = A - avir_off[Ga];
                                        for(c=0; c < bvirtpi[Gc]; c++) {
                                            C = bvir_off[Gc] + c;
                                            cb = WmAfE->params->colidx[C][B];
                                            W2[Ga][a][cb] = W1[Gab][ab][c];
                                        }
                                    }
                                }

                                ji = SIJAB->params->rowidx[J][I];

                                for(Gd=0; Gd < nirreps; Gd++) {
                                    Gkd = Gk ^ Gd;
                                    Gcb = Gkd ^ GH;
                                    Ga = Gd ^ Gij ^ GS;

                                    kd = WmAfE->row_offset[Gkd][K];
                                    WmAfE->matrix[Gkd] = dpd_block_matrix(avirtpi[Gd], WmAfE->params->coltot[Gkd^GH]);
                                    buf4_mat_irrep_rd_block(WmAfE, Gkd, kd, avirtpi[Gd]);
                                    Z = dpd_block_matrix(avirtpi[Ga], avirtpi[Gd]);

                                    nrows = avirtpi[Ga];
                                    ncols = avirtpi[Gd];
                                    nlinks = WmAfE->params->coltot[Gkd^GH];

                                    if(nrows && ncols && nlinks)
                                        C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, W2[Ga][0], nlinks,
                                                WmAfE->matrix[Gkd][0], nlinks, 0.0, Z[0], ncols);

                                    for(a=0; a < avirtpi[Ga]; a++) {
                                        A = avir_off[Ga] + a;
                                        for(d=0; d < avirtpi[Gd]; d++) {
                                            DD = avir_off[Gd] + d;
                                            ad = SIJAB->params->colidx[A][DD];
                                            da = SIJAB->params->colidx[DD][A];
                                            SIJAB->matrix[Gij][ji][ad] += -Z[a][d];
                                            SIJAB->matrix[Gij][ji][da] += Z[a][d];
                                        }
                                    }

                                    free_dpd_block(WmAfE->matrix[Gkd], avirtpi[Gd], WmAfE->params->coltot[Gkd^GH]);
                                    free_dpd_block(Z, avirtpi[Ga], avirtpi[Gd]);
                                }

                                /* S_JkDc <-- 1/2 t_IJkABc W_IDAB */

                                jk = SIjAb->params->rowidx[J][K];

                                for(Gd=0; Gd < nirreps; Gd++) {
                                    Gid = Gi ^ Gd;
                                    Gab = Gid ^ GH;
                                    Gc = Gab ^ Gijk ^ GX3;

                                    id = WMAFE->row_offset[Gid][I];
                                    WMAFE->matrix[Gid] = dpd_block_matrix(avirtpi[Gd], WMAFE->params->coltot[Gid^GH]);
                                    buf4_mat_irrep_rd_block(WMAFE, Gid, id, avirtpi[Gd]);
                                    Z = dpd_block_matrix(bvirtpi[Gc], avirtpi[Gd]);

                                    nrows = bvirtpi[Gc];
                                    ncols = avirtpi[Gd];
                                    nlinks = WMAFE->params->coltot[Gid^GH];

                                    if(nrows && ncols && nlinks)
                                        C_DGEMM('t', 't', nrows, ncols, nlinks, 0.5, W1[Gab][0], nrows,
                                                WMAFE->matrix[Gid][0], nlinks, 0.0, Z[0], ncols);

                                    for(c=0; c < bvirtpi[Gc]; c++) {
                                        C = bvir_off[Gc] + c;
                                        for(d=0; d < avirtpi[Gd]; d++) {
                                            DD = avir_off[Gd] + d;
                                            dc = SIjAb->params->colidx[DD][C];
                                            SIjAb->matrix[Gjk][jk][dc] += Z[c][d];
                                        }
                                    }

                                    free_dpd_block(Z, bvirtpi[Gc], avirtpi[Gd]);
                                    free_dpd_block(WMAFE->matrix[Gid], avirtpi[Gd], WMAFE->params->coltot[Gid^GH]);
                                }

                                /* t_JkBd <-- t_IJkABc W_IdAc */
                                /* sort W1(AB,c) to W3(B,Ac) */
                                for(Gab=0; Gab < nirreps; Gab++) {
                                    Gc = Gab ^ Gijk ^ GX3;
                                    for(ab=0; ab < FAA->params->coltot[Gab]; ab++) {
                                        A = FAA->params->colorb[Gab][ab][0];
                                        B = FAA->params->colorb[Gab][ab][1];
                                        Gb = FAA->params->ssym[B];
                                        b = B - avir_off[Gb];
                                        for(c=0; c < bvirtpi[Gc]; c++) {
                                            C = bvir_off[Gc] + c;
                                            ac = WMaFe->params->colidx[A][C];
                                            W3[Gb][b][ac] = W1[Gab][ab][c];
                                        }
                                    }
                                }

                                jk = SIjAb->params->rowidx[J][K];

                                for(Gd=0; Gd < nirreps; Gd++) {
                                    Gid = Gi ^ Gd;
                                    Gac = Gid ^ GH;
                                    Gb = Gac ^ Gijk ^ GX3;

                                    id = WMaFe->row_offset[Gid][I];
                                    WMaFe->matrix[Gid] = dpd_block_matrix(bvirtpi[Gd], WMaFe->params->coltot[Gid^GH]);
                                    buf4_mat_irrep_rd_block(WMaFe, Gid, id, bvirtpi[Gd]);
                                    Z = dpd_block_matrix(avirtpi[Gb], bvirtpi[Gd]);

                                    nrows = avirtpi[Gb];
                                    ncols = bvirtpi[Gd];
                                    nlinks = WMaFe->params->coltot[Gid^GH];

                                    if(nrows && ncols && nlinks)
                                        C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, W3[Gb][0], nlinks,
                                                WMaFe->matrix[Gid][0], nlinks, 0.0, Z[0], ncols);

                                    for(b=0; b < avirtpi[Gb]; b++) {
                                        B = avir_off[Gb] + b;
                                        for(d=0; d < bvirtpi[Gd]; d++) {
                                            DD = bvir_off[Gd] + d;
                                            bd = SIjAb->params->colidx[B][DD];
                                            SIjAb->matrix[Gjk][jk][bd] += Z[b][d];
                                        }
                                    }

                                    free_dpd_block(WMaFe->matrix[Gid], bvirtpi[Gd], WMaFe->params->coltot[Gid^GH]);
                                    free_dpd_block(Z, avirtpi[Gb], bvirtpi[Gd]);
                                }

                                /* S_MIAB <--- +t_IJkABc W_JkMc */
                                /* S_IMAB <--- -t_IJkABc W_JkMc */

                                jk = WMnIe->params->rowidx[J][K];

                                for(Gm=0; Gm < nirreps; Gm++) {
                                    Gmi = Gm ^ Gi;
                                    Gab = Gmi ^ GS;
                                    Gc = Gab ^ Gijk ^ GX3;

                                    mc = WMnIe->col_offset[Gjk][Gm];

                                    nrows = FAA->params->coltot[Gab];
                                    ncols = aoccpi[Gm];
                                    nlinks = bvirtpi[Gc];

                                    Z = dpd_block_matrix(nrows, ncols);

                                    if(nrows && ncols && nlinks)
                                        C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, W1[Gab][0], nlinks,
                                                &(WMnIe->matrix[Gjk][jk][mc]), nlinks, 0.0, Z[0], ncols);

                                    for(m=0; m < ncols; m++) {
                                        M = aocc_off[Gm] + m;
                                        mi = SIJAB->params->rowidx[M][I];
                                        im = SIJAB->params->rowidx[I][M];
                                        for(ab=0; ab < nrows; ab++) {
                                            SIJAB->matrix[Gmi][mi][ab] += Z[ab][m];
                                            SIJAB->matrix[Gmi][im][ab] -= Z[ab][m];
                                        }
                                    }

                                    free_dpd_block(Z, nrows, ncols);
                                }

                                /* t_MkBc <-- 1/2 t_IJkABc W_IJMA */
                                /* sort W(AB,c) to W(A,Bc) */
                                for(Gab=0; Gab < nirreps; Gab++) {
                                    Gc = Gab ^ Gijk ^ GX3;
                                    for(ab=0; ab < FAA->params->coltot[Gab]; ab++ ){
                                        A = FAA->params->colorb[Gab][ab][0];
                                        B = FAA->params->colorb[Gab][ab][1];
                                        Ga = FAA->params->rsym[A];
                                        a = A - avir_off[Ga];
                                        for(c=0; c < bvirtpi[Gc]; c++) {
                                            C = bvir_off[Gc] + c;
                                            bc = SIjAb->params->colidx[B][C];
                                            W3[Ga][a][bc] = W1[Gab][ab][c];
                                        }
                                    }
                                }

                                ij = WMNIE->params->rowidx[I][J];

                                for(Gm=0; Gm < nirreps; Gm++) {
                                    Gmk = Gm ^ Gk;
                                    Gbc = Gmk ^ GS;
                                    Ga = Gbc ^ Gijk ^ GX3;

                                    ma = WMNIE->col_offset[Gij][Gm];

                                    nrows = SIjAb->params->coltot[Gmk^GS];
                                    ncols = aoccpi[Gm];
                                    nlinks = avirtpi[Ga];

                                    Z = dpd_block_matrix(nrows, ncols);

                                    if(nrows && ncols && nlinks)
                                        C_DGEMM('t', 't', nrows, ncols, nlinks, 0.5, W3[Ga][0], nrows,
                                                &(WMNIE->matrix[Gij][ij][ma]), nlinks, 0.0, Z[0], ncols);

                                    for(m=0; m < aoccpi[Gm]; m++) {
                                        M = aocc_off[Gm] + m;
                                        mk = SIjAb->params->rowidx[M][K];
                                        for(bc=0; bc < nrows; bc++) {
                                            SIjAb->matrix[Gmk][mk][bc] += Z[bc][m];
                                        }
                                    }

                                    free_dpd_block(Z, nrows, ncols);
                                }

                                /* S_ImBc <-- t_IJkABc W_kJmA */
                                /* sort W(AB,c) to W(A,Bc) */
                                for(Gab=0; Gab < nirreps; Gab++) {
                                    Gc = Gab ^ Gijk ^ GX3;  /* assumes totally symmetric */
                                    for(ab=0; ab < FAA->params->coltot[Gab]; ab++ ){
                                        A = FAA->params->colorb[Gab][ab][0];
                                        B = FAA->params->colorb[Gab][ab][1];
                                        Ga = FAA->params->rsym[A];
                                        a = A - avir_off[Ga];
                                        for(c=0; c < bvirtpi[Gc]; c++) {
                                            C = bvir_off[Gc] + c;
                                            bc = SIjAb->params->colidx[B][C];
                                            W3[Ga][a][bc] = W1[Gab][ab][c];
                                        }
                                    }
                                }

                                kj = WmNiE->params->rowidx[K][J];

                                for(Gm=0; Gm < nirreps; Gm++) {
                                    Gim = Gi ^ Gm;
                                    Gbc = Gim ^ GS;
                                    Ga = Gbc ^ Gijk ^ GX3;

                                    ma = WmNiE->col_offset[Gjk][Gm];

                                    nrows = SIjAb->params->coltot[Gim^GS];
                                    ncols = boccpi[Gm];
                                    nlinks = avirtpi[Ga];

                                    Z = dpd_block_matrix(nrows, ncols);

                                    if(nrows && ncols && nlinks)
                                        C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0, W3[Ga][0], nrows,
                                                &(WmNiE->matrix[Gjk][kj][ma]), nlinks, 0.0, Z[0], ncols);

                                    for(m=0; m < boccpi[Gm]; m++) {
                                        M = bocc_off[Gm] + m;
                                        im = SIjAb->params->rowidx[I][M];
                                        for(bc=0; bc < nrows; bc++) {
                                            SIjAb->matrix[Gim][im][bc] += Z[bc][m];
                                        }
                                    }

                                    free_dpd_block(Z, nrows, ncols);
                                }
                            }

                        } /* k */
                    } /* j */
                } /* i */

                for(Gab=0; Gab < nirreps; Gab++) {
                    /* This will need to change for non-totally-symmetric cases */
                    Gc = Gab ^ Gijk ^ GX3;
                    free_dpd_block(W1[Gab], FAA->params->coltot[Gab], bvirtpi[Gc]);
                }
                for(Ga=0; Ga < nirreps; Ga++) {
                    Gcb = Ga ^ Gijk ^ GX3; /* assumes totally symmetric */
                    free_dpd_block(W2[Ga], avirtpi[Ga], WmAfE->params->coltot[Gcb]);
                }
                for(Gb=0; Gb < nirreps; Gb++) {
                    Gac = Gb ^ Gijk ^ GX3; /* assumes totally symmetric */
                    free_dpd_block(W3[Gb], avirtpi[Gb], WMaFe->params->coltot[Gac]);
                }

            } /* Gk */
        } /* Gj */
    } /* Gi */

    free(W1);
    free(W2);
    free(W3);
    free_int_matrix(W_offset);

    file2_close(&fIJ);
    file2_close(&fAB);
    file2_close(&fij);
    file2_close(&fab);

    if (do_singles) {

        file2_mat_wrt(SIA);
        file2_mat_close(SIA);
        file2_mat_wrt(Sia);
        file2_mat_close(Sia);

        for(h=0; h < nirreps; h++) {
            buf4_mat_irrep_close(DAA, h);
            buf4_mat_irrep_close(DAB, h);
        }
    }

    if (do_doubles) {

        file2_mat_close(FME);
        file2_mat_close(Fme);

        for(h=0; h < nirreps; h++) {
            buf4_mat_irrep_close(WMNIE, h);
            buf4_mat_irrep_close(WMnIe, h);
            buf4_mat_irrep_close(WmNiE, h);
        }
        for(h=0; h < nirreps; h++) {
            buf4_mat_irrep_wrt(SIJAB, h);
            buf4_mat_irrep_close(SIJAB, h);
        }
        for(h=0; h < nirreps; h++) {
            buf4_mat_irrep_wrt(SIjAb, h);
            buf4_mat_irrep_close(SIjAb, h);
        }
    }
}

void DPD::cc3_sigma_UHF_BBA(dpdbuf4 *C2BB, dpdbuf4 *C2AB, dpdbuf4 *C2BA,
                            dpdbuf4 *FBB, dpdbuf4 *FAB, dpdbuf4 *FBA,
                            dpdbuf4 *EBB, dpdbuf4 *EAB, dpdbuf4 *EBA,
                            int do_singles, dpdbuf4 *DBB, dpdbuf4 *DBA, dpdfile2 *SIA, dpdfile2 *Sia,
                            int do_doubles, dpdfile2 *FME, dpdfile2 *Fme,
                            dpdbuf4 *Wmafe, dpdbuf4 *WMaFe, dpdbuf4 *WmAfE,
                            dpdbuf4 *Wmnie, dpdbuf4 *WMnIe, dpdbuf4 *WmNiE,
                            dpdbuf4 *Sijab, dpdbuf4 *SIjAb, int *aoccpi, int *aocc_off, int *boccpi,
                            int *bocc_off, int *avirtpi, int *avir_off, int *bvirtpi, int *bvir_off,
                            double omega, std::string out)
{
   std::shared_ptr<psi::PsiOutStream> printer=(out=="outfile"?outfile:
            std::shared_ptr<OutFile>(new OutFile(out)));
   int h, nirreps, S_irr;
    int Gi, Gj, Gk, Gijk;
    int Ga, Gb, Gc, Gab;
    int Gij, ij, ji, Gjk, jk, Gbc, bc, Gcb, cb;
    int Gd, Gkd, kd, d, DD, ad, da, dc, Gid, id, Gac, ac, bd;
    int i, j, k, I, J, K;
    int a, b, c, A, B, C;
    int ab;
    double ***W1, ***W2, ***W3;
    dpdfile2 fIJ, fAB, fij, fab;
    int nrows, ncols, nlinks;
    int **W_offset, offset;
    dpdbuf4 SiJaB;
    double **Z;
    int Gm, m, Gmi, mi, im, mc, M;
    int Gmk, mk, ma, kj, Gim;
    int GS, GH, GU, GC, GX3;

    GC = C2BB->file.my_irrep;
    GU = EBB->file.my_irrep;

    if (do_singles)
        GH = DBB->file.my_irrep;
    else if (do_doubles)
        GH = Wmafe->file.my_irrep;

    GS = Sijab->file.my_irrep;
    GX3 = GU^GC;

    nirreps = C2BB->params->nirreps;

    W_offset = init_int_matrix(nirreps, nirreps);
    for(Gab=0; Gab < nirreps; Gab++) {
        for(Ga=0,offset=0; Ga < nirreps; Ga++) {
            Gb = Ga ^ Gab;
            W_offset[Gab][Ga] = offset;
            offset += bvirtpi[Ga] * bvirtpi[Gb];
        }
    }

    if (do_singles) {

        file2_mat_init(SIA);
        file2_mat_rd(SIA);
        file2_mat_init(Sia);
        file2_mat_rd(Sia);

        for(h=0; h < nirreps; h++) {
            buf4_mat_irrep_init(DBB, h);
            buf4_mat_irrep_rd(DBB, h);
            buf4_mat_irrep_init(DBA, h);
            buf4_mat_irrep_rd(DBA, h);
        }
    }

    if (do_doubles) {

        file2_mat_init(FME);
        file2_mat_rd(FME);
        file2_mat_init(Fme);
        file2_mat_rd(Fme);

        for(h=0; h < nirreps; h++) {
            buf4_mat_irrep_init(WmNiE, h);
            buf4_mat_irrep_rd(WmNiE, h);

            buf4_mat_irrep_init(Wmnie, h);
            buf4_mat_irrep_rd(Wmnie, h);

            buf4_mat_irrep_init(WMnIe, h);
            buf4_mat_irrep_rd(WMnIe, h);
        }
        for(h=0; h < nirreps; h++) {
            buf4_mat_irrep_init(Sijab, h);
            buf4_mat_irrep_rd(Sijab, h);
        }

        /* put result here until end of function */
        S_irr = SIjAb->file.my_irrep;
        buf4_init(&SiJaB, PSIF_CC_TMP0, S_irr, 23, 29, 23, 29, 0, "CC3 SiJaB");
        for(h=0; h < nirreps; h++) {
            buf4_mat_irrep_init(&SiJaB, h);
        }
    }

    file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
    file2_init(&fij, PSIF_CC_OEI, 0, 2, 2, "fij");
    file2_init(&fab, PSIF_CC_OEI, 0, 3, 3, "fab");

    /* target T3 amplitudes go in here */
    W1 = (double ***) malloc(nirreps * sizeof(double **));
    W2 = (double ***) malloc(nirreps * sizeof(double **));
    W3 = (double ***) malloc(nirreps * sizeof(double **));

    for(Gi=0; Gi < nirreps; Gi++) {
        for(Gj=0; Gj < nirreps; Gj++) {
            Gij = Gi ^ Gj;
            for(Gk=0; Gk < nirreps; Gk++) {
                Gijk = Gi ^ Gj ^ Gk;
                Gjk = Gj ^ Gk;

                for(Gab=0; Gab < nirreps; Gab++) {
                    Gc = Gab ^ Gijk ^ GX3;
                    W1[Gab] = dpd_block_matrix(FBB->params->coltot[Gab], avirtpi[Gc]);
                }
                for(Ga=0; Ga < nirreps; Ga++) {
                    Gcb = Ga ^ Gijk ^ GX3;
                    W2[Ga] = dpd_block_matrix(bvirtpi[Ga], WMaFe->params->coltot[Gcb]);
                }
                for(Gb=0; Gb < nirreps; Gb++) {
                    Gac = Gb ^ Gijk ^ GX3;
                    W3[Gb] = dpd_block_matrix(bvirtpi[Gb], WmAfE->params->coltot[Gac]);
                }

                for(i=0; i < boccpi[Gi]; i++) {
                    I = bocc_off[Gi] + i;
                    for(j=0; j < boccpi[Gj]; j++) {
                        J = bocc_off[Gj] + j;
                        for(k=0; k < aoccpi[Gk]; k++) {
                            K = aocc_off[Gk] + k;

                            T3_AAB(W1, nirreps, I, Gi, J, Gj, K, Gk, C2BB, C2BA, C2AB,
                                   FBB, FBA, FAB, EBB, EBA, EAB, &fij, &fIJ, &fab, &fAB,
                                   boccpi, bocc_off, aoccpi, aocc_off, bvirtpi, bvir_off, avirtpi, avir_off, omega);

                            if (do_singles) {
                                /* S_KC <-- 1/4 t_ijKabC <ij||ab> */

                                Gc = Gk ^ GS;
                                Gab = Gij ^ GH;

                                ij = DBB->params->rowidx[I][J];

                                nrows = DBB->params->coltot[Gij^GH];
                                ncols = avirtpi[Gc];

                                if(nrows && ncols)
                                    C_DGEMV('t', nrows, ncols, 0.25, W1[Gab][0], ncols, DBB->matrix[Gij][ij], 1,
                                            1.0, SIA->matrix[Gk][k], 1);

                                /* S_ia <-- t_ijKabC <jK|bC> */

                                Ga = Gi ^ GS;
                                Gbc = Gjk ^ GH;

                                jk = DBA->params->rowidx[J][K];

                                for(Gab=0; Gab < nirreps; Gab++) {
                                    Gb = Ga ^ Gab;
                                    Gc = Gb ^ Gbc;

                                    ab = W_offset[Gab][Ga];
                                    bc = DBA->col_offset[Gjk][Gb];

                                    nrows = bvirtpi[Ga];
                                    ncols = bvirtpi[Gb] * avirtpi[Gc];

                                    if(nrows && ncols)
                                        C_DGEMV('n', nrows, ncols, 1.0, W1[Gab][ab], ncols, &(DBA->matrix[Gjk][jk][bc]), 1,
                                                1.0, Sia->matrix[Gi][i], 1);
                                }
                            } /* end do_singles */

                            if (do_doubles ) {
                                /* S_ijab <-- t_ijKabC F_KC */
                                Gc = Gk ^ GH;
                                Gab = Gij ^ GS;

                                nrows = Sijab->params->coltot[Gij^GS];
                                ncols = avirtpi[Gc];
                                ij = Sijab->params->rowidx[I][J];

                                if(nrows && ncols)
                                    C_DGEMV('n', nrows, ncols, 1.0, W1[Gab][0], ncols, FME->matrix[Gk][k], 1,
                                            1.0, Sijab->matrix[Gij][ij], 1);

                                /* S_jKbC <-- t_ijKabC F_ia */
                                Ga = Gi ^ GH;
                                Gbc = Gjk ^ GS;

                                jk = C2BA->params->rowidx[J][K];

                                for(Gab=0; Gab < nirreps; Gab++) {
                                    Gb = Ga ^ Gab;
                                    Gc = Gb ^ Gbc;

                                    ab = W_offset[Gab][Ga];
                                    bc = SiJaB.col_offset[Gjk][Gb];

                                    nrows = bvirtpi[Ga];
                                    ncols = bvirtpi[Gb] * avirtpi[Gc];

                                    if(nrows && ncols)
                                        C_DGEMV('t', nrows, ncols, 1.0, W1[Gab][ab], ncols, Fme->matrix[Gi][i], 1,
                                                1.0, &(SiJaB.matrix[Gjk][jk][bc]), 1);
                                }

                                /* S_jida <-- t_ijKabC W_KdCb */
                                /* sort W1(ab,C) to W2(a,Cb) */
                                for(Gab=0; Gab < nirreps; Gab++) {
                                    Gc = Gab ^ Gijk ^ GX3;
                                    for(ab=0; ab < FBB->params->coltot[Gab]; ab++) {
                                        A = FBB->params->colorb[Gab][ab][0];
                                        B = FBB->params->colorb[Gab][ab][1];
                                        Ga = FBB->params->rsym[A];
                                        a = A - bvir_off[Ga];
                                        for(c=0; c < avirtpi[Gc]; c++) {
                                            C = avir_off[Gc] + c;
                                            cb = WMaFe->params->colidx[C][B];
                                            W2[Ga][a][cb] = W1[Gab][ab][c];
                                        }
                                    }
                                }

                                ji = Sijab->params->rowidx[J][I];

                                for(Gd=0; Gd < nirreps; Gd++) {
                                    Gkd = Gk ^ Gd;
                                    Gcb = Gkd ^ GH;
                                    Ga = Gd ^ Gij ^ GS;

                                    kd = WMaFe->row_offset[Gkd][K];
                                    WMaFe->matrix[Gkd] = dpd_block_matrix(bvirtpi[Gd], WMaFe->params->coltot[Gkd^GH]);
                                    buf4_mat_irrep_rd_block(WMaFe, Gkd, kd, bvirtpi[Gd]);
                                    Z = dpd_block_matrix(bvirtpi[Ga], bvirtpi[Gd]);

                                    nrows = bvirtpi[Ga];
                                    ncols = bvirtpi[Gd];
                                    nlinks = WMaFe->params->coltot[Gkd^GH];

                                    if(nrows && ncols && nlinks)
                                        C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, W2[Ga][0], nlinks,
                                                WMaFe->matrix[Gkd][0], nlinks, 0.0, Z[0], ncols);

                                    for(a=0; a < bvirtpi[Ga]; a++) {
                                        A = bvir_off[Ga] + a;
                                        for(d=0; d < bvirtpi[Gd]; d++) {
                                            DD = bvir_off[Gd] + d;
                                            ad = Sijab->params->colidx[A][DD];
                                            da = Sijab->params->colidx[DD][A];
                                            Sijab->matrix[Gij][ji][ad] += -Z[a][d];
                                            Sijab->matrix[Gij][ji][da] += Z[a][d];
                                        }
                                    }

                                    free_dpd_block(WMaFe->matrix[Gkd], bvirtpi[Gd], WMaFe->params->coltot[Gkd^GH]);
                                    free_dpd_block(Z, bvirtpi[Ga], bvirtpi[Gd]);
                                }

                                /* t_jKcD <-- 1/2 t_ijKabC W_idab */

                                jk = SiJaB.params->rowidx[J][K];

                                for(Gd=0; Gd < nirreps; Gd++) {
                                    Gid = Gi ^ Gd;
                                    Gab = Gid ^ GH;
                                    Gc = Gab ^ Gijk ^ GX3;

                                    id = Wmafe->row_offset[Gid][I];
                                    Wmafe->matrix[Gid] = dpd_block_matrix(bvirtpi[Gd], Wmafe->params->coltot[Gid^GH]);
                                    buf4_mat_irrep_rd_block(Wmafe, Gid, id, bvirtpi[Gd]);
                                    Z = dpd_block_matrix(avirtpi[Gc], bvirtpi[Gd]);

                                    nrows = avirtpi[Gc];
                                    ncols = bvirtpi[Gd];
                                    nlinks = Wmafe->params->coltot[Gid^GH];

                                    if(nrows && ncols && nlinks)
                                        C_DGEMM('t', 't', nrows, ncols, nlinks, 0.5, W1[Gab][0], nrows,
                                                Wmafe->matrix[Gid][0], nlinks, 0.0, Z[0], ncols);

                                    for(c=0; c < avirtpi[Gc]; c++) {
                                        C = avir_off[Gc] + c;
                                        for(d=0; d < bvirtpi[Gd]; d++) {
                                            DD = bvir_off[Gd] + d;
                                            dc = SiJaB.params->colidx[DD][C];
                                            SiJaB.matrix[Gjk][jk][dc] += Z[c][d];
                                        }
                                    }

                                    free_dpd_block(Z, avirtpi[Gc], bvirtpi[Gd]);
                                    free_dpd_block(Wmafe->matrix[Gid], bvirtpi[Gd], Wmafe->params->coltot[Gid^GH]);
                                }

                                /* t_jKbD <-- t_ijKabC W_iDaC */
                                /* sort W1(ab,C) to W3(b,aC) */
                                for(Gab=0; Gab < nirreps; Gab++) {
                                    Gc = Gab ^ Gijk ^ GX3;
                                    for(ab=0; ab < FBB->params->coltot[Gab]; ab++) {
                                        A = FBB->params->colorb[Gab][ab][0];
                                        B = FBB->params->colorb[Gab][ab][1];
                                        Gb = FBB->params->ssym[B];
                                        b = B - bvir_off[Gb];
                                        for(c=0; c < avirtpi[Gc]; c++) {
                                            C = avir_off[Gc] + c;
                                            ac = WmAfE->params->colidx[A][C];
                                            W3[Gb][b][ac] = W1[Gab][ab][c];
                                        }
                                    }
                                }

                                jk = SiJaB.params->rowidx[J][K];

                                for(Gd=0; Gd < nirreps; Gd++) {
                                    Gid = Gi ^ Gd;
                                    Gac = Gid ^ GH;
                                    Gb = Gac ^ Gijk ^ GX3;

                                    id = WmAfE->row_offset[Gid][I];
                                    WmAfE->matrix[Gid] = dpd_block_matrix(avirtpi[Gd], WmAfE->params->coltot[Gid^GH]);
                                    buf4_mat_irrep_rd_block(WmAfE, Gid, id, avirtpi[Gd]);
                                    Z = dpd_block_matrix(bvirtpi[Gb], avirtpi[Gd]);

                                    nrows = bvirtpi[Gb];
                                    ncols = avirtpi[Gd];
                                    nlinks = WmAfE->params->coltot[Gid^GH];

                                    if(nrows && ncols && nlinks)
                                        C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, W3[Gb][0], nlinks,
                                                WmAfE->matrix[Gid][0], nlinks, 0.0, Z[0], ncols);

                                    for(b=0; b < bvirtpi[Gb]; b++) {
                                        B = bvir_off[Gb] + b;
                                        for(d=0; d < avirtpi[Gd]; d++) {
                                            DD = avir_off[Gd] + d;
                                            bd = SiJaB.params->colidx[B][DD];
                                            SiJaB.matrix[Gjk][jk][bd] += Z[b][d];
                                        }
                                    }

                                    free_dpd_block(WmAfE->matrix[Gid], avirtpi[Gd], WmAfE->params->coltot[Gid^GH]);
                                    free_dpd_block(Z, bvirtpi[Gb], avirtpi[Gd]);
                                }

                                /* S_miab <--- +t_ijKabC W_jKmC */
                                /* S_imab <--- -t_ijKabC W_jKmC */

                                jk = WmNiE->params->rowidx[J][K];

                                for(Gm=0; Gm < nirreps; Gm++) {
                                    Gmi = Gm ^ Gi;
                                    Gab = Gmi ^ GS;
                                    Gc = Gab ^ Gijk ^ GX3;

                                    mc = WmNiE->col_offset[Gjk][Gm];

                                    nrows = Sijab->params->coltot[Gab];
                                    ncols = boccpi[Gm];
                                    nlinks = avirtpi[Gc];

                                    Z = dpd_block_matrix(nrows, ncols);

                                    if(nrows && ncols && nlinks)
                                        C_DGEMM('n', 't', nrows, ncols, nlinks, 1.0, W1[Gab][0], nlinks,
                                                &(WmNiE->matrix[Gjk][jk][mc]), nlinks, 0.0, Z[0], ncols);

                                    for(m=0; m < ncols; m++) {
                                        M = bocc_off[Gm] + m;
                                        mi = Sijab->params->rowidx[M][I];
                                        im = Sijab->params->rowidx[I][M];
                                        for(ab=0; ab < nrows; ab++) {
                                            Sijab->matrix[Gmi][mi][ab] += Z[ab][m];
                                            Sijab->matrix[Gmi][im][ab] -= Z[ab][m];
                                        }
                                    }

                                    free_dpd_block(Z, nrows, ncols);
                                }

                                /* t_mKbC <-- 1/2 t_ijKabC W_ijma */
                                /* sort W(ab,C) to W(a,bC) */
                                for(Gab=0; Gab < nirreps; Gab++) {
                                    Gc = Gab ^ Gijk ^ GX3;
                                    for(ab=0; ab < FBB->params->coltot[Gab]; ab++ ){
                                        A = FBB->params->colorb[Gab][ab][0];
                                        B = FBB->params->colorb[Gab][ab][1];
                                        Ga = FBB->params->rsym[A];
                                        a = A - bvir_off[Ga];
                                        for(c=0; c < avirtpi[Gc]; c++) {
                                            C = avir_off[Gc] + c;
                                            bc = SiJaB.params->colidx[B][C];
                                            W3[Ga][a][bc] = W1[Gab][ab][c];
                                        }
                                    }
                                }

                                ij = Wmnie->params->rowidx[I][J];

                                for(Gm=0; Gm < nirreps; Gm++) {
                                    Gmk = Gm ^ Gk;
                                    Gbc = Gmk ^ GS;
                                    Ga = Gbc ^ Gijk ^ GX3;

                                    ma = Wmnie->col_offset[Gij][Gm];

                                    nrows = SiJaB.params->coltot[Gmk^GS];
                                    ncols = boccpi[Gm];
                                    nlinks = bvirtpi[Ga];

                                    Z = dpd_block_matrix(nrows, ncols);

                                    if(nrows && ncols && nlinks)
                                        C_DGEMM('t', 't', nrows, ncols, nlinks, 0.5, W3[Ga][0], nrows,
                                                &(Wmnie->matrix[Gij][ij][ma]), nlinks, 0.0, Z[0], ncols);

                                    for(m=0; m < boccpi[Gm]; m++) {
                                        M = bocc_off[Gm] + m;
                                        mk = SiJaB.params->rowidx[M][K];
                                        for(bc=0; bc < nrows; bc++) {
                                            SiJaB.matrix[Gmk][mk][bc] += Z[bc][m];
                                        }
                                    }

                                    free_dpd_block(Z, nrows, ncols);
                                }

                                /* S_iMbC <-- t_ijKabC W_KjMa */
                                /* sort W(ab,C) to W(a,bC) */
                                for(Gab=0; Gab < nirreps; Gab++) {
                                    Gc = Gab ^ Gijk ^ GX3;
                                    for(ab=0; ab < FBB->params->coltot[Gab]; ab++ ){
                                        A = FBB->params->colorb[Gab][ab][0];
                                        B = FBB->params->colorb[Gab][ab][1];
                                        Ga = FBB->params->rsym[A];
                                        a = A - bvir_off[Ga];
                                        for(c=0; c < avirtpi[Gc]; c++) {
                                            C = avir_off[Gc] + c;
                                            bc = SiJaB.params->colidx[B][C];
                                            W3[Ga][a][bc] = W1[Gab][ab][c];
                                        }
                                    }
                                }

                                kj = WMnIe->params->rowidx[K][J];

                                for(Gm=0; Gm < nirreps; Gm++) {
                                    Gim = Gi ^ Gm;
                                    Gbc = Gim ^ GS;
                                    Ga = Gbc ^ Gijk ^ GX3;

                                    ma = WMnIe->col_offset[Gjk][Gm];

                                    nrows = SiJaB.params->coltot[Gim^GS];
                                    ncols = aoccpi[Gm];
                                    nlinks = bvirtpi[Ga];

                                    Z = dpd_block_matrix(nrows, ncols);

                                    if(nrows && ncols && nlinks)
                                        C_DGEMM('t', 't', nrows, ncols, nlinks, 1.0, W3[Ga][0], nrows,
                                                &(WMnIe->matrix[Gjk][kj][ma]), nlinks, 0.0, Z[0], ncols);

                                    for(m=0; m < aoccpi[Gm]; m++) {
                                        M = aocc_off[Gm] + m;
                                        im = SiJaB.params->rowidx[I][M];
                                        for(bc=0; bc < nrows; bc++) {
                                            SiJaB.matrix[Gim][im][bc] += Z[bc][m];
                                        }
                                    }

                                    free_dpd_block(Z, nrows, ncols);
                                }
                            } /* end do_doubles */

                        } /* k */
                    } /* j */
                } /* i */

                for(Gab=0; Gab < nirreps; Gab++) {
                    Gc = Gab ^ Gijk ^ GX3;
                    free_dpd_block(W1[Gab], FBB->params->coltot[Gab], avirtpi[Gc]);
                }
                for(Ga=0; Ga < nirreps; Ga++) {
                    Gcb = Ga ^ Gijk ^ GX3;
                    free_dpd_block(W2[Ga], bvirtpi[Ga], WMaFe->params->coltot[Gcb]);
                }
                for(Gb=0; Gb < nirreps; Gb++) {
                    Gac = Gb ^ Gijk ^ GX3;
                    free_dpd_block(W3[Gb], bvirtpi[Gb], WmAfE->params->coltot[Gac]);
                }
            } /* Gk */
        } /* Gj */
    } /* Gi */

    free(W1);
    free(W2);
    free(W3);
    free_int_matrix(W_offset);

    file2_close(&fIJ);
    file2_close(&fAB);
    file2_close(&fij);
    file2_close(&fab);

    if (do_singles) {

        file2_mat_wrt(SIA);
        file2_mat_close(SIA);
        file2_mat_wrt(Sia);
        file2_mat_close(Sia);

        for(h=0; h < nirreps; h++) {
            buf4_mat_irrep_close(DBB, h);
            buf4_mat_irrep_close(DBA, h);
        }
    }

    if (do_doubles) {

        file2_mat_close(FME);
        file2_mat_close(Fme);

        for(h=0; h < nirreps; h++) {
            buf4_mat_irrep_close(WmNiE, h);
            buf4_mat_irrep_close(Wmnie, h);
            buf4_mat_irrep_close(WMnIe, h);
        }

        for(h=0; h < nirreps; h++) {
            buf4_mat_irrep_wrt(Sijab, h);
            buf4_mat_irrep_close(Sijab, h);
        }
        for(h=0; h < nirreps; h++) {
            buf4_mat_irrep_wrt(&SiJaB, h);
            buf4_mat_irrep_close(&SiJaB, h);
        }
        buf4_sort(&SiJaB, PSIF_CC_TMP0, qpsr, 22, 28, "CC3 SIjAb");
        buf4_close(&SiJaB);
        buf4_init(&SiJaB, PSIF_CC_TMP0, S_irr, 22, 28, 22, 28, 0, "CC3 SIjAb");
        buf4_axpy(&SiJaB, SIjAb, 1.0);
        buf4_close(&SiJaB);
    }

}

}
