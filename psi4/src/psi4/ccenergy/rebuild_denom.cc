/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
 */

/*! \file
    \ingroup CCENERGY
    \brief Rebuilds denominators in case.
    \warning This is only used for PCM-CC with PTED.
*/

#include "Params.h"
#include "MOInfo.h"
#include "ccwave.h"

#include <cstdio>
#include <cstdlib>

#include "psi4/libdpd/dpd.h"

namespace psi {
namespace ccenergy {

void rebuild_denom_rhf(const MOInfo &);
void rebuild_denom_uhf(const MOInfo &);

void CCEnergyWavefunction::rebuild_denom() {
    if (params_.ref == 2) {
        rebuild_denom_uhf(moinfo_);
    } else { /* RHF and ROHF */
        rebuild_denom_rhf(moinfo_);
    }
}

void rebuild_denom_rhf(const MOInfo &moinfo_) {
    int i, j, a, b, ij, ab;
    int I, J, A, B;
    int isym, jsym, asym, bsym;
    double fii, fjj, faa, fbb;

    int nirreps = moinfo_.nirreps;
    int *occpi = moinfo_.occpi;
    int *virtpi = moinfo_.virtpi;
    int *openpi = moinfo_.openpi;
    int *occ_off = moinfo_.occ_off;
    int *vir_off = moinfo_.vir_off;

    dpdfile2 fIJ;
    global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    global_dpd_->file2_mat_init(&fIJ);
    global_dpd_->file2_mat_rd(&fIJ);

    dpdfile2 fij;
    global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 0, 0, "fij");
    global_dpd_->file2_mat_init(&fij);
    global_dpd_->file2_mat_rd(&fij);

    dpdfile2 fAB;
    global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
    global_dpd_->file2_mat_init(&fAB);
    global_dpd_->file2_mat_rd(&fAB);

    dpdfile2 fab;
    global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 1, 1, "fab");
    global_dpd_->file2_mat_init(&fab);
    global_dpd_->file2_mat_rd(&fab);

    /* Alpha one-electron denominator */
    dpdfile2 dIA;
    global_dpd_->file2_init(&dIA, PSIF_CC_OEI, 0, 0, 1, "dIA");
    global_dpd_->file2_mat_init(&dIA);

    for (int h = 0; h < nirreps; h++) {
        for (i = 0; i < occpi[h]; i++) {
            fii = fIJ.matrix[h][i][i];

            for (a = 0; a < (virtpi[h] - openpi[h]); a++) {
                faa = fAB.matrix[h][a][a];

                dIA.matrix[h][i][a] = 1.0 / (fii - faa);
            }
        }
    }

    global_dpd_->file2_mat_wrt(&dIA);
    global_dpd_->file2_mat_close(&dIA);
    global_dpd_->file2_close(&dIA);

    /* Beta one-electron denominator */
    dpdfile2 dia;
    global_dpd_->file2_init(&dia, PSIF_CC_OEI, 0, 0, 1, "dia");
    global_dpd_->file2_mat_init(&dia);

    for (int h = 0; h < nirreps; h++) {
        for (i = 0; i < (occpi[h] - openpi[h]); i++) {
            fii = fij.matrix[h][i][i];

            for (a = 0; a < virtpi[h]; a++) {
                faa = fab.matrix[h][a][a];

                dia.matrix[h][i][a] = 1.0 / (fii - faa);
            }
        }
    }

    global_dpd_->file2_mat_wrt(&dia);
    global_dpd_->file2_mat_close(&dia);
    global_dpd_->file2_close(&dia);

    /* Alpha-alpha two-electron denominator */
    dpdfile4 dIJAB;
    global_dpd_->file4_init(&dIJAB, PSIF_CC_DENOM, 0, 1, 6, "dIJAB");

    for (int h = 0; h < nirreps; h++) {
        global_dpd_->file4_mat_irrep_init(&dIJAB, h);

        /* Loop over the rows */
        for (ij = 0; ij < dIJAB.params->rowtot[h]; ij++) {
            i = dIJAB.params->roworb[h][ij][0];
            j = dIJAB.params->roworb[h][ij][1];
            isym = dIJAB.params->psym[i];
            jsym = dIJAB.params->qsym[j];

            /* Convert to relative orbital index */
            I = i - occ_off[isym];
            J = j - occ_off[jsym];

            fii = fIJ.matrix[isym][I][I];
            fjj = fIJ.matrix[jsym][J][J];

            /* Loop over the columns */
            for (ab = 0; ab < dIJAB.params->coltot[h]; ab++) {
                a = dIJAB.params->colorb[h][ab][0];
                b = dIJAB.params->colorb[h][ab][1];
                asym = dIJAB.params->rsym[a];
                bsym = dIJAB.params->ssym[b];

                /* Convert to relative orbital index */
                A = a - vir_off[asym];
                B = b - vir_off[bsym];

                faa = fAB.matrix[asym][A][A];
                fbb = fAB.matrix[bsym][B][B];

                dIJAB.matrix[h][ij][ab] = ((A >= (virtpi[asym] - openpi[asym])) || (B >= (virtpi[bsym] - openpi[bsym]))
                                               ? 0.0
                                               : 1.0 / (fii + fjj - faa - fbb));
            }
        }

        global_dpd_->file4_mat_irrep_wrt(&dIJAB, h);
        global_dpd_->file4_mat_irrep_close(&dIJAB, h);
    }

    global_dpd_->file4_close(&dIJAB);

    /* Beta-beta two-electron denominator */
    dpdfile4 dijab;
    global_dpd_->file4_init(&dijab, PSIF_CC_DENOM, 0, 1, 6, "dijab");

    for (int h = 0; h < nirreps; h++) {
        global_dpd_->file4_mat_irrep_init(&dijab, h);

        /* Loop over the rows */
        for (ij = 0; ij < dijab.params->rowtot[h]; ij++) {
            i = dijab.params->roworb[h][ij][0];
            j = dijab.params->roworb[h][ij][1];
            isym = dijab.params->psym[i];
            jsym = dijab.params->qsym[j];

            /* Convert to relative orbital index */
            I = i - occ_off[isym];
            J = j - occ_off[jsym];

            fii = fij.matrix[isym][I][I];
            fjj = fij.matrix[jsym][J][J];

            /* Loop over the columns */
            for (ab = 0; ab < dijab.params->coltot[h]; ab++) {
                a = dijab.params->colorb[h][ab][0];
                b = dijab.params->colorb[h][ab][1];
                asym = dijab.params->rsym[a];
                bsym = dijab.params->ssym[b];

                /* Convert to relative orbital index */
                A = a - vir_off[asym];
                B = b - vir_off[bsym];

                faa = fab.matrix[asym][A][A];
                fbb = fab.matrix[bsym][B][B];

                dijab.matrix[h][ij][ab] = ((I >= (occpi[isym] - openpi[isym])) || (J >= (occpi[jsym] - openpi[jsym]))
                                               ? 0.0
                                               : 1.0 / (fii + fjj - faa - fbb));
            }
        }

        global_dpd_->file4_mat_irrep_wrt(&dijab, h);
        global_dpd_->file4_mat_irrep_close(&dijab, h);
    }

    global_dpd_->file4_close(&dijab);

    /* Alpha-beta two-electron denominator */
    dpdfile4 dIjAb;
    global_dpd_->file4_init(&dIjAb, PSIF_CC_DENOM, 0, 0, 5, "dIjAb");

    for (int h = 0; h < nirreps; h++) {
        global_dpd_->file4_mat_irrep_init(&dIjAb, h);

        /* Loop over the rows */
        for (ij = 0; ij < dIjAb.params->rowtot[h]; ij++) {
            i = dIjAb.params->roworb[h][ij][0];
            j = dIjAb.params->roworb[h][ij][1];
            isym = dIjAb.params->psym[i];
            jsym = dIjAb.params->qsym[j];

            /* Convert to relative orbital index */
            I = i - occ_off[isym];
            J = j - occ_off[jsym];

            fii = fIJ.matrix[isym][I][I];
            fjj = fij.matrix[jsym][J][J];

            /* Loop over the columns */
            for (ab = 0; ab < dIjAb.params->coltot[h]; ab++) {
                a = dIjAb.params->colorb[h][ab][0];
                b = dIjAb.params->colorb[h][ab][1];
                asym = dIjAb.params->rsym[a];
                bsym = dIjAb.params->ssym[b];

                /* Convert to relative orbital index */
                A = a - vir_off[asym];
                B = b - vir_off[bsym];

                faa = fAB.matrix[asym][A][A];
                fbb = fab.matrix[bsym][B][B];

                dIjAb.matrix[h][ij][ab] = ((A >= (virtpi[asym] - openpi[asym])) || (J >= (occpi[jsym] - openpi[jsym]))
                                               ? 0.0
                                               : 1.0 / (fii + fjj - faa - fbb));
            }
        }

        global_dpd_->file4_mat_irrep_wrt(&dIjAb, h);
        global_dpd_->file4_mat_irrep_close(&dIjAb, h);
    }

    global_dpd_->file4_close(&dIjAb);

    global_dpd_->file2_mat_close(&fIJ);
    global_dpd_->file2_mat_close(&fij);
    global_dpd_->file2_mat_close(&fAB);
    global_dpd_->file2_mat_close(&fab);
    global_dpd_->file2_close(&fIJ);
    global_dpd_->file2_close(&fij);
    global_dpd_->file2_close(&fAB);
    global_dpd_->file2_close(&fab);
}

void rebuild_denom_uhf(const MOInfo &moinfo_) {
    int i, j, a, b, ij, ab, I, J, A, B, isym, jsym, asym, bsym;
    double fii, fjj, faa, fbb;

    int nirreps = moinfo_.nirreps;
    int *aoccpi = moinfo_.aoccpi;
    int *avirtpi = moinfo_.avirtpi;
    int *boccpi = moinfo_.boccpi;
    int *bvirtpi = moinfo_.bvirtpi;
    int *aocc_off = moinfo_.aocc_off;
    int *bocc_off = moinfo_.bocc_off;
    int *avir_off = moinfo_.avir_off;
    int *bvir_off = moinfo_.bvir_off;

    dpdfile2 fIJ;
    global_dpd_->file2_init(&fIJ, PSIF_CC_OEI, 0, 0, 0, "fIJ");
    global_dpd_->file2_mat_init(&fIJ);
    global_dpd_->file2_mat_rd(&fIJ);

    dpdfile2 fij;
    global_dpd_->file2_init(&fij, PSIF_CC_OEI, 0, 2, 2, "fij");
    global_dpd_->file2_mat_init(&fij);
    global_dpd_->file2_mat_rd(&fij);

    dpdfile2 fAB;
    global_dpd_->file2_init(&fAB, PSIF_CC_OEI, 0, 1, 1, "fAB");
    global_dpd_->file2_mat_init(&fAB);
    global_dpd_->file2_mat_rd(&fAB);

    dpdfile2 fab;
    global_dpd_->file2_init(&fab, PSIF_CC_OEI, 0, 3, 3, "fab");
    global_dpd_->file2_mat_init(&fab);
    global_dpd_->file2_mat_rd(&fab);

    dpdfile2 dIA;
    global_dpd_->file2_init(&dIA, PSIF_CC_OEI, 0, 0, 1, "dIA");
    global_dpd_->file2_mat_init(&dIA);
    for (int h = 0; h < nirreps; h++) {
        for (int i = 0; i < aoccpi[h]; i++) {
            fii = fIJ.matrix[h][i][i];
            for (int a = 0; a < avirtpi[h]; a++) {
                faa = fAB.matrix[h][a][a];
                dIA.matrix[h][i][a] = 1.0 / (fii - faa);
            }
        }
    }
    global_dpd_->file2_mat_wrt(&dIA);
    global_dpd_->file2_mat_close(&dIA);
    global_dpd_->file2_close(&dIA);

    global_dpd_->file2_init(&dIA, PSIF_CC_OEI, 0, 2, 3, "dia");
    global_dpd_->file2_mat_init(&dIA);
    for (int h = 0; h < nirreps; h++) {
        for (i = 0; i < boccpi[h]; i++) {
            fii = fij.matrix[h][i][i];
            for (a = 0; a < bvirtpi[h]; a++) {
                faa = fab.matrix[h][a][a];
                dIA.matrix[h][i][a] = 1.0 / (fii - faa);
            }
        }
    }
    global_dpd_->file2_mat_wrt(&dIA);
    global_dpd_->file2_mat_close(&dIA);
    global_dpd_->file2_close(&dIA);

    dpdfile4 dIJAB;
    global_dpd_->file4_init(&dIJAB, PSIF_CC_DENOM, 0, 1, 6, "dIJAB");
    for (int h = 0; h < nirreps; h++) {
        global_dpd_->file4_mat_irrep_init(&dIJAB, h);
        for (ij = 0; ij < dIJAB.params->rowtot[h]; ij++) {
            i = dIJAB.params->roworb[h][ij][0];
            j = dIJAB.params->roworb[h][ij][1];
            isym = dIJAB.params->psym[i];
            jsym = dIJAB.params->qsym[j];
            I = i - aocc_off[isym];
            J = j - aocc_off[jsym];
            fii = fIJ.matrix[isym][I][I];
            fjj = fIJ.matrix[jsym][J][J];

            for (ab = 0; ab < dIJAB.params->coltot[h]; ab++) {
                a = dIJAB.params->colorb[h][ab][0];
                b = dIJAB.params->colorb[h][ab][1];
                asym = dIJAB.params->rsym[a];
                bsym = dIJAB.params->ssym[b];
                A = a - avir_off[asym];
                B = b - avir_off[bsym];
                faa = fAB.matrix[asym][A][A];
                fbb = fAB.matrix[bsym][B][B];

                dIJAB.matrix[h][ij][ab] = 1.0 / (fii + fjj - faa - fbb);
            }
        }
        global_dpd_->file4_mat_irrep_wrt(&dIJAB, h);
        global_dpd_->file4_mat_irrep_close(&dIJAB, h);
    }
    global_dpd_->file4_close(&dIJAB);

    global_dpd_->file4_init(&dIJAB, PSIF_CC_DENOM, 0, 11, 16, "dijab");
    for (int h = 0; h < nirreps; h++) {
        global_dpd_->file4_mat_irrep_init(&dIJAB, h);
        for (ij = 0; ij < dIJAB.params->rowtot[h]; ij++) {
            i = dIJAB.params->roworb[h][ij][0];
            j = dIJAB.params->roworb[h][ij][1];
            isym = dIJAB.params->psym[i];
            jsym = dIJAB.params->qsym[j];
            I = i - bocc_off[isym];
            J = j - bocc_off[jsym];
            fii = fij.matrix[isym][I][I];
            fjj = fij.matrix[jsym][J][J];

            for (ab = 0; ab < dIJAB.params->coltot[h]; ab++) {
                a = dIJAB.params->colorb[h][ab][0];
                b = dIJAB.params->colorb[h][ab][1];
                asym = dIJAB.params->rsym[a];
                bsym = dIJAB.params->ssym[b];
                A = a - bvir_off[asym];
                B = b - bvir_off[bsym];
                faa = fab.matrix[asym][A][A];
                fbb = fab.matrix[bsym][B][B];

                dIJAB.matrix[h][ij][ab] = 1.0 / (fii + fjj - faa - fbb);
            }
        }
        global_dpd_->file4_mat_irrep_wrt(&dIJAB, h);
        global_dpd_->file4_mat_irrep_close(&dIJAB, h);
    }
    global_dpd_->file4_close(&dIJAB);

    global_dpd_->file4_init(&dIJAB, PSIF_CC_DENOM, 0, 22, 28, "dIjAb");
    for (int h = 0; h < nirreps; h++) {
        global_dpd_->file4_mat_irrep_init(&dIJAB, h);
        for (ij = 0; ij < dIJAB.params->rowtot[h]; ij++) {
            i = dIJAB.params->roworb[h][ij][0];
            j = dIJAB.params->roworb[h][ij][1];
            isym = dIJAB.params->psym[i];
            jsym = dIJAB.params->qsym[j];
            I = i - aocc_off[isym];
            J = j - bocc_off[jsym];
            fii = fIJ.matrix[isym][I][I];
            fjj = fij.matrix[jsym][J][J];

            for (ab = 0; ab < dIJAB.params->coltot[h]; ab++) {
                a = dIJAB.params->colorb[h][ab][0];
                b = dIJAB.params->colorb[h][ab][1];
                asym = dIJAB.params->rsym[a];
                bsym = dIJAB.params->ssym[b];
                A = a - avir_off[asym];
                B = b - bvir_off[bsym];
                faa = fAB.matrix[asym][A][A];
                fbb = fab.matrix[bsym][B][B];

                dIJAB.matrix[h][ij][ab] = 1.0 / (fii + fjj - faa - fbb);
            }
        }
        global_dpd_->file4_mat_irrep_wrt(&dIJAB, h);
        global_dpd_->file4_mat_irrep_close(&dIJAB, h);
    }
    global_dpd_->file4_close(&dIJAB);

    global_dpd_->file2_mat_close(&fIJ);
    global_dpd_->file2_mat_close(&fij);
    global_dpd_->file2_mat_close(&fAB);
    global_dpd_->file2_mat_close(&fab);
    global_dpd_->file2_close(&fIJ);
    global_dpd_->file2_close(&fij);
    global_dpd_->file2_close(&fAB);
    global_dpd_->file2_close(&fab);
}
}
}  // namespace psi::ccenergy
