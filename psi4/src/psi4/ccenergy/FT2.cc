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
    \ingroup CCENERGY
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <cstdlib>
#include "psi4/libdpd/dpd.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libqt/qt.h"
#include "Params.h"
#include "MOInfo.h"
#include "ccwave.h"

namespace psi { namespace ccenergy {

void CCEnergyWavefunction::FT2(void)
{
    dpdfile2 tIA, tia, t1;
    dpdbuf4 newtIJAB, newtijab, newtIjAb, t2, t2a, t2b;
    dpdbuf4 F_anti, F;
    dpdbuf4 Z, X;
    int Gie, Gij, Gab, nrows, ncols, nlinks, Gi, Ge, Gj, i, I;

    if(params_.ref == 0) { /** RHF **/

        /*     dpd_buf4_init(&newtIjAb, CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb"); */

        /*     dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z1(ij,ab)"); */
        /*     dpd_buf4_init(&F, CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>"); */
        /*     dpd_file2_init(&tIA, CC_OEI, 0, 0, 1, "tIA"); */
        /*     dpd_contract244(&tIA, &F, &Z, 1, 0, 0, 1, 0); */
        /*     dpd_file2_close(&tIA);  */
        /*     dpd_buf4_close(&F); */

        /*     dpd_buf4_sort(&Z, CC_TMP0, qpsr, 0, 5, "Z2(ji,ba)"); */
        /*     dpd_buf4_axpy(&Z, &newtIjAb, 1.0); */
        /*     dpd_buf4_close(&Z); */
        /*     dpd_buf4_init(&Z, CC_TMP0, 0, 0, 5, 0, 5, 0, "Z2(ji,ba)"); */
        /*     dpd_buf4_axpy(&Z, &newtIjAb, 1.0); */
        /*     dpd_buf4_close(&Z); */

        /*     dpd_buf4_close(&newtIjAb); */

        /* t(ij,ab) <-- t(j,e) * <ie|ab> + t(i,e) * <je|ba> */
        /* OOC code added 3/23/05, TDC */
        if(params_.df){
            dpdbuf4 OV, VV, Tov, Tovov, Toovv;
            // (Q|ia) = (Q|ea) t_ie
            dpd_set_default(1);
            global_dpd_->buf4_init(&VV, PSIF_CC_OEI, 0, 10, 43, 13, 43, 0, "B(VV|Q)");
            global_dpd_->buf4_init(&Tov, PSIF_CC_TMP0, 0, 27, 43, 27, 43, 0, "T(OV|Q)");
            dpd_set_default(0);
            global_dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "tIA");
            global_dpd_->contract244(&t1, &VV, &Tov, 1, 0, 0, 1.0, 0.0);
            global_dpd_->file2_close(&t1);
            global_dpd_->buf4_close(&VV);

            // t_jbia = (Q|jb) (Q|ia)
            dpd_set_default(1);
            global_dpd_->buf4_init(&OV, PSIF_CC_OEI, 0, 43, 27, 43, 27, 0, "B(Q|OV)");
            global_dpd_->buf4_init(&Tovov, PSIF_CC_TMP0, 0, 27, 27, 27, 27, 0, "T(OV|OV)");
            dpd_set_default(0);
            global_dpd_->contract444(&Tov, &OV, &Tovov, 0, 1, 1.0, 0.0);
            global_dpd_->buf4_sort(&Tovov, PSIF_CC_TMP0, prqs, 0, 5, "T<OO|VV>");
            global_dpd_->buf4_close(&OV);
            global_dpd_->buf4_close(&Tovov);

            global_dpd_->buf4_init(&Toovv, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "T<OO|VV>");
            global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
            // T_IjAb -> T2 IjAb
            global_dpd_->buf4_axpy(&Toovv, &t2, 1.0);
            // T_IjAb -> T2 jIBa
            global_dpd_->buf4_sort_axpy(&Toovv, PSIF_CC_TAMPS, qpsr, 0, 5, "New tIjAb", 1.0);
            global_dpd_->buf4_close(&t2);
            global_dpd_->buf4_close(&Toovv);
        }else{
            global_dpd_->buf4_init(&X, PSIF_CC_TMP0, 0, 0, 5, 0, 5, 0, "X(Ij,Ab)");
            global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
            global_dpd_->file2_init(&t1, PSIF_CC_OEI, 0, 0, 1, "tIA");
            global_dpd_->file2_mat_init(&t1);
            global_dpd_->file2_mat_rd(&t1);
            for(Gie=0; Gie < moinfo_.nirreps; Gie++) {
                Gab = Gie; /* F is totally symmetric */
                Gij = Gab; /* T2 is totally symmetric */
                global_dpd_->buf4_mat_irrep_init(&X, Gij);
                ncols = F.params->coltot[Gie];

                for(Gi=0; Gi < moinfo_.nirreps; Gi++) {
                    Gj = Ge = Gi^Gie; /* T1 is totally symmetric */

                    nlinks = moinfo_.virtpi[Ge];
                    nrows = moinfo_.occpi[Gj];

                    global_dpd_->buf4_mat_irrep_init_block(&F, Gie, nlinks);

                    for(i=0; i < moinfo_.occpi[Gi]; i++) {
                        I = F.params->poff[Gi] + i;
                        global_dpd_->buf4_mat_irrep_rd_block(&F, Gie, F.row_offset[Gie][I], nlinks);

                        if(nrows && ncols && nlinks)
                            C_DGEMM('n','n',nrows,ncols,nlinks,1.0,t1.matrix[Gj][0],nlinks,F.matrix[Gie][0],ncols,
                                    0.0,X.matrix[Gij][X.row_offset[Gij][I]],ncols);
                    }

                    global_dpd_->buf4_mat_irrep_close_block(&F, Gie, nlinks);
                }

                global_dpd_->buf4_mat_irrep_wrt(&X, Gij);
                global_dpd_->buf4_mat_irrep_close(&X, Gij);
            }
            global_dpd_->file2_mat_close(&t1);
            global_dpd_->file2_close(&t1);
            global_dpd_->buf4_close(&F);
            global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");
            global_dpd_->buf4_axpy(&X, &t2, 1);
            global_dpd_->buf4_close(&t2);
            global_dpd_->buf4_sort_axpy(&X, PSIF_CC_TAMPS, qpsr, 0, 5, "New tIjAb", 1);
            global_dpd_->buf4_close(&X);
        }
    }
    else if(params_.ref == 1) { /** ROHF **/

        global_dpd_->buf4_init(&newtIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "New tIJAB");
        global_dpd_->buf4_init(&newtijab, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "New tijab");
        global_dpd_->buf4_init(&newtIjAb, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

        global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
        global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 0, 1, "tia");

        /*** AA ***/

        global_dpd_->buf4_init(&F_anti, PSIF_CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
        global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
        global_dpd_->contract424(&F_anti, &tIA, &t2, 1, 1, 1, 1, 0);
        global_dpd_->buf4_sort(&t2, PSIF_CC_TMP0, qprs, 0, 7, "T (JI,A>B)");
        global_dpd_->buf4_close(&t2);
        global_dpd_->buf4_init(&t2a, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
        global_dpd_->buf4_init(&t2b, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (JI,A>B)");
        global_dpd_->buf4_axpy(&t2b, &t2a, -1);
        global_dpd_->buf4_axpy(&t2a, &newtIJAB, 1);
        global_dpd_->buf4_close(&t2b);
        global_dpd_->buf4_close(&t2a);
        global_dpd_->buf4_close(&F_anti);

        /*** BB ***/

        global_dpd_->buf4_init(&F_anti, PSIF_CC_FINTS, 0, 10, 7, 10, 5, 1, "F <ia|bc>");
        global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
        global_dpd_->contract424(&F_anti, &tia, &t2, 1, 1, 1, 1, 0);
        global_dpd_->buf4_sort(&t2, PSIF_CC_TMP0, qprs, 0, 7, "T (JI,A>B)");
        global_dpd_->buf4_close(&t2);
        global_dpd_->buf4_init(&t2a, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
        global_dpd_->buf4_init(&t2b, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (JI,A>B)");
        global_dpd_->buf4_axpy(&t2b, &t2a, -1);
        global_dpd_->buf4_axpy(&t2a, &newtijab, 1);
        global_dpd_->buf4_close(&t2b);
        global_dpd_->buf4_close(&t2a);
        global_dpd_->buf4_close(&F_anti);

        /*** AB ***/

        global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 10, 5, 10, 5, 0, "F <ia|bc>");
        global_dpd_->contract424(&F, &tia, &newtIjAb, 1, 1, 1, 1, 1);
        global_dpd_->buf4_close(&F);
        global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 11, 5, 11, 5, 0, "F <ai|bc>");
        global_dpd_->contract244(&tIA, &F, &newtIjAb, 1, 0, 0, 1, 1);
        global_dpd_->buf4_close(&F);

        global_dpd_->file2_close(&tIA);
        global_dpd_->file2_close(&tia);

        global_dpd_->buf4_close(&newtIJAB);
        global_dpd_->buf4_close(&newtijab);
        global_dpd_->buf4_close(&newtIjAb);
    }
    else if(params_.ref == 2) { /*** UHF ***/

        global_dpd_->buf4_init(&newtIJAB, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "New tIJAB");
        global_dpd_->buf4_init(&newtijab, PSIF_CC_TAMPS, 0, 10, 17, 12, 17, 0, "New tijab");
        global_dpd_->buf4_init(&newtIjAb, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");

        global_dpd_->file2_init(&tIA, PSIF_CC_OEI, 0, 0, 1, "tIA");
        global_dpd_->file2_init(&tia, PSIF_CC_OEI, 0, 2, 3, "tia");

        /*** AA ***/

        global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 20, 7, 20, 5, 1, "F <IA|BC>");
        global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
        global_dpd_->contract424(&F, &tIA, &t2, 1, 1, 1, 1, 0);
        global_dpd_->buf4_sort(&t2, PSIF_CC_TMP0, qprs, 0, 7, "T (JI,A>B)");
        global_dpd_->buf4_close(&t2);
        global_dpd_->buf4_init(&t2a, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (IJ,A>B)");
        global_dpd_->buf4_init(&t2b, PSIF_CC_TMP0, 0, 0, 7, 0, 7, 0, "T (JI,A>B)");
        global_dpd_->buf4_axpy(&t2b, &t2a, -1);
        global_dpd_->buf4_axpy(&t2a, &newtIJAB, 1);
        global_dpd_->buf4_close(&t2b);
        global_dpd_->buf4_close(&t2a);
        global_dpd_->buf4_close(&F);

        /*** BB ***/

        global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 30, 17, 30, 15, 1, "F <ia|bc>");
        global_dpd_->buf4_init(&t2, PSIF_CC_TMP0, 0, 10, 17, 10, 17, 0, "T (ij,a>b)");
        global_dpd_->contract424(&F, &tia, &t2, 1, 1, 1, 1, 0);
        global_dpd_->buf4_sort(&t2, PSIF_CC_TMP0, qprs, 10, 17, "T (ji,a>b)");
        global_dpd_->buf4_close(&t2);
        global_dpd_->buf4_init(&t2a, PSIF_CC_TMP0, 0, 10, 17, 10, 17, 0, "T (ij,a>b)");
        global_dpd_->buf4_init(&t2b, PSIF_CC_TMP0, 0, 10, 17, 10, 17, 0, "T (ji,a>b)");
        global_dpd_->buf4_axpy(&t2b, &t2a, -1);
        global_dpd_->buf4_axpy(&t2a, &newtijab, 1);
        global_dpd_->buf4_close(&t2b);
        global_dpd_->buf4_close(&t2a);
        global_dpd_->buf4_close(&F);

        /*** AB ***/

        global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 24, 28, 24, 28, 0, "F <Ia|Bc>");
        global_dpd_->contract424(&F, &tia, &newtIjAb, 1, 1, 1, 1, 1);
        global_dpd_->buf4_close(&F);
        global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 28, 26, 28, 26, 0, "F <Ab|Ci>");
        global_dpd_->contract244(&tIA, &F, &newtIjAb, 1, 2, 0, 1, 1);
        global_dpd_->buf4_close(&F);

        global_dpd_->file2_close(&tIA);
        global_dpd_->file2_close(&tia);

        global_dpd_->buf4_close(&newtIJAB);
        global_dpd_->buf4_close(&newtijab);
        global_dpd_->buf4_close(&newtIjAb);

    }

}
}} // namespace psi::ccenergy
