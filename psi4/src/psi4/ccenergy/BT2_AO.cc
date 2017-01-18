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
#include <cstring>
#include <cmath>
#include "psi4/libciomr/libciomr.h"
#include "psi4/libpsio/psio.h"
#include "psi4/libiwl/iwl.h"
#include "psi4/libdpd/dpd.h"
#include "psi4/libqt/qt.h"
#include "psi4/psifiles.h"
#include "Params.h"
#include "MOInfo.h"
#include "ccwave.h"

namespace psi { namespace ccenergy {


void CCEnergyWavefunction::BT2_AO(void)
{
    int h, nirreps, i, Gc, Gd, Ga, Gb, ij;
    double ***C, **X;
    double ***Ca, ***Cb;
    int *sopi, *virtpi;
    int *avirtpi, *bvirtpi;
    int **T2_cd_row_start, **T2_pq_row_start, offset, cd, pq;
    int **T2_CD_row_start, **T2_Cd_row_start;
    dpdbuf4 tau, t2, tau1_AO, tau2_AO;
    dpdfile4 T;
    psio_address next;
    struct iwlbuf InBuf;
    int lastbuf;
    double tolerance=1e-14;
    double **integrals;
    int **tau1_cols, **tau2_cols, *num_ints;
    int counter=0, counterAA=0, counterBB=0, counterAB=0;

    nirreps = moinfo_.nirreps;
    sopi = moinfo_.sopi;

    T2_pq_row_start = init_int_matrix(nirreps,nirreps);
    for(h=0; h < nirreps; h++) {
        for(Gc=0,offset=0; Gc < nirreps; Gc++) {
            Gd = Gc ^ h;
            T2_pq_row_start[h][Gc] = offset;
            offset += sopi[Gc] * sopi[Gd];
        }
    }

    if(params_.ref == 0 || params_.ref == 1) { /** RHF or ROHF **/
        virtpi = moinfo_.virtpi;
        C = moinfo_.Cv;

        T2_cd_row_start = init_int_matrix(nirreps,nirreps);
        for(h=0; h < nirreps; h++) {
            for(Gc=0,offset=0; Gc < nirreps; Gc++) {
                Gd = Gc ^ h;
                T2_cd_row_start[h][Gc] = offset;
                offset += virtpi[Gc] * virtpi[Gd];
            }
        }
    }
    else if(params_.ref == 2) {  /** UHF **/
        avirtpi = moinfo_.avirtpi;
        bvirtpi = moinfo_.bvirtpi;
        Ca = moinfo_.Cav;
        Cb = moinfo_.Cbv;

        T2_CD_row_start = init_int_matrix(nirreps,nirreps);
        for(h=0; h < nirreps; h++) {
            for(Gc=0,offset=0; Gc < nirreps; Gc++) {
                Gd = Gc ^ h;
                T2_CD_row_start[h][Gc] = offset;
                offset += avirtpi[Gc] * avirtpi[Gd];
            }
        }
        T2_cd_row_start = init_int_matrix(nirreps,nirreps);
        for(h=0; h < nirreps; h++) {
            for(Gc=0,offset=0; Gc < nirreps; Gc++) {
                Gd = Gc ^ h;
                T2_cd_row_start[h][Gc] = offset;
                offset += bvirtpi[Gc] * bvirtpi[Gd];
            }
        }
        T2_Cd_row_start = init_int_matrix(nirreps,nirreps);
        for(h=0; h < nirreps; h++) {
            for(Gc=0,offset=0; Gc < nirreps; Gc++) {
                Gd = Gc ^ h;
                T2_Cd_row_start[h][Gc] = offset;
                offset += avirtpi[Gc] * bvirtpi[Gd];
            }
        }

    }

    if(params_.ref == 0) { /** RHF **/

        if(params_.aobasis == "DISK") {

            dpd_set_default(1);
            global_dpd_->buf4_init(&tau1_AO, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjPq (1)");
            global_dpd_->buf4_scm(&tau1_AO, 0.0);

            dpd_set_default(0);
            global_dpd_->buf4_init(&tau, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");

            halftrans(&tau, 0, &tau1_AO, 1, C, C, nirreps, T2_cd_row_start, T2_pq_row_start,
                      virtpi, virtpi, sopi, 0, 1.0, 0.0);

            global_dpd_->buf4_close(&tau);
            global_dpd_->buf4_close(&tau1_AO);

            /* Transpose tau1_AO for better memory access patterns */
            dpd_set_default(1);
            global_dpd_->buf4_init(&tau1_AO, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjPq (1)");
            global_dpd_->buf4_sort(&tau1_AO, PSIF_CC_TMP0, rspq, 5, 0, "tauPqIj (1)");
            global_dpd_->buf4_close(&tau1_AO);


            global_dpd_->buf4_init(&tau1_AO, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "tauPqIj (1)");
            global_dpd_->buf4_init(&tau2_AO, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "tauPqIj (2)");
            global_dpd_->buf4_scm(&tau2_AO, 0.0);

            if(params_.df){
                dpdbuf4 B;
                // 5 = unpacked. eventually use perm sym and pair number 8
                global_dpd_->buf4_init(&B, PSIF_CC_OEI, 0, 5, 43, 8, 43, 0, "B(pq|Q)");
                global_dpd_->contract444_df(&B, &tau1_AO, &tau2_AO, 1.0, 0.0);
                global_dpd_->buf4_close(&B);
            }else{
                for(h=0; h < nirreps; h++) {
                    global_dpd_->buf4_mat_irrep_init(&tau1_AO, h);
                    global_dpd_->buf4_mat_irrep_rd(&tau1_AO, h);
                    global_dpd_->buf4_mat_irrep_init(&tau2_AO, h);
                }

                iwl_buf_init(&InBuf, PSIF_SO_TEI, tolerance, 1, 1);

                lastbuf = InBuf.lastbuf;

                counter += AO_contribute(&InBuf, &tau1_AO, &tau2_AO);

                while(!lastbuf) {
                    iwl_buf_fetch(&InBuf);
                    lastbuf = InBuf.lastbuf;

                    counter += AO_contribute(&InBuf, &tau1_AO, &tau2_AO);
                }

                iwl_buf_close(&InBuf, 1);

                if(params_.print & 2) outfile->Printf( "     *** Processed %d SO integrals for <ab||cd> --> T2\n", counter);

                for(h=0; h < nirreps; h++) {
                    global_dpd_->buf4_mat_irrep_wrt(&tau2_AO, h);
                    global_dpd_->buf4_mat_irrep_close(&tau2_AO, h);
                    global_dpd_->buf4_mat_irrep_close(&tau1_AO, h);
                }
            }
            global_dpd_->buf4_close(&tau1_AO);
//            global_dpd_->buf4_print(&tau2_AO, outfile, 1);
//            exit(1);
            global_dpd_->buf4_close(&tau2_AO);

            /* Transpose tau2_AO for the half-backtransformation */
            dpd_set_default(1);
            global_dpd_->buf4_init(&tau2_AO, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "tauPqIj (2)");
            global_dpd_->buf4_sort(&tau2_AO, PSIF_CC_TAMPS, rspq, 0, 5, "tauIjPq (2)");
            global_dpd_->buf4_close(&tau2_AO);

            global_dpd_->buf4_init(&tau2_AO, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjPq (2)");

            dpd_set_default(0);
            global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

            halftrans(&t2, 0, &tau2_AO, 1, C, C, nirreps, T2_cd_row_start, T2_pq_row_start,
                      virtpi, virtpi, sopi, 1, 1.0, 1.0);

            global_dpd_->buf4_close(&t2);
            global_dpd_->buf4_close(&tau2_AO);

        }
        else if(params_.aobasis == "DIRECT") {
            global_dpd_->file4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, "tauIjAb");
            global_dpd_->file4_cache_del(&T);
            global_dpd_->file4_close(&T);

            global_dpd_->file4_init(&T, PSIF_CC_TAMPS, 0, 0, 5, "New tIjAb");
            global_dpd_->file4_cache_del(&T);
            global_dpd_->file4_close(&T);

            /* close the CC_TAMPS file for cints to use it */
            psio_close(PSIF_CC_TAMPS, 1);

            int statusvalue=system("cints --cc_bt2");

            /* re-open CCC_TAMPS for remaining terms */
            psio_open(PSIF_CC_TAMPS, PSIO_OPEN_OLD);
        }

    }
    else if(params_.ref == 1) { /** ROHF **/

        /************************************* AA *****************************************/

        dpd_set_default(1);
        global_dpd_->buf4_init(&tau1_AO, PSIF_CC_TAMPS, 0, 2, 5, 2, 5, 0, "tauIJPQ (1)");
        global_dpd_->buf4_scm(&tau1_AO, 0.0);

        dpd_set_default(0);
        global_dpd_->buf4_init(&tau, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tauIJAB");

        halftrans(&tau, 0, &tau1_AO, 1, C, C, nirreps, T2_cd_row_start, T2_pq_row_start,
                  virtpi, virtpi, sopi, 0, 1.0, 0.0);

        global_dpd_->buf4_close(&tau);
        global_dpd_->buf4_close(&tau1_AO);

        /* Transpose tau1_AO for better memory access patterns */
        dpd_set_default(1);
        global_dpd_->buf4_init(&tau1_AO, PSIF_CC_TAMPS, 0, 2, 5, 2, 5, 1, "tauIJPQ (1)");
        global_dpd_->buf4_sort(&tau1_AO, PSIF_CC_TMP0, rspq, 5, 2, "tauPQIJ (1)");
        global_dpd_->buf4_close(&tau1_AO);

        global_dpd_->buf4_init(&tau1_AO, PSIF_CC_TMP0, 0, 5, 2, 5, 2, 0, "tauPQIJ (1)");
        global_dpd_->buf4_init(&tau2_AO, PSIF_CC_TMP0, 0, 5, 2, 5, 2, 0, "tauPQIJ (2)");
        global_dpd_->buf4_scm(&tau2_AO, 0.0);

        for(h=0; h < nirreps; h++) {
            global_dpd_->buf4_mat_irrep_init(&tau1_AO, h);
            global_dpd_->buf4_mat_irrep_rd(&tau1_AO, h);
            global_dpd_->buf4_mat_irrep_init(&tau2_AO, h);
        }

        iwl_buf_init(&InBuf, PSIF_SO_TEI, tolerance, 1, 1);

        lastbuf = InBuf.lastbuf;

        counterAA += AO_contribute(&InBuf, &tau1_AO, &tau2_AO);

        while(!lastbuf) {
            iwl_buf_fetch(&InBuf);
            lastbuf = InBuf.lastbuf;

            counterAA += AO_contribute(&InBuf, &tau1_AO, &tau2_AO);
        }

        iwl_buf_close(&InBuf, 1);

        if(params_.print & 2) outfile->Printf( "     *** Processed %d SO integrals for <AB||CD> --> T2\n", counterAA);

        for(h=0; h < nirreps; h++) {
            global_dpd_->buf4_mat_irrep_wrt(&tau2_AO, h);
            global_dpd_->buf4_mat_irrep_close(&tau2_AO, h);
            global_dpd_->buf4_mat_irrep_close(&tau1_AO, h);
        }
        global_dpd_->buf4_close(&tau1_AO);
        global_dpd_->buf4_close(&tau2_AO);


        /* Transpose tau2_AO for the half-backtransformation */
        dpd_set_default(1);
        global_dpd_->buf4_init(&tau2_AO, PSIF_CC_TMP0, 0, 5, 2, 5, 2, 0, "tauPQIJ (2)");
        global_dpd_->buf4_sort(&tau2_AO, PSIF_CC_TAMPS, rspq, 2, 5, "tauIJPQ (2)");
        global_dpd_->buf4_close(&tau2_AO);

        global_dpd_->buf4_init(&tau2_AO, PSIF_CC_TAMPS, 0, 2, 5, 2, 5, 0, "tauIJPQ (2)");

        dpd_set_default(0);
        global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tIJAB");

        halftrans(&t2, 0, &tau2_AO, 1, C, C, nirreps, T2_cd_row_start, T2_pq_row_start,
                  virtpi, virtpi, sopi, 1, 0.5, 1.0);

        global_dpd_->buf4_close(&t2);
        global_dpd_->buf4_close(&tau2_AO);

        /************************************* BB *****************************************/

        dpd_set_default(1);
        global_dpd_->buf4_init(&tau1_AO, PSIF_CC_TAMPS, 0, 2, 5, 2, 5, 0, "tauijpq (1)");
        global_dpd_->buf4_scm(&tau1_AO, 0.0);

        dpd_set_default(0);
        global_dpd_->buf4_init(&tau, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tauijab");

        halftrans(&tau, 0, &tau1_AO, 1, C, C, nirreps, T2_cd_row_start, T2_pq_row_start,
                  virtpi, virtpi, sopi, 0, 1.0, 0.0);

        global_dpd_->buf4_close(&tau);
        global_dpd_->buf4_close(&tau1_AO);

        /* Transpose tau1_AO for better memory access patterns */
        dpd_set_default(1);
        global_dpd_->buf4_init(&tau1_AO, PSIF_CC_TAMPS, 0, 2, 5, 2, 5, 1, "tauijpq (1)");
        global_dpd_->buf4_sort(&tau1_AO, PSIF_CC_TMP0, rspq, 5, 2, "taupqij (1)");
        global_dpd_->buf4_close(&tau1_AO);

        global_dpd_->buf4_init(&tau1_AO, PSIF_CC_TMP0, 0, 5, 2, 5, 2, 0, "taupqij (1)");
        global_dpd_->buf4_init(&tau2_AO, PSIF_CC_TMP0, 0, 5, 2, 5, 2, 0, "taupqij (2)");
        global_dpd_->buf4_scm(&tau2_AO, 0.0);

        for(h=0; h < nirreps; h++) {
            global_dpd_->buf4_mat_irrep_init(&tau1_AO, h);
            global_dpd_->buf4_mat_irrep_rd(&tau1_AO, h);
            global_dpd_->buf4_mat_irrep_init(&tau2_AO, h);
        }

        iwl_buf_init(&InBuf, PSIF_SO_TEI, tolerance, 1, 1);

        lastbuf = InBuf.lastbuf;

        counterBB += AO_contribute(&InBuf, &tau1_AO, &tau2_AO);

        while(!lastbuf) {
            iwl_buf_fetch(&InBuf);
            lastbuf = InBuf.lastbuf;

            counterBB += AO_contribute(&InBuf, &tau1_AO, &tau2_AO);
        }

        iwl_buf_close(&InBuf, 1);

        if(params_.print & 2) outfile->Printf( "     *** Processed %d SO integrals for <ab||cd> --> T2\n", counterBB);

        for(h=0; h < nirreps; h++) {
            global_dpd_->buf4_mat_irrep_wrt(&tau2_AO, h);
            global_dpd_->buf4_mat_irrep_close(&tau2_AO, h);
            global_dpd_->buf4_mat_irrep_close(&tau1_AO, h);
        }
        global_dpd_->buf4_close(&tau1_AO);
        global_dpd_->buf4_close(&tau2_AO);


        /* Transpose tau2_AO for the half-backtransformation */
        dpd_set_default(1);
        global_dpd_->buf4_init(&tau2_AO, PSIF_CC_TMP0, 0, 5, 2, 5, 2, 0, "taupqij (2)");
        global_dpd_->buf4_sort(&tau2_AO, PSIF_CC_TAMPS, rspq, 2, 5, "tauijpq (2)");
        global_dpd_->buf4_close(&tau2_AO);

        global_dpd_->buf4_init(&tau2_AO, PSIF_CC_TAMPS, 0, 2, 5, 2, 5, 0, "tauijpq (2)");

        dpd_set_default(0);
        global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tijab");

        halftrans(&t2, 0, &tau2_AO, 1, C, C, nirreps, T2_cd_row_start, T2_pq_row_start,
                  virtpi, virtpi, sopi, 1, 0.5, 1.0);

        global_dpd_->buf4_close(&t2);
        global_dpd_->buf4_close(&tau2_AO);

        /************************************* AB *****************************************/

        dpd_set_default(1);
        global_dpd_->buf4_init(&tau1_AO, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjPq (1)");
        global_dpd_->buf4_scm(&tau1_AO, 0.0);

        dpd_set_default(0);
        global_dpd_->buf4_init(&tau, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjAb");

        halftrans(&tau, 0, &tau1_AO, 1, C, C, nirreps, T2_cd_row_start, T2_pq_row_start,
                  virtpi, virtpi, sopi, 0, 1.0, 0.0);

        global_dpd_->buf4_close(&tau);
        global_dpd_->buf4_close(&tau1_AO);


        /* Transpose tau1_AO for better memory access patterns */
        dpd_set_default(1);
        global_dpd_->buf4_init(&tau1_AO, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjPq (1)");
        global_dpd_->buf4_sort(&tau1_AO, PSIF_CC_TMP0, rspq, 5, 0, "tauPqIj (1)");
        global_dpd_->buf4_close(&tau1_AO);

        global_dpd_->buf4_init(&tau1_AO, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "tauPqIj (1)");
        global_dpd_->buf4_init(&tau2_AO, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "tauPqIj (2)");
        global_dpd_->buf4_scm(&tau2_AO, 0.0);

        for(h=0; h < nirreps; h++) {
            global_dpd_->buf4_mat_irrep_init(&tau1_AO, h);
            global_dpd_->buf4_mat_irrep_rd(&tau1_AO, h);
            global_dpd_->buf4_mat_irrep_init(&tau2_AO, h);
        }

        iwl_buf_init(&InBuf, PSIF_SO_TEI, tolerance, 1, 1);

        lastbuf = InBuf.lastbuf;

        counterAB += AO_contribute(&InBuf, &tau1_AO, &tau2_AO);

        while(!lastbuf) {
            iwl_buf_fetch(&InBuf);
            lastbuf = InBuf.lastbuf;

            counterAB += AO_contribute(&InBuf, &tau1_AO, &tau2_AO);
        }

        iwl_buf_close(&InBuf, 1);

        if(params_.print & 2) outfile->Printf( "     *** Processed %d SO integrals for <Ab|Cd> --> T2\n", counterAB);

        for(h=0; h < nirreps; h++) {
            global_dpd_->buf4_mat_irrep_wrt(&tau2_AO, h);
            global_dpd_->buf4_mat_irrep_close(&tau2_AO, h);
            global_dpd_->buf4_mat_irrep_close(&tau1_AO, h);
        }
        global_dpd_->buf4_close(&tau1_AO);
        global_dpd_->buf4_close(&tau2_AO);

        /* Transpose tau2_AO for the half-backtransformation */
        dpd_set_default(1);
        global_dpd_->buf4_init(&tau2_AO, PSIF_CC_TMP0, 0, 5, 0, 5, 0, 0, "tauPqIj (2)");
        global_dpd_->buf4_sort(&tau2_AO, PSIF_CC_TAMPS, rspq, 0, 5, "tauIjPq (2)");
        global_dpd_->buf4_close(&tau2_AO);

        global_dpd_->buf4_init(&tau2_AO, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "tauIjPq (2)");

        dpd_set_default(0);
        global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 0, 5, 0, 5, 0, "New tIjAb");

        halftrans(&t2, 0, &tau2_AO, 1, C, C, nirreps, T2_cd_row_start, T2_pq_row_start,
                  virtpi, virtpi, sopi, 1, 1.0, 1.0);

        global_dpd_->buf4_close(&t2);
        global_dpd_->buf4_close(&tau2_AO);

    }  /** ROHF **/
    else if(params_.ref == 2) { /** UHF **/

        /************************************* AA *****************************************/

        dpd_set_default(1);
        global_dpd_->buf4_init(&tau1_AO, PSIF_CC_TAMPS, 0, 2, 5, 2, 5, 0, "tauIJPQ (1)");
        global_dpd_->buf4_scm(&tau1_AO, 0.0);

        dpd_set_default(0);
        global_dpd_->buf4_init(&tau, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "tauIJAB");

        halftrans(&tau, 0, &tau1_AO, 1, Ca, Ca, nirreps, T2_CD_row_start, T2_pq_row_start,
                  avirtpi, avirtpi, sopi, 0, 1.0, 0.0);

        global_dpd_->buf4_close(&tau);
        global_dpd_->buf4_close(&tau1_AO);

        /* Transpose tau1_AO for better memory access patterns */
        dpd_set_default(1);
        global_dpd_->buf4_init(&tau1_AO, PSIF_CC_TAMPS, 0, 2, 5, 2, 5, 1, "tauIJPQ (1)");
        global_dpd_->buf4_sort(&tau1_AO, PSIF_CC_TMP0, rspq, 5, 2, "tauPQIJ (1)");
        global_dpd_->buf4_close(&tau1_AO);

        global_dpd_->buf4_init(&tau1_AO, PSIF_CC_TMP0, 0, 5, 2, 5, 2, 0, "tauPQIJ (1)");
        global_dpd_->buf4_init(&tau2_AO, PSIF_CC_TMP0, 0, 5, 2, 5, 2, 0, "tauPQIJ (2)");
        global_dpd_->buf4_scm(&tau2_AO, 0.0);

        for(h=0; h < nirreps; h++) {
            global_dpd_->buf4_mat_irrep_init(&tau1_AO, h);
            global_dpd_->buf4_mat_irrep_rd(&tau1_AO, h);
            global_dpd_->buf4_mat_irrep_init(&tau2_AO, h);
        }

        iwl_buf_init(&InBuf, PSIF_SO_TEI, tolerance, 1, 1);

        lastbuf = InBuf.lastbuf;

        counterAA += AO_contribute(&InBuf, &tau1_AO, &tau2_AO);

        while(!lastbuf) {
            iwl_buf_fetch(&InBuf);
            lastbuf = InBuf.lastbuf;

            counterAA += AO_contribute(&InBuf, &tau1_AO, &tau2_AO);
        }

        iwl_buf_close(&InBuf, 1);

        if(params_.print & 2) outfile->Printf( "     *** Processed %d SO integrals for <AB||CD> --> T2\n", counterAA);

        for(h=0; h < nirreps; h++) {
            global_dpd_->buf4_mat_irrep_wrt(&tau2_AO, h);
            global_dpd_->buf4_mat_irrep_close(&tau2_AO, h);
            global_dpd_->buf4_mat_irrep_close(&tau1_AO, h);
        }
        global_dpd_->buf4_close(&tau1_AO);
        global_dpd_->buf4_close(&tau2_AO);


        /* Transpose tau2_AO for the half-backtransformation */
        dpd_set_default(1);
        global_dpd_->buf4_init(&tau2_AO, PSIF_CC_TMP0, 0, 5, 2, 5, 2, 0, "tauPQIJ (2)");
        global_dpd_->buf4_sort(&tau2_AO, PSIF_CC_TAMPS, rspq, 2, 5, "tauIJPQ (2)");
        global_dpd_->buf4_close(&tau2_AO);

        global_dpd_->buf4_init(&tau2_AO, PSIF_CC_TAMPS, 0, 2, 5, 2, 5, 0, "tauIJPQ (2)");

        dpd_set_default(0);
        global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 2, 5, 2, 7, 0, "New tIJAB");

        halftrans(&t2, 0, &tau2_AO, 1, Ca, Ca, nirreps, T2_CD_row_start, T2_pq_row_start,
                  avirtpi, avirtpi, sopi, 1, 0.5, 1.0);

        global_dpd_->buf4_close(&t2);
        global_dpd_->buf4_close(&tau2_AO);

        /************************************* BB *****************************************/

        dpd_set_default(1);
        global_dpd_->buf4_init(&tau1_AO, PSIF_CC_TAMPS, 0, 12, 15, 12, 15, 0, "tauijpq (1)");
        global_dpd_->buf4_scm(&tau1_AO, 0.0);

        dpd_set_default(0);
        global_dpd_->buf4_init(&tau, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "tauijab");

        halftrans(&tau, 0, &tau1_AO, 1, Cb, Cb, nirreps, T2_cd_row_start, T2_pq_row_start,
                  bvirtpi, bvirtpi, sopi, 0, 1.0, 0.0);

        global_dpd_->buf4_close(&tau);
        global_dpd_->buf4_close(&tau1_AO);

        /* Transpose tau1_AO for better memory access patterns */
        dpd_set_default(1);
        global_dpd_->buf4_init(&tau1_AO, PSIF_CC_TAMPS, 0, 12, 15, 12, 15, 1, "tauijpq (1)");
        global_dpd_->buf4_sort(&tau1_AO, PSIF_CC_TMP0, rspq, 15, 12, "taupqij (1)");
        global_dpd_->buf4_close(&tau1_AO);

        global_dpd_->buf4_init(&tau1_AO, PSIF_CC_TMP0, 0, 15, 12, 15, 12, 0, "taupqij (1)");
        global_dpd_->buf4_init(&tau2_AO, PSIF_CC_TMP0, 0, 15, 12, 15, 12, 0, "taupqij (2)");
        global_dpd_->buf4_scm(&tau2_AO, 0.0);

        for(h=0; h < nirreps; h++) {
            global_dpd_->buf4_mat_irrep_init(&tau1_AO, h);
            global_dpd_->buf4_mat_irrep_rd(&tau1_AO, h);
            global_dpd_->buf4_mat_irrep_init(&tau2_AO, h);
        }

        iwl_buf_init(&InBuf, PSIF_SO_TEI, tolerance, 1, 1);

        lastbuf = InBuf.lastbuf;

        counterBB += AO_contribute(&InBuf, &tau1_AO, &tau2_AO);

        while(!lastbuf) {
            iwl_buf_fetch(&InBuf);
            lastbuf = InBuf.lastbuf;

            counterBB += AO_contribute(&InBuf, &tau1_AO, &tau2_AO);
        }

        iwl_buf_close(&InBuf, 1);

        if(params_.print & 2) outfile->Printf( "     *** Processed %d SO integrals for <ab||cd> --> T2\n", counterBB);

        for(h=0; h < nirreps; h++) {
            global_dpd_->buf4_mat_irrep_wrt(&tau2_AO, h);
            global_dpd_->buf4_mat_irrep_close(&tau2_AO, h);
            global_dpd_->buf4_mat_irrep_close(&tau1_AO, h);
        }
        global_dpd_->buf4_close(&tau1_AO);
        global_dpd_->buf4_close(&tau2_AO);


        /* Transpose tau2_AO for the half-backtransformation */
        dpd_set_default(1);
        global_dpd_->buf4_init(&tau2_AO, PSIF_CC_TMP0, 0, 15, 12, 15, 12, 0, "taupqij (2)");
        global_dpd_->buf4_sort(&tau2_AO, PSIF_CC_TAMPS, rspq, 12, 15, "tauijpq (2)");
        global_dpd_->buf4_close(&tau2_AO);

        global_dpd_->buf4_init(&tau2_AO, PSIF_CC_TAMPS, 0, 12, 15, 12, 15, 0, "tauijpq (2)");

        dpd_set_default(0);
        global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 12, 15, 12, 17, 0, "New tijab");

        halftrans(&t2, 0, &tau2_AO, 1, Cb, Cb, nirreps, T2_cd_row_start, T2_pq_row_start,
                  bvirtpi, bvirtpi, sopi, 1, 0.5, 1.0);

        global_dpd_->buf4_close(&t2);
        global_dpd_->buf4_close(&tau2_AO);

        /************************************* AB *****************************************/

        dpd_set_default(1);
        global_dpd_->buf4_init(&tau1_AO, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjPq (1)");
        global_dpd_->buf4_scm(&tau1_AO, 0.0);

        dpd_set_default(0);
        global_dpd_->buf4_init(&tau, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjAb");

        halftrans(&tau, 0, &tau1_AO, 1, Ca, Cb, nirreps, T2_Cd_row_start, T2_pq_row_start,
                  avirtpi, bvirtpi, sopi, 0, 1.0, 0.0);

        global_dpd_->buf4_close(&tau);
        global_dpd_->buf4_close(&tau1_AO);


        /* Transpose tau1_AO for better memory access patterns */
        dpd_set_default(1);
        global_dpd_->buf4_init(&tau1_AO, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjPq (1)");
        global_dpd_->buf4_sort(&tau1_AO, PSIF_CC_TMP0, rspq, 28, 22, "tauPqIj (1)");
        global_dpd_->buf4_close(&tau1_AO);

        global_dpd_->buf4_init(&tau1_AO, PSIF_CC_TMP0, 0, 28, 22, 28, 22, 0, "tauPqIj (1)");
        global_dpd_->buf4_init(&tau2_AO, PSIF_CC_TMP0, 0, 28, 22, 28, 22, 0, "tauPqIj (2)");
        global_dpd_->buf4_scm(&tau2_AO, 0.0);

        for(h=0; h < nirreps; h++) {
            global_dpd_->buf4_mat_irrep_init(&tau1_AO, h);
            global_dpd_->buf4_mat_irrep_rd(&tau1_AO, h);
            global_dpd_->buf4_mat_irrep_init(&tau2_AO, h);
        }

        iwl_buf_init(&InBuf, PSIF_SO_TEI, tolerance, 1, 1);

        lastbuf = InBuf.lastbuf;

        counterAB += AO_contribute(&InBuf, &tau1_AO, &tau2_AO);

        while(!lastbuf) {
            iwl_buf_fetch(&InBuf);
            lastbuf = InBuf.lastbuf;

            counterAB += AO_contribute(&InBuf, &tau1_AO, &tau2_AO);
        }

        iwl_buf_close(&InBuf, 1);

        if(params_.print & 2) outfile->Printf( "     *** Processed %d SO integrals for <Ab|Cd> --> T2\n", counterAB);

        for(h=0; h < nirreps; h++) {
            global_dpd_->buf4_mat_irrep_wrt(&tau2_AO, h);
            global_dpd_->buf4_mat_irrep_close(&tau2_AO, h);
            global_dpd_->buf4_mat_irrep_close(&tau1_AO, h);
        }
        global_dpd_->buf4_close(&tau1_AO);
        global_dpd_->buf4_close(&tau2_AO);

        /* Transpose tau2_AO for the half-backtransformation */
        dpd_set_default(1);
        global_dpd_->buf4_init(&tau2_AO, PSIF_CC_TMP0, 0, 28, 22, 28, 22, 0, "tauPqIj (2)");
        global_dpd_->buf4_sort(&tau2_AO, PSIF_CC_TAMPS, rspq, 22, 28, "tauIjPq (2)");
        global_dpd_->buf4_close(&tau2_AO);

        global_dpd_->buf4_init(&tau2_AO, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "tauIjPq (2)");

        dpd_set_default(0);
        global_dpd_->buf4_init(&t2, PSIF_CC_TAMPS, 0, 22, 28, 22, 28, 0, "New tIjAb");

        halftrans(&t2, 0, &tau2_AO, 1, Ca, Cb, nirreps, T2_Cd_row_start, T2_pq_row_start,
                  avirtpi, bvirtpi, sopi, 1, 1.0, 1.0);

        global_dpd_->buf4_close(&t2);
        global_dpd_->buf4_close(&tau2_AO);

    }  /** UHF **/

    if(params_.ref == 0 || params_.ref == 1)
        free_int_matrix(T2_cd_row_start);
    else if(params_.ref ==2) {
        free_int_matrix(T2_CD_row_start);
        free_int_matrix(T2_cd_row_start);
        free_int_matrix(T2_Cd_row_start);
    }

    free_int_matrix(T2_pq_row_start);

    /* Reset the default dpd back to 0 --- this stuff gets really ugly */
    dpd_set_default(0);

}

}} // namespace psi::ccenergy
