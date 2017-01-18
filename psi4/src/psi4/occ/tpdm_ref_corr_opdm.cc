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

/* This code includes reference and correlation opdm contributions. */
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libpsio/psio.hpp"
#include "psi4/libmints/matrix.h"
#include "occwave.h"
#include "defines.h"


using namespace std;

namespace psi{ namespace occwave{

/*=======================*/
/*  tpdm_ref()           */
/*=======================*/
void OCCWave::tpdm_ref()
{
    dpdbuf4 G;

    psio_->open(PSIF_OCC_DENSITY, PSIO_OPEN_OLD);

 if (reference_ == "RESTRICTED") {
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
      for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
	if (wfn_type_ != "OMP2") global_dpd_->buf4_mat_irrep_rd(&G, h);
        #pragma omp parallel for
        for(int ij = 0; ij < G.params->rowtot[h]; ++ij){
            int i = G.params->roworb[h][ij][0];
            int j = G.params->roworb[h][ij][1];
            for(int kl = 0; kl < G.params->coltot[h]; ++kl){
                int k = G.params->colorb[h][kl][0];
                int l = G.params->colorb[h][kl][1];
		if (i == k && j == l) G.matrix[h][ij][kl] += 1.0;
		if (i == l && j == k) G.matrix[h][ij][kl] -= 0.25;
		if (i == j && k == l) G.matrix[h][ij][kl] -= 0.25;
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }

    global_dpd_->buf4_close(&G);

 }// end RHF

 else if (reference_ == "UNRESTRICTED") {
    // Alpha-Alpha spin case
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
      for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
	if (wfn_type_ != "OMP2") global_dpd_->buf4_mat_irrep_rd(&G, h);
        #pragma omp parallel for
        for(int ij = 0; ij < G.params->rowtot[h]; ++ij){
            int i = G.params->roworb[h][ij][0];
            int j = G.params->roworb[h][ij][1];
            for(int kl = 0; kl < G.params->coltot[h]; ++kl){
                int k = G.params->colorb[h][kl][0];
                int l = G.params->colorb[h][kl][1];
		if (i == k && j == l) G.matrix[h][ij][kl] += 0.25;
		if (i == l && j == k) G.matrix[h][ij][kl] -= 0.25;
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);


    // Beta-Beta spin case
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "TPDM <oo|oo>");
      for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
	if (wfn_type_ != "OMP2") global_dpd_->buf4_mat_irrep_rd(&G, h);
        #pragma omp parallel for
        for(int ij = 0; ij < G.params->rowtot[h]; ++ij){
            int i = G.params->roworb[h][ij][0];
            int j = G.params->roworb[h][ij][1];
            for(int kl = 0; kl < G.params->coltot[h]; ++kl){
                int k = G.params->colorb[h][kl][0];
                int l = G.params->colorb[h][kl][1];
		if (i == k && j == l) G.matrix[h][ij][kl] += 0.25;
		if (i == l && j == k) G.matrix[h][ij][kl] -= 0.25;
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);


    // Alpha-Beta spin case
     global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                 ID("[O,o]"), ID("[O,o]"), 0, "TPDM <Oo|Oo>");
      for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
	if (wfn_type_ != "OMP2") global_dpd_->buf4_mat_irrep_rd(&G, h);
        #pragma omp parallel for
        for(int ij = 0; ij < G.params->rowtot[h]; ++ij){
            int i = G.params->roworb[h][ij][0];
            int j = G.params->roworb[h][ij][1];
            for(int kl = 0; kl < G.params->coltot[h]; ++kl){
                int k = G.params->colorb[h][kl][0];
                int l = G.params->colorb[h][kl][1];
		if (i == k && j == l) G.matrix[h][ij][kl] += 0.25;
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

 }// end UHF

    psio_->close(PSIF_OCC_DENSITY, 1);

} // end of twopdm_ref


/*=======================*/
/*  tpdm_corr_opdm()     */
/*=======================*/
void OCCWave::tpdm_corr_opdm()
{

    dpdbuf4 G;

    psio_->open(PSIF_OCC_DENSITY, PSIO_OPEN_OLD);

 if (reference_ == "RESTRICTED") {
    // TPDM <OO|OO>
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
      for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        #pragma omp parallel for
        for(int ij = 0; ij < G.params->rowtot[h]; ++ij){
            int i = G.params->roworb[h][ij][0];
            int j = G.params->roworb[h][ij][1];
            for(int kl = 0; kl < G.params->coltot[h]; ++kl){
                int k = G.params->colorb[h][kl][0];
                int l = G.params->colorb[h][kl][1];

		int hi = G.params->psym[i];
		int hj = G.params->qsym[j];
		int hk = G.params->rsym[k];
		int hl = G.params->ssym[l];

		int ii = i - G.params->poff[hi];
		int jj = j - G.params->qoff[hj];
		int kk = k - G.params->roff[hk];
		int ll = l - G.params->soff[hl];

		if (i == k && hj == hl) G.matrix[h][ij][kl] += 0.5 * gamma1corr->get(hj,jj,ll);
		if (j == l && hi == hk) G.matrix[h][ij][kl] += 0.5 * gamma1corr->get(hi,ii,kk);
		if (i == l && hj == hk) G.matrix[h][ij][kl] -= 0.125 * gamma1corr->get(hj,jj,kk);
		if (j == k && hi == hl) G.matrix[h][ij][kl] -= 0.125 * gamma1corr->get(hi,ii,ll);
		if (i == j && hk == hl) G.matrix[h][ij][kl] -= 0.125 * gamma1corr->get(hk,kk,ll);
		if (k == l && hi == hj) G.matrix[h][ij][kl] -= 0.125 * gamma1corr->get(hi,ii,jj);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

    //Print
    if (print_ > 3) {
      global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
      global_dpd_->buf4_print(&G, "outfile", 1);
      global_dpd_->buf4_close(&G);
    }


    // TPDM <OO|VV>
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
      for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        #pragma omp parallel for
        for(int ij = 0; ij < G.params->rowtot[h]; ++ij){
            int i = G.params->roworb[h][ij][0];
            int j = G.params->roworb[h][ij][1];
            for(int ab = 0; ab < G.params->coltot[h]; ++ab){
                int a = G.params->colorb[h][ab][0];
                int b = G.params->colorb[h][ab][1];
		int ha = G.params->rsym[a];
		int hb = G.params->ssym[b];
		int aa = a - G.params->roff[ha] + occpiA[ha];
		int bb = b - G.params->soff[hb] + occpiA[hb];
                if (i == j && ha == hb) G.matrix[h][ij][ab] -= 0.125 * gamma1corr->get(ha, aa, bb);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

    //Print
    if (print_ > 3) {
      global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
      global_dpd_->buf4_print(&G, "outfile", 1);
      global_dpd_->buf4_close(&G);
    }

    // TPDM <OV|OV>
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
      for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        if (wfn_type_ != "OMP2") global_dpd_->buf4_mat_irrep_rd(&G, h);
        #pragma omp parallel for
        for(int ia = 0; ia < G.params->rowtot[h]; ++ia){
            int i = G.params->roworb[h][ia][0];
            int a = G.params->roworb[h][ia][1];
            for(int jb = 0; jb < G.params->coltot[h]; ++jb){
                int j = G.params->colorb[h][jb][0];
                int b = G.params->colorb[h][jb][1];

		int ha = G.params->qsym[a];
		int hb = G.params->ssym[b];
		int aa = a - G.params->qoff[ha] + occpiA[ha];
		int bb = b - G.params->soff[hb] + occpiA[hb];

		if (i == j && ha == hb) {
                    if (wfn_type_ == "OMP2") G.matrix[h][ia][jb] = 0.5 * gamma1corr->get(ha,aa,bb);
                    else G.matrix[h][ia][jb] += 0.5 * gamma1corr->get(ha,aa,bb);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

    //Print
    if (print_ > 3) {
      global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
      global_dpd_->buf4_print(&G, "outfile", 1);
      global_dpd_->buf4_close(&G);
    }

 }// end if (reference_ == "RESTRICTED")

 else if (reference_ == "UNRESTRICTED") {
    // TPDM OOOO-Block
    // Alpha-Alpha spin case
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
      for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        #pragma omp parallel for
        for(int ij = 0; ij < G.params->rowtot[h]; ++ij){
            int i = G.params->roworb[h][ij][0];
            int j = G.params->roworb[h][ij][1];
            for(int kl = 0; kl < G.params->coltot[h]; ++kl){
                int k = G.params->colorb[h][kl][0];
                int l = G.params->colorb[h][kl][1];

		int hi = G.params->psym[i];
		int hj = G.params->qsym[j];
		int hk = G.params->rsym[k];
		int hl = G.params->ssym[l];

		int ii = i - G.params->poff[hi];
		int jj = j - G.params->qoff[hj];
		int kk = k - G.params->roff[hk];
		int ll = l - G.params->soff[hl];

		if (i == k && hj == hl) G.matrix[h][ij][kl] += 0.25 * gamma1corrA->get(hj,jj,ll);
		if (j == l && hi == hk) G.matrix[h][ij][kl] += 0.25 * gamma1corrA->get(hi,ii,kk);
		if (i == l && hj == hk) G.matrix[h][ij][kl] -= 0.25 * gamma1corrA->get(hj,jj,kk);
		if (j == k && hi == hl) G.matrix[h][ij][kl] -= 0.25 * gamma1corrA->get(hi,ii,ll);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);


    // Beta-Beta spin case
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "TPDM <oo|oo>");
      for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        #pragma omp parallel for
        for(int ij = 0; ij < G.params->rowtot[h]; ++ij){
            int i = G.params->roworb[h][ij][0];
            int j = G.params->roworb[h][ij][1];
            for(int kl = 0; kl < G.params->coltot[h]; ++kl){
                int k = G.params->colorb[h][kl][0];
                int l = G.params->colorb[h][kl][1];

		int hi = G.params->psym[i];
		int hj = G.params->qsym[j];
		int hk = G.params->rsym[k];
		int hl = G.params->ssym[l];

		int ii = i - G.params->poff[hi];
		int jj = j - G.params->qoff[hj];
		int kk = k - G.params->roff[hk];
		int ll = l - G.params->soff[hl];

		if (i == k && hj == hl) G.matrix[h][ij][kl] += 0.25 * gamma1corrB->get(hj,jj,ll);
		if (j == l && hi == hk) G.matrix[h][ij][kl] += 0.25 * gamma1corrB->get(hi,ii,kk);
		if (i == l && hj == hk) G.matrix[h][ij][kl] -= 0.25 * gamma1corrB->get(hj,jj,kk);
		if (j == k && hi == hl) G.matrix[h][ij][kl] -= 0.25 * gamma1corrB->get(hi,ii,ll);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);


    // Alpha-Beta spin case
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                 ID("[O,o]"), ID("[O,o]"), 0, "TPDM <Oo|Oo>");
      for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        global_dpd_->buf4_mat_irrep_rd(&G, h);
        #pragma omp parallel for
        for(int ij = 0; ij < G.params->rowtot[h]; ++ij){
            int i = G.params->roworb[h][ij][0];
            int j = G.params->roworb[h][ij][1];
            for(int kl = 0; kl < G.params->coltot[h]; ++kl){
                int k = G.params->colorb[h][kl][0];
                int l = G.params->colorb[h][kl][1];

		int hi = G.params->psym[i];
		int hj = G.params->qsym[j];
		int hk = G.params->rsym[k];
		int hl = G.params->ssym[l];

		int ii = i - G.params->poff[hi];
		int jj = j - G.params->qoff[hj];
		int kk = k - G.params->roff[hk];
		int ll = l - G.params->soff[hl];

		if (i == k && hj == hl) G.matrix[h][ij][kl] += 0.25 * gamma1corrB->get(hj,jj,ll);
		if (j == l && hi == hk) G.matrix[h][ij][kl] += 0.25 * gamma1corrA->get(hi,ii,kk);
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

    //Print
    if (print_ > 3) {
      global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
      global_dpd_->buf4_print(&G, "outfile", 1);
      global_dpd_->buf4_close(&G);

      global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "TPDM <oo|oo>");
      global_dpd_->buf4_print(&G, "outfile", 1);
      global_dpd_->buf4_close(&G);

      global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "TPDM <Oo|Oo>");
      global_dpd_->buf4_print(&G, "outfile", 1);
      global_dpd_->buf4_close(&G);
    }

    // TPDM <OV|OV>
    // Alpha-Alpha spin case
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
      for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        if (wfn_type_ != "OMP2") global_dpd_->buf4_mat_irrep_rd(&G, h);
        #pragma omp parallel for
        for(int ia = 0; ia < G.params->rowtot[h]; ++ia){
            int i = G.params->roworb[h][ia][0];
            int a = G.params->roworb[h][ia][1];
            for(int jb = 0; jb < G.params->coltot[h]; ++jb){
                int j = G.params->colorb[h][jb][0];
                int b = G.params->colorb[h][jb][1];
		int ha = G.params->qsym[a];
		int hb = G.params->ssym[b];
		int aa = a - G.params->qoff[ha] + occpiA[ha];
		int bb = b - G.params->soff[hb] + occpiA[hb];
		if (i == j && ha == hb) {
                    if (wfn_type_ == "OMP2") G.matrix[h][ia][jb] = 0.25 * gamma1corrA->get(ha,aa,bb);
                    else G.matrix[h][ia][jb] += 0.25 * gamma1corrA->get(ha,aa,bb);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);


    // Beta-Beta spin case
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "TPDM <ov|ov>");
      for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        if (wfn_type_ != "OMP2") global_dpd_->buf4_mat_irrep_rd(&G, h);
        #pragma omp parallel for
        for(int ia = 0; ia < G.params->rowtot[h]; ++ia){
            int i = G.params->roworb[h][ia][0];
            int a = G.params->roworb[h][ia][1];
            for(int jb = 0; jb < G.params->coltot[h]; ++jb){
                int j = G.params->colorb[h][jb][0];
                int b = G.params->colorb[h][jb][1];
		int ha = G.params->qsym[a];
		int hb = G.params->ssym[b];
		int aa = a - G.params->qoff[ha] + occpiB[ha];
		int bb = b - G.params->soff[hb] + occpiB[hb];
		if (i == j && ha == hb) {
                    if (wfn_type_ == "OMP2") G.matrix[h][ia][jb] = 0.25 * gamma1corrB->get(ha,aa,bb);
                    else G.matrix[h][ia][jb] += 0.25 * gamma1corrB->get(ha,aa,bb);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);


    // Alpha-Beta spin case
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "TPDM <Ov|Ov>");
      for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
        if (wfn_type_ != "OMP2") global_dpd_->buf4_mat_irrep_rd(&G, h);
        #pragma omp parallel for
        for(int ia = 0; ia < G.params->rowtot[h]; ++ia){
            int i = G.params->roworb[h][ia][0];
            int a = G.params->roworb[h][ia][1];
            for(int jb = 0; jb < G.params->coltot[h]; ++jb){
                int j = G.params->colorb[h][jb][0];
                int b = G.params->colorb[h][jb][1];
		int ha = G.params->qsym[a];
		int hb = G.params->ssym[b];
		int aa = a - G.params->qoff[ha] + occpiB[ha];
		int bb = b - G.params->soff[hb] + occpiB[hb];
		if (i == j && ha == hb) {
                    if (wfn_type_ == "OMP2") G.matrix[h][ia][jb] = 0.25 * gamma1corrB->get(ha,aa,bb);
                    else G.matrix[h][ia][jb] += 0.25 * gamma1corrB->get(ha,aa,bb);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

    //Print
    if (print_ > 3) {
      global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
      global_dpd_->buf4_print(&G, "outfile", 1);
      global_dpd_->buf4_close(&G);

      global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "TPDM <ov|ov>");
      global_dpd_->buf4_print(&G, "outfile", 1);
      global_dpd_->buf4_close(&G);

      global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "TPDM <Ov|Ov>");
      global_dpd_->buf4_print(&G, "outfile", 1);
      global_dpd_->buf4_close(&G);
    }


    // TPDM <Vo|Vo>
    // Alpha-Beta spin case
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,o]"), ID("[V,o]"),
                  ID("[V,o]"), ID("[V,o]"), 0, "TPDM <Vo|Vo>");
      for(int h = 0; h < nirrep_; ++h){
        global_dpd_->buf4_mat_irrep_init(&G, h);
	if (wfn_type_ != "OMP2") global_dpd_->buf4_mat_irrep_rd(&G, h);
        #pragma omp parallel for
        for(int ai = 0; ai < G.params->rowtot[h]; ++ai){
            int a = G.params->roworb[h][ai][0];
            int i = G.params->roworb[h][ai][1];
            for(int bj = 0; bj < G.params->coltot[h]; ++bj){
                int b = G.params->colorb[h][bj][0];
                int j = G.params->colorb[h][bj][1];
		int ha = G.params->psym[a];
		int hb = G.params->rsym[b];
		int aa = a - G.params->poff[ha] + occpiA[ha];
		int bb = b - G.params->roff[hb] + occpiA[hb];
		if (i == j && ha == hb) {
                    if (wfn_type_ == "OMP2") G.matrix[h][ai][bj] = 0.25 * gamma1corrA->get(ha,aa,bb);
                    else G.matrix[h][ai][bj] += 0.25 * gamma1corrA->get(ha,aa,bb);
                }
            }
        }
        global_dpd_->buf4_mat_irrep_wrt(&G, h);
        global_dpd_->buf4_mat_irrep_close(&G, h);
    }
    global_dpd_->buf4_close(&G);

    //Print
    if (print_ > 3) {
      global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,o]"), ID("[V,o]"),
                  ID("[V,o]"), ID("[V,o]"), 0, "TPDM <Vo|Vo>");
      global_dpd_->buf4_print(&G, "outfile", 1);
      global_dpd_->buf4_close(&G);
    }


    // G_IABJ = -GIAJB so I do not need to OVVO block,
    //  however for proper contraction (to avoid construction of <OV||VV>)  in the GFock.cc I need it.
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
    global_dpd_->buf4_sort(&G, PSIF_OCC_DENSITY , pqsr, ID("[O,V]"), ID("[V,O]"), "TPDM <OV|VO>");
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[V,O]"),
                  ID("[O,V]"), ID("[V,O]"), 0, "TPDM <OV|VO>");
    global_dpd_->buf4_scm(&G, -1.0);
    global_dpd_->buf4_close(&G);

    /*
    dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[V,O]"),
                  ID("[O,V]"), ID("[V,O]"), 0, "TPDM <OV|VO>");
    dpd_buf4_print(&G, outfile, 1);
    dpd_buf4_close(&G);
    */

    // G_iabj = -Giajb so I do not need to ovvo block,
    //  however for proper contraction (to avoid construction of <ov||vv>)  in the GFock.cc I need it.
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "TPDM <ov|ov>");
    global_dpd_->buf4_sort(&G, PSIF_OCC_DENSITY , pqsr, ID("[o,v]"), ID("[v,o]"), "TPDM <ov|vo>");
    global_dpd_->buf4_close(&G);
    global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,v]"), ID("[v,o]"),
                  ID("[o,v]"), ID("[v,o]"), 0, "TPDM <ov|vo>");
    global_dpd_->buf4_scm(&G, -1.0);
    global_dpd_->buf4_close(&G);

    /*
    dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,v]"), ID("[v,o]"),
                  ID("[o,v]"), ID("[v,o]"), 0, "TPDM <ov|vo>");
    dpd_buf4_print(&G, outfile, 1);
    dpd_buf4_close(&G);
    */

 }// end if (reference_ == "UNRESTRICTED")
    psio_->close(PSIF_OCC_DENSITY, 1);

}
}} // End Namespaces
