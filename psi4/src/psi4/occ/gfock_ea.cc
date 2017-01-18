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

#include "psi4/libiwl/iwl.hpp"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libmints/matrix.h"
#include "occwave.h"
#include "defines.h"
#include "dpd.h"


using namespace psi;
using namespace std;

namespace psi{ namespace occwave{

void OCCWave::gfock_ea()
{

//outfile->Printf("\n gfock_ea is starting... \n");
//===========================================================================================
//========================= RHF =============================================================
//===========================================================================================
if (reference_ == "RESTRICTED") {

	// Initialize
	Ftilde = std::shared_ptr<Matrix>(new Matrix("MO-basis GFM-EA", nirrep_, nmopi_, nmopi_));
	Ftilde->zero();
	Ftilde->add(HmoA);
        Ftilde->scale(2.0);
	Ftilde->subtract(GFock);

        // 2e-part
	dpdbuf4 G, K, X, T, Y;
	dpdfile2 GF, F, G1;

	psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
	psio_->open(PSIF_OCC_DENSITY, PSIO_OPEN_OLD);

/************************************************************************************************/
/*********************************** Write G1PDM to DPD files ***********************************/
/************************************************************************************************/
    // OO Block
    global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('O'), ID('O'), "OPDM <O|O>");
    global_dpd_->file2_mat_init(&G1);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < occpiA[h]; ++i){
            for(int j = 0 ; j < occpiA[h]; ++j){
		G1.matrix[h][i][j] = g1symm->get(h, i, j);
            }
        }
    }
    global_dpd_->file2_mat_wrt(&G1);
    global_dpd_->file2_close(&G1);

    // VV Block
    global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "OPDM <V|V>");
    global_dpd_->file2_mat_init(&G1);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < virtpiA[h]; ++i){
            for(int j = 0 ; j < virtpiA[h]; ++j){
                G1.matrix[h][i][j] = g1symm->get(h, i + occpiA[h], j + occpiA[h]);
            }
        }
    }
    global_dpd_->file2_mat_wrt(&G1);
    global_dpd_->file2_close(&G1);

    // OV Block
    global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('O'), ID('V'), "OPDM <O|V>");
    global_dpd_->file2_mat_init(&G1);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < occpiA[h]; ++i){
            for(int j = 0 ; j < virtpiA[h]; ++j){
		G1.matrix[h][i][j] = g1symm->get(h, i, j + occpiA[h]);
            }
        }
    }
    global_dpd_->file2_mat_wrt(&G1);
    global_dpd_->file2_close(&G1);

    // VO Block
    global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('V'), ID('O'), "OPDM <V|O>");
    global_dpd_->file2_mat_init(&G1);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < virtpiA[h]; ++i){
            for(int j = 0 ; j < occpiA[h]; ++j){
                G1.matrix[h][i][j] = g1symm->get(h, i + occpiA[h], j);
            }
        }
    }
    global_dpd_->file2_mat_wrt(&G1);
    global_dpd_->file2_close(&G1);

/************************************************************************************************/
/*********************************** Build Fij **************************************************/
/************************************************************************************************/
	// Build Fij
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('O'), ID('O'), "Ftilde <O|O>");

        // Fij += 2* \sum{m,n} (ij|mn) * G_mn
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O>=O]+"), ID("[O>=O]+"), 0, "MO Ints (OO|OO)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('O'), ID('O'), "OPDM <O|O>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 2.0, 0.0);
	global_dpd_->buf4_close(&K);

        // Fij += -1 * \sum{m,n} <ij|mn> * G_mn
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,O]"), ints->DPD_ID("[O,O]"),
                  ints->DPD_ID("[O,O]"), ints->DPD_ID("[O,O]"), 0, "MO Ints <OO|OO>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, -1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

	// Fij += 4 * \sum{m,e} (ij|me) * G_me
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O>=O]+"), ID("[O,V]"), 0, "MO Ints (OO|OV)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('O'), ID('V'), "OPDM <O|V>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 4.0, 1.0);
	global_dpd_->buf4_close(&K);

	// Fij += -1 * \sum{m,e} <ij|me> * G_me
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,O]"), ints->DPD_ID("[O,V]"),
                  ints->DPD_ID("[O,O]"), ints->DPD_ID("[O,V]"), 0, "MO Ints <OO|OV>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, -1.0, 1.0);
	// Fij += -1 * \sum{m,e} <ji|me> * G_me
	global_dpd_->contract422(&K, &G1, &GF, 0, 1, -1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

        // Fij += 2 * \sum{e,f} (ij|ef) * G_ef
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>=O]+"), ID("[V>=V]+"), 0, "MO Ints (OO|VV)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "OPDM <V|V>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&K);

        // Fij += -1 * \sum{e,f} <ij|ef> * G_ef
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO|VV>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, -1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

        // Close GF File
        global_dpd_->file2_close(&GF);

/************************************************************************************************/
/*********************************** Build Fia **************************************************/
/************************************************************************************************/
	// Build Fia
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('O'), ID('V'), "Ftilde <O|V>");

        // Fia += 2* \sum{m,n} (ia|mn) * G_mn
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O>=O]+"), ID("[O,V]"), 0, "MO Ints (OO|OV)");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , rspq, ID("[O,V]"), ID("[O,O]"), "MO Ints (OV|OO)");
        global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,O]"),
                  ID("[O,V]"), ID("[O,O]"), 0, "MO Ints (OV|OO)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('O'), ID('O'), "OPDM <O|O>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 2.0, 0.0);
	global_dpd_->buf4_close(&K);

	// Fia += -1 * \sum{m,n} <ia|mn> * G_mn
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,O]"), ints->DPD_ID("[O,V]"),
                  ints->DPD_ID("[O,O]"), ints->DPD_ID("[O,V]"), 0, "MO Ints <OO|OV>");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , rspq, ID("[O,V]"), ID("[O,O]"), "MO Ints <OV|OO>");
        global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,O]"),
                  ID("[O,V]"), ID("[O,O]"), 0, "MO Ints <OV|OO>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, -1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

	// Fia += 4 * \sum{m,e} (ia|me) * G_me
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV|OV)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('O'), ID('V'), "OPDM <O|V>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 4.0, 1.0);
	global_dpd_->buf4_close(&K);

	// Fia += -1 * \sum{m,e} <ia|me> * G_me
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,V]"), ints->DPD_ID("[O,V]"),
                  ints->DPD_ID("[O,V]"), ints->DPD_ID("[O,V]"), 0, "MO Ints <OV|OV>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, -1.0, 1.0);
	// Fia += -1 * \sum{m,e} <ai|me> * G_me
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , qprs, ID("[V,O]"), ID("[O,V]"), "MO Ints <VO|OV>");
        global_dpd_->buf4_close(&K);
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[V,O]"), ints->DPD_ID("[V,O]"),
                  ints->DPD_ID("[V,O]"), ints->DPD_ID("[O,V]"), 0, "MO Ints <VO|OV>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 1, -1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

        // Fia += 2 * \sum{e,f} (ia|ef) * G_ef
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V>=V]+"), 0, "MO Ints (OV|VV)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "OPDM <V|V>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&K);

        // Fia += -1 * \sum{e,f} <ia|ef> * G_ef
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV|VV>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, -1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

        // Close GF File
        global_dpd_->file2_close(&GF);

/************************************************************************************************/
/*********************************** Build Fab **************************************************/
/************************************************************************************************/
	// Build Fab
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "Ftilde <V|V>");

        // Fab += 2* \sum{m,n} (ab|mn) * G_mn
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O>=O]+"), ID("[V>=V]+"), 0, "MO Ints (OO|VV)");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , rspq, ID("[V,V]"), ID("[O,O]"), "MO Ints (VV|OO)");
        global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[O,O]"),
                  ID("[V,V]"), ID("[O,O]"), 0, "MO Ints (VV|OO)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('O'), ID('O'), "OPDM <O|O>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 2.0, 0.0);
	global_dpd_->buf4_close(&K);

	// Fab += -1 * \sum{m,n} <ab|mn> * G_mn
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,O]"), ints->DPD_ID("[V,V]"),
                  ints->DPD_ID("[O,O]"), ints->DPD_ID("[V,V]"), 0, "MO Ints <OO|VV>");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , rspq, ID("[V,V]"), ID("[O,O]"), "MO Ints <VV|OO>");
        global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[O,O]"),
                  ID("[V,V]"), ID("[O,O]"), 0, "MO Ints <VV|OO>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, -1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

	// Fab += 4 * \sum{m,e} (ab|me) * G_me
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V>=V]+"), 0, "MO Ints (OV|VV)");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , rspq, ID("[V,V]"), ID("[O,V]"), "MO Ints (VV|OV)");
        global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[O,V]"),
                  ID("[V,V]"), ID("[O,V]"), 0, "MO Ints (VV|OV)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('O'), ID('V'), "OPDM <O|V>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 4.0, 1.0);
	global_dpd_->buf4_close(&K);

	// Fab += -1 * \sum{m,e} <ab|me> * G_me
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,V]"), ints->DPD_ID("[V,V]"),
                  ints->DPD_ID("[O,V]"), ints->DPD_ID("[V,V]"), 0, "MO Ints <OV|VV>");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , rspq, ID("[V,V]"), ID("[O,V]"), "MO Ints <VV|OV>");
        global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[O,V]"),
                  ID("[V,V]"), ID("[O,V]"), 0, "MO Ints <VV|OV>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, -1.0, 1.0);
	// Fab += -1 * \sum{m,e} <ba|me> * G_me
	global_dpd_->contract422(&K, &G1, &GF, 0, 1, -1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

        // Fab += 2 * \sum{e,f} (ab|ef) * G_ef
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                 ID("[V>=V]+"), ID("[V>=V]+"), 0, "MO Ints (VV|VV)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "OPDM <V|V>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&K);

        // Fab += -1 * \sum{e,f} <ab|ef> * G_ef
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "MO Ints <VV|VV>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, -1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

        // Close GF File
        global_dpd_->file2_close(&GF);

        // Close Integral DPD Files
	psio_->close(PSIF_LIBTRANS_DPD, 1);

/************************************************************************************************/
/*********************************** Load *******************************************************/
/************************************************************************************************/
	// Load Fij
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('O'), ID('O'), "Ftilde <O|O>");
	global_dpd_->file2_mat_init(&GF);
	global_dpd_->file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < occpiA[h]; ++i){
            for(int j = 0 ; j < occpiA[h]; ++j){
                Ftilde->add(h, i, j, GF.matrix[h][i][j]);
            }
	  }
	}
	global_dpd_->file2_close(&GF);

	// Load Fia & Fai
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('O'), ID('V'), "Ftilde <O|V>");
	global_dpd_->file2_mat_init(&GF);
	global_dpd_->file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < occpiA[h]; ++i){
            for(int a = 0 ; a < virtpiA[h]; ++a){
                Ftilde->add(h, i, a + occpiA[h], GF.matrix[h][i][a]);
                Ftilde->add(h, a + occpiA[h], i, GF.matrix[h][i][a]);
            }
	  }
	}
	global_dpd_->file2_close(&GF);

	// Load Fab
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "Ftilde <V|V>");
	global_dpd_->file2_mat_init(&GF);
	global_dpd_->file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
          for(int a = 0 ; a < virtpiA[h]; ++a){
            for(int b = 0 ; b < virtpiA[h]; ++b){
                Ftilde->add(h, a + occpiA[h], b + occpiA[h], GF.matrix[h][a][b]);
            }
	  }
	}
	global_dpd_->file2_close(&GF);


	psio_->close(PSIF_OCC_DENSITY, 1);
	if (print_ > 2) Ftilde->print();

}// end if (reference_ == "RESTRICTED")



//===========================================================================================
//========================= UHF =============================================================
//===========================================================================================
else if (reference_ == "UNRESTRICTED") {

	// Initialize
	FtildeA = std::shared_ptr<Matrix>(new Matrix("MO-basis Alpha GFM-EA", nirrep_, nmopi_, nmopi_));
	FtildeB = std::shared_ptr<Matrix>(new Matrix("MO-basis Beta GFM-EA", nirrep_, nmopi_, nmopi_));
	FtildeA->zero();
	FtildeB->zero();
	FtildeA->add(HmoA);
	FtildeB->add(HmoB);
	FtildeA->subtract(GFockA);
	FtildeB->subtract(GFockB);

        // 2e-part
	dpdbuf4 G, K, X, T, Y;
	dpdfile2 GF, F, G1;

	psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
	psio_->open(PSIF_OCC_DENSITY, PSIO_OPEN_OLD);

/************************************************************************************************/
/*********************************** Write G1PDM to DPD files ***********************************/
/************************************************************************************************/
    // OO Block
    // AA-Spin
    global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('O'), ID('O'), "OPDM <O|O>");
    global_dpd_->file2_mat_init(&G1);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < occpiA[h]; ++i){
            for(int j = 0 ; j < occpiA[h]; ++j){
		G1.matrix[h][i][j] = g1symmA->get(h, i, j);
            }
        }
    }
    global_dpd_->file2_mat_wrt(&G1);
    global_dpd_->file2_close(&G1);

    // BB-Spin
    global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('o'), ID('o'), "OPDM <o|o>");
    global_dpd_->file2_mat_init(&G1);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < occpiB[h]; ++i){
            for(int j = 0 ; j < occpiB[h]; ++j){
		G1.matrix[h][i][j] = g1symmB->get(h, i, j);
            }
        }
    }
    global_dpd_->file2_mat_wrt(&G1);
    global_dpd_->file2_close(&G1);


    // VV Block
    // AA-Spin
    global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "OPDM <V|V>");
    global_dpd_->file2_mat_init(&G1);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < virtpiA[h]; ++i){
            for(int j = 0 ; j < virtpiA[h]; ++j){
                G1.matrix[h][i][j] = g1symmA->get(h, i + occpiA[h], j + occpiA[h]);
            }
        }
    }
    global_dpd_->file2_mat_wrt(&G1);
    global_dpd_->file2_close(&G1);

    // BB-Spin
    global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('v'), ID('v'), "OPDM <v|v>");
    global_dpd_->file2_mat_init(&G1);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < virtpiB[h]; ++i){
            for(int j = 0 ; j < virtpiB[h]; ++j){
                G1.matrix[h][i][j] = g1symmB->get(h, i + occpiB[h], j + occpiB[h]);
            }
        }
    }
    global_dpd_->file2_mat_wrt(&G1);
    global_dpd_->file2_close(&G1);


    // OV Block
    // AA-Spin
    global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('O'), ID('V'), "OPDM <O|V>");
    global_dpd_->file2_mat_init(&G1);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < occpiA[h]; ++i){
            for(int j = 0 ; j < virtpiA[h]; ++j){
		G1.matrix[h][i][j] = g1symmA->get(h, i, j + occpiA[h]);
            }
        }
    }
    global_dpd_->file2_mat_wrt(&G1);
    global_dpd_->file2_close(&G1);

    // BB-Spin
    global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('o'), ID('v'), "OPDM <o|v>");
    global_dpd_->file2_mat_init(&G1);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < occpiB[h]; ++i){
            for(int j = 0 ; j < virtpiB[h]; ++j){
		G1.matrix[h][i][j] = g1symmB->get(h, i, j + occpiB[h]);
            }
        }
    }
    global_dpd_->file2_mat_wrt(&G1);
    global_dpd_->file2_close(&G1);


    // VO Block
    // AA-Spin
    global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('V'), ID('O'), "OPDM <V|O>");
    global_dpd_->file2_mat_init(&G1);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < virtpiA[h]; ++i){
            for(int j = 0 ; j < occpiA[h]; ++j){
                G1.matrix[h][i][j] = g1symmA->get(h, i + occpiA[h], j);
            }
        }
    }
    global_dpd_->file2_mat_wrt(&G1);
    global_dpd_->file2_close(&G1);

    // BB-Spin
    global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('v'), ID('o'), "OPDM <v|o>");
    global_dpd_->file2_mat_init(&G1);
    for(int h = 0; h < nirrep_; ++h){
        for(int i = 0 ; i < virtpiB[h]; ++i){
            for(int j = 0 ; j < occpiB[h]; ++j){
                G1.matrix[h][i][j] = g1symmB->get(h, i + occpiB[h], j);
            }
        }
    }
    global_dpd_->file2_mat_wrt(&G1);
    global_dpd_->file2_close(&G1);

/************************************************************************************************/
/*********************************** Build FIJ **************************************************/
/************************************************************************************************/
	// Build FIJ
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('O'), ID('O'), "Ftilde <O|O>");

        // FIJ += \sum{M,N} (IJ||MN) * G_MN
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "MO Ints <OO||OO>");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,O]"), ID("[O,O]"), "MO Ints (OO||OO)");
	global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "MO Ints (OO||OO)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('O'), ID('O'), "OPDM <O|O>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 0.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

	// FIJ += \sum{M,E} (IJ||ME) * G_ME
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "MO Ints <OO||OV>");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,O]"), ID("[O,V]"), "MO Ints (OO||OV)");
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "MO Ints (OO||OV)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('O'), ID('V'), "OPDM <O|V>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

	// FIJ += \sum{E,M} (IJ||EM) * G_EM
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "MO Ints (OO||OV)");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , pqsr, ID("[O,O]"), ID("[V,O]"), "MO Ints (OO||VO)");
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,O]"),
                  ID("[O,O]"), ID("[V,O]"), 0, "MO Ints (OO||VO)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('V'), ID('O'), "OPDM <V|O>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

        // FIJ += \sum{E,F} (IJ||EF) * G_EF
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV||OV>");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,O]"), ID("[V,V]"), "MO Ints (OO||VV)");
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints (OO||VV)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "OPDM <V|V>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

        // FIJ += \sum{m,n} (IJ|mn) * G_mn
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[o,o]"),
                  ID("[O>=O]+"), ID("[o>=o]+"), 0, "MO Ints (OO|oo)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('o'), ID('o'), "OPDM <o|o>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

	// FIJ += 2 * \sum{m,e} (IJ|me) * G_me
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[o,v]"),
                  ID("[O>=O]+"), ID("[o,v]"), 0, "MO Ints (OO|ov)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('o'), ID('v'), "OPDM <o|v>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

        // FIJ += \sum{e,f} (IJ|ef) * G_ef
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[v,v]"),
                  ID("[O>=O]+"), ID("[v>=v]+"), 0, "MO Ints (OO|vv)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('v'), ID('v'), "OPDM <v|v>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

        // Close GF File
        global_dpd_->file2_close(&GF);

/************************************************************************************************/
/*********************************** Build Fij **************************************************/
/************************************************************************************************/
	// Build Fij
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('o'), ID('o'), "Ftilde <o|o>");

        // Fij += \sum{m,n} (ij||mn) * G_mn
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "MO Ints <oo||oo>");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[o,o]"), ID("[o,o]"), "MO Ints (oo||oo)");
	global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "MO Ints (oo||oo)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('o'), ID('o'), "OPDM <o|o>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 0.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

	// Fij += \sum{m,e} (ij||me) * G_me
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,v]"),
                  ID("[o,o]"), ID("[o,v]"), 0, "MO Ints <oo||ov>");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[o,o]"), ID("[o,v]"), "MO Ints (oo||ov)");
        global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,v]"),
                  ID("[o,o]"), ID("[o,v]"), 0, "MO Ints (oo||ov)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('o'), ID('v'), "OPDM <o|v>");
        global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 1.0);
        global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

	// Fij += \sum{e,m} (ij||em) * G_em
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,v]"),
                  ID("[o,o]"), ID("[o,v]"), 0, "MO Ints (oo||ov)");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , pqsr, ID("[o,o]"), ID("[v,o]"), "MO Ints (oo||vo)");
        global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,o]"),
                  ID("[o,o]"), ID("[v,o]"), 0, "MO Ints (oo||vo)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('v'), ID('o'), "OPDM <v|o>");
        global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 1.0);
        global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

        // Fij += \sum{e,f} (ij||ef) * G_ef
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov||ov>");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[o,o]"), ID("[v,v]"), "MO Ints (oo||vv)");
        global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints (oo||vv)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('v'), ID('v'), "OPDM <v|v>");
        global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 1.0);
        global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

        // Fij += \sum{M,N} (MN|ij) * G_MN
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[o,o]"),
                  ID("[O>=O]+"), ID("[o>=o]+"), 0, "MO Ints (OO|oo)");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , rspq, ID("[o,o]"), ID("[O,O]"), "MO Ints (oo|OO)");
	global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[O,O]"),
                  ID("[o,o]"), ID("[O,O]"), 0, "MO Ints (oo|OO)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('O'), ID('O'), "OPDM <O|O>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

	// Fij += 2 * \sum{M,E} (ME|ij) * G_ME
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,o]"),
                  ID("[O,V]"), ID("[o>=o]+"), 0, "MO Ints (OV|oo)");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , rspq, ID("[o,o]"), ID("[O,V]"), "MO Ints (oo|OV)");
	global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[O,V]"),
                  ID("[o,o]"), ID("[O,V]"), 0, "MO Ints (oo|OV)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('O'), ID('V'), "OPDM <O|V>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

        // Fij += \sum{E,F} (EF|ij) * G_EF
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[o,o]"),
                  ID("[V>=V]+"), ID("[o>=o]+"), 0, "MO Ints (VV|oo)");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , rspq, ID("[o,o]"), ID("[V,V]"), "MO Ints (oo|VV)");
	global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[V,V]"),
                  ID("[o,o]"), ID("[V,V]"), 0, "MO Ints (oo|VV)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "OPDM <V|V>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

        // Close GF File
        global_dpd_->file2_close(&GF);

/************************************************************************************************/
/*********************************** Build FIA **************************************************/
/************************************************************************************************/
	// Build FIA
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('O'), ID('V'), "Ftilde <O|V>");

        // FIA += \sum{M,N} (IA||MN) * G_MN
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "MO Ints <OO||OV>");
        // (IA||MN) = <MI||NA>
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , qspr, ID("[O,V]"), ID("[O,O]"), "MO Ints (OV||OO)");
	global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,O]"),
                  ID("[O,V]"), ID("[O,O]"), 0, "MO Ints (OV||OO)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('O'), ID('O'), "OPDM <O|O>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 0.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

	// FIA += \sum{M,E} (IA||ME) * G_ME
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[O,V]"), ID("[O,V]"), "MO Ints (OV||OV)");
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV||OV)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('O'), ID('V'), "OPDM <O|V>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

	// FIA += \sum{E,M} (IA||EM) * G_EM
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints (OV||OV)");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , pqsr, ID("[O,V]"), ID("[V,O]"), "MO Ints (OV||VO)");
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,O]"),
                  ID("[O,V]"), ID("[V,O]"), 0, "MO Ints (OV||VO)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('V'), ID('O'), "OPDM <V|O>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

        // FIA += \sum{E,F} (IA||EF) * G_EF = \sum{E,F} [(IA|EF)-<IA|EF>] * G_EF
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V>=V]+"), 0, "MO Ints (OV|VV)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "OPDM <V|V>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV|VV>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, -1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

        // FIA += \sum{m,n} (IA|mn) * G_mn
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,o]"),
                  ID("[O,V]"), ID("[o>=o]+"), 0, "MO Ints (OV|oo)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('o'), ID('o'), "OPDM <o|o>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

	// FIA += 2 * \sum{m,e} (IA|me) * G_me
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('o'), ID('v'), "OPDM <o|v>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

        // FIA += \sum{e,f} (IA|ef) * G_ef
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[v,v]"),
                  ID("[O,V]"), ID("[v>=v]+"), 0, "MO Ints (OV|vv)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('v'), ID('v'), "OPDM <v|v>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

        // Close GF File
        global_dpd_->file2_close(&GF);

/************************************************************************************************/
/*********************************** Build Fia **************************************************/
/************************************************************************************************/
	// Build Fia
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('o'), ID('v'), "Ftilde <o|v>");

        // Fia += \sum{m,n} (ia||mn) * G_mn
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,v]"),
                  ID("[o,o]"), ID("[o,v]"), 0, "MO Ints <oo||ov>");
        // (ia||mn) = <mi||na>
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , qspr, ID("[o,v]"), ID("[o,o]"), "MO Ints (ov||oo)");
	global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,o]"),
                  ID("[o,v]"), ID("[o,o]"), 0, "MO Ints (ov||oo)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('o'), ID('o'), "OPDM <o|o>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 0.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

	// Fia += \sum{m,e} (ia||me) * G_me
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , prqs, ID("[o,v]"), ID("[o,v]"), "MO Ints (ov||ov)");
        global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints (ov||ov)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('o'), ID('v'), "OPDM <o|v>");
        global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 1.0);
        global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

	// Fia += \sum{e,m} (ia||em) * G_em
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints (ov||ov)");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , pqsr, ID("[o,v]"), ID("[v,o]"), "MO Ints (ov||vo)");
        global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,o]"),
                  ID("[o,v]"), ID("[v,o]"), 0, "MO Ints (ov||vo)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('v'), ID('o'), "OPDM <v|o>");
        global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 1.0);
        global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

        // Fia += \sum{e,f} (ia||ef) * G_ef = \sum{e,f} [(ia|ef)-<ia|ef>] * G_ef
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v>=v]+"), 0, "MO Ints (ov|vv)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('v'), ID('v'), "OPDM <v|v>");
        global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 1.0);
        global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <ov|vv>");
        global_dpd_->contract422(&K, &G1, &GF, 0, 0, -1.0, 1.0);
        global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

        // Fia += \sum{M,N} (MN|ia) * G_MN
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[o,v]"),
                  ID("[O>=O]+"), ID("[o,v]"), 0, "MO Ints (OO|ov)");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , rspq, ID("[o,v]"), ID("[O,O]"), "MO Ints (ov|OO)");
	global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[O,O]"),
                  ID("[o,v]"), ID("[O,O]"), 0, "MO Ints (ov|OO)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('O'), ID('O'), "OPDM <O|O>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

	// Fia += 2 * \sum{M,E} (ME|ia) * G_ME
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[o,v]"),
                  ID("[O,V]"), ID("[o,v]"), 0, "MO Ints (OV|ov)");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , rspq, ID("[o,v]"), ID("[O,V]"), "MO Ints (ov|OV)");
	global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[O,V]"),
                  ID("[o,v]"), ID("[O,V]"), 0, "MO Ints (ov|OV)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('O'), ID('V'), "OPDM <O|V>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

        // Fia += \sum{E,F} (EF|ia) * G_EF
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[o,v]"),
                  ID("[V>=V]+"), ID("[o,v]"), 0, "MO Ints (VV|ov)");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , rspq, ID("[o,v]"), ID("[V,V]"), "MO Ints (ov|VV)");
	global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[V,V]"),
                  ID("[o,v]"), ID("[V,V]"), 0, "MO Ints (ov|VV)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "OPDM <V|V>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

        // Close GF File
        global_dpd_->file2_close(&GF);

/************************************************************************************************/
/*********************************** Build FAB **************************************************/
/************************************************************************************************/
	// Build FAB
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "Ftilde <V|V>");

        // FAB += \sum{M,N} (AB||MN) * G_MN
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV||OV>");
        // (AB||MN) = <MA||NB>
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , qspr, ID("[V,V]"), ID("[O,O]"), "MO Ints (VV||OO)");
	global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[O,O]"),
                  ID("[V,V]"), ID("[O,O]"), 0, "MO Ints (VV||OO)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('O'), ID('O'), "OPDM <O|O>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 0.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

	// FAB += \sum{M,E} (AB||ME) * G_ME = \sum{M,E} [(AB|ME)-<AB|ME>] * G_ME
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[O,V]"),
                  ID("[V>=V]+"), ID("[O,V]"), 0, "MO Ints (VV|OV)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('O'), ID('V'), "OPDM <O|V>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV|VV>");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , rspq, ID("[V,V]"), ID("[O,V]"), "MO Ints <VV|OV>");
	global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[O,V]"),
                  ID("[V,V]"), ID("[O,V]"), 0, "MO Ints <VV|OV>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, -1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

	// FAB += \sum{E,M} (AB||EM) * G_EM = \sum{E,M} [(AB|EM)-<AB|EM>] * G_EM
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('V'), ID('O'), "OPDM <V|O>");
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[O,V]"),
                  ID("[V>=V]+"), ID("[O,V]"), 0, "MO Ints (VV|OV)");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , pqsr, ID("[V,V]"), ID("[V,O]"), "MO Ints (VV|VO)");
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,O]"),
                  ID("[V,V]"), ID("[V,O]"), 0, "MO Ints (VV|VO)");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[O,V]"),
                  ID("[V,V]"), ID("[O,V]"), 0, "MO Ints <VV|OV>");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , pqsr, ID("[V,V]"), ID("[V,O]"), "MO Ints <VV|VO>");
	global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,O]"),
                  ID("[V,V]"), ID("[V,O]"), 0, "MO Ints <VV|VO>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, -1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

        // FAB += \sum{E,F} (AB||EF) * G_EF = \sum{E,F} [(AB|EF)-<AB|EF>] * G_EF
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V>=V]+"), ID("[V>=V]+"), 0, "MO Ints (VV|VV)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "OPDM <V|V>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "MO Ints <VV|VV>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, -1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

        // FAB += \sum{m,n} (AB|mn) * G_mn
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[o,o]"),
                  ID("[V>=V]+"), ID("[o>=o]+"), 0, "MO Ints (VV|oo)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('o'), ID('o'), "OPDM <o|o>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

	// FAB += 2 * \sum{m,e} (AB|me) * G_me
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[o,v]"),
                  ID("[V>=V]+"), ID("[o,v]"), 0, "MO Ints (VV|ov)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('o'), ID('v'), "OPDM <o|v>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

        // FAB += \sum{e,f} (AB|ef) * G_ef
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[v,v]"),
                  ID("[V>=V]+"), ID("[v>=v]+"), 0, "MO Ints (VV|vv)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('v'), ID('v'), "OPDM <v|v>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

        // Close GF File
        global_dpd_->file2_close(&GF);

/************************************************************************************************/
/*********************************** Build Fab **************************************************/
/************************************************************************************************/
	// Build Fab
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('v'), ID('v'), "Ftilde <v|v>");

        // Fab += \sum{m,n} (ab||mn) * G_mn
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov||ov>");
        // (ab||mn) = <ma||nb>
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , qspr, ID("[v,v]"), ID("[o,o]"), "MO Ints (vv||oo)");
        global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[o,o]"),
                  ID("[v,v]"), ID("[o,o]"), 0, "MO Ints (vv||oo)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('o'), ID('o'), "OPDM <o|o>");
        global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 0.0);
        global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

	// Fab += \sum{m,e} (ab||me) * G_me = \sum{m,e} [(ab|me)-<ab|me>] * G_me
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[o,v]"),
                  ID("[v>=v]+"), ID("[o,v]"), 0, "MO Ints (vv|ov)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('o'), ID('v'), "OPDM <o|v>");
        global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 1.0);
        global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <ov|vv>");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , rspq, ID("[v,v]"), ID("[o,v]"), "MO Ints <vv|ov>");
        global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[o,v]"),
                  ID("[v,v]"), ID("[o,v]"), 0, "MO Ints <vv|ov>");
        global_dpd_->contract422(&K, &G1, &GF, 0, 0, -1.0, 1.0);
        global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

	// Fab += \sum{e,m} (ab||em) * G_em = \sum{e,m} [(ab|em)-<ab|em>] * G_em
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('v'), ID('o'), "OPDM <v|o>");
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[o,v]"),
                  ID("[v>=v]+"), ID("[o,v]"), 0, "MO Ints (vv|ov)");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , pqsr, ID("[v,v]"), ID("[v,o]"), "MO Ints (vv|vo)");
        global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,o]"),
                  ID("[v,v]"), ID("[v,o]"), 0, "MO Ints (vv|vo)");
        global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 1.0);
        global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[o,v]"),
                  ID("[v,v]"), ID("[o,v]"), 0, "MO Ints <vv|ov>");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , pqsr, ID("[v,v]"), ID("[v,o]"), "MO Ints <vv|vo>");
        global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,o]"),
                  ID("[v,v]"), ID("[v,o]"), 0, "MO Ints <vv|vo>");
        global_dpd_->contract422(&K, &G1, &GF, 0, 0, -1.0, 1.0);
        global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

        // Fab += \sum{e,f} (ab||ef) * G_ef = \sum{e,f} [(ab|ef)-<ab|ef>] * G_ef
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v>=v]+"), ID("[v>=v]+"), 0, "MO Ints (vv|vv)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('v'), ID('v'), "OPDM <v|v>");
        global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 1.0);
        global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "MO Ints <vv|vv>");
        global_dpd_->contract422(&K, &G1, &GF, 0, 0, -1.0, 1.0);
        global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

        // Fab += \sum{M,N} (ab|MN) * G_MN
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[v,v]"),
                  ID("[O>=O]+"), ID("[v>=v]+"), 0, "MO Ints (OO|vv)");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , rspq, ID("[v,v]"), ID("[O,O]"), "MO Ints (vv|OO)");
	global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[O,O]"),
                  ID("[v,v]"), ID("[O,O]"), 0, "MO Ints (vv|OO)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('O'), ID('O'), "OPDM <O|O>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

	// Fab += 2 * \sum{M,E} (ab|ME) * G_ME
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[v,v]"),
                  ID("[O,V]"), ID("[v>=v]+"), 0, "MO Ints (OV|vv)");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , rspq, ID("[v,v]"), ID("[O,V]"), "MO Ints (vv|OV)");
	global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[O,V]"),
                  ID("[v,v]"), ID("[O,V]"), 0, "MO Ints (vv|OV)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('O'), ID('V'), "OPDM <O|V>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

        // Fab += \sum{E,F} (ab|EF) * G_EF
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[v,v]"),
                  ID("[V>=V]+"), ID("[v>=v]+"), 0, "MO Ints (VV|vv)");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , rspq, ID("[v,v]"), ID("[V,V]"), "MO Ints (vv|VV)");
	global_dpd_->buf4_close(&K);
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[V,V]"),
                  ID("[v,v]"), ID("[V,V]"), 0, "MO Ints (vv|VV)");
        global_dpd_->file2_init(&G1, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "OPDM <V|V>");
	global_dpd_->contract422(&K, &G1, &GF, 0, 0, 1.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G1);

        // Close GF File
        global_dpd_->file2_close(&GF);

/************************************************************************************************/
/*********************************** Load *******************************************************/
/************************************************************************************************/
	// Load FIJ
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('O'), ID('O'), "Ftilde <O|O>");
	global_dpd_->file2_mat_init(&GF);
	global_dpd_->file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < occpiA[h]; ++i){
            for(int j = 0 ; j < occpiA[h]; ++j){
                FtildeA->add(h, i, j, GF.matrix[h][i][j]);
            }
	  }
	}
	global_dpd_->file2_close(&GF);

	// Load Fij
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('o'), ID('o'), "Ftilde <o|o>");
	global_dpd_->file2_mat_init(&GF);
	global_dpd_->file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < occpiB[h]; ++i){
            for(int j = 0 ; j < occpiB[h]; ++j){
                FtildeB->add(h, i, j, GF.matrix[h][i][j]);
            }
	  }
	}
	global_dpd_->file2_close(&GF);

	// Load FIA & FAI
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('O'), ID('V'), "Ftilde <O|V>");
	global_dpd_->file2_mat_init(&GF);
	global_dpd_->file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < occpiA[h]; ++i){
            for(int a = 0 ; a < virtpiA[h]; ++a){
                FtildeA->add(h, i, a + occpiA[h], GF.matrix[h][i][a]);
                FtildeA->add(h, a + occpiA[h], i, GF.matrix[h][i][a]);
            }
	  }
	}
	global_dpd_->file2_close(&GF);

	// Load Fia & Fai
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('o'), ID('v'), "Ftilde <o|v>");
	global_dpd_->file2_mat_init(&GF);
	global_dpd_->file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < occpiB[h]; ++i){
            for(int a = 0 ; a < virtpiB[h]; ++a){
                FtildeB->add(h, i, a + occpiB[h], GF.matrix[h][i][a]);
                FtildeB->add(h, a + occpiB[h], i, GF.matrix[h][i][a]);
            }
	  }
	}
	global_dpd_->file2_close(&GF);

	// Load Fab
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "Ftilde <V|V>");
	global_dpd_->file2_mat_init(&GF);
	global_dpd_->file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
          for(int a = 0 ; a < virtpiA[h]; ++a){
            for(int b = 0 ; b < virtpiA[h]; ++b){
                FtildeA->add(h, a + occpiA[h], b + occpiA[h], GF.matrix[h][a][b]);
            }
	  }
	}
	global_dpd_->file2_close(&GF);

	// Load Fab
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('v'), ID('v'), "Ftilde <v|v>");
	global_dpd_->file2_mat_init(&GF);
	global_dpd_->file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
          for(int a = 0 ; a < virtpiB[h]; ++a){
            for(int b = 0 ; b < virtpiB[h]; ++b){
                FtildeB->add(h, a + occpiB[h], b + occpiB[h], GF.matrix[h][a][b]);
            }
	  }
	}
	global_dpd_->file2_close(&GF);

	psio_->close(PSIF_OCC_DENSITY, 1);
	if (print_ > 2) {
            FtildeA->print();
            FtildeB->print();
        }

}// end if (reference_ == "UNRESTRICTED")
//outfile->Printf("\n gfock_ea done. \n");

} // End main
}} // End Namespaces
