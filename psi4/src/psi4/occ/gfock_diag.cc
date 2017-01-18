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

#include "psi4/libtrans/integraltransform.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsio/psio.hpp"
#include "occwave.h"
#include "defines.h"


using namespace std;

namespace psi{ namespace occwave{

void OCCWave::gfock_diag()
{
//outfile->Printf("\n gfock_diag is starting... \n");
//===========================================================================================
//========================= RHF =============================================================
//===========================================================================================
if (reference_ == "RESTRICTED") {
        // 2e-part
	dpdbuf4 G, K, X, T, Y;
	dpdfile2 GF;

	psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
	psio_->open(PSIF_OCC_DENSITY, PSIO_OPEN_OLD);
        psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);

	// Build Fij
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('O'), ID('O'), "GF <O|O>");

	// Fij += 4 * \sum{m,n,k} <km|ni> * G_kmnj
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,O]"), ints->DPD_ID("[O,O]"),
                  ints->DPD_ID("[O,O]"), ints->DPD_ID("[O,O]"), 0, "MO Ints <OO|OO>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
	global_dpd_->contract442(&K, &G, &GF, 3, 3, 4.0, 0.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

	// Fij += 4 * \sum{e,m,f} <me|if> * G_mejf
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,V]"), ints->DPD_ID("[O,V]"),
                  ints->DPD_ID("[O,V]"), ints->DPD_ID("[O,V]"), 0, "MO Ints <OV|OV>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
	global_dpd_->contract442(&K, &G, &GF, 2, 2, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

	// Fij += 8 * \sum{e,m,f} <mi|ef> * G_mjef
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,O]"), ints->DPD_ID("[V,V]"),
                  ints->DPD_ID("[O,O]"), ints->DPD_ID("[V,V]"), 0, "MO Ints <OO|VV>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
	global_dpd_->contract442(&K, &G, &GF, 1, 1, 8.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);
	global_dpd_->file2_close(&GF);

	// Build Fab
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "GF <V|V>");

// Form X intermediate
if (wfn_type_ != "OMP2") {
   // Build X intermediate
   if (twopdm_abcd_type == "DIRECT" ) {
        // With this algorithm cost changes to v5 => o2v4 + o2v3, roughly v/o times faster
 	// X_MNFA = 2\sum{E,C} [2t_MN^CE(1) - t_MN^EC(1)] * <FA|CE>
        global_dpd_->buf4_init(&X, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[O,O]"),
                  ID("[V,V]"), ID("[O,O]"), 0, "X <VV|OO>");
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "MO Ints <VV|VV>");
        if (wfn_type_ == "OMP3" || wfn_type_ == "OMP2.5") {
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Tau_1 <OO|VV>");
        }
        else if (wfn_type_ == "OCEPA") {
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Tau <OO|VV>");
        }
        global_dpd_->contract444(&K, &T, &X, 0, 0, 2.0, 0.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&T);
        global_dpd_->buf4_sort(&X, PSIF_OCC_DENSITY, rspq, ID("[O,O]"), ID("[V,V]"), "X <OO|VV>");
	global_dpd_->buf4_close(&X);

        // OMP2.5
        if (wfn_type_ == "OMP2.5") {
        global_dpd_->buf4_init(&X, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "X <OO|VV>");
        global_dpd_->buf4_scm(&X, 0.5);
	global_dpd_->buf4_close(&X);
        }
    }
}// end if (wfn_type_ != "OMP2")


	// Fab += 4 * \sum{m,e,n} <ma|ne> * G_mbne
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,V]"), ints->DPD_ID("[O,V]"),
                  ints->DPD_ID("[O,V]"), ints->DPD_ID("[O,V]"), 0, "MO Ints <OV|OV>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
	global_dpd_->contract442(&K, &G, &GF, 1, 1, 4.0, 0.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

	// Fab += 8 * \sum{m,e,n} <nm|ae> * G_nmbe
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,O]"), ints->DPD_ID("[V,V]"),
                  ints->DPD_ID("[O,O]"), ints->DPD_ID("[V,V]"), 0, "MO Ints <OO|VV>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
	global_dpd_->contract442(&K, &G, &GF, 2, 2, 8.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

if (wfn_type_ != "OMP2") {
   if (twopdm_abcd_type == "DIRECT" ) {
       	// Fab += \sum{m,n,f} X_mnfa * t_mn^fb(1)
        global_dpd_->buf4_init(&X, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "X <OO|VV>");
        if (wfn_type_ == "OMP3" || wfn_type_ == "OMP2.5") {
           global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
        }
        else if (wfn_type_ == "OCEPA") {
           global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
        }
	global_dpd_->contract442(&X, &T, &GF, 3, 3, 1.0, 1.0);
	global_dpd_->buf4_close(&X);
	global_dpd_->buf4_close(&T);
	global_dpd_->file2_close(&GF);
   }

   else if (twopdm_abcd_type == "COMPUTE" ) {
	// Fab += 4 * \sum{c,e,f} <ce|fa> * G_cefb
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[V,V]"), ints->DPD_ID("[V,V]"),
                  ints->DPD_ID("[V,V]"), ints->DPD_ID("[V,V]"), 0, "MO Ints <VV|VV>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");
	global_dpd_->contract442(&K, &G, &GF, 3, 3, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);
	global_dpd_->file2_close(&GF);
   }
}// end if (wfn_type_ != "OMP2")

        // close the integral file
	psio_->close(PSIF_LIBTRANS_DPD, 1);
        psio_->close(PSIF_OCC_DPD, 1);

	// Load dpd_file2 to SharedMatrix (GFock)
	// Load Fij
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('O'), ID('O'), "GF <O|O>");
	global_dpd_->file2_mat_init(&GF);
	global_dpd_->file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < occpiA[h]; ++i){
            for(int j = 0 ; j < occpiA[h]; ++j){
                GFock->add(h, i, j, GF.matrix[h][i][j]);
            }
	  }
	}
	global_dpd_->file2_close(&GF);

        // Load Fab
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "GF <V|V>");
	global_dpd_->file2_mat_init(&GF);
	global_dpd_->file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
	  for(int a = 0 ; a < virtpiA[h]; ++a){
            for(int b = 0 ; b < virtpiA[h]; ++b){
                GFock->add(h, a + occpiA[h], b + occpiA[h], GF.matrix[h][a][b]);
            }
	  }
	}
	global_dpd_->file2_close(&GF);

	psio_->close(PSIF_OCC_DENSITY, 1);
	if (print_ > 2) GFock->print();

}// end if (reference_ == "RESTRICTED")


//===========================================================================================
//========================= UHF =============================================================
//===========================================================================================
else if (reference_ == "UNRESTRICTED") {

/********************************************************************************************/
/************************** 2e-part *********************************************************/
/********************************************************************************************/
	dpdbuf4 G, K, T, L, X;
	dpdfile2 GF;

	psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
	psio_->open(PSIF_OCC_DENSITY, PSIO_OPEN_OLD);
        psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);

/********************************************************************************************/
/************************** Build X intermediates *******************************************/
/********************************************************************************************/
if (wfn_type_ != "OMP2") {
   if (twopdm_abcd_type == "DIRECT" ) {
 	// X_MNAC = 1/4 \sum{E,F} t_MN^EF * <AC||EF> = 1/2 \sum{E,F} t_MN^EF * <AC|EF>
        global_dpd_->buf4_init(&X, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[O,O]"),
                  ID("[V,V]"), ID("[O,O]"), 0, "X <VV|OO>");
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "MO Ints <VV|VV>");
        if (wfn_type_ == "OMP3" || wfn_type_ == "OMP2.5") {
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
        }
        else if (wfn_type_ == "OCEPA") {
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
        }
        global_dpd_->contract444(&K, &T, &X, 0, 0, 0.5, 0.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&T);
        global_dpd_->buf4_sort(&X, PSIF_OCC_DENSITY, rspq, ID("[O,O]"), ID("[V,V]"), "X <OO|VV>");
	global_dpd_->buf4_close(&X);

        // OMP2.5
        if (wfn_type_ == "OMP2.5") {
        global_dpd_->buf4_init(&X, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "X <OO|VV>");
        global_dpd_->buf4_scm(&X, 0.5);
	global_dpd_->buf4_close(&X);
        }

	// X_mnac = 1/4 \sum{e,f} t_mn^ef * <ac||ef> = 1/2 \sum{e,f} t_mn^ef * <ac|ef>
        global_dpd_->buf4_init(&X, PSIF_OCC_DENSITY, 0, ID("[v,v]"), ID("[o,o]"),
                  ID("[v,v]"), ID("[o,o]"), 0, "X <vv|oo>");
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "MO Ints <vv|vv>");
        if (wfn_type_ == "OMP3" || wfn_type_ == "OMP2.5") {
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
        }
        else if (wfn_type_ == "OCEPA") {
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2 <oo|vv>");
        }
        global_dpd_->contract444(&K, &T, &X, 0, 0, 0.5, 0.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&T);
        global_dpd_->buf4_sort(&X, PSIF_OCC_DENSITY, rspq, ID("[o,o]"), ID("[v,v]"), "X <oo|vv>");
	global_dpd_->buf4_close(&X);

        // OMP2.5
        if (wfn_type_ == "OMP2.5") {
        global_dpd_->buf4_init(&X, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "X <oo|vv>");
        global_dpd_->buf4_scm(&X, 0.5);
	global_dpd_->buf4_close(&X);
        }

        // X_MnAc = \sum{E,f} t_Mn^Ef * <Ac|Ef>
        global_dpd_->buf4_init(&X, PSIF_OCC_DENSITY, 0, ID("[V,v]"), ID("[O,o]"),
                  ID("[V,v]"), ID("[O,o]"), 0, "X <Vv|Oo>");
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "MO Ints <Vv|Vv>");
        if (wfn_type_ == "OMP3" || wfn_type_ == "OMP2.5") {
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
        }
        else if (wfn_type_ == "OCEPA") {
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");
        }
        global_dpd_->contract444(&K, &T, &X, 0, 0, 1.0, 0.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&T);
        global_dpd_->buf4_sort(&X, PSIF_OCC_DENSITY, rspq, ID("[O,o]"), ID("[V,v]"), "X <Oo|Vv>");
	global_dpd_->buf4_close(&X);

        // OMP2.5
        if (wfn_type_ == "OMP2.5") {
        global_dpd_->buf4_init(&X, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "X <Oo|Vv>");
        global_dpd_->buf4_scm(&X, 0.5);
	global_dpd_->buf4_close(&X);
        }

  }// end main if for X
} // end if (wfn_type_ != "OMP2")


/********************************************************************************************/
/************************** OO-Block ********************************************************/
/********************************************************************************************/
	// Build FIJ
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('O'), ID('O'), "GF <O|O>");

	// FIJ = 2 * \sum{M,N,K} <MN||KI> * G_MNKJ
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "MO Ints <OO||OO>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
	global_dpd_->contract442(&K, &G, &GF, 3, 3, 2.0, 0.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

	// FIJ += 2 * \sum{E,F,M} <EF||MI> * G_MJEF = 2 * \sum{E,F,M} <MI||EF> * G_MJEF
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
	global_dpd_->contract442(&K, &G, &GF, 1, 1, 2.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

	// FIJ += 4 * \sum{E,F,M} <EM||FI> * G_MEJF = 4 * \sum{E,F,M} <ME||IF> * G_MEJF
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV||OV>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
	global_dpd_->contract442(&K, &G, &GF, 2, 2, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

	// FIJ += 4 * \sum{m,N,k} <Nm|Ik> * G_NmJk
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "MO Ints <Oo|Oo>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "TPDM <Oo|Oo>");
	global_dpd_->contract442(&K, &G, &GF, 2, 2, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

	// FIJ += 4 * \sum{e,F,m} <Fe|Im> * G_JmFe = 4 * \sum{e,F,m} <Im|Fe> * G_JmFe
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "TPDM <Oo|Vv>");
	global_dpd_->contract442(&K, &G, &GF, 0, 0, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

	// FIJ += 4 * \sum{e,f,M} <Me|If> * G_MeJf
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "MO Ints <Ov|Ov>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "TPDM <Ov|Ov>");
	global_dpd_->contract442(&K, &G, &GF, 2, 2, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

if (wfn_type_ != "OMP2") {
        // FIJ += 4 * \sum{E,F,m} <Em|If> * G_EmJf = 4 * \sum{E,F,m} <If|Em> * G_JfEm => new
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "MO Ints <Ov|Vo>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "TPDM <Ov|Vo>");
	global_dpd_->contract442(&K, &G, &GF, 0, 0, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);
}// end if (wfn_type_ != "OMP2")

	// Close
	global_dpd_->file2_close(&GF);


	// Build Fij
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('o'), ID('o'), "GF <o|o>");

	// Fij = 2 * \sum{m,n,k} <mn||ki> * G_mnkj
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "MO Ints <oo||oo>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "TPDM <oo|oo>");
	global_dpd_->contract442(&K, &G, &GF, 3, 3, 2.0, 0.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

	// Fij += 2 * \sum{e,f,m} <ef||mi> * G_mjef = 2 * \sum{e,f,m} <mi||ef> * G_mjef
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "TPDM <oo|vv>");
	global_dpd_->contract442(&K, &G, &GF, 1, 1, 2.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

	// Fij += 4 * \sum{e,f,m} <em||fi> * G_mejf = 4 * \sum{e,f,m} <me||if> * G_mejf
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov||ov>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "TPDM <ov|ov>");
	global_dpd_->contract442(&K, &G, &GF, 2, 2, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

	// Fij += 4 * \sum{M,n,K} <Mn|Ki> * G_MnKj
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "MO Ints <Oo|Oo>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "TPDM <Oo|Oo>");
	global_dpd_->contract442(&K, &G, &GF, 3, 3, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

	// Fij += 4 * \sum{E,f,M} <Ef|Mi> * G_MjEf = 4 * \sum{E,f,M} <Mi|Ef> * G_MjEf
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "TPDM <Oo|Vv>");
	global_dpd_->contract442(&K, &G, &GF, 1, 1, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

	// Fij += 4 * \sum{E,F,m} <Em|Fi> * G_EmFj
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[V,o]"),
                  ID("[V,o]"), ID("[V,o]"), 0, "MO Ints <Vo|Vo>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,o]"), ID("[V,o]"),
                  ID("[V,o]"), ID("[V,o]"), 0, "TPDM <Vo|Vo>");
	global_dpd_->contract442(&K, &G, &GF, 3, 3, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

if (wfn_type_ != "OMP2") {
	// Fij += 4 * \sum{e,F,M} <Me|Fi> * G_MeFj => new
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "MO Ints <Ov|Vo>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "TPDM <Ov|Vo>");
	global_dpd_->contract442(&K, &G, &GF, 3, 3, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);
}// end if (wfn_type_ != "OMP2")

	// Close
	global_dpd_->file2_close(&GF);

/********************************************************************************************/
/************************** VV-Block ********************************************************/
/********************************************************************************************/
	// Build FAB
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "GF <V|V>");

	// FAB = 2 * \sum{M,N,E} <MN||EA> * G_MNEB
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
	global_dpd_->contract442(&K, &G, &GF, 3, 3, 2.0, 0.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);


if (wfn_type_ != "OMP2") {
   if (twopdm_abcd_type == "DIRECT" ) {
       	// FAB += 2 * \sum{E,F,C} <EF||CA> * G_EFCB = \sum{M,N,C} X_MNAC * t_MN^BC
        global_dpd_->buf4_init(&X, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "X <OO|VV>");
        if (wfn_type_ == "OMP3" || wfn_type_ == "OMP2.5") {
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
        }
        else if (wfn_type_ == "OCEPA") {
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
        }
	global_dpd_->contract442(&X, &T, &GF, 2, 2, 1.0, 1.0);
	global_dpd_->buf4_close(&X);
	global_dpd_->buf4_close(&T);
   }

   else if (twopdm_abcd_type == "COMPUTE" ) {
        // FAB = 2 * \sum{E,F,C} <EF||CA> * G_EFCB => new
	// FAB = 4 * \sum{E,F,C} <EF|CA> * G_EFCB => new
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[V,V]"), ints->DPD_ID("[V,V]"),
                  ints->DPD_ID("[V,V]"), ints->DPD_ID("[V,V]"), 0, "MO Ints <VV|VV>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");
	global_dpd_->contract442(&K, &G, &GF, 3, 3, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);
   }
}// end if (wfn_type_ != "OMP2")


	// FAB += 4 * \sum{M,N,E} <ME||NA> * G_MENB
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV||OV>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
	global_dpd_->contract442(&K, &G, &GF, 3, 3, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

	// FAB = 4 * \sum{m,N,e} <Nm|Ae> * G_NmBe
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "TPDM <Oo|Vv>");
	global_dpd_->contract442(&K, &G, &GF, 2, 2, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

	// FAB = 4 * \sum{m,n,E} <Em|An> * G_EmBn
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[V,o]"),
                  ID("[V,o]"), ID("[V,o]"), 0, "MO Ints <Vo|Vo>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,o]"), ID("[V,o]"),
                  ID("[V,o]"), ID("[V,o]"), 0, "TPDM <Vo|Vo>");
	global_dpd_->contract442(&K, &G, &GF, 2, 2, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

if (wfn_type_ != "OMP2") {
	// FAB = 4 * \sum{M,n,e} <Me|An> * G_MeBn => new
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "MO Ints <Ov|Vo>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "TPDM <Ov|Vo>");
	global_dpd_->contract442(&K, &G, &GF, 2, 2, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);
}// end if (wfn_type_ != "OMP2")

if (wfn_type_ != "OMP2") {
   if (twopdm_abcd_type == "DIRECT" ) {
       	// FAB += 4 * \sum{e,F,c} <Fe|Ac> * G_FeBc =  \sum{M,n,C} X_MnAc * t_Mn^Bc
        global_dpd_->buf4_init(&X, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "X <Oo|Vv>");
        if (wfn_type_ == "OMP3" || wfn_type_ == "OMP2.5") {
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
        }
        else if (wfn_type_ == "OCEPA") {
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");
        }
	global_dpd_->contract442(&X, &T, &GF, 2, 2, 1.0, 1.0);
	global_dpd_->buf4_close(&X);
	global_dpd_->buf4_close(&T);
   }

   else if (twopdm_abcd_type == "COMPUTE" ) {
	// FAB = 4 * \sum{e,F,c} <Fe|Ac> * G_FeBc => new
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[V,v]"), ints->DPD_ID("[V,v]"),
                  ints->DPD_ID("[V,v]"), ints->DPD_ID("[V,v]"), 0, "MO Ints <Vv|Vv>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "TPDM <Vv|Vv>");
	global_dpd_->contract442(&K, &G, &GF, 2, 2, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);
   }
}// end if (wfn_type_ != "OMP2")
	// Close
	global_dpd_->file2_close(&GF);


	// Build Fab
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('v'), ID('v'), "GF <v|v>");

	// Fab = 2 * \sum{m,n,e} <mn||ea> * G_mneb
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "TPDM <oo|vv>");
	global_dpd_->contract442(&K, &G, &GF, 3, 3, 2.0, 0.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

if (wfn_type_ != "OMP2") {

}// end if (wfn_type_ != "OMP2")

if (wfn_type_ != "OMP2") {
   if (twopdm_abcd_type == "DIRECT" ) {
       	// Fab += 2 * \sum{e,f,c} <ef||ca> * G_efcb = \sum{m,n,c} X_mnac * t_mn^bc
        global_dpd_->buf4_init(&X, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "X <oo|vv>");
        if (wfn_type_ == "OMP3" || wfn_type_ == "OMP2.5") {
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
        }
        else if (wfn_type_ == "OCEPA") {
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2 <oo|vv>");
        }
	global_dpd_->contract442(&X, &T, &GF, 2, 2, 1.0, 1.0);
	global_dpd_->buf4_close(&X);
	global_dpd_->buf4_close(&T);
   }

   else if (twopdm_abcd_type == "COMPUTE" ) {
        // Fab = 2 * \sum{efc} <ef||ca> * G_efcb => new
	// Fab = 4 * \sum{efc} <ef|ca> * G_efcb => new
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[v,v]"), ints->DPD_ID("[v,v]"),
                  ints->DPD_ID("[v,v]"), ints->DPD_ID("[v,v]"), 0, "MO Ints <vv|vv>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "TPDM <vv|vv>");
	global_dpd_->contract442(&K, &G, &GF, 3, 3, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);
   }
}// end if (wfn_type_ != "OMP2")


	// Fab += 4 * \sum{m,n,e} <me||na> * G_menb
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov||ov>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "TPDM <ov|ov>");
	global_dpd_->contract442(&K, &G, &GF, 3, 3, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

	// Fab = 4 * \sum{M,n,E} <Mn|Ea> * G_MnEb
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "TPDM <Oo|Vv>");
	global_dpd_->contract442(&K, &G, &GF, 3, 3, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

	// Fab = 4 * \sum{M,N,e} <Me|Na> * G_MeNb
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "MO Ints <Ov|Ov>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "TPDM <Ov|Ov>");
	global_dpd_->contract442(&K, &G, &GF, 3, 3, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

if (wfn_type_ != "OMP2") {
	// Fab = 4 * \sum{m,N,E} <Em|Na> * G_EmNb = 4 * \sum{M,n,e} <Na|Em> * G_NbEm => new
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "MO Ints <Ov|Vo>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "TPDM <Ov|Vo>");
	global_dpd_->contract442(&K, &G, &GF, 1, 1, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);
}// end if (wfn_type_ != "OMP2")

if (wfn_type_ != "OMP2") {
   if (twopdm_abcd_type == "DIRECT" ) {
       	// Fab += 4 * \sum{E,f,C} <Ef|Ca> * G_EfCb =  \sum{M,n,C} X_MnCa * t_Mn^Cb
        global_dpd_->buf4_init(&X, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "X <Oo|Vv>");
        if (wfn_type_ == "OMP3" || wfn_type_ == "OMP2.5") {
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
        }
        else if (wfn_type_ == "OCEPA") {
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");
        }
	global_dpd_->contract442(&X, &T, &GF, 3, 3, 1.0, 1.0);
	global_dpd_->buf4_close(&X);
	global_dpd_->buf4_close(&T);
   }

   else if (twopdm_abcd_type == "COMPUTE" ) {
	// Fab = 4 * \sum{e,F,c} <Ef|Ca> * G_EfCb => new
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[V,v]"), ints->DPD_ID("[V,v]"),
                  ints->DPD_ID("[V,v]"), ints->DPD_ID("[V,v]"), 0, "MO Ints <Vv|Vv>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "TPDM <Vv|Vv>");
	global_dpd_->contract442(&K, &G, &GF, 3, 3, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

   }
}// end if (wfn_type_ != "OMP2")

	// Close
	global_dpd_->file2_close(&GF);
	psio_->close(PSIF_LIBTRANS_DPD, 1);
        psio_->close(PSIF_OCC_DPD, 1);

/********************************************************************************************/
/************************** Load dpd_file2 to SharedMatrix (GFock) **************************/
/********************************************************************************************/
	// Load FIJ
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('O'), ID('O'), "GF <O|O>");
	global_dpd_->file2_mat_init(&GF);
	global_dpd_->file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < occpiA[h]; ++i){
            for(int j = 0 ; j < occpiA[h]; ++j){
                GFockA->add(h, i, j, GF.matrix[h][i][j]);
            }
	  }
	}
	global_dpd_->file2_close(&GF);

	// Load Fij
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('o'), ID('o'), "GF <o|o>");
	global_dpd_->file2_mat_init(&GF);
	global_dpd_->file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < occpiB[h]; ++i){
            for(int j = 0 ; j < occpiB[h]; ++j){
                GFockB->add(h, i, j, GF.matrix[h][i][j]);
            }
	  }
	}
	global_dpd_->file2_close(&GF);

	// Load FAB
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "GF <V|V>");
	global_dpd_->file2_mat_init(&GF);
	global_dpd_->file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
	  for(int a = 0 ; a < virtpiA[h]; ++a){
            for(int b = 0 ; b < virtpiA[h]; ++b){
                GFockA->add(h, a + occpiA[h], b + occpiA[h], GF.matrix[h][a][b]);
            }
	  }
	}
	global_dpd_->file2_close(&GF);

	// Load Fab
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('v'), ID('v'), "GF <v|v>");
	global_dpd_->file2_mat_init(&GF);
	global_dpd_->file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
	  for(int a = 0 ; a < virtpiB[h]; ++a){
            for(int b = 0 ; b < virtpiB[h]; ++b){
                GFockB->add(h, a + occpiB[h], b + occpiB[h], GF.matrix[h][a][b]);
            }
	  }
	}
	global_dpd_->file2_close(&GF);
	psio_->close(PSIF_OCC_DENSITY, 1);

        // Print
	if (print_ > 1) {
	  GFockA->print();
	  GFockB->print();
	}

}// end if (reference_ == "UNRESTRICTED")
//outfile->Printf("\n gfock_diag done. \n");

}// end gfock_diag
}} // End Namespaces
