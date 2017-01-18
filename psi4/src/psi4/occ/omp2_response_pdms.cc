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

#include "psi4/libqt/qt.h"
#include "psi4/libtrans/integraltransform.h"
#include "psi4/libmints/matrix.h"
#include "psi4/libpsio/psio.hpp"
#include "occwave.h"
#include "defines.h"


using namespace std;

namespace psi{ namespace occwave{

void OCCWave::omp2_response_pdms()
{
 //outfile->Printf("\n omp2_response_pdms is starting... \n");

 if (reference_ == "RESTRICTED") {
        // initialize
	gamma1corr->zero();
	g1symm->zero();

        // Build G intermediates
        timer_on("G int");
	omp2_g_int();
        timer_off("G int");

        // Build OPDM
        timer_on("OPDM");
        // OO-block
	#pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < aoccpiA[h]; ++i){
            for(int j = 0 ; j < aoccpiA[h]; ++j){
                g1symm->set(h, i, j, GooA->get(h, i, j) + GooA->get(h, j, i));
            }
	  }
	}

        // VV-Block
	#pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int a = 0 ; a < avirtpiA[h]; ++a){
            for(int b = 0 ; b < avirtpiA[h]; ++b){
                int aa = a + occpiA[h];
                int bb = b + occpiA[h];
                g1symm->set(h, aa, bb, GvvA->get(h, a, b) + GvvA->get(h, b, a));
            }
	  }
	}

	g1symm->scale(-1.0);
	gamma1corr->copy(g1symm); // correlation opdm

        // REF contribution
	for(int h=0; h<nirrep_; h++) {
	  if (occpiA[h] != 0) {
	    for (int i=0; i<occpiA[h];i++) {
	      g1symm->add(h,i,i,2.0);
	    }
	  }
	}
        timer_off("OPDM");

        // print
        if (print_ > 2) g1symm->print();

        // Build TPDM
        timer_on("TPDM OOVV");
	omp2_tpdm_oovv();
        timer_off("TPDM OOVV");
        timer_on("TPDM REF");
	tpdm_ref();
        timer_off("TPDM REF");
        timer_on("TPDM CORR OPDM");
	tpdm_corr_opdm();
        timer_off("TPDM CORR OPDM");

 }// end if (reference_ == "RESTRICTED")

 else if (reference_ == "UNRESTRICTED") {
        // Initialize
	gamma1corrA->zero();
	gamma1corrB->zero();
	g1symmA->zero();
	g1symmB->zero();

        // Build G intermediates
        timer_on("G int");
	omp2_g_int();
        timer_off("G int");

        // OPDM
        timer_on("OPDM");
	// OO-block alpha contrb.
	#pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < aoccpiA[h]; ++i){
            for(int j = 0 ; j < aoccpiA[h]; ++j){
                g1symmA->set(h, i, j, GooA->get(h, i, j) + GooA->get(h, j, i));
            }
	  }
	}

	// OO-block beta contrb.
	#pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < aoccpiB[h]; ++i){
            for(int j = 0 ; j < aoccpiB[h]; ++j){
                g1symmB->set(h, i, j, GooB->get(h, i, j) + GooB->get(h, j, i));
            }
	  }
	}

	// VV-block alpha contrb.
        #pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int a = 0 ; a < avirtpiA[h]; ++a){
            for(int b = 0 ; b < avirtpiA[h]; ++b){
                int aa = a + occpiA[h];
                int bb = b + occpiA[h];
                g1symmA->set(h, aa, bb, GvvA->get(h, a, b) + GvvA->get(h, b, a));
            }
	  }
	}

        // VV-block beta contrb.
        #pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int a = 0 ; a < avirtpiB[h]; ++a){
            for(int b = 0 ; b < avirtpiB[h]; ++b){
                int aa = a + occpiB[h];
                int bb = b + occpiB[h];
                g1symmB->set(h, aa, bb, GvvB->get(h, a, b) + GvvB->get(h, b, a));
            }
	  }
	}


	g1symmA->scale(-0.5);
	g1symmB->scale(-0.5);
	gamma1corrA->copy(g1symmA); // correlation opdm
	gamma1corrB->copy(g1symmB); // correlation opdm

        // REF contribution
	// alpha contrb.
        #pragma omp parallel for
	for(int h=0; h<nirrep_; h++) {
	  if (occpiA[h] != 0) {
	    for (int i=0; i<occpiA[h];i++) {
	      g1symmA->add(h,i,i,1.0);
	    }
	  }
	}

	// beta contrb.
        #pragma omp parallel for
	for(int h=0; h<nirrep_; h++) {
	  if (occpiB[h] != 0) {
	    for (int i=0; i<occpiB[h];i++) {
	      g1symmB->add(h,i,i,1.0);
	    }
	  }
	}
        timer_off("OPDM");

        // print
        if (print_ > 2) {
	   g1symmA->print();
	   g1symmB->print();

        }

        // TPDM
        timer_on("TPDM OOVV");
	omp2_tpdm_oovv();
        timer_off("TPDM OOVV");
        timer_on("TPDM REF");
	tpdm_ref();
        timer_off("TPDM REF");
        timer_on("TPDM CORR OPDM");
	tpdm_corr_opdm();
        timer_off("TPDM CORR OPDM");


 }// end if (reference_ == "UNRESTRICTED")

  //outfile->Printf("\n omp2_response_pdms done... \n");

} // end of omp2_response_pdms


void OCCWave::omp2_g_int()
{
        //outfile->Printf("\n omp2_g_int is starting... \n");

 if (reference_ == "RESTRICTED") {
	GooA->zero();
	GvvA->zero();


	dpdbuf4 T, Tau;
	dpdfile2 Go,Gv;

	psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);
        psio_->open(PSIF_OCC_DENSITY, PSIO_OPEN_OLD);

	global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T <OO|VV>");
	global_dpd_->buf4_init(&Tau, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Tau <OO|VV>");

	// G_mi = \sum{n,e,f} t_mn^ef * tau_in^ef
	global_dpd_->file2_init(&Go, PSIF_OCC_DENSITY, 0, ID('O'), ID('O'), "G <O|O>");
	global_dpd_->contract442(&T, &Tau, &Go, 0, 0, 1.0, 0.0);
	global_dpd_->file2_close(&Go);

	// G_ae = -\sum{m,n,f} t_mn^ef * tau_mn^af
	global_dpd_->file2_init(&Gv, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "G <V|V>");
	global_dpd_->contract442(&Tau, &T, &Gv, 2, 2, -1.0, 0.0);
	global_dpd_->file2_close(&Gv);

	global_dpd_->buf4_close(&T);
	global_dpd_->buf4_close(&Tau);

	// Load dpd_file2 to Matrix (Goo)
	global_dpd_->file2_init(&Go, PSIF_OCC_DENSITY, 0, ID('O'), ID('O'), "G <O|O>");
	global_dpd_->file2_mat_init(&Go);
	global_dpd_->file2_mat_rd(&Go);
        #pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < aoccpiA[h]; ++i){
            for(int j = 0 ; j < aoccpiA[h]; ++j){
                GooA->set(h, i, j, Go.matrix[h][i][j]);
            }
	  }
	}
	global_dpd_->file2_close(&Go);


	// Load dpd_file2 to Matrix (Gvv)
	global_dpd_->file2_init(&Gv, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "G <V|V>");
	global_dpd_->file2_mat_init(&Gv);
	global_dpd_->file2_mat_rd(&Gv);
        #pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < avirtpiA[h]; ++i){
            for(int j = 0 ; j < avirtpiA[h]; ++j){
                GvvA->set(h, i, j, Gv.matrix[h][i][j]);
            }
	  }
	}
	global_dpd_->file2_close(&Gv);

	psio_->close(PSIF_OCC_DPD, 1);
        psio_->close(PSIF_OCC_DENSITY, 1);


	if (print_ > 3) {
	  GooA->print();
	  GvvA->print();
	}

 }// end if (reference_ == "RESTRICTED")

 else if (reference_ == "UNRESTRICTED") {
	GooA->zero();
	GooB->zero();
	GvvA->zero();
	GvvB->zero();

	dpdbuf4 TAA, TAB, TBB;
        dpdbuf4 TAA2, TBB2;
        dpdbuf4 TAB2_;
	//dpdbuf4 TAA, TAB, TBB, TAA2, TAB2_, TBB2;
	dpdfile2 G;

	psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);
        psio_->open(PSIF_OCC_DENSITY, PSIO_OPEN_OLD);

	// Open amplitude files
	global_dpd_->buf4_init(&TAA, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
	global_dpd_->buf4_init(&TBB, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
	global_dpd_->buf4_init(&TAB, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
	global_dpd_->buf4_init(&TAA2, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
	global_dpd_->buf4_init(&TBB2, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
	global_dpd_->buf4_init(&TAB2_, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");

	// Occupied-Occupied block
	// Alpha-Alpha spin case
	// G_IM = 1/2 \sum{N,E,F} t_IN^EF * l_EF^MN = 1/2 \sum{N,E,F} t_IN^EF * t_MN^EF
	global_dpd_->file2_init(&G, PSIF_OCC_DENSITY, 0, ID('O'), ID('O'), "G <O|O>");
	global_dpd_->contract442(&TAA, &TAA2, &G, 0, 0, 0.5, 0.0);

	// G_IM += \sum{n,E,f} t_In^Ef * l_Ef^Mn = \sum{N,E,F} t_In^Ef * t_Mn^Ef
	global_dpd_->contract442(&TAB, &TAB2_, &G, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&G);


	// Beta-Beta spin case
	// G_im = 1/2 \sum{n,e,f} t_in^ef * l_ef^mn = 1/2 \sum{n,e,f} t_in^ef * t_mn^ef
	global_dpd_->file2_init(&G, PSIF_OCC_DENSITY, 0, ID('o'), ID('o'), "G <o|o>");
	global_dpd_->contract442(&TBB, &TBB2, &G, 0, 0, 0.5, 0.0);

	// G_im  += \sum{N,e,F} t_Ni^Fe * l_Fe^Nm = \sum{N,e,F} t_Ni^Fe * t_Nm^Fe
	global_dpd_->contract442(&TAB, &TAB2_, &G, 1, 1, 1.0, 1.0);
	global_dpd_->file2_close(&G);



	// Virtual-Virtual block
	// Alpha-Alpha spin case
	// G_EA = -1/2 \sum{M,N,F} t_MN^AF * l_EF^MN = -1/2 \sum{M,N,F} t_MN^AF * t_MN^EF
	global_dpd_->file2_init(&G, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "G <V|V>");
	global_dpd_->contract442(&TAA, &TAA2, &G, 2, 2, -0.5, 0.0);

	// G_EA += - \sum{M,n,f} t_Mn^Af * l_Ef^Mn = - \sum{M,n,f} t_Mn^Af * t_Mn^Ef
	global_dpd_->contract442(&TAB, &TAB2_, &G, 2, 2, -1.0, 1.0);
	global_dpd_->file2_close(&G);

	// Beta-Beta spin case
	// G_ea = -1/2 \sum{m,n,f} t_mn^af * l_ef^mn = -1/2 \sum{m,n,f} t_mn^af * t_mn^ef
	global_dpd_->file2_init(&G, PSIF_OCC_DENSITY, 0, ID('v'), ID('v'), "G <v|v>");
	global_dpd_->contract442(&TBB, &TBB2, &G, 2, 2, -0.5, 0.0);

	// G_ea += - \sum{M,n,F} t_Mn^Fa * l_Fe^Mn = - \sum{M,n,F} t_Mn^Fa * t_Mn^Fe
	global_dpd_->contract442(&TAB, &TAB2_, &G, 3, 3, -1.0, 1.0);
	global_dpd_->file2_close(&G);

	// Close amplitude files
	global_dpd_->buf4_close(&TAA);
	global_dpd_->buf4_close(&TBB);
	global_dpd_->buf4_close(&TAB);
	global_dpd_->buf4_close(&TAA2);
	global_dpd_->buf4_close(&TBB2);
	global_dpd_->buf4_close(&TAB2_);


	// Load dpd_file2 to Matrix (Goo)
	// Alpha-Alpha spin case
	global_dpd_->file2_init(&G, PSIF_OCC_DENSITY, 0, ID('O'), ID('O'), "G <O|O>");
	global_dpd_->file2_mat_init(&G);
	global_dpd_->file2_mat_rd(&G);
        #pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < aoccpiA[h]; ++i){
            for(int j = 0 ; j < aoccpiA[h]; ++j){
                GooA->set(h, i, j, G.matrix[h][i][j]);
            }
	  }
	}
	global_dpd_->file2_close(&G);

	// Beta-Beta spin case
	global_dpd_->file2_init(&G, PSIF_OCC_DENSITY, 0, ID('o'), ID('o'), "G <o|o>");
	global_dpd_->file2_mat_init(&G);
	global_dpd_->file2_mat_rd(&G);
        #pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < aoccpiB[h]; ++i){
            for(int j = 0 ; j < aoccpiB[h]; ++j){
                GooB->set(h, i, j, G.matrix[h][i][j]);
            }
	  }
	}
	global_dpd_->file2_close(&G);



	// Load dpd_file2 to Matrix (Gvv)
	// Alpha-Alpha spin case
	global_dpd_->file2_init(&G, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "G <V|V>");
	global_dpd_->file2_mat_init(&G);
	global_dpd_->file2_mat_rd(&G);
        #pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < avirtpiA[h]; ++i){
            for(int j = 0 ; j < avirtpiA[h]; ++j){
                GvvA->set(h, i, j, G.matrix[h][i][j]);
            }
	  }
	}
	global_dpd_->file2_close(&G);

	// Beta-Beta spin case
	global_dpd_->file2_init(&G, PSIF_OCC_DENSITY, 0, ID('v'), ID('v'), "G <v|v>");
	global_dpd_->file2_mat_init(&G);
	global_dpd_->file2_mat_rd(&G);
        #pragma omp parallel for
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < avirtpiB[h]; ++i){
            for(int j = 0 ; j < avirtpiB[h]; ++j){
                GvvB->set(h, i, j, G.matrix[h][i][j]);
            }
	  }
	}
	global_dpd_->file2_close(&G);

	psio_->close(PSIF_OCC_DPD, 1);
        psio_->close(PSIF_OCC_DENSITY, 1);

	if (print_ > 3) {
	  GooA->print();
	  GooB->print();
	  GvvA->print();
	  GvvB->print();
	}

 }// end if (reference_ == "UNRESTRICTED")

  //outfile->Printf("\n omp2_g_int done... \n");

} // end of omp2_g_int
}} // End Namespaces
