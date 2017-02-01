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

void OCCWave::ocepa_g_int()
{
 //outfile->Printf("\n G_int is starting... \n");
 if (reference_ == "RESTRICTED") {
       // initialize
	GooA->zero();
	GvvA->zero();

	dpdbuf4 Tau, T;
	dpdfile2 G;

	psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);
        psio_->open(PSIF_OCC_DENSITY, PSIO_OPEN_OLD);


	// Open amplitude files
	global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
	global_dpd_->buf4_init(&Tau, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Tau <OO|VV>");

	// Occupied-Occupied block
	// G_IJ = \sum{M,E,F} t_IM^EF(1) * (2l_EF^JM - l_FE^JM)
	// G_IJ = \sum{M,E,F} t_IM^EF(1) * (2t_JM^EF - t_MJ^EF)
	global_dpd_->file2_init(&G, PSIF_OCC_DENSITY, 0, ID('O'), ID('O'), "G <O|O>");
	global_dpd_->contract442(&T, &Tau, &G, 0, 0, 1.0, 0.0);
	global_dpd_->file2_close(&G);


	// Virtual-Virtual block
	// G_AB = -\sum{M,N,E} t_MN^BE(1) * (2l_AE^MN - l_EA^MN)
	// G_AB = -\sum{M,N,E} t_MN^BE(1) * (2t_MN^AE - t_MN^EA)
	global_dpd_->file2_init(&G, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "G <V|V>");
	global_dpd_->contract442(&Tau, &T, &G, 2, 2, -1.0, 0.0);
	global_dpd_->file2_close(&G);


	// Close amplitude files
	global_dpd_->buf4_close(&T);
	global_dpd_->buf4_close(&Tau);


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

	psio_->close(PSIF_OCC_DPD, 1);
        psio_->close(PSIF_OCC_DENSITY, 1);

	if (print_ > 1) {
	  GooA->print();
	  GvvA->print();
	}

 }// end if (reference_ == "RESTRICTED")

 else if (reference_ == "UNRESTRICTED") {
        // initialize
	GooA->zero();
	GooB->zero();
	GvvA->zero();
	GvvB->zero();

	dpdbuf4 TAA, TAB, TBB, LAA, LAB, LBB;
	dpdfile2 G;

	psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);
        psio_->open(PSIF_OCC_DENSITY, PSIO_OPEN_OLD);

	// Open amplitude files
	global_dpd_->buf4_init(&TAA, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
	global_dpd_->buf4_init(&TBB, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2 <oo|vv>");
	global_dpd_->buf4_init(&TAB, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");
	global_dpd_->buf4_init(&LAA, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
	global_dpd_->buf4_init(&LBB, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2 <oo|vv>");
	global_dpd_->buf4_init(&LAB, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");


	// Occupied-Occupied block
	// Alpha-Alpha spin case
	// G_IM = 1/2 \sum{N,E,F} t_IN^EF(1) * l_EF^MN(1) = 1/2 \sum{N,E,F} t_IN^EF(1) * t_MN^EF(1)
	global_dpd_->file2_init(&G, PSIF_OCC_DENSITY, 0, ID('O'), ID('O'), "G <O|O>");
	global_dpd_->contract442(&TAA, &LAA, &G, 0, 0, 0.5, 0.0);

	// G_IM += \sum{n,E,f} t_In^Ef(1) * l_Ef^Mn(1) = \sum{N,E,F} t_In^Ef(1) * t_Mn^Ef(1)
	global_dpd_->contract442(&TAB, &LAB, &G, 0, 0, 1.0, 1.0);
	global_dpd_->file2_close(&G);


	// Beta-Beta spin case
	// G_im = 1/2 \sum{n,e,f} t_in^ef(1) * l_ef^mn(1) = 1/2 \sum{n,e,f} t_in^ef(1) * t_mn^ef(1)
	global_dpd_->file2_init(&G, PSIF_OCC_DENSITY, 0, ID('o'), ID('o'), "G <o|o>");
	global_dpd_->contract442(&TBB, &LBB, &G, 0, 0, 0.5, 0.0);

	// G_im  += \sum{N,e,F} t_Ni^Fe(1) * l_Fe^Nm(1) = \sum{N,e,F} t_Ni^Fe(1) * t_Nm^Fe(1)
	global_dpd_->contract442(&TAB, &LAB, &G, 1, 1, 1.0, 1.0);
	global_dpd_->file2_close(&G);



	// Virtual-Virtual block
	// Alpha-Alpha spin case
	// G_EA = -1/2 \sum{M,N,F} t_MN^AF(1) * l_EF^MN(1) = -1/2 \sum{M,N,F} t_MN^AF(1) * t_MN^EF(1)
	global_dpd_->file2_init(&G, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "G <V|V>");
	global_dpd_->contract442(&TAA, &LAA, &G, 2, 2, -0.5, 0.0);

	// G_EA += - \sum{M,n,f} t_Mn^Af(1) * l_Ef^Mn(1) = - \sum{M,n,f} t_Mn^Af(1) * t_Mn^Ef(1)
	global_dpd_->contract442(&TAB, &LAB, &G, 2, 2, -1.0, 1.0);
	global_dpd_->file2_close(&G);


	// Beta-Beta spin case
	// G_ea = -1/2 \sum{m,n,f} t_mn^af(1) * l_ef^mn(1) = -1/2 \sum{m,n,f} t_mn^af(1) * t_mn^ef(1)
	global_dpd_->file2_init(&G, PSIF_OCC_DENSITY, 0, ID('v'), ID('v'), "G <v|v>");
	global_dpd_->contract442(&TBB, &LBB, &G, 2, 2, -0.5, 0.0);

	// G_ea += - \sum{M,n,F} t_Mn^Fa(1) * l_Fe^Mn(1) = - \sum{M,n,F} t_Mn^Fa(1) * t_Mn^Fe(1)
	global_dpd_->contract442(&TAB, &LAB, &G, 3, 3, -1.0, 1.0);
	global_dpd_->file2_close(&G);


	// Close amplitude files
	global_dpd_->buf4_close(&TAA);
	global_dpd_->buf4_close(&TBB);
	global_dpd_->buf4_close(&TAB);
	global_dpd_->buf4_close(&LAA);
	global_dpd_->buf4_close(&LBB);
	global_dpd_->buf4_close(&LAB);


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

	if (print_ > 1) {
	  GooA->print();
	  GooB->print();
	  GvvA->print();
	  GvvB->print();
	}
 }// end if (reference_ == "UNRESTRICTED")

  //outfile->Printf("\n G_int done... \n");

} // end of G_int
}} // End Namespaces
