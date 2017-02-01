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


using namespace std;

namespace psi{ namespace occwave{

void OCCWave::gfock()
{

//outfile->Printf("\n omp3_gfock is starting... \n");
//===========================================================================================
//========================= RHF =============================================================
//===========================================================================================
if (reference_ == "RESTRICTED") {

	// Initialize
	GFock->zero();

        // 1e-part
	HG1->zero();
	HG1->gemm(false, false, 1.0, HmoA, g1symm, 0.0);
	GFock->add(HG1);
	Ecc_rdm = HG1->trace() + Enuc; // One-electron contribution to MP2L
        Eopdm = Ecc_rdm;
	//outfile->Printf("\tOPDM energy (a.u.)          : %12.14f\n", Ecc_rdm);

        // 2e-part
	dpdbuf4 G, K, X, T, Y;
	dpdfile2 GF;

	psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
	psio_->open(PSIF_OCC_DENSITY, PSIO_OPEN_OLD);
        psio_->open(PSIF_OCC_DPD, PSIO_OPEN_OLD);


if (wfn_type_ != "OMP2") {
   // Build X intermediate
   if (twopdm_abcd_type == "DIRECT" ) {
        // With this algorithm cost changes to v5 => o2v4 + o2v3, roughly v/o times faster
 	// X_MNIC = 2\sum{E,F} t_MN^EF(1) * <IC|EF>
        global_dpd_->buf4_init(&X, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "X <OO|OV>");
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV|VV>");
        if (wfn_type_ == "OMP3" || wfn_type_ == "OMP2.5") {
             global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
        }
        else if (wfn_type_ == "OCEPA") {
             global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
        }
        global_dpd_->contract444(&T, &K, &X, 0, 0, 2.0, 0.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&T);
	global_dpd_->buf4_close(&X);

        // OMP2.5
        if (wfn_type_ == "OMP2.5") {
        global_dpd_->buf4_init(&X, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "X <OO|OV>");
        global_dpd_->buf4_scm(&X, 0.5);
	global_dpd_->buf4_close(&X);
        }
    }
}// end if (wfn_type_ != "OMP2")


/************************************************************************************************/
/*********************************** Build Fai **************************************************/
/************************************************************************************************/
	// Build Fai
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('V'), ID('O'), "GF <V|O>");

	// Fai += 4 * \sum{m,n,k} <km|na> * G_kmni
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,O]"), ints->DPD_ID("[O,V]"),
                  ints->DPD_ID("[O,O]"), ints->DPD_ID("[O,V]"), 0, "MO Ints <OO|OV>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
	global_dpd_->contract442(&K, &G, &GF, 3, 3, 4.0, 0.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

if (wfn_type_ == "OMP2" && incore_iabc_ == 1) {
        dpdfile2 G;

        // Build Virtual-Virtual block of correlation OPDM
        global_dpd_->file2_init(&G, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "CORR OPDM <V|V>");
        global_dpd_->file2_mat_init(&G);
        for(int h = 0; h < nirrep_; ++h){
            for(int i = 0 ; i < avirtpiA[h]; ++i){
                for(int j = 0 ; j < avirtpiA[h]; ++j){
                    G.matrix[h][i][j] = gamma1corr->get(h, i + occpiA[h], j + occpiA[h]);
                }
            }
        }
        global_dpd_->file2_mat_wrt(&G);
        global_dpd_->file2_close(&G);

        // Sort
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                      ID("[O,V]"), ID("[V>=V]+"), 0, "MO Ints (OV|VV)");
        global_dpd_->buf4_sort(&K, PSIF_LIBTRANS_DPD , qprs, ID("[V,O]"), ID("[V,V]"), "MO Ints (VO|VV)");
	global_dpd_->buf4_close(&K);

	// Fai += 2 * \sum{e,f} (ai|ef) * Gcorr(e,f)
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[V,O]"), ints->DPD_ID("[V,V]"),
                  ints->DPD_ID("[V,O]"), ints->DPD_ID("[V,V]"), 0, "MO Ints (VO|VV)");
        global_dpd_->file2_init(&G, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "CORR OPDM <V|V>");
        global_dpd_->contract422(&K, &G, &GF, 0, 0, 2.0, 1.0);
	global_dpd_->buf4_close(&K);
        global_dpd_->file2_close(&G);
}// end if (wfn_type_ == "OMP2" && incore_iabc_ == 1)

else if (wfn_type_ == "OMP2" && incore_iabc_ == 0) {
      	IWL ERIIN(psio_.get(), PSIF_OCC_IABC, 0.0, 1, 1);
	int ilsti,nbuf,index,fi;
	double value = 0.0;
        double summ = 0.0;

 do
 {
        ilsti = ERIIN.last_buffer();
        nbuf = ERIIN.buffer_count();

   fi = 0;
   for (int idx=0; idx < nbuf; idx++ )
   {

        int i = ERIIN.labels()[fi];
            i = abs(i);
        int e = ERIIN.labels()[fi+1];
        int a = ERIIN.labels()[fi+2];
        int f = ERIIN.labels()[fi+3];
        value = ERIIN.values()[idx];
        fi += 4;

        int i_pitzer = qt2pitzerA[i];
        int e_pitzer = qt2pitzerA[e];
        int a_pitzer = qt2pitzerA[a];
        int f_pitzer = qt2pitzerA[f];

        int hi = mosym[i_pitzer];
        int ha = mosym[a_pitzer];
        int he = mosym[e_pitzer];
        int hf = mosym[f_pitzer];

        if (hi == ha && he == hf) {
            int ii = pitzer2symblk[i_pitzer];
            int ee = pitzer2symblk[e_pitzer];
            int aa = pitzer2symblk[a_pitzer];
            int ff = pitzer2symblk[f_pitzer];
            summ = 2.0 * value * gamma1corr->get(he,ee,ff);
            GFock->add(ha, aa, ii, summ);
        }

   }
        if(!ilsti)
	  ERIIN.fetch();

 } while(!ilsti);
}// end else if (wfn_type_ == "OMP2" && incore_iabc_ == 0)

else if (wfn_type_ != "OMP2") {
	// Fai += 4 * \sum{e,m,f} <me|af> * G_meif
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,V]"), ints->DPD_ID("[V,V]"),
                  ints->DPD_ID("[O,V]"), ints->DPD_ID("[V,V]"), 0, "MO Ints <OV|VV>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
	global_dpd_->contract442(&K, &G, &GF, 2, 2, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);
}// end else if (wfn_type_ != "OMP2")

if (wfn_type_ == "OMP2" && incore_iabc_ == 0) {
      // Fai += 8 * \sum{e,m,f} <ma|ef> * G_mief
      	IWL ERIIN(psio_.get(), PSIF_OCC_IABC, 0.0, 1, 1);
	int ilsti,nbuf,index,fi;
	double value = 0;

       SymBlockMatrix *Goovv = new SymBlockMatrix("TPDM <OO|VV>", nirrep_, oo_pairpiAA, vv_pairpiAA);
       Goovv->zero();
       global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                     ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
       Goovv->set(G);
       global_dpd_->buf4_close(&G);

 do
 {
        ilsti = ERIIN.last_buffer();
        nbuf = ERIIN.buffer_count();

   fi = 0;
   for (int idx=0; idx < nbuf; idx++ )
   {

        int m = ERIIN.labels()[fi];
            m = abs(m);
        int a = ERIIN.labels()[fi+1];
        int e = ERIIN.labels()[fi+2];
        int f = ERIIN.labels()[fi+3];
        value = ERIIN.values()[idx];
        fi += 4;

        int m_pitzer = qt2pitzerA[m];
        int a_pitzer = qt2pitzerA[a];
        int e_pitzer = qt2pitzerA[e];
        int f_pitzer = qt2pitzerA[f];

        int hm = mosym[m_pitzer];
        int ha = mosym[a_pitzer];
        int he = mosym[e_pitzer];
        int hf = mosym[f_pitzer];

        int hma = ha^hm;
        int hef = he^hf;

        int E = e - nooA;// convert to vir_qt
        int F = f - nooA;

        if (hma == hef) {
            double summ = 0.0;
            for (int i=0; i < occpiA[ha]; i++) {
                 int I = i + occ_offA[ha];
                 int mi = oo_pairidxAA->get(hma, m, I);
                 int ef = vv_pairidxAA->get(hef, E, F);
                 summ = 8.0 * value * Goovv->get(hma, mi, ef);
                 int aa = pitzer2symblk[a_pitzer];
                 GFock->add(ha, aa, i, summ);
            }
        }

   }
        if(!ilsti)
	  ERIIN.fetch();

 } while(!ilsti);
       delete Goovv;

} // end if (wfn_type_ == "OMP2" && incore_iabc_ == 0)

else {
        // Fai += 8 * \sum{e,m,f} <ma|ef> * G_mief
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,V]"), ints->DPD_ID("[V,V]"),
                  ints->DPD_ID("[O,V]"), ints->DPD_ID("[V,V]"), 0, "MO Ints <OV|VV>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
	global_dpd_->contract442(&K, &G, &GF, 1, 1, 8.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);
	global_dpd_->file2_close(&GF);
}// end else

/************************************************************************************************/
/*********************************** Build Fia **************************************************/
/************************************************************************************************/
	// Build Fia
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('O'), ID('V'), "GF <O|V>");

	// Fia += 4 * \sum{m,e,n} <mi|ne> * G_mane
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,O]"), ints->DPD_ID("[O,V]"),
                  ints->DPD_ID("[O,O]"), ints->DPD_ID("[O,V]"), 0, "MO Ints <OO|OV>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
	global_dpd_->contract442(&K, &G, &GF, 1, 1, 4.0, 0.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

	// Fia += 8 * \sum{m,e,n} <nm|ie> * G_nmae
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,O]"), ints->DPD_ID("[O,V]"),
                  ints->DPD_ID("[O,O]"), ints->DPD_ID("[O,V]"), 0, "MO Ints <OO|OV>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
	global_dpd_->contract442(&K, &G, &GF, 2, 2, 8.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

if (wfn_type_ != "OMP2") {
   if (twopdm_abcd_type == "DIRECT" ) {
       	// Fia += \sum{m,n,c} X_mnic * (2t_mn^ac(1) - tmn^ca(1))
        global_dpd_->buf4_init(&X, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "X <OO|OV>");
        if (wfn_type_ == "OMP3" || wfn_type_ == "OMP2.5") {
           global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Tau_1 <OO|VV>");
        }
        else if (wfn_type_ == "OCEPA") {
           global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "Tau <OO|VV>");
        }
	global_dpd_->contract442(&X, &T, &GF, 2, 2, 1.0, 1.0);
	global_dpd_->buf4_close(&X);
	global_dpd_->buf4_close(&T);
	global_dpd_->file2_close(&GF);
   }

   else if (twopdm_abcd_type == "COMPUTE" ) {
       	// Fia += 4 * \sum{e,f,c} <ce|fi> * G_cefa = 4 * \sum{e,f,c} <if|ce> * G_afce
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV|VV>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");
	global_dpd_->contract442(&K, &G, &GF, 0, 0, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);
	global_dpd_->file2_close(&GF);
   }
}// end if (wfn_type_ != "OMP2")

	psio_->close(PSIF_LIBTRANS_DPD, 1);
        psio_->close(PSIF_OCC_DPD, 1);

/************************************************************************************************/
/*********************************** Load *******************************************************/
/************************************************************************************************/
	// Load Fai
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('V'), ID('O'), "GF <V|O>");
	global_dpd_->file2_mat_init(&GF);
	global_dpd_->file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
	  for(int a = 0 ; a < virtpiA[h]; ++a){
            for(int i = 0 ; i < occpiA[h]; ++i){
                GFock->add(h, a + occpiA[h], i, GF.matrix[h][a][i]);
            }
	  }
	}
	global_dpd_->file2_close(&GF);

	// Load Fia
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('O'), ID('V'), "GF <O|V>");
	global_dpd_->file2_mat_init(&GF);
	global_dpd_->file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < occpiA[h]; ++i){
            for(int a = 0 ; a < virtpiA[h]; ++a){
                GFock->add(h, i, a + occpiA[h], GF.matrix[h][i][a]);
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
/************************** Initialize ******************************************************/
/********************************************************************************************/
	GFockA->zero();
	GFockB->zero();

/********************************************************************************************/
/************************** 1e-part *********************************************************/
/********************************************************************************************/
	HG1A->zero();
	HG1B->zero();
	HG1A->gemm(false, false, 1.0, HmoA, g1symmA, 0.0);
	HG1B->gemm(false, false, 1.0, HmoB, g1symmB, 0.0);
	GFockA->add(HG1A);
	GFockB->add(HG1B);
	Ecc_rdm = HG1A->trace() + HG1B->trace() + Enuc; // One-electron contribution to MP2L
        Eopdm = Ecc_rdm;

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
        // With this algorithm cost changes to 3*o2v4 + 4*ov4 => 4*(o3v3 + o3v2), ~5-times faster
 	// X_MNIC = \sum{E,F} t_MN^EF * <IC||EF> = 2\sum{E,F} t_MN^EF * <IC|EF>
        global_dpd_->buf4_init(&X, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "X <OO|OV>");
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV|VV>");
        if (wfn_type_ == "OMP3" || wfn_type_ == "OMP2.5") {
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
        }
        else if (wfn_type_ == "OCEPA") {
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
        }
        global_dpd_->contract444(&T, &K, &X, 0, 0, 2.0, 0.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&T);
	global_dpd_->buf4_close(&X);

        // OMP2.5
        if (wfn_type_ == "OMP2.5") {
        global_dpd_->buf4_init(&X, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "X <OO|OV>");
        global_dpd_->buf4_scm(&X, 0.5);
	global_dpd_->buf4_close(&X);
        }

	// X_mnic = \sum{e,f} t_mn^ef * <ic||ef> = 2\sum{e,f} t_mn^ef * <ic|ef>
        global_dpd_->buf4_init(&X, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[o,v]"),
                  ID("[o,o]"), ID("[o,v]"), 0, "X <oo|ov>");
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <ov|vv>");
        if (wfn_type_ == "OMP3" || wfn_type_ == "OMP2.5") {
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
        }
        else if (wfn_type_ == "OCEPA") {
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2 <oo|vv>");
        }
        global_dpd_->contract444(&T, &K, &X, 0, 0, 2.0, 0.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&T);
	global_dpd_->buf4_close(&X);

        // OMP2.5
        if (wfn_type_ == "OMP2.5") {
        global_dpd_->buf4_init(&X, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[o,v]"),
                  ID("[o,o]"), ID("[o,v]"), 0, "X <oo|ov>");
        global_dpd_->buf4_scm(&X, 0.5);
	global_dpd_->buf4_close(&X);
        }

        // X_MnIc = \sum{E,f} t_Mn^Ef * <Ic|Ef>
        global_dpd_->buf4_init(&X, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[O,v]"),
                  ID("[O,o]"), ID("[O,v]"), 0, "X <Oo|Ov>");
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <Ov|Vv>");
        if (wfn_type_ == "OMP3" || wfn_type_ == "OMP2.5") {
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
        }
        else if (wfn_type_ == "OCEPA") {
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");
        }
        global_dpd_->contract444(&T, &K, &X, 0, 0, 1.0, 0.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&T);
	global_dpd_->buf4_close(&X);

        // OMP2.5
        if (wfn_type_ == "OMP2.5") {
        global_dpd_->buf4_init(&X, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[O,v]"),
                  ID("[O,o]"), ID("[O,v]"), 0, "X <Oo|Ov>");
        global_dpd_->buf4_scm(&X, 0.5);
	global_dpd_->buf4_close(&X);
        }

        // X_MnCi = \sum{E,f} t_Mn^Ef * <Ci|Ef>
        global_dpd_->buf4_init(&X, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[V,o]"),
                  ID("[O,o]"), ID("[V,o]"), 0, "X <Oo|Vo>");
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[V,v]"),
                  ID("[V,o]"), ID("[V,v]"), 0, "MO Ints <Vo|Vv>");
        if (wfn_type_ == "OMP3" || wfn_type_ == "OMP2.5") {
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
        }
        else if (wfn_type_ == "OCEPA") {
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2 <Oo|Vv>");
        }
        global_dpd_->contract444(&T, &K, &X, 0, 0, 1.0, 0.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&T);
	global_dpd_->buf4_close(&X);

        // OMP2.5
        if (wfn_type_ == "OMP2.5") {
        global_dpd_->buf4_init(&X, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[V,o]"),
                  ID("[O,o]"), ID("[V,o]"), 0, "X <Oo|Vo>");
        global_dpd_->buf4_scm(&X, 0.5);
	global_dpd_->buf4_close(&X);
        }

  }// end main if for X
} // end if (wfn_type_ != "OMP2") {

/********************************************************************************************/
/************************** VO-Block ********************************************************/
/********************************************************************************************/
	// Build FAJ
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('V'), ID('O'), "GF <V|O>");

	// FAJ = 2 * \sum{M,N,K} <MN||KA> * G_MNKJ
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "MO Ints <OO||OV>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
	global_dpd_->contract442(&K, &G, &GF, 3, 3, 2.0, 0.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

        /*
	// FAJ += 2 * \sum{E,M,F} <EF||MA> * G_MJEF = 2 * \sum{E,M,F} <MA||EF> * G_MJEF
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV||VV>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
	dpd_contract442(&K, &G, &GF, 1, 1, 2.0, 1.0);
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
        */

        // FAJ += 2 * \sum{E,M,F} <EF||MA> * G_MJEF = 2 * \sum{E,M,F} <MA||EF> * G_MJEF
        //      = 4 * \sum{E,M,F} <MA|EF> * G_MJEF
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV|VV>");
        global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
        global_dpd_->contract442(&K, &G, &GF, 1, 1, 4.0, 1.0);
        global_dpd_->buf4_close(&G);

        /*
	// FAJ += 4 * \sum{E,M,F} <EM||FA> * G_MEJF = 4 * \sum{E,M,F} <ME||AF> * G_MEJF
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV||VV>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
	dpd_contract442(&K, &G, &GF, 2, 2, 4.0, 1.0);
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
        */

     	// FAJ += 4 * \sum{E,M,F} <EM||FA> * G_MEJF = 4 * \sum{E,M,F} <ME||AF> * G_MEJF
     	//      = 4 * \sum{E,M,F} <ME|AF> * G_MEJF - 4 * \sum{E,M,F} <ME|FA> * G_MEJF
        //      = 4 * \sum{E,M,F} <ME|AF> * G_MEJF + 4 * \sum{E,M,F} <ME|FA> * G_MEFJ => 1st term
        global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
        global_dpd_->contract442(&K, &G, &GF, 2, 2, 4.0, 1.0);
        global_dpd_->buf4_close(&G);

	// FAJ += 4 * \sum{E,M,F} <EM||FA> * G_MEJF = 4 * \sum{E,M,F} <ME||AF> * G_MEJF
     	//      = 4 * \sum{E,M,F} <ME|AF> * G_MEJF - 4 * \sum{E,M,F} <ME|FA> * G_MEJF
        //      = 4 * \sum{E,M,F} <ME|AF> * G_MEJF + 4 * \sum{E,M,F} <ME|FA> * G_MEFJ => 2nd term
        global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[V,O]"),
                  ID("[O,V]"), ID("[V,O]"), 0, "TPDM <OV|VO>");
        global_dpd_->contract442(&K, &G, &GF, 3, 3, 4.0, 1.0);
        global_dpd_->buf4_close(&K);
        global_dpd_->buf4_close(&G);

	// FAJ += 4 * \sum{m,N,k} <Nm|Ak> * G_NmJk
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,o]"),
                  ID("[O,o]"), ID("[V,o]"), 0, "MO Ints <Oo|Vo>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "TPDM <Oo|Oo>");
	global_dpd_->contract442(&K, &G, &GF, 2, 2, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

	// FAJ += 4 * \sum{e,F,m} <Fe|Am> * G_JmFe = 4 * \sum{e,F,m} <Am|Fe> * G_JmFe
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[V,v]"),
                  ID("[V,o]"), ID("[V,v]"), 0, "MO Ints <Vo|Vv>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "TPDM <Oo|Vv>");
	global_dpd_->contract442(&K, &G, &GF, 0, 0, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

if (wfn_type_ != "OMP2") {
	// FAJ += 4 * \sum{E,f,m} <Em|Af> * G_EmJf  => this new
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[V,v]"),
                  ID("[V,o]"), ID("[V,v]"), 0, "MO Ints <Vo|Vv>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,o]"), ID("[O,v]"),
                  ID("[V,o]"), ID("[O,v]"), 0, "TPDM <Vo|Ov>");
	global_dpd_->contract442(&K, &G, &GF, 2, 2, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);
}// end if (wfn_type_ != "OMP2")

	// FAJ += 4 * \sum{e,f,M} <Me|Af> * G_MeJf
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <Ov|Vv>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "TPDM <Ov|Ov>");
	global_dpd_->contract442(&K, &G, &GF, 2, 2, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

	// Close
	global_dpd_->file2_close(&GF);


	// Build Faj
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('v'), ID('o'), "GF <v|o>");

	// Faj = 2 * \sum{m,n,k} <mn||ka> * G_mnkj
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,v]"),
                  ID("[o,o]"), ID("[o,v]"), 0, "MO Ints <oo||ov>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "TPDM <oo|oo>");
	global_dpd_->contract442(&K, &G, &GF, 3, 3, 2.0, 0.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

        /*
        // Faj += 2 * \sum{e,m,f} <ef||ma> * G_mjef = 2 * \sum{e,m,f} <ma||ef> * G_mjef
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <ov||vv>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "TPDM <oo|vv>");
	dpd_contract442(&K, &G, &GF, 1, 1, 2.0, 1.0);
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
        */

        // Faj += 2 * \sum{e,m,f} <ef||ma> * G_mjef = 2 * \sum{e,m,f} <ma||ef> * G_mjef
	//      = 4 * \sum{e,m,f} <ma|ef> * G_mjef
        global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <ov|vv>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "TPDM <oo|vv>");
	global_dpd_->contract442(&K, &G, &GF, 1, 1, 4.0, 1.0);
	global_dpd_->buf4_close(&G);

	/*
        // Faj += 4 * \sum{e,m,f} <em||fa> * G_mejf = 4 * \sum{e,m,f} <me||af> * G_mejf
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <ov||vv>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "TPDM <ov|ov>");
	dpd_contract442(&K, &G, &GF, 2, 2, 4.0, 1.0);
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	*/

        // Faj += 4 * \sum{e,m,f} <em||fa> * G_mejf = 4 * \sum{e,m,f} <me||af> * G_mejf
        //      = 4 * \sum{e,m,f} <me|af> * G_mejf - 4 * \sum{e,m,f} <me|fa> * G_mejf
        //      = 4 * \sum{e,m,f} <me|af> * G_mejf + 4 * \sum{e,m,f} <me|fa> * G_mefj => 1st term
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "TPDM <ov|ov>");
	global_dpd_->contract442(&K, &G, &GF, 2, 2, 4.0, 1.0);
	global_dpd_->buf4_close(&G);

        // Faj += 4 * \sum{e,m,f} <em||fa> * G_mejf = 4 * \sum{e,m,f} <me||af> * G_mejf
        //      = 4 * \sum{e,m,f} <me|af> * G_mejf - 4 * \sum{e,m,f} <me|fa> * G_mejf
        //      = 4 * \sum{e,m,f} <me|af> * G_mejf + 4 * \sum{e,m,f} <me|fa> * G_mefj => 2nd term
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,v]"), ID("[v,o]"),
                  ID("[o,v]"), ID("[v,o]"), 0, "TPDM <ov|vo>");
	global_dpd_->contract442(&K, &G, &GF, 3, 3, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

	// Faj += 4 * \sum{M,n,K} <Mn|Ka> * G_MnKj
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,v]"),
                  ID("[O,o]"), ID("[O,v]"), 0, "MO Ints <Oo|Ov>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "TPDM <Oo|Oo>");
	global_dpd_->contract442(&K, &G, &GF, 3, 3, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

	// Faj += 4 * \sum{E,f,M} <Ef|Ma> * G_MjEf = 4 * \sum{E,f,M} <Ma|Ef> * G_MjEf
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <Ov|Vv>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "TPDM <Oo|Vv>");
	global_dpd_->contract442(&K, &G, &GF, 1, 1, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

if (wfn_type_ != "OMP2") {
	// Faj += 4 * \sum{E,f,M} <Me|Fa> * G_MeFj => this new
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <Ov|Vv>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "TPDM <Ov|Vo>");
	global_dpd_->contract442(&K, &G, &GF, 3, 3, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);
}// end if (wfn_type_ != "OMP2")

	// Faj += 4 * \sum{E,F,m} <Em|Fa> * G_EmFj
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[V,v]"),
                  ID("[V,o]"), ID("[V,v]"), 0, "MO Ints <Vo|Vv>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,o]"), ID("[V,o]"),
                  ID("[V,o]"), ID("[V,o]"), 0, "TPDM <Vo|Vo>");
	global_dpd_->contract442(&K, &G, &GF, 3, 3, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

	// Close
	global_dpd_->file2_close(&GF);


/********************************************************************************************/
/************************** OV-Block ********************************************************/
/********************************************************************************************/
	// Build FIB
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('O'), ID('V'), "GF <O|V>");

	// FIB = 2 * \sum{M,N,E} <MN||EI> * G_MNEB = 2 * \sum{M,N,E} <MN||IE> * G_MNBE
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "MO Ints <OO||OV>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
	global_dpd_->contract442(&K, &G, &GF, 2, 2, 2.0, 0.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

        /*
	// FIB += 2 * \sum{E,F,C} <EF||CI> * G_EFCB = 2 * \sum{E,F,C} <IC||FE> * G_BCFE
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV||VV>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");
	dpd_contract442(&K, &G, &GF, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
        */

if (wfn_type_ != "OMP2") {
   if (twopdm_abcd_type == "DIRECT" ) {
       	// FIB += 2 * \sum{E,F,C} <EF||CI> * G_EFCB = 1/4  \sum{M,N,C} X_MNIC * t_MN^BC
        global_dpd_->buf4_init(&X, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "X <OO|OV>");
        if (wfn_type_ == "OMP3" || wfn_type_ == "OMP2.5") {
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
        }
        else if (wfn_type_ == "OCEPA") {
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2 <OO|VV>");
        }
	global_dpd_->contract442(&X, &T, &GF, 2, 2, 0.25, 1.0);
	global_dpd_->buf4_close(&X);
	global_dpd_->buf4_close(&T);
   }

   else if (twopdm_abcd_type == "COMPUTE" ) {
       	// FIB += 2 * \sum{E,F,C} <EF||CI> * G_EFCB = 2 * \sum{E,F,C} <IC||FE> * G_BCFE
       	//      = 4 * \sum{E,F,C} <IC|FE> * G_BCFE
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV|VV>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");
	global_dpd_->contract442(&K, &G, &GF, 0, 0, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);
   }
}// end if (wfn_type_ != "OMP2")

	// FIB += 4 * \sum{M,N,E} <ME||NI> * G_MENB = 4 * \sum{M,N,E} <NI||ME> * G_NBME
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "MO Ints <OO||OV>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
	global_dpd_->contract442(&K, &G, &GF, 1, 1, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);


if (wfn_type_ != "OMP2") {
   if (twopdm_abcd_type == "DIRECT" ) {
       	// FIB += 4 * \sum{e,F,c} <Fe|Ic> * G_FeBc =  \sum{M,n,C} X_MnIc * t_Mn^Bc
        global_dpd_->buf4_init(&X, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[O,v]"),
                  ID("[O,o]"), ID("[O,v]"), 0, "X <Oo|Ov>");
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
	// FIB += 4 * \sum{e,F,c} <Fe|Ic> * G_FeBc = 4 * \sum{e,F,c} <Ic|Fe> * G_BcFe  => this is new
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <Ov|Vv>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "TPDM <Vv|Vv>");
	global_dpd_->contract442(&K, &G, &GF, 0, 0, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);
   }
}// end if (wfn_type_ != "OMP2")

	// FIB += 4 * \sum{m,N,e} <Nm|Ie> * G_NmBe
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,v]"),
                  ID("[O,o]"), ID("[O,v]"), 0, "MO Ints <Oo|Ov>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "TPDM <Oo|Vv>");
	global_dpd_->contract442(&K, &G, &GF, 2, 2, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

	// FIB += 4 * \sum{m,n,E} <Em|In> * G_EmBn = 4 * \sum{m,n,E} <In|Em> * G_BnEm
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,o]"),
                 ID("[O,o]"), ID("[V,o]"), 0, "MO Ints <Oo|Vo>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,o]"), ID("[V,o]"),
                  ID("[V,o]"), ID("[V,o]"), 0, "TPDM <Vo|Vo>");
	global_dpd_->contract442(&K, &G, &GF, 0, 0, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

if (wfn_type_ != "OMP2") {
	// FIB += 4 * \sum{M,n,e} <Me|In> * G_MeBn = 4 * \sum{M,n,e} <In|Me> * G_BnMe  => this is new
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,v]"),
                  ID("[O,o]"), ID("[O,v]"), 0, "MO Ints <Oo|Ov>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,o]"), ID("[O,v]"),
                  ID("[V,o]"), ID("[O,v]"), 0, "TPDM <Vo|Ov>");
	global_dpd_->contract442(&K, &G, &GF, 0, 0, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);
}// end if (wfn_type_ != "OMP2")

	// Close
	global_dpd_->file2_close(&GF);



	// Build Fib
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('o'), ID('v'), "GF <o|v>");

	// Fib = 2 * \sum{m,n,e} <mn||ei> * G_mneb = 2 * \sum{m,n,e} <mn||ie> * G_mnbe
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,v]"),
                  ID("[o,o]"), ID("[o,v]"), 0, "MO Ints <oo||ov>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "TPDM <oo|vv>");
	global_dpd_->contract442(&K, &G, &GF, 2, 2, 2.0, 0.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

        /*
	// Fib += 2 * \sum{e,f,c} <ef||ci> * G_efcb = 2 * \sum{e,f,c} <ic||fe> * G_bcfe
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <ov||vv>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "TPDM <vv|vv>");
	dpd_contract442(&K, &G, &GF, 0, 0, 2.0, 1.0);
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
        */

if (wfn_type_ != "OMP2") {
   if (twopdm_abcd_type == "DIRECT" ) {
       	// Fib += 2 * \sum{e,f,c} <ef||ci> * G_efcb = 1/4  \sum{m,n,c} X_mnic * t_mn^bc
        global_dpd_->buf4_init(&X, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[o,v]"),
                  ID("[o,o]"), ID("[o,v]"), 0, "X <oo|ov>");
        if (wfn_type_ == "OMP3" || wfn_type_ == "OMP2.5") {
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
        }
        else if (wfn_type_ == "OCEPA") {
            global_dpd_->buf4_init(&T, PSIF_OCC_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2 <oo|vv>");
        }
	global_dpd_->contract442(&X, &T, &GF, 2, 2, 0.25, 1.0);
	global_dpd_->buf4_close(&X);
	global_dpd_->buf4_close(&T);
   }

   else if (twopdm_abcd_type == "COMPUTE" ) {
        // Fib += 2 * \sum{e,f,c} <ef||ci> * G_efcb = 2 * \sum{e,f,c} <ic||fe> * G_bcfe
        //      = 4 * \sum{e,f,c} <ic|fe> * G_bcfe
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <ov|vv>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "TPDM <vv|vv>");
	global_dpd_->contract442(&K, &G, &GF, 0, 0, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);
   }
}// end if (wfn_type_ != "OMP2")


	// Fib += 4 * \sum{m,n,e} <me||ni> * G_menb = 4 * \sum{m,n,e} <ni||me> * G_nbme
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,v]"),
                  ID("[o,o]"), ID("[o,v]"), 0, "MO Ints <oo||ov>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "TPDM <ov|ov>");
	global_dpd_->contract442(&K, &G, &GF, 1, 1, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

if (wfn_type_ != "OMP2") {
   if (twopdm_abcd_type == "DIRECT" ) {
       	// Fib += 4 * \sum{E,f,C} <Ef|Ci> * G_EfCb =  \sum{M,n,C} X_MnCi * t_Mn^Cb
        global_dpd_->buf4_init(&X, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[V,o]"),
                  ID("[O,o]"), ID("[V,o]"), 0, "X <Oo|Vo>");
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
	// Fib += 4 * \sum{E,f,C} <Ef|Ci> * G_EfCb = 4 * \sum{E,f,C} <Ci|Ef> * G_CbEf  => this is new
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[V,v]"),
                  ID("[V,o]"), ID("[V,v]"), 0, "MO Ints <Vo|Vv>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "TPDM <Vv|Vv>");
	global_dpd_->contract442(&K, &G, &GF, 1, 1, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);
   }
}// end if (wfn_type_ != "OMP2")

	// Fib += 4 * \sum{M,n,E} <Mn|Ei> * G_MnEb
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,o]"),
                  ID("[O,o]"), ID("[V,o]"), 0, "MO Ints <Oo|Vo>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "TPDM <Oo|Vv>");
	global_dpd_->contract442(&K, &G, &GF, 3, 3, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

	// Fib += 4 * \sum{M,N,e} <Me|Ni> * G_MeNb = 4 * \sum{M,N,e} <Ni|Me> * G_NbMe
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,v]"),
                 ID("[O,o]"), ID("[O,v]"), 0, "MO Ints <Oo|Ov>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "TPDM <Ov|Ov>");
	global_dpd_->contract442(&K, &G, &GF, 1, 1, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);

if (wfn_type_ != "OMP2") {
	// Fib += 4 * \sum{m,N,E} <Em|Ni> * G_EmNb = 4 * \sum{m,N,E} <Ni|Em> * G_NbEm   => this is new
	global_dpd_->buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,o]"),
                   ID("[O,o]"), ID("[V,o]"), 0, "MO Ints <Oo|Vo>");
	global_dpd_->buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "TPDM <Ov|Vo>");
	global_dpd_->contract442(&K, &G, &GF, 1, 1, 4.0, 1.0);
	global_dpd_->buf4_close(&K);
	global_dpd_->buf4_close(&G);
}// end if (wfn_type_ != "OMP2")

	// Close
	global_dpd_->file2_close(&GF);
	psio_->close(PSIF_LIBTRANS_DPD, 1);
        psio_->close(PSIF_OCC_DPD, 1);

/********************************************************************************************/
/************************** Load dpd_file2 to SharedMatrix (GFock) **************************/
/********************************************************************************************/
	// Load FAI
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('V'), ID('O'), "GF <V|O>");
	global_dpd_->file2_mat_init(&GF);
	global_dpd_->file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
	  for(int a = 0 ; a < virtpiA[h]; ++a){
            for(int i = 0 ; i < occpiA[h]; ++i){
                GFockA->add(h, a + occpiA[h], i, GF.matrix[h][a][i]);
            }
	  }
	}
	global_dpd_->file2_close(&GF);

	// Load Fai
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('v'), ID('o'), "GF <v|o>");
	global_dpd_->file2_mat_init(&GF);
	global_dpd_->file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
	  for(int a = 0 ; a < virtpiB[h]; ++a){
            for(int i = 0 ; i < occpiB[h]; ++i){
                GFockB->add(h, a + occpiB[h], i, GF.matrix[h][a][i]);
            }
	  }
	}
	global_dpd_->file2_close(&GF);

	// Load FIA
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('O'), ID('V'), "GF <O|V>");
	global_dpd_->file2_mat_init(&GF);
	global_dpd_->file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < occpiA[h]; ++i){
            for(int a = 0 ; a < virtpiA[h]; ++a){
                GFockA->add(h, i, a + occpiA[h], GF.matrix[h][i][a]);
            }
	  }
	}
	global_dpd_->file2_close(&GF);

	// Load Fia
	global_dpd_->file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('o'), ID('v'), "GF <o|v>");
	global_dpd_->file2_mat_init(&GF);
	global_dpd_->file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < occpiB[h]; ++i){
            for(int a = 0 ; a < virtpiB[h]; ++a){
                GFockB->add(h, i, a + occpiB[h], GF.matrix[h][i][a]);
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
//outfile->Printf("\n omp3_gfock done. \n");

} // End main
}} // End Namespaces
