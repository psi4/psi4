/** Standard library includes */
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>

 
/** Required PSI4 includes */
#include <psifiles.h>
#include <libciomr/libciomr.h>
#include <libpsio/psio.h> 
#include <libchkpt/chkpt.h>
#include <libpsio/psio.hpp>
#include <libchkpt/chkpt.hpp> 
#include <libiwl/iwl.h>
#include <libqt/qt.h>
#include <libtrans/mospace.h>
#include <libtrans/integraltransform.h>


/** Required libmints includes */
#include <libmints/mints.h>
#include <libmints/factory.h>
#include <libmints/wavefunction.h>

#include "occwave.h"
#include "defines.h"

using namespace boost;
using namespace psi;
using namespace std;

namespace psi{ namespace occwave{

void OCCWave::gfock_diag()
{
//fprintf(outfile,"\n gfock_diag is starting... \n"); fflush(outfile);
//===========================================================================================
//========================= RHF =============================================================
//===========================================================================================
if (reference_ == "RESTRICTED") {
        // 2e-part
	dpdbuf4 G, K;
	dpdfile2 GF;
	
	psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
	psio_->open(PSIF_OCC_DENSITY, PSIO_OPEN_OLD);

	// Build Fij
	dpd_file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('O'), ID('O'), "GF <O|O>");  
	
	// Fij += 4 * \sum{m,n,k} <km|ni> * G_kmnj
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,O]"), ints->DPD_ID("[O,O]"),
                  ints->DPD_ID("[O,O]"), ints->DPD_ID("[O,O]"), 0, "MO Ints <OO|OO>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
	dpd_contract442(&K, &G, &GF, 3, 3, 4.0, 0.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// Fij += 4 * \sum{e,m,f} <me|if> * G_mejf 
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,V]"), ints->DPD_ID("[O,V]"),
                  ints->DPD_ID("[O,V]"), ints->DPD_ID("[O,V]"), 0, "MO Ints <OV|OV>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
	dpd_contract442(&K, &G, &GF, 2, 2, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// Fij += 8 * \sum{e,m,f} <mi|ef> * G_mjef 
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,O]"), ints->DPD_ID("[V,V]"),
                  ints->DPD_ID("[O,O]"), ints->DPD_ID("[V,V]"), 0, "MO Ints <OO|VV>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
	dpd_contract442(&K, &G, &GF, 1, 1, 8.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	dpd_file2_close(&GF);

	// Build Fab
	dpd_file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "GF <V|V>"); 
	
	// Fab += 4 * \sum{m,e,n} <ma|ne> * G_mbne
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,V]"), ints->DPD_ID("[O,V]"),
                  ints->DPD_ID("[O,V]"), ints->DPD_ID("[O,V]"), 0, "MO Ints <OV|OV>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
	dpd_contract442(&K, &G, &GF, 1, 1, 4.0, 0.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// Fab += 8 * \sum{m,e,n} <nm|ae> * G_nmbe
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,O]"), ints->DPD_ID("[V,V]"),
                  ints->DPD_ID("[O,O]"), ints->DPD_ID("[V,V]"), 0, "MO Ints <OO|VV>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
	dpd_contract442(&K, &G, &GF, 2, 2, 8.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);

if (wfn_type_ == "OMP3" || wfn_type_ == "OCEPA") { 
	// Fab += 4 * \sum{c,e,f} <ce|fa> * G_cefb
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[V,V]"), ints->DPD_ID("[V,V]"),
                  ints->DPD_ID("[V,V]"), ints->DPD_ID("[V,V]"), 0, "MO Ints <VV|VV>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");
	dpd_contract442(&K, &G, &GF, 3, 3, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	dpd_file2_close(&GF);
}// end if (wfn_type_ == "OMP3" || wfn_type_ == "OCEPA") { 

        // close the integral file
	psio_->close(PSIF_LIBTRANS_DPD, 1);

	// Load dpd_file2 to SharedMatrix (GFock)
	// Load Fij
	dpd_file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('O'), ID('O'), "GF <O|O>");  
	dpd_file2_mat_init(&GF);
	dpd_file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < occpiA[h]; ++i){
            for(int j = 0 ; j < occpiA[h]; ++j){
                GFock->add(h, i, j, GF.matrix[h][i][j]);
            }
	  }
	}
	dpd_file2_close(&GF);

        // Load Fab
	dpd_file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "GF <V|V>");  
	dpd_file2_mat_init(&GF);
	dpd_file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
	  for(int a = 0 ; a < virtpiA[h]; ++a){
            for(int b = 0 ; b < virtpiA[h]; ++b){
                GFock->add(h, a + occpiA[h], b + occpiA[h], GF.matrix[h][a][b]);
            }
	  }
	}
	dpd_file2_close(&GF);

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
	dpdbuf4 G, K;
	dpdfile2 GF;
	
	psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
	psio_->open(PSIF_OCC_DENSITY, PSIO_OPEN_OLD);
	
/********************************************************************************************/
/************************** OO-Block ********************************************************/
/********************************************************************************************/ 
	// Build FIJ
	dpd_file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('O'), ID('O'), "GF <O|O>");  
	
	// FIJ = 2 * \sum{M,N,K} <MN||KI> * G_MNKJ
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "MO Ints <OO||OO>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
	dpd_contract442(&K, &G, &GF, 3, 3, 2.0, 0.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// FIJ += 2 * \sum{E,F,M} <EF||MI> * G_MJEF = 2 * \sum{E,F,M} <MI||EF> * G_MJEF 
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
	dpd_contract442(&K, &G, &GF, 1, 1, 2.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// FIJ += 4 * \sum{E,F,M} <EM||FI> * G_MEJF = 4 * \sum{E,F,M} <ME||IF> * G_MEJF 
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV||OV>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
	dpd_contract442(&K, &G, &GF, 2, 2, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// FIJ += 4 * \sum{m,N,k} <Nm|Ik> * G_NmJk  
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "MO Ints <Oo|Oo>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "TPDM <Oo|Oo>");
	dpd_contract442(&K, &G, &GF, 2, 2, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// FIJ += 4 * \sum{e,F,m} <Fe|Im> * G_JmFe = 4 * \sum{e,F,m} <Im|Fe> * G_JmFe
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "TPDM <Oo|Vv>");
	dpd_contract442(&K, &G, &GF, 0, 0, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// FIJ += 4 * \sum{e,f,M} <Me|If> * G_MeJf 
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "MO Ints <Ov|Ov>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "TPDM <Ov|Ov>");
	dpd_contract442(&K, &G, &GF, 2, 2, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
if (wfn_type_ == "OMP3" || wfn_type_ == "OCEPA") { 
        // FIJ += 4 * \sum{E,F,m} <Em|If> * G_EmJf = 4 * \sum{E,F,m} <If|Em> * G_JfEm => new 
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "MO Ints <Ov|Vo>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "TPDM <Ov|Vo>");
	dpd_contract442(&K, &G, &GF, 0, 0, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
}// end if (wfn_type_ == "OMP3" || wfn_type_ == "OCEPA") { 
	
	// Close 
	dpd_file2_close(&GF);


	// Build Fij
	dpd_file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('o'), ID('o'), "GF <o|o>");  
	
	// Fij = 2 * \sum{m,n,k} <mn||ki> * G_mnkj
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "MO Ints <oo||oo>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "TPDM <oo|oo>");
	dpd_contract442(&K, &G, &GF, 3, 3, 2.0, 0.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// Fij += 2 * \sum{e,f,m} <ef||mi> * G_mjef = 2 * \sum{e,f,m} <mi||ef> * G_mjef 
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "TPDM <oo|vv>");
	dpd_contract442(&K, &G, &GF, 1, 1, 2.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// Fij += 4 * \sum{e,f,m} <em||fi> * G_mejf = 4 * \sum{e,f,m} <me||if> * G_mejf 
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov||ov>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "TPDM <ov|ov>");
	dpd_contract442(&K, &G, &GF, 2, 2, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// Fij += 4 * \sum{M,n,K} <Mn|Ki> * G_MnKj 
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "MO Ints <Oo|Oo>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "TPDM <Oo|Oo>");
	dpd_contract442(&K, &G, &GF, 3, 3, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// Fij += 4 * \sum{E,f,M} <Ef|Mi> * G_MjEf = 4 * \sum{E,f,M} <Mi|Ef> * G_MjEf
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "TPDM <Oo|Vv>");
	dpd_contract442(&K, &G, &GF, 1, 1, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// Fij += 4 * \sum{E,F,m} <Em|Fi> * G_EmFj 
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[V,o]"),
                  ID("[V,o]"), ID("[V,o]"), 0, "MO Ints <Vo|Vo>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,o]"), ID("[V,o]"),
                  ID("[V,o]"), ID("[V,o]"), 0, "TPDM <Vo|Vo>");
	dpd_contract442(&K, &G, &GF, 3, 3, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
if (wfn_type_ == "OMP3" || wfn_type_ == "OCEPA") { 
	// Fij += 4 * \sum{e,F,M} <Me|Fi> * G_MeFj => new 
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "MO Ints <Ov|Vo>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "TPDM <Ov|Vo>");
	dpd_contract442(&K, &G, &GF, 3, 3, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
}// end if (wfn_type_ == "OMP3" || wfn_type_ == "OCEPA") { 

	// Close 
	dpd_file2_close(&GF);
	
/********************************************************************************************/
/************************** VV-Block ********************************************************/
/********************************************************************************************/ 
	// Build FAB
	dpd_file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "GF <V|V>"); 
	
	// FAB = 2 * \sum{M,N,E} <MN||EA> * G_MNEB
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "MO Ints <OO||VV>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
	dpd_contract442(&K, &G, &GF, 3, 3, 2.0, 0.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);

if (wfn_type_ == "OMP3" || wfn_type_ == "OCEPA") { 
	// FAB = 2 * \sum{E,F,C} <EF||CA> * G_EFCB => new
	// FAB = 4 * \sum{E,F,C} <EF|CA> * G_EFCB => new
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[V,V]"), ints->DPD_ID("[V,V]"),
                  ints->DPD_ID("[V,V]"), ints->DPD_ID("[V,V]"), 0, "MO Ints <VV|VV>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,V]"), ID("[V,V]"),
                  ID("[V,V]"), ID("[V,V]"), 0, "TPDM <VV|VV>");
	dpd_contract442(&K, &G, &GF, 3, 3, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
}// end if (wfn_type_ == "OMP3" || wfn_type_ == "OCEPA") { 
	
	// FAB += 4 * \sum{M,N,E} <ME||NA> * G_MENB
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "MO Ints <OV||OV>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
	dpd_contract442(&K, &G, &GF, 3, 3, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// FAB = 4 * \sum{m,N,e} <Nm|Ae> * G_NmBe
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "TPDM <Oo|Vv>");
	dpd_contract442(&K, &G, &GF, 2, 2, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// FAB = 4 * \sum{m,n,E} <Em|An> * G_EmBn
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[V,o]"),
                  ID("[V,o]"), ID("[V,o]"), 0, "MO Ints <Vo|Vo>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,o]"), ID("[V,o]"),
                  ID("[V,o]"), ID("[V,o]"), 0, "TPDM <Vo|Vo>");
	dpd_contract442(&K, &G, &GF, 2, 2, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
if (wfn_type_ == "OMP3" || wfn_type_ == "OCEPA") { 
	// FAB = 4 * \sum{M,n,e} <Me|An> * G_MeBn => new
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "MO Ints <Ov|Vo>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "TPDM <Ov|Vo>");
	dpd_contract442(&K, &G, &GF, 2, 2, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// FAB = 4 * \sum{e,F,c} <Fe|Ac> * G_FeBc => new
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[V,v]"), ints->DPD_ID("[V,v]"),
                  ints->DPD_ID("[V,v]"), ints->DPD_ID("[V,v]"), 0, "MO Ints <Vv|Vv>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "TPDM <Vv|Vv>");
	dpd_contract442(&K, &G, &GF, 2, 2, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
}// end if (wfn_type_ == "OMP3" || wfn_type_ == "OCEPA") { 

	// Close
	dpd_file2_close(&GF);
	
	
	// Build Fab
	dpd_file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('v'), ID('v'), "GF <v|v>"); 
	
	// Fab = 2 * \sum{m,n,e} <mn||ea> * G_mneb
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "MO Ints <oo||vv>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "TPDM <oo|vv>");
	dpd_contract442(&K, &G, &GF, 3, 3, 2.0, 0.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);

if (wfn_type_ == "OMP3" || wfn_type_ == "OCEPA") { 
	// Fab = 2 * \sum{efc} <ef||ca> * G_efcb => new
	// Fab = 4 * \sum{efc} <ef|ca> * G_efcb => new
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[v,v]"), ints->DPD_ID("[v,v]"),
                  ints->DPD_ID("[v,v]"), ints->DPD_ID("[v,v]"), 0, "MO Ints <vv|vv>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[v,v]"), ID("[v,v]"),
                  ID("[v,v]"), ID("[v,v]"), 0, "TPDM <vv|vv>");
	dpd_contract442(&K, &G, &GF, 3, 3, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
}// end if (wfn_type_ == "OMP3" || wfn_type_ == "OCEPA") { 
	
	// Fab += 4 * \sum{m,n,e} <me||na> * G_menb
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "MO Ints <ov||ov>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "TPDM <ov|ov>");
	dpd_contract442(&K, &G, &GF, 3, 3, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// Fab = 4 * \sum{M,n,E} <Mn|Ea> * G_MnEb
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "MO Ints <Oo|Vv>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "TPDM <Oo|Vv>");
	dpd_contract442(&K, &G, &GF, 3, 3, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// Fab = 4 * \sum{M,N,e} <Me|Na> * G_MeNb
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "MO Ints <Ov|Ov>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "TPDM <Ov|Ov>");
	dpd_contract442(&K, &G, &GF, 3, 3, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
if (wfn_type_ == "OMP3" || wfn_type_ == "OCEPA") { 
	// Fab = 4 * \sum{m,N,E} <Em|Na> * G_EmNb = 4 * \sum{M,n,e} <Na|Em> * G_NbEm => new
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "MO Ints <Ov|Vo>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[O,v]"), ID("[V,o]"),
                  ID("[O,v]"), ID("[V,o]"), 0, "TPDM <Ov|Vo>");
	dpd_contract442(&K, &G, &GF, 1, 1, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);

	// Fab = 4 * \sum{e,F,c} <Ef|Ca> * G_EfCb => new
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[V,v]"), ints->DPD_ID("[V,v]"),
                  ints->DPD_ID("[V,v]"), ints->DPD_ID("[V,v]"), 0, "MO Ints <Vv|Vv>");
	dpd_buf4_init(&G, PSIF_OCC_DENSITY, 0, ID("[V,v]"), ID("[V,v]"),
                  ID("[V,v]"), ID("[V,v]"), 0, "TPDM <Vv|Vv>");
	dpd_contract442(&K, &G, &GF, 3, 3, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
}// end if (wfn_type_ == "OMP3" || wfn_type_ == "OCEPA") { 

	// Close
	dpd_file2_close(&GF);
	psio_->close(PSIF_LIBTRANS_DPD, 1);

/********************************************************************************************/
/************************** Load dpd_file2 to SharedMatrix (GFock) **************************/
/********************************************************************************************/ 
	// Load FIJ
	dpd_file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('O'), ID('O'), "GF <O|O>");  
	dpd_file2_mat_init(&GF);
	dpd_file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < occpiA[h]; ++i){
            for(int j = 0 ; j < occpiA[h]; ++j){
                GFockA->add(h, i, j, GF.matrix[h][i][j]);
            }
	  }
	}
	dpd_file2_close(&GF);
	
	// Load Fij
	dpd_file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('o'), ID('o'), "GF <o|o>");  
	dpd_file2_mat_init(&GF);
	dpd_file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < occpiB[h]; ++i){
            for(int j = 0 ; j < occpiB[h]; ++j){
                GFockB->add(h, i, j, GF.matrix[h][i][j]);
            }
	  }
	}
	dpd_file2_close(&GF);
	
	// Load FAB
	dpd_file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('V'), ID('V'), "GF <V|V>");  
	dpd_file2_mat_init(&GF);
	dpd_file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
	  for(int a = 0 ; a < virtpiA[h]; ++a){
            for(int b = 0 ; b < virtpiA[h]; ++b){
                GFockA->add(h, a + occpiA[h], b + occpiA[h], GF.matrix[h][a][b]);
            }
	  }
	}
	dpd_file2_close(&GF);
	
	// Load Fab
	dpd_file2_init(&GF, PSIF_OCC_DENSITY, 0, ID('v'), ID('v'), "GF <v|v>");  
	dpd_file2_mat_init(&GF);
	dpd_file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
	  for(int a = 0 ; a < virtpiB[h]; ++a){
            for(int b = 0 ; b < virtpiB[h]; ++b){
                GFockB->add(h, a + occpiB[h], b + occpiB[h], GF.matrix[h][a][b]);
            }
	  }
	}
	dpd_file2_close(&GF);
	psio_->close(PSIF_OCC_DENSITY, 1);

        // Print
	if (print_ > 1) {
	  GFockA->print();
	  GFockB->print();
	}
	
}// end if (reference_ == "UNRESTRICTED") 
//fprintf(outfile,"\n gfock_diag done. \n"); fflush(outfile);

}// end gfock_diag 
}} // End Namespaces

