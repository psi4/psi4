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

#include "omp2wave.h"
#include "defines.h"

using namespace boost;
using namespace psi;
using namespace std;

namespace psi{ namespace omp2wave{

void OMP2Wave::GFockmo()
{
//fprintf(outfile,"\n GFockmo is starting... \n"); fflush(outfile);
//===========================================================================================
//========================= RHF =============================================================
//===========================================================================================
if (reference == "RHF") {

	// Initialize
	GFock->zero();

        // 1e-part
	HG1->zero();
	HG1->gemm(false, false, 1.0, HmoA, g1symm, 0.0);
	GFock->add(HG1);
	Emp2_rdm=HG1->trace() + Enuc; // One-electron contribution to MP2L

        // 2e-part
	dpdbuf4 G, K;
	dpdfile2 GF;
	
	psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
	psio_->open(PSIF_OMP2_DENSITY, PSIO_OPEN_OLD);
	
	// Build Fai
	dpd_file2_init(&GF, PSIF_OMP2_DENSITY, 0, ID('V'), ID('O'), "GF <V|O>"); 
	
	// Fai += 4 * \sum{m,n,k} <km|na> * G_kmni
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,O]"), ints->DPD_ID("[O,V]"),
                  ints->DPD_ID("[O,O]"), ints->DPD_ID("[O,V]"), 0, "MO Ints <OO|OV>");
	dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
	dpd_contract442(&K, &G, &GF, 3, 3, 4.0, 0.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// Fai += 4 * \sum{e,m,f} <me|af> * G_meif 
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,V]"), ints->DPD_ID("[V,V]"),
                  ints->DPD_ID("[O,V]"), ints->DPD_ID("[V,V]"), 0, "MO Ints <OV|VV>");
	dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
	dpd_contract442(&K, &G, &GF, 2, 2, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// Fai += 8 * \sum{e,m,f} <ma|ef> * G_mief 
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,V]"), ints->DPD_ID("[V,V]"),
                  ints->DPD_ID("[O,V]"), ints->DPD_ID("[V,V]"), 0, "MO Ints <OV|VV>");
	dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
	dpd_contract442(&K, &G, &GF, 1, 1, 8.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	dpd_file2_close(&GF);
	
	
	// Build Fia
	dpd_file2_init(&GF, PSIF_OMP2_DENSITY, 0, ID('O'), ID('V'), "GF <O|V>");  
	
	// Fia += 4 * \sum{m,e,n} <mi|ne> * G_mane
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,O]"), ints->DPD_ID("[O,V]"),
                  ints->DPD_ID("[O,O]"), ints->DPD_ID("[O,V]"), 0, "MO Ints <OO|OV>");
	dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
	dpd_contract442(&K, &G, &GF, 1, 1, 4.0, 0.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// Fia += 8 * \sum{m,e,n} <nm|ie> * G_nmae
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ints->DPD_ID("[O,O]"), ints->DPD_ID("[O,V]"),
                  ints->DPD_ID("[O,O]"), ints->DPD_ID("[O,V]"), 0, "MO Ints <OO|OV>");
	dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
	dpd_contract442(&K, &G, &GF, 2, 2, 8.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	dpd_file2_close(&GF);

        // close the integral file
	psio_->close(PSIF_LIBTRANS_DPD, 1);


	// Load dpd_file2 to SharedMatrix (GFock)
	// Load Fai
	dpd_file2_init(&GF, PSIF_OMP2_DENSITY, 0, ID('V'), ID('O'), "GF <V|O>");  
	dpd_file2_mat_init(&GF);
	dpd_file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
	  for(int a = 0 ; a < virtpiA[h]; ++a){
            for(int i = 0 ; i < occpiA[h]; ++i){
                GFock->add(h, a + occpiA[h], i, GF.matrix[h][a][i]);
            }
	  }
	}
	dpd_file2_close(&GF);
	
	// Load Fia
	dpd_file2_init(&GF, PSIF_OMP2_DENSITY, 0, ID('O'), ID('V'), "GF <O|V>");  
	dpd_file2_mat_init(&GF);
	dpd_file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < occpiA[h]; ++i){
            for(int a = 0 ; a < virtpiA[h]; ++a){
                GFock->add(h, i, a + occpiA[h], GF.matrix[h][i][a]);
            }
	  }
	}
	dpd_file2_close(&GF);

	psio_->close(PSIF_OMP2_DENSITY, 1);
	if (print_ > 2) GFock->print();
	
}// end if (reference == "RHF") 



//===========================================================================================
//========================= UHF =============================================================
//===========================================================================================
else if (reference == "UHF") {

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
	Emp2_rdm = HG1A->trace() + HG1B->trace() + Enuc; // One-electron contribution to MP2L

/********************************************************************************************/
/************************** 2e-part *********************************************************/
/********************************************************************************************/ 
	dpdbuf4 G, K;
	dpdfile2 GF;  
	
	psio_->open(PSIF_LIBTRANS_DPD, PSIO_OPEN_OLD);
	psio_->open(PSIF_OMP2_DENSITY, PSIO_OPEN_OLD);
	
/********************************************************************************************/
/************************** VO-Block ********************************************************/
/********************************************************************************************/ 
	// Build FAJ
	dpd_file2_init(&GF, PSIF_OMP2_DENSITY, 0, ID('V'), ID('O'), "GF <V|O>"); 
	
	// FAJ = 2 * \sum{M,N,K} <MN||KA> * G_MNKJ
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "MO Ints <OO||OV>");
	dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
	dpd_contract442(&K, &G, &GF, 3, 3, 2.0, 0.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);

        /*	
	// FAJ += 2 * \sum{E,M,F} <EF||MA> * G_MJEF = 2 * \sum{E,M,F} <MA||EF> * G_MJEF 
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV||VV>");
	dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
	dpd_contract442(&K, &G, &GF, 1, 1, 2.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
        */

        // FAJ += 2 * \sum{E,M,F} <EF||MA> * G_MJEF = 2 * \sum{E,M,F} <MA||EF> * G_MJEF 
        //      = 4 * \sum{E,M,F} <MA|EF> * G_MJEF         
        dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV|VV>");
        dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
        dpd_contract442(&K, &G, &GF, 1, 1, 4.0, 1.0);
        dpd_buf4_close(&G);
 
	
        /*
	// FAJ += 4 * \sum{E,M,F} <EM||FA> * G_MEJF = 4 * \sum{E,M,F} <ME||AF> * G_MEJF
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,V]"), ID("[V,V]"),
                  ID("[O,V]"), ID("[V,V]"), 0, "MO Ints <OV||VV>");
	dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
	dpd_contract442(&K, &G, &GF, 2, 2, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
        */

      	// FAJ += 4 * \sum{E,M,F} <EM||FA> * G_MEJF = 4 * \sum{E,M,F} <ME||AF> * G_MEJF
     	//      = 4 * \sum{E,M,F} <ME|AF> * G_MEJF - 4 * \sum{E,M,F} <ME|FA> * G_MEJF 
        //      = 4 * \sum{E,M,F} <ME|AF> * G_MEJF + 4 * \sum{E,M,F} <ME|FA> * G_MEFJ => 1st term
        dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
        dpd_contract442(&K, &G, &GF, 2, 2, 4.0, 1.0);
        dpd_buf4_close(&G);
 
	// FAJ += 4 * \sum{E,M,F} <EM||FA> * G_MEJF = 4 * \sum{E,M,F} <ME||AF> * G_MEJF
     	//      = 4 * \sum{E,M,F} <ME|AF> * G_MEJF - 4 * \sum{E,M,F} <ME|FA> * G_MEJF 
        //      = 4 * \sum{E,M,F} <ME|AF> * G_MEJF + 4 * \sum{E,M,F} <ME|FA> * G_MEFJ => 2nd term
        dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,V]"), ID("[V,O]"),
                  ID("[O,V]"), ID("[V,O]"), 0, "TPDM <OV|VO>");
        dpd_contract442(&K, &G, &GF, 3, 3, 4.0, 1.0);
        dpd_buf4_close(&K);
        dpd_buf4_close(&G);

	
	// FAJ += 4 * \sum{m,N,k} <Nm|Ak> * G_NmJk 
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,o]"),
                  ID("[O,o]"), ID("[V,o]"), 0, "MO Ints <Oo|Vo>");
	dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "TPDM <Oo|Oo>");
	dpd_contract442(&K, &G, &GF, 2, 2, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// FAJ += 4 * \sum{e,F,m} <Fe|Am> * G_JmFe = 4 * \sum{e,F,m} <Am|Fe> * G_JmFe
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[V,v]"),
                  ID("[V,o]"), ID("[V,v]"), 0, "MO Ints <Vo|Vv>");
	dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "TPDM <Oo|Vv>");
	dpd_contract442(&K, &G, &GF, 0, 0, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// FAJ += 4 * \sum{e,f,M} <Me|Af> * G_MeJf
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <Ov|Vv>");
	dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "TPDM <Ov|Ov>");
	dpd_contract442(&K, &G, &GF, 2, 2, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);

	// Close
	dpd_file2_close(&GF);
	
	
	// Build Faj
	dpd_file2_init(&GF, PSIF_OMP2_DENSITY, 0, ID('v'), ID('o'), "GF <v|o>"); 
	
	// Faj = 2 * \sum{m,n,k} <mn||ka> * G_mnkj
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,v]"),
                  ID("[o,o]"), ID("[o,v]"), 0, "MO Ints <oo||ov>");
	dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "TPDM <oo|oo>");
	dpd_contract442(&K, &G, &GF, 3, 3, 2.0, 0.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);

        /*	
	// Faj += 2 * \sum{e,m,f} <ef||ma> * G_mjef = 2 * \sum{e,m,f} <ma||ef> * G_mjef 
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <ov||vv>");
	dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "TPDM <oo|vv>");
	dpd_contract442(&K, &G, &GF, 1, 1, 2.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
        */

        // Faj += 2 * \sum{e,m,f} <ef||ma> * G_mjef = 2 * \sum{e,m,f} <ma||ef> * G_mjef 
	//      = 4 * \sum{e,m,f} <ma|ef> * G_mjef 
        dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <ov|vv>");
	dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "TPDM <oo|vv>");
	dpd_contract442(&K, &G, &GF, 1, 1, 4.0, 1.0); 
	dpd_buf4_close(&G);

        /*	
	// Faj += 4 * \sum{e,m,f} <em||fa> * G_mejf = 4 * \sum{e,m,f} <me||af> * G_mejf
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,v]"), ID("[v,v]"),
                  ID("[o,v]"), ID("[v,v]"), 0, "MO Ints <ov||vv>");
	dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "TPDM <ov|ov>");
	dpd_contract442(&K, &G, &GF, 2, 2, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
        */

        // Faj += 4 * \sum{e,m,f} <em||fa> * G_mejf = 4 * \sum{e,m,f} <me||af> * G_mejf
        //      = 4 * \sum{e,m,f} <me|af> * G_mejf - 4 * \sum{e,m,f} <me|fa> * G_mejf 
        //      = 4 * \sum{e,m,f} <me|af> * G_mejf + 4 * \sum{e,m,f} <me|fa> * G_mefj => 1st term
	dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "TPDM <ov|ov>");
	dpd_contract442(&K, &G, &GF, 2, 2, 4.0, 1.0); 
	dpd_buf4_close(&G);
 
        // Faj += 4 * \sum{e,m,f} <em||fa> * G_mejf = 4 * \sum{e,m,f} <me||af> * G_mejf
        //      = 4 * \sum{e,m,f} <me|af> * G_mejf - 4 * \sum{e,m,f} <me|fa> * G_mejf 
        //      = 4 * \sum{e,m,f} <me|af> * G_mejf + 4 * \sum{e,m,f} <me|fa> * G_mefj => 2nd term
	dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[o,v]"), ID("[v,o]"),
                  ID("[o,v]"), ID("[v,o]"), 0, "TPDM <ov|vo>");
	dpd_contract442(&K, &G, &GF, 3, 3, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);

	
	// Faj += 4 * \sum{M,n,K} <Mn|Ka> * G_MnKj 
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,v]"),
                  ID("[O,o]"), ID("[O,v]"), 0, "MO Ints <Oo|Ov>");
	dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                  ID("[O,o]"), ID("[O,o]"), 0, "TPDM <Oo|Oo>");
	dpd_contract442(&K, &G, &GF, 3, 3, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// Faj += 4 * \sum{E,f,M} <Ef|Ma> * G_MjEf = 4 * \sum{E,f,M} <Ma|Ef> * G_MjEf 
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,v]"), ID("[V,v]"),
                  ID("[O,v]"), ID("[V,v]"), 0, "MO Ints <Ov|Vv>");
	dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "TPDM <Oo|Vv>");
	dpd_contract442(&K, &G, &GF, 1, 1, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// Faj += 4 * \sum{E,F,m} <Em|Fa> * G_EmFj
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[V,o]"), ID("[V,v]"),
                  ID("[V,o]"), ID("[V,v]"), 0, "MO Ints <Vo|Vv>");
	dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[V,o]"), ID("[V,o]"),
                  ID("[V,o]"), ID("[V,o]"), 0, "TPDM <Vo|Vo>");
	dpd_contract442(&K, &G, &GF, 3, 3, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);

	// Close
	dpd_file2_close(&GF);
	
	
/********************************************************************************************/
/************************** OV-Block ********************************************************/
/********************************************************************************************/ 
	// Build FIB
	dpd_file2_init(&GF, PSIF_OMP2_DENSITY, 0, ID('O'), ID('V'), "GF <O|V>");  
	
	// FIB = 2 * \sum{M,N,E} <MN||EI> * G_MNEB = 2 * \sum{M,N,E} <MN||IE> * G_MNBE
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "MO Ints <OO||OV>");
	dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
	dpd_contract442(&K, &G, &GF, 2, 2, 2.0, 0.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// FIB += 4 * \sum{M,N,E} <ME||NI> * G_MENB = 4 * \sum{M,N,E} <NI||ME> * G_NBME
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,O]"), ID("[O,V]"),
                  ID("[O,O]"), ID("[O,V]"), 0, "MO Ints <OO||OV>");
	dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,V]"), ID("[O,V]"),
                  ID("[O,V]"), ID("[O,V]"), 0, "TPDM <OV|OV>");
	dpd_contract442(&K, &G, &GF, 1, 1, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// FIB += 4 * \sum{m,N,e} <Nm|Ie> * G_NmBe 
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,v]"),
                  ID("[O,o]"), ID("[O,v]"), 0, "MO Ints <Oo|Ov>");
	dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "TPDM <Oo|Vv>");
	dpd_contract442(&K, &G, &GF, 2, 2, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// FIB += 4 * \sum{m,n,E} <Em|In> * G_EmBn = 4 * \sum{m,n,E} <In|Em> * G_BnEm
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,o]"),
                 ID("[O,o]"), ID("[V,o]"), 0, "MO Ints <Oo|Vo>");
	dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[V,o]"), ID("[V,o]"),
                  ID("[V,o]"), ID("[V,o]"), 0, "TPDM <Vo|Vo>");
	dpd_contract442(&K, &G, &GF, 0, 0, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// Close
	dpd_file2_close(&GF);
	
	
	// Build Fib
	dpd_file2_init(&GF, PSIF_OMP2_DENSITY, 0, ID('o'), ID('v'), "GF <o|v>");  
	
	// Fib = 2 * \sum{m,n,e} <mn||ei> * G_mneb = 2 * \sum{m,n,e} <mn||ie> * G_mnbe
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,v]"),
                  ID("[o,o]"), ID("[o,v]"), 0, "MO Ints <oo||ov>");
	dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "TPDM <oo|vv>");
	dpd_contract442(&K, &G, &GF, 2, 2, 2.0, 0.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// Fib += 4 * \sum{m,n,e} <me||ni> * G_menb = 4 * \sum{m,n,e} <ni||me> * G_nbme
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[o,o]"), ID("[o,v]"),
                  ID("[o,o]"), ID("[o,v]"), 0, "MO Ints <oo||ov>");
	dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[o,v]"), ID("[o,v]"),
                  ID("[o,v]"), ID("[o,v]"), 0, "TPDM <ov|ov>");
	dpd_contract442(&K, &G, &GF, 1, 1, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// Fib += 4 * \sum{M,n,E} <Mn|Ei> * G_MnEb 
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[V,o]"),
                  ID("[O,o]"), ID("[V,o]"), 0, "MO Ints <Oo|Vo>");
	dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "TPDM <Oo|Vv>");
	dpd_contract442(&K, &G, &GF, 3, 3, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// Fib += 4 * \sum{M,N,e} <Me|Ni> * G_MeNb = 4 * \sum{M,N,e} <Ni|Me> * G_NbMe
	dpd_buf4_init(&K, PSIF_LIBTRANS_DPD, 0, ID("[O,o]"), ID("[O,v]"),
                 ID("[O,o]"), ID("[O,v]"), 0, "MO Ints <Oo|Ov>");
	dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,v]"), ID("[O,v]"),
                  ID("[O,v]"), ID("[O,v]"), 0, "TPDM <Ov|Ov>");
	dpd_contract442(&K, &G, &GF, 1, 1, 4.0, 1.0); 
	dpd_buf4_close(&K);
	dpd_buf4_close(&G);
	
	// Close
	dpd_file2_close(&GF);

        // close the integral file
	psio_->close(PSIF_LIBTRANS_DPD, 1);
	
/********************************************************************************************/
/************************** Load dpd_file2 to SharedMatrix (GFock) **************************/
/********************************************************************************************/ 
	// Load FAI
	dpd_file2_init(&GF, PSIF_OMP2_DENSITY, 0, ID('V'), ID('O'), "GF <V|O>");  
	dpd_file2_mat_init(&GF);
	dpd_file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
	  for(int a = 0 ; a < virtpiA[h]; ++a){
            for(int i = 0 ; i < occpiA[h]; ++i){
                GFockA->add(h, a + occpiA[h], i, GF.matrix[h][a][i]);
            }
	  }
	}
	dpd_file2_close(&GF);
	
	// Load Fai
	dpd_file2_init(&GF, PSIF_OMP2_DENSITY, 0, ID('v'), ID('o'), "GF <v|o>");  
	dpd_file2_mat_init(&GF);
	dpd_file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
	  for(int a = 0 ; a < virtpiB[h]; ++a){
            for(int i = 0 ; i < occpiB[h]; ++i){
                GFockB->add(h, a + occpiB[h], i, GF.matrix[h][a][i]);
            }
	  }
	}
	dpd_file2_close(&GF);
	
	// Load FIA
	dpd_file2_init(&GF, PSIF_OMP2_DENSITY, 0, ID('O'), ID('V'), "GF <O|V>");  
	dpd_file2_mat_init(&GF);
	dpd_file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < occpiA[h]; ++i){
            for(int a = 0 ; a < virtpiA[h]; ++a){
                GFockA->add(h, i, a + occpiA[h], GF.matrix[h][i][a]);
            }
	  }
	}
	dpd_file2_close(&GF);
	
	// Load Fia
	dpd_file2_init(&GF, PSIF_OMP2_DENSITY, 0, ID('o'), ID('v'), "GF <o|v>");  
	dpd_file2_mat_init(&GF);
	dpd_file2_mat_rd(&GF);
	for(int h = 0; h < nirrep_; ++h){
	  for(int i = 0 ; i < occpiB[h]; ++i){
            for(int a = 0 ; a < virtpiB[h]; ++a){
                GFockB->add(h, i, a + occpiB[h], GF.matrix[h][i][a]);
            }
	  }
	}
	dpd_file2_close(&GF);
	
	psio_->close(PSIF_OMP2_DENSITY, 1);

        // Print
	if (print_ > 2) {
	  GFockA->print();
	  GFockB->print();
	}
	
}// end if (reference == "UHF") 
//fprintf(outfile,"\n GFockmo done. \n"); fflush(outfile);

}// end GFockmo
}} // End Namespaces

