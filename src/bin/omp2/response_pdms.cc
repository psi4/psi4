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
#include <libdpd/dpd.h>
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

void OMP2Wave::response_pdms()
{   
     //fprintf(outfile,"\n response_pdms is starting... \n"); fflush(outfile);
/********************************************************************************************/
/************************** Initialize ******************************************************/
/********************************************************************************************/ 
	gamma1corrA->zero();
	gamma1corrB->zero();
	g1symmA->zero();
	g1symmB->zero();

/********************************************************************************************/
/************************** Build G intermediates *******************************************/
/********************************************************************************************/  
        timer_on("G int");
	G_int(); 
        timer_off("G int");
   
/********************************************************************************************/
/************************** (OO) Block whole ************************************************/
/********************************************************************************************/ 
        timer_on("OPDM");
	// alpha contrb.
	for (int i=nfrzc; i<nooA;i++) {
	  for (int j=nfrzc; j<nooA;j++) { 
	    int i2 = c1topitzerA[i];
	    int j2 = c1topitzerA[j];
	    int i_sym = pitzer2symblk[i2];
	    int j_sym = pitzer2symblk[j2];
	    int h1=mosym[i2];
	    int h2=mosym[j2];
	    if (h1 == h2){
	      int i3 = i_sym - frzcpi[h1];
	      int j3 = j_sym - frzcpi[h1];
	      g1symmA->set(h1,i_sym,j_sym,GooA->get(h1,i3,j3) + GooA->get(h1,j3,i3));
	    }
	  }
	}
	
	// beta contrb.
	for (int i=nfrzc; i<nooB;i++) {
	  for (int j=nfrzc; j<nooB;j++) { 
	    int i2 = c1topitzerB[i];
	    int j2 = c1topitzerB[j];
	    int i_sym = pitzer2symblk[i2];
	    int j_sym = pitzer2symblk[j2];
	    int h1=mosym[i2];
	    int h2=mosym[j2];
	    if (h1 == h2){
	      int i3 = i_sym - frzcpi[h1];
	      int j3 = j_sym - frzcpi[h1];
	      g1symmB->set(h1,i_sym,j_sym,GooB->get(h1,i3,j3) + GooB->get(h1,j3,i3));
	    }
	  }
	}
	
/********************************************************************************************/
/************************** (VV) Block whole ************************************************/
/********************************************************************************************/ 
	 // alpha contrb.
	 for (int a=nooA; a<npop;a++) {
	  for (int b=nooA; b<npop;b++) { 
	    int a2 = c1topitzerA[a];
	    int b2 = c1topitzerA[b];
	    int a_sym = pitzer2symblk[a2];
	    int b_sym = pitzer2symblk[b2];
	    int h1=mosym[a2];
	    int h2=mosym[b2];
	    if (h1 == h2){
	      int a3 = a_sym - occpiA[h1];
	      int b3 = b_sym - occpiA[h1];
	      g1symmA->set(h1,a_sym,b_sym,GvvA->get(h1,a3,b3) + GvvA->get(h1,b3,a3));
	    }
	  }
	 }
	 
	 // beta contrb.
	 for (int a=nooB; a<npop;a++) {
	  for (int b=nooB; b<npop;b++) { 
	    int a2 = c1topitzerB[a];
	    int b2 = c1topitzerB[b];
	    int a_sym = pitzer2symblk[a2];
	    int b_sym = pitzer2symblk[b2];
	    int h1=mosym[a2];
	    int h2=mosym[b2];
	    if (h1 == h2){
	      int a3 = a_sym - occpiB[h1];
	      int b3 = b_sym - occpiB[h1];
	      g1symmB->set(h1,a_sym,b_sym,GvvB->get(h1,a3,b3) + GvvB->get(h1,b3,a3));
	    }
	  }
	 }
	 
	g1symmA->scale(-0.5);
	g1symmB->scale(-0.5);
	gamma1corrA->copy(g1symmA); // correlation opdm
	gamma1corrB->copy(g1symmB); // correlation opdm
   
/********************************************************************************************/
/************************** HF Contribution *************************************************/
/********************************************************************************************/
	// alpha contrb.
	for(int h=0; h<nirreps; h++) {
	  if (occpiA[h] != 0) {
	    for (int i=0; i<occpiA[h];i++) {
	      g1symmA->add(h,i,i,1.0);
	    }
	  }
	}
	
	// beta contrb.
	for(int h=0; h<nirreps; h++) {
	  if (occpiB[h] != 0) {
	    for (int i=0; i<occpiB[h];i++) {
	      g1symmB->add(h,i,i,1.0);
	    }
	  }
	}
        timer_off("OPDM");

/********************************************************************************************/
/************************** TPDMs ***********************************************************/
/********************************************************************************************/   
        timer_on("TPDM OOVV");
	twopdm_oovv();
        timer_off("TPDM OOVV");
        timer_on("TPDM REF");
	twopdm_ref(); 
        timer_off("TPDM REF");
        timer_on("TPDM CORR OPDM");
	twopdm_corr_opdm();
        timer_off("TPDM CORR OPDM");

/********************************************************************************************/
/************************** Print ***********************************************************/
/********************************************************************************************/  
      if (print_ > 1) {
	g1symmA->print();
	g1symmB->print();
      }
      //fprintf(outfile,"\n response_pdms done... \n"); fflush(outfile);

} // end of response_pdms


void OMP2Wave::G_int()
{  
        //fprintf(outfile,"\n G_int is starting... \n"); fflush(outfile);
	GooA->zero();
	GooB->zero();
	GvvA->zero();
	GvvB->zero();

	dpdbuf4 TAA, TAB, TBB;
   dpdbuf4 TAA2, TBB2;
   dpdbuf4 TAB2_;
	//dpdbuf4 TAA, TAB, TBB, TAA2, TAB2_, TBB2;
	dpdfile2 G;
	
	psio_->open(PSIF_OMP2_DPD, PSIO_OPEN_OLD);  
        psio_->open(PSIF_OMP2_DENSITY, PSIO_OPEN_OLD);
	
	// Open amplitude files
	dpd_buf4_init(&TAA, PSIF_OMP2_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
	dpd_buf4_init(&TBB, PSIF_OMP2_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
	dpd_buf4_init(&TAB, PSIF_OMP2_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
	dpd_buf4_init(&TAA2, PSIF_OMP2_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
	dpd_buf4_init(&TBB2, PSIF_OMP2_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
	dpd_buf4_init(&TAB2_, PSIF_OMP2_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                  ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
	
	// Occupied-Occupied block
	// Alpha-Alpha spin case
	// G_IM = 1/2 \sum{N,E,F} t_IN^EF * l_EF^MN = 1/2 \sum{N,E,F} t_IN^EF * t_MN^EF
	dpd_file2_init(&G, PSIF_OMP2_DENSITY, 0, ID('O'), ID('O'), "G <O|O>");  
	dpd_contract442(&TAA, &TAA2, &G, 0, 0, 0.5, 0.0);
	
	// G_IM += \sum{n,E,f} t_In^Ef * l_Ef^Mn = \sum{N,E,F} t_In^Ef * t_Mn^Ef
	dpd_contract442(&TAB, &TAB2_, &G, 0, 0, 1.0, 1.0);
	dpd_file2_close(&G);
	
	
	// Beta-Beta spin case
	// G_im = 1/2 \sum{n,e,f} t_in^ef * l_ef^mn = 1/2 \sum{n,e,f} t_in^ef * t_mn^ef
	dpd_file2_init(&G, PSIF_OMP2_DENSITY, 0, ID('o'), ID('o'), "G <o|o>");  
	dpd_contract442(&TBB, &TBB2, &G, 0, 0, 0.5, 0.0);
	
	// G_im  += \sum{N,e,F} t_Ni^Fe * l_Fe^Nm = \sum{N,e,F} t_Ni^Fe * t_Nm^Fe
	dpd_contract442(&TAB, &TAB2_, &G, 1, 1, 1.0, 1.0);
	dpd_file2_close(&G);
	
	
	
	// Virtual-Virtual block
	// Alpha-Alpha spin case
	// G_EA = -1/2 \sum{M,N,F} t_MN^AF * l_EF^MN = -1/2 \sum{M,N,F} t_MN^AF * t_MN^EF
	dpd_file2_init(&G, PSIF_OMP2_DENSITY, 0, ID('V'), ID('V'), "G <V|V>");  
	dpd_contract442(&TAA, &TAA2, &G, 2, 2, -0.5, 0.0); 
	
	// G_EA += - \sum{M,n,f} t_Mn^Af * l_Ef^Mn = - \sum{M,n,f} t_Mn^Af * t_Mn^Ef
	dpd_contract442(&TAB, &TAB2_, &G, 2, 2, -1.0, 1.0); 
	dpd_file2_close(&G);
	
	// Beta-Beta spin case
	// G_ea = -1/2 \sum{m,n,f} t_mn^af * l_ef^mn = -1/2 \sum{m,n,f} t_mn^af * t_mn^ef
	dpd_file2_init(&G, PSIF_OMP2_DENSITY, 0, ID('v'), ID('v'), "G <v|v>");  
	dpd_contract442(&TBB, &TBB2, &G, 2, 2, -0.5, 0.0); 
	
	// G_ea += - \sum{M,n,F} t_Mn^Fa * l_Fe^Mn = - \sum{M,n,F} t_Mn^Fa * t_Mn^Fe
	dpd_contract442(&TAB, &TAB2_, &G, 3, 3, -1.0, 1.0); 
	dpd_file2_close(&G);
	
	
	// Close amplitude files
	dpd_buf4_close(&TAA);
	dpd_buf4_close(&TBB);
	dpd_buf4_close(&TAB);
	dpd_buf4_close(&TAA2);
	dpd_buf4_close(&TBB2);
	dpd_buf4_close(&TAB2_);
	
	
	
	// Load dpd_file2 to Matrix (Goo)
	// Alpha-Alpha spin case
	dpd_file2_init(&G, PSIF_OMP2_DENSITY, 0, ID('O'), ID('O'), "G <O|O>");  
	dpd_file2_mat_init(&G);
	dpd_file2_mat_rd(&G);
	for(int h = 0; h < nirreps; ++h){
	  for(int i = 0 ; i < aoccpiA[h]; ++i){
            for(int j = 0 ; j < aoccpiA[h]; ++j){
                GooA->set(h, i, j, G.matrix[h][i][j]);
            }
	  }
	}
	dpd_file2_close(&G);
	
	// Beta-Beta spin case
	dpd_file2_init(&G, PSIF_OMP2_DENSITY, 0, ID('o'), ID('o'), "G <o|o>");  
	dpd_file2_mat_init(&G);
	dpd_file2_mat_rd(&G);
	for(int h = 0; h < nirreps; ++h){
	  for(int i = 0 ; i < aoccpiB[h]; ++i){
            for(int j = 0 ; j < aoccpiB[h]; ++j){
                GooB->set(h, i, j, G.matrix[h][i][j]);
            }
	  }
	}
	dpd_file2_close(&G);
	
	
	
	// Load dpd_file2 to Matrix (Gvv)
	// Alpha-Alpha spin case
	dpd_file2_init(&G, PSIF_OMP2_DENSITY, 0, ID('V'), ID('V'), "G <V|V>"); 
	dpd_file2_mat_init(&G);
	dpd_file2_mat_rd(&G);
	for(int h = 0; h < nirreps; ++h){
	  for(int i = 0 ; i < avirtpiA[h]; ++i){
            for(int j = 0 ; j < avirtpiA[h]; ++j){
                GvvA->set(h, i, j, G.matrix[h][i][j]);
            }
	  }
	}
	dpd_file2_close(&G);
	
	// Beta-Beta spin case
	dpd_file2_init(&G, PSIF_OMP2_DENSITY, 0, ID('v'), ID('v'), "G <v|v>");  
	dpd_file2_mat_init(&G);
	dpd_file2_mat_rd(&G);
	for(int h = 0; h < nirreps; ++h){
	  for(int i = 0 ; i < avirtpiB[h]; ++i){
            for(int j = 0 ; j < avirtpiB[h]; ++j){
                GvvB->set(h, i, j, G.matrix[h][i][j]);
            }
	  }
	}
	dpd_file2_close(&G);
	
	psio_->close(PSIF_OMP2_DPD, 1);  
        psio_->close(PSIF_OMP2_DENSITY, 1);

	if (print_ > 1) {
	  GooA->print();
	  GooB->print();
	  GvvA->print();
	  GvvB->print();
	}
	
	//printf(outfile,"\n G_int done... \n"); fflush(outfile);

} // end of G_int



void OMP2Wave::twopdm_oovv()
{      
    dpdbuf4 G, T;
    
    psio_->open(PSIF_OMP2_DPD, PSIO_OPEN_OLD);  
    psio_->open(PSIF_OMP2_DENSITY, PSIO_OPEN_OLD);
    
    // Alpha-Alpha spin case
    dpd_buf4_init(&T, PSIF_OMP2_DPD, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "T2_1 <OO|VV>");
    dpd_buf4_copy(&T, PSIF_OMP2_DENSITY, "TPDM <OO|VV>");
    dpd_buf4_close(&T);
    
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,O]"), ID("[V,V]"),
                  ID("[O,O]"), ID("[V,V]"), 0, "TPDM <OO|VV>");
    dpd_buf4_scm(&G, 0.25);
    dpd_buf4_close(&G);
    
    // Beta-Beta spin case
    dpd_buf4_init(&T, PSIF_OMP2_DPD, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "T2_1 <oo|vv>");
    dpd_buf4_copy(&T, PSIF_OMP2_DENSITY, "TPDM <oo|vv>");
    dpd_buf4_close(&T);
    
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[o,o]"), ID("[v,v]"),
                  ID("[o,o]"), ID("[v,v]"), 0, "TPDM <oo|vv>");
    dpd_buf4_scm(&G, 0.25);
    dpd_buf4_close(&G);
    
    // Alpha-Beta spin case
    dpd_buf4_init(&T, PSIF_OMP2_DPD, 0, ID("[O,o]"), ID("[V,v]"),
                 ID("[O,o]"), ID("[V,v]"), 0, "T2_1 <Oo|Vv>");
    dpd_buf4_copy(&T, PSIF_OMP2_DENSITY, "TPDM <Oo|Vv>");
    dpd_buf4_close(&T);
    
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,o]"), ID("[V,v]"),
                 ID("[O,o]"), ID("[V,v]"), 0, "TPDM <Oo|Vv>");
    dpd_buf4_scm(&G, 0.25);
    dpd_buf4_close(&G);
    
    
    psio_->close(PSIF_OMP2_DPD, 1);  
    psio_->close(PSIF_OMP2_DENSITY, 1);

} // end of twopdm_oovv



void OMP2Wave::twopdm_ref()
{
    dpdbuf4 G;
    
    psio_->open(PSIF_OMP2_DENSITY, PSIO_OPEN_OLD);
    
    // Alpha-Alpha spin case
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,O]"), ID("[O,O]"),
                  ID("[O,O]"), ID("[O,O]"), 0, "TPDM <OO|OO>");
      for(int h = 0; h < nirreps; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
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
        dpd_buf4_mat_irrep_wrt(&G, h);
        dpd_buf4_mat_irrep_close(&G, h);
    }    
    dpd_buf4_close(&G);
    
    
    // Beta-Beta spin case
    dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[o,o]"), ID("[o,o]"),
                  ID("[o,o]"), ID("[o,o]"), 0, "TPDM <oo|oo>");    
      for(int h = 0; h < nirreps; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
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
        dpd_buf4_mat_irrep_wrt(&G, h);
        dpd_buf4_mat_irrep_close(&G, h);
    }    
    dpd_buf4_close(&G);
    
    
    // Alpha-Beta spin case
     dpd_buf4_init(&G, PSIF_OMP2_DENSITY, 0, ID("[O,o]"), ID("[O,o]"),
                 ID("[O,o]"), ID("[O,o]"), 0, "TPDM <Oo|Oo>"); 
      for(int h = 0; h < nirreps; ++h){
        dpd_buf4_mat_irrep_init(&G, h);
        for(int ij = 0; ij < G.params->rowtot[h]; ++ij){
            int i = G.params->roworb[h][ij][0];
            int j = G.params->roworb[h][ij][1];
            for(int kl = 0; kl < G.params->coltot[h]; ++kl){
                int k = G.params->colorb[h][kl][0];
                int l = G.params->colorb[h][kl][1];
		if (i == k && j == l) G.matrix[h][ij][kl] += 0.25;
            }
        }
        dpd_buf4_mat_irrep_wrt(&G, h);
        dpd_buf4_mat_irrep_close(&G, h);
    }    
    dpd_buf4_close(&G);
    
    psio_->close(PSIF_OMP2_DENSITY, 1);
    
} // end of twopdm_ref



}} // End Namespaces

